function coordTable = computePerCoordinateSEM_VGF_v001(dbFile)
% computePerCoordinateSEM_VGF_v001  Per-coordinate SEM and bias for VGF recovery.
%
% VGF-outcome sibling of computePerCoordinateSEM_v2_001.m. Same v058
% database, same (pipeline, betaGen, VGF, fs, alpha, sigma) coordinate
% grid, same SQL aggregation strategy -- the only change is that bias is
% computed against VGF_gen (p.vgf_value) instead of beta_gen
% (p.generated_beta). The two AVG() terms needed for this were already
% being pulled by the beta-version query (as meanVGFrec, unused
% downstream); this version turns them into the actual outcome instead
% of a side reading.
%
% Built to answer the pre-reg's own VGF gating condition (Sections 3.2,
% 4.2): "Assessment of VGF consistency (low SEM) ... conducted
% conditional on adequate beta parameter recovery (SEM < 0.011)." This
% script produces the VGF side of that comparison; checkVGFSEMCentroid_v001.m
% reads this alongside the existing perCoordinateSEM_v2_001.mat (beta) and
% applies the gate.
%
% No beta-calibrated MDC/SEM-stratum thresholds are applied here.
% docs/LMM_VGF_Desiderata_v003.md (D6) already found that the beta
% threshold (MDC/2.77 = 0.0108) is a unit mismatch when applied literally
% to a mm/s outcome -- VGF has no established clinical-significance
% threshold of its own. This script reports raw meanBias/sem in native
% VGF units (mm/s) and leaves stratification/interpretation to the
% consuming script or the manuscript text, rather than hardcoding a
% threshold that would silently misapply beta's own MDC to a different
% outcome.
%
% USAGE:
%   T = computePerCoordinateSEM_VGF_v001()                  % default: v058 DB
%   T = computePerCoordinateSEM_VGF_v001('path/to/db.db')
%
% OUTPUT:
%   coordTable - table with columns:
%     pipeline, filterType, regressType, betaGen, VGF, fs, alpha, sigma,
%     meanBias, sem, nReps, meanVGFrec
%   (meanBias/sem are VGF bias/SEM here, i.e. VGF_gen - VGF_rec; the
%   column names are kept identical to the beta version so downstream
%   centroid-lookup code can be reused unchanged, just pointed at a
%   different .mat file.)
%
% PREREG REF: Section 3.2 (VGF Measurement Error Considerations),
%             Section 4.2 (Secondary Validation: VGF Parameter Recovery)
% Author: Fraser, D.S. (2026)

    arguments
        dbFile (1,1) string = defaultDBpath()
    end

    %% Connect and pull data
    fprintf("=== Per-Coordinate SEM Computation, VGF outcome (v058) ===\n");
    fprintf("  Database: %s\n", dbFile);
    if ~isfile(dbFile)
        error("computePerCoordinateSEM_VGF:DBnotFound", "Database not found: %s", dbFile);
    end

    conn = sqlite(dbFile);
    cleanupConn = onCleanup(@() close(conn));

    % Same dedup/aggregation strategy as computePerCoordinateSEM_v2_001.m
    % (see that file's own comments for the SEM-from-SQL-aggregates trick
    % and the crash-restart dedup rationale) -- only the bias expression
    % changes, from (generated_beta - r.beta) to (vgf_value - r.vgf).
    sql = strjoin([
        "WITH dedup AS ("
        "  SELECT MIN(rowid) AS rid, config_id"
        "  FROM results"
        "  WHERE success = 1 AND vgf IS NOT NULL"
        "  GROUP BY config_id"
        ")"
        "SELECT p.filter_type, p.regress_type,"
        "  p.generated_beta, p.vgf_value, p.sampling_rate,"
        "  CAST(p.noise_type AS REAL) AS noise_type, p.noise_magnitude,"
        "  AVG(p.vgf_value - r.vgf) AS meanBias,"
        "  AVG((p.vgf_value - r.vgf)*(p.vgf_value - r.vgf)) AS meanBiasSq,"
        "  COUNT(*) AS nReps,"
        "  AVG(r.vgf) AS meanVGFrec"
        "FROM param_configs p"
        "INNER JOIN dedup d ON p.config_id = d.config_id"
        "INNER JOIN results r ON r.rowid = d.rid"
        "GROUP BY p.filter_type, p.regress_type,"
        "  p.generated_beta, p.vgf_value, p.sampling_rate,"
        "  p.noise_type, p.noise_magnitude"
    ], " ");

    fprintf("  Fetching aggregated data (single SQL GROUP BY)...");
    tic;
    coordTable = fetch(conn, sql);
    elapsed = toc;
    fprintf(" %d coordinate groups in %.1f s\n", height(coordTable), elapsed);

    if isempty(coordTable)
        error("computePerCoordinateSEM_VGF:NoData", "Query returned 0 rows. Is the database populated?");
    end

    % Guard: dedup CTE should guarantee nReps <= repeatTrial (5).
    nOver = sum(double(coordTable.nReps) > 5);
    if nOver > 0
        error("computePerCoordinateSEM_VGF:ExcessReps", "%s", ...
            sprintf("%d coordinates have nReps > 5 after dedup. " + ...
            "The CTE deduplication did not work as expected.", nOver));
    end

    % Sample SD from E[X^2] and E[X]^2, Bessel-corrected -- identical
    % arithmetic to the beta version, applied to the VGF bias instead.
    popVar = coordTable.meanBiasSq - coordTable.meanBias.^2;
    popVar = max(popVar, 0);  % guard floating-point rounding below zero
    n = double(coordTable.nReps);
    besselFactor = n ./ max(n - 1, 1);
    coordTable.sem = sqrt(besselFactor .* popVar);
    coordTable.meanBiasSq = [];

    %% Build pipeline labels (identical scheme to the beta version)
    ftLabels = strings(height(coordTable), 1);
    ftLabels(coordTable.filter_type == 2) = "BWFD";
    ftLabels(coordTable.filter_type == 6) = "SG";
    ftLabels(ftLabels == "") = "F" + string(coordTable.filter_type(ftLabels == ""));

    rtLabels = strings(height(coordTable), 1);
    rtLabels(coordTable.regress_type == 3) = "OLS";
    rtLabels(coordTable.regress_type == 4) = "LMLS";
    rtLabels(coordTable.regress_type == 5) = "IRLS";
    rtLabels(rtLabels == "") = "R" + string(coordTable.regress_type(rtLabels == ""));

    coordTable.pipeline = categorical(ftLabels + "-" + rtLabels);

    coordTable = renamevars(coordTable, "generated_beta", "betaGen");
    coordTable = renamevars(coordTable, "vgf_value", "VGF");
    coordTable = renamevars(coordTable, "sampling_rate", "fs");
    coordTable = renamevars(coordTable, "noise_type", "alpha");
    coordTable = renamevars(coordTable, "noise_magnitude", "sigma");
    coordTable = renamevars(coordTable, "filter_type", "filterType");
    coordTable = renamevars(coordTable, "regress_type", "regressType");

    %% Flag incomplete coordinates
    nIncomplete = sum(coordTable.nReps < 5);
    if nIncomplete > 0
        fprintf("  WARNING: %d coordinates have < 5 reps (SEM unreliable)\n", nIncomplete);
    end
    nNanSEM = sum(isnan(coordTable.sem));
    if nNanSEM > 0
        fprintf("  WARNING: %d coordinates have NaN SEM (single rep?)\n", nNanSEM);
    end

    %% Summary (no adequacy stratum -- see header note on why not)
    nTotal = height(coordTable);
    fprintf("\n  === VGF BIAS/SEM SUMMARY (native VGF units, mm/s) ===\n");
    fprintf("  Total coordinates: %d\n", nTotal);
    fprintf("  meanBias: median=%.3f  [%.3f, %.3f] (5th-95th pctile)\n", ...
        median(coordTable.meanBias, "omitnan"), ...
        prctile(coordTable.meanBias, 5), prctile(coordTable.meanBias, 95));
    fprintf("  sem:      median=%.3f  [%.3f, %.3f] (5th-95th pctile)\n", ...
        median(coordTable.sem, "omitnan"), ...
        prctile(coordTable.sem, 5), prctile(coordTable.sem, 95));

    fprintf("\n  === PER-PIPELINE MEDIAN BIAS/SEM ===\n");
    pipes = categories(coordTable.pipeline);
    for pIdx = 1:numel(pipes)
        mask = coordTable.pipeline == pipes{pIdx};
        fprintf("  %-12s  meanBias=%8.3f  sem=%8.3f\n", pipes{pIdx}, ...
            median(coordTable.meanBias(mask), "omitnan"), ...
            median(coordTable.sem(mask), "omitnan"));
    end

    %% Save
    outFile = "perCoordinateSEM_VGF_v001.mat";
    save(outFile, "coordTable", "-v7.3");
    fprintf("\n  Saved: %s (%d coordinates)\n", outFile, nTotal);

end

function p = defaultDBpath()
% Resolve default database path relative to this file's location.
    thisDir = fileparts(mfilename("fullpath"));
    p = fullfile(thisDir, "..", "results", "powerlaw_debug_v058.db");
end
