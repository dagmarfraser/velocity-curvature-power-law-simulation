function out = checkGridTempoAtFold_v004(dbFile)
% checkGridTempoAtFold_v004  V0 gate of SPEC_TempoFoldSweep_v001.
%
% Does the v058 grid's zero-noise fold sit at a constant REALISED tempo?
% Reads realised trajectory duration from the v058 DB (results.duration =
% Toolchain_func_v032 genDuration, 10 orbits: Toolchain_caller_v058.m L203),
% so f0 = ORBITS/duration needs no assumption about grid geometry or units.
%
% For each (pipeline, fs, VGF) at sigma = 0, alpha = 0: the beta_gen -> beta_rec
% curve, its peak, the drop to the last node, and f0 at the peak.
% P1 predicts f0 at the peak is ~constant across VGF (spec section 3).
%
% Run on BlueBEAR next to the DB (aggregation is pushed into SQLite; the
% result is ~5,500 rows). Writes results/gridTempoAtFold_v004.mat.
%
% USAGE:  out = checkGridTempoAtFold_v004();          % default v058 DB path
%         out = checkGridTempoAtFold_v004("/path/powerlaw_debug_v058.db");
% v002: default DB is powerlaw_debug_v058.db (the v058 production DB; SESSION_LOG
%       Part 2 naming note), multiverse name as fallback; REAL aggregates wrapped in
%       COALESCE because MATLAB sqlite fetch throws on NULL in REAL columns.
% v003: tempo looked up per (beta_gen, VGF, fs) coordinate from any pipeline, not per
%       curve (v002 L96 errored where IRLS lacked the node at sigma = 0); missing
%       coordinates reported by pipeline x fs (length guard vs pipeline-only);
%       per-curve nNodes and lastBeta so truncated curves are visible.
% v004: node match tolerance 1e-4 (DB stores beta_gen to 5 significant figures; v003
%       printed NaN for the 1/3 and 2/3 tempo lines); analytic check of realised tempo
%       against the grid ellipse, f0 = VGF / closed-integral(kappa^beta ds); fixed-beta
%       counterfactual for the peak-tempo constancy (P1).
% Fraser, D.S. (2026)  v004

    arguments
        dbFile (1,1) string = defaultDB()
    end
    ORBITS  = 10;         % Toolchain_caller_v058.m L203 (cfg.orbitCount)
    DUR_TOL = 1e-9;       % duration must not depend on pipeline or repetition
    GRID_A  = 235;  GRID_B = 109;   % least-squares ellipse fit to the v058 nu = 2 shape, canvas px (extents 240 x 109.85; Finding #236, checkGridShapeUnits_v001)
    if ~isfile(dbFile)
        error("gridTempo:noDB", "%s", "DB not found: " + dbFile);
    end
    conn = sqlite(dbFile, "readonly");
    cleanupConn = onCleanup(@() close(conn));

    % Dedup exactly as computePerCoordinateSEM_v2_001 (MIN(rowid) per config_id),
    % restricted first to sigma = 0, alpha = 0 so the CTE never scans the full grid.
    sql = strjoin([
        "WITH sel AS ("
        "  SELECT config_id FROM param_configs"
        "  WHERE noise_magnitude = 0 AND CAST(noise_type AS REAL) = 0"
        "), dedup AS ("
        "  SELECT MIN(r.rowid) AS rid, r.config_id"
        "  FROM results r INNER JOIN sel s ON r.config_id = s.config_id"
        "  WHERE r.success = 1 AND r.beta IS NOT NULL"
        "  GROUP BY r.config_id"
        ")"
        "SELECT p.filter_type, p.regress_type, p.generated_beta, p.vgf_value,"
        "  p.sampling_rate,"
        "  COALESCE(AVG(r.beta), -999) AS betaRec, COALESCE(AVG(r.duration), -999) AS dur,"
        "  COALESCE(MIN(r.duration), -999) AS durMin, COALESCE(MAX(r.duration), -999) AS durMax,"
        "  COUNT(*) AS nReps"
        "FROM param_configs p"
        "INNER JOIN dedup d ON p.config_id = d.config_id"
        "INNER JOIN results r ON r.rowid = d.rid"
        "GROUP BY p.filter_type, p.regress_type, p.generated_beta, p.vgf_value, p.sampling_rate"
    ], " ");
    fprintf("Querying %s (sigma = 0, alpha = 0) ... ", dbFile); tic;
    T = fetch(conn, sql);
    fprintf("%d rows in %.1f s\n", height(T), toc);
    if isempty(T), error("gridTempo:noRows", "Query returned 0 rows."); end
    for v = ["filter_type" "regress_type" "generated_beta" "vgf_value" "sampling_rate" ...
             "betaRec" "dur" "durMin" "durMax" "nReps"]
        T.(v) = double(T.(v));
    end
    for v = ["betaRec" "dur" "durMin" "durMax"]
        T.(v)(T.(v) == -999) = NaN;             % COALESCE sentinel back to NaN
    end
    nBad = sum(isnan(T.betaRec) | isnan(T.dur) | isnan(T.durMin) | isnan(T.durMax));
    if nBad > 0
        error("gridTempo:nullValues", "%s", sprintf("%d aggregated rows have NULL beta or duration.", nBad));
    end
    T.pipeline = pipelineLabel(T.filter_type, T.regress_type);

    % Guard: realised duration is a property of (beta_gen, VGF, fs) only.
    G = groupsummary(T, ["generated_beta" "vgf_value" "sampling_rate"], "range", "dur");
    badRep  = max(T.durMax - T.durMin);
    badPipe = max(G.range_dur);
    if badRep > DUR_TOL || badPipe > DUR_TOL
        error("gridTempo:durationVaries", "%s", sprintf( ...
            "duration varies across reps (max %.3g s) or pipelines (max %.3g s); " + ...
            "f0 = %d/duration is not a property of the coordinate.", badRep, badPipe, ORBITS));
    end
    T.f0 = ORBITS ./ T.dur;

    % Coordinate-level tempo: pipeline-free (guard above), so any pipeline present
    % at (beta_gen, VGF, fs) supplies it.
    F = groupsummary(T, ["generated_beta" "vgf_value" "sampling_rate"], "mean", "f0");
    F = renamevars(F, "mean_f0", "f0");

    % Missing coordinates: full crossing of the observed levels vs what came back.
    bgN = unique(T.generated_beta); vgN = unique(T.vgf_value); fsN = unique(T.sampling_rate);
    pipes = unique(T.pipeline);
    fprintf("\nRows: %d of %d expected (%d beta x %d VGF x %d fs x %d pipelines)\n", height(T), ...
        numel(bgN)*numel(vgN)*numel(fsN)*numel(pipes), numel(bgN), numel(vgN), numel(fsN), numel(pipes));
    fprintf("Missing (beta, VGF) cells: coord-absent = no pipeline has it (length guard);\n");
    fprintf("pipeline-only = other pipelines have it (NULL beta, e.g. IRLS at sigma = 0)\n");
    fprintf("%-10s %4s %13s %14s\n", "pipeline", "fs", "coord-absent", "pipeline-only");
    for p = pipes'
        for fsv = fsN'
            nAll  = sum(F.sampling_rate == fsv);
            nHave = sum(T.pipeline == p & T.sampling_rate == fsv);
            fprintf("%-10s %4d %13d %14d\n", p, fsv, numel(bgN)*numel(vgN) - nAll, nAll - nHave);
        end
    end

    % Per (pipeline, fs, VGF): peak of the curve and realised tempo there.
    keys = unique(T(:, ["pipeline" "sampling_rate" "vgf_value"]), "rows");
    n = height(keys);
    [pkBeta, pkRec, drop, f0Pk, f0At13, f0At23, lastBeta, nNodes] = deal(NaN(n,1));
    for k = 1:n
        c = T(T.pipeline == keys.pipeline(k) & T.sampling_rate == keys.sampling_rate(k) & ...
              T.vgf_value == keys.vgf_value(k), :);
        c = sortrows(c, "generated_beta");
        nNodes(k)   = height(c);
        lastBeta(k) = c.generated_beta(end);
        [pkRec(k), ip] = max(c.betaRec);
        pkBeta(k) = c.generated_beta(ip);
        drop(k)   = pkRec(k) - c.betaRec(end);    % drop > 0 implies an interior peak
        f0Pk(k)   = c.f0(ip);
        f0At13(k) = coordTempo(F, 1/3, keys.vgf_value(k), keys.sampling_rate(k));
        f0At23(k) = coordTempo(F, 2/3, keys.vgf_value(k), keys.sampling_rate(k));
    end
    S = [keys, table(nNodes, lastBeta, pkBeta, pkRec, drop, f0Pk, f0At13, f0At23)];

    % P1 readout: spread of peak location in beta vs in realised tempo, across VGF.
    fprintf("\nPer pipeline x fs, across VGF nodes (folds: interior peak, drop > 0.01):\n");
    fprintf("%-10s %4s  %-16s %-17s %-10s %s\n", "pipeline", "fs", "peak beta_gen", "f0 at peak (Hz)", "nFold/nVGF", "nTruncated");
    grp = unique(S(:, ["pipeline" "sampling_rate"]), "rows");
    for g = 1:height(grp)
        s = S(S.pipeline == grp.pipeline(g) & S.sampling_rate == grp.sampling_rate(g), :);
        nTr = sum(s.lastBeta < max(bgN) - 1e-9);
        f = s(s.drop > 0.01, :);
        if isempty(f)
            fprintf("%-10s %4d  no fold at any VGF node%31s%d\n", grp.pipeline(g), grp.sampling_rate(g), "", nTr);
            continue
        end
        fprintf("%-10s %4d  %5.3f - %5.3f    %5.2f - %5.2f      %2d/%-7d %d\n", grp.pipeline(g), ...
            grp.sampling_rate(g), min(f.pkBeta), max(f.pkBeta), min(f.f0Pk), max(f.f0Pk), ...
            height(f), height(s), nTr);
    end
    fprintf("\nRealised tempo across VGF (coordinate-level; NaN = coordinate absent):\n");
    fprintf("  beta_gen = 1/3: %.2f - %.2f Hz (%d absent);  beta_gen = 2/3: %.2f - %.2f Hz (%d absent)\n", ...
        min(S.f0At13), max(S.f0At13), sum(isnan(S.f0At13)), min(S.f0At23), max(S.f0At23), sum(isnan(S.f0At23)));
    fprintf("  (spec section 1 predicts ~0.49-1.78 and ~2.74-10.0 Hz if the stated grid geometry is right)\n");

    % Analytic check: realised tempo vs the grid ellipse (GRID_A, GRID_B in canvas units,
    % the same units as VGF; confirmed by this comparison, not assumed).
    F.f0Analytic = arrayfun(@(v, b) ellipseTempo(v, b, GRID_A, GRID_B), F.vgf_value, F.generated_beta);
    relErr = abs(F.f0 ./ F.f0Analytic - 1);
    fprintf("\nAnalytic tempo (ellipse a = %g, b = %g): max |realised/analytic - 1| = %.3g over %d coordinates (beta_gen > 0)\n", ...
        GRID_A, GRID_B, max(relErr(F.generated_beta > 0)), sum(F.generated_beta > 0));

    % P1 counterfactual: if the fold sat at a fixed beta_gen, tempo at the peak would span
    % the full VGF range at that beta. Compare with the observed spread of tempo at the peak.
    fsRef = 120;  bRef = F.generated_beta(abs(F.generated_beta - 0.5333) < 1e-3 & F.sampling_rate == fsRef);
    fx = F.f0(F.sampling_rate == fsRef & abs(F.generated_beta - bRef(1)) < 1e-9);
    s  = S(S.pipeline == "SG-IRLS" & S.sampling_rate == fsRef & S.drop > 0.01, :);
    fprintf("P1 counterfactual (fs %d): tempo at fixed beta_gen %.4f spans %.2f-%.2f Hz (ratio %.2f);\n", ...
        fsRef, bRef(1), min(fx), max(fx), max(fx)/min(fx));
    fprintf("  observed SG-IRLS tempo at the peak spans %.2f-%.2f Hz (ratio %.2f) while peak beta_gen spans %.3f-%.3f\n", ...
        min(s.f0Pk), max(s.f0Pk), max(s.f0Pk)/min(s.f0Pk), min(s.pkBeta), max(s.pkBeta));

    out = struct("perCoordinate", T, "perCurve", S, "perNode", F, "gridAxes", [GRID_A GRID_B], "dbFile", dbFile, "orbits", ORBITS, ...
                 "runDate", string(datetime("now")));
    outFile = fullfile(projectRoot(), "results", "gridTempoAtFold_v004.mat");
    if ~isfolder(fileparts(outFile))
        error("gridTempo:outDir", "%s", "Missing folder: " + fileparts(outFile));
    end
    save(outFile, "out", "-v7.3");
    fprintf("\nSaved: %s\n", outFile);
end

function v = coordTempo(F, b0, vgf, fsv)
    % Tempo at one (beta_gen, VGF, fs) coordinate; NaN (counted by the caller) if absent.
    i = find(abs(F.generated_beta - b0) < 1e-4 & F.vgf_value == vgf & F.sampling_rate == fsv, 1);
    if isempty(i), v = NaN; else, v = F.f0(i); end
end

function lab = pipelineLabel(ft, rt)
    % DB codes as computePerCoordinateSEM_v2_001 L109-117; unknown codes fail loud.
    fMap = containers.Map([2 6], {'BWFD','SG'});
    rMap = containers.Map([3 4 5], {'OLS','LMLS','IRLS'});
    uf = setdiff(unique(ft), cell2mat(keys(fMap)));  ur = setdiff(unique(rt), cell2mat(keys(rMap)));
    if ~isempty(uf) || ~isempty(ur)
        error("gridTempo:codes", "%s", "Unknown filter/regress codes: " + mat2str(uf) + " / " + mat2str(ur));
    end
    lab = strings(numel(ft), 1);
    for k = 1:numel(ft), lab(k) = string(fMap(ft(k))) + "-" + string(rMap(rt(k))); end
end

function r = projectRoot()
    r = fileparts(fileparts(mfilename("fullpath")));   % src/.. = project root
end

function p = defaultDB()
    % v058 production DB is powerlaw_debug_v058.db despite Toolchain_caller_v058.m L98
    % (SESSION_LOG, extractCompressionData naming note); try both, fail loud if neither.
    cands = fullfile(projectRoot(), "results", ["powerlaw_debug_v058.db" "powerlaw_multiverse_v058.db"]);
    i = find(isfile(cands), 1);
    if isempty(i)
        error("gridTempo:noDefaultDB", "%s", "Neither default DB found: " + strjoin(cands, " | "));
    end
    p = string(cands(i));
end

function f0 = ellipseTempo(vgf, beta, a, b)
    % Orbital frequency of an ellipse traced at speed v = vgf * kappa^-beta.
    kap = @(p) a*b ./ (a^2*sin(p).^2 + b^2*cos(p).^2).^1.5;
    dsp = @(p) sqrt(a^2*sin(p).^2 + b^2*cos(p).^2);
    f0  = vgf / integral(@(p) kap(p).^beta .* dsp(p), 0, 2*pi, "RelTol", 1e-10);
end
