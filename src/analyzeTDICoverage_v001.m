function results = analyzeTDICoverage_v001(opts)
% analyzeTDICoverage_v001  Prereg S7.5 candidates TDI and Coverage Probability.
%
% Prereg (prereg_v101.docx, Notation table + S7.5 body text, verbatim):
%   TDI  "Total Deviation Index (Lin, 2000); the boundary within which a
%        specified proportion (e.g., 90%) of prediction/observation
%        differences fall. Considered as a supplementary agreement metric
%        to CCC where distributional assumptions are uncertain."
%   CP   "Coverage Probability; proportion of predictions falling within
%        MDC (0.03) of observed values. Directly aligned with the clinical
%        detection framework."
% Both hedged ("potentially including"), same status as Bland-Altman, ICC,
% pattern preservation. Previously declined in D6/CCC_MEASURES_REFERENCE.md
% measure 7's note as "distributional assumptions unavailable per trial" --
% that objection applies to a PER-TRIAL construction (n=15 pairs), not to
% the pooled (trial x pipeline-pair) construction Bland-Altman (Finding
% #126) already used for the same deltaBetaPred/deltaBetaObs population;
% pooling sidesteps it, so both metrics are built here on that same pooled
% population, not attempted per-trial.
%
% CP(MDC) turns out to already be sitting in `blandAltman_v001.mat` as
% `pctWithinMDC` -- the prereg's CP is literally the same quantity Finding
% #126 already computed, just not labelled or cross-referenced as CP. This
% script recomputes it independently from `constellationMetrics_v002.mat`
% (not by reading Finding #126's saved output) and hard-errors if the two
% do not agree, before adding TDI and a second CP boundary (2xMDC, for
% context only, not a prereg-specified figure).
%
% TDI is reported two ways given the prereg's own "distributional
% assumptions uncertain" hedge, tested rather than assumed:
%   TDI90_param = Lin (2000)'s normal-theory approximation,
%                 z_(0.95) * sqrt(bias^2 + sdDiff^2)
%   TDI90_emp   = the empirical 90th percentile of |diff| (no distributional
%                 assumption; the honest fallback where Bland-Altman
%                 (Finding #126) already showed real proportional bias,
%                 i.e. non-normal differences, for Hickman PLAC/HALO)
% The two are reported side by side, not reconciled to one number, because
% their divergence is itself the answer to whether the normal-theory
% shortcut was safe to use.
%
% INPUTS (name-value)
%   MetricsMat      path to constellationMetrics_v002.mat
%   BlandAltmanMat  path to blandAltman_v001.mat (cross-check only)
%   NoiseModel      default "shaped_xu"
%   Proportion      TDI/CP target proportion, default 0.90 (prereg's own example)
%   SaveMat/SaveTxt/SavePng   default true
%
% OUTPUT
%   results.perDataset(:).{dataset,n,bias,sdDiff,cpMDC,cp2xMDC,
%                           tdiParam,tdiEmp}
%
% Fraser, D.S. (2026)  v001
    arguments
        opts.MetricsMat (1,1) string = fullfile(fileparts(mfilename("fullpath")), "constellationMetrics_v002.mat")
        opts.BlandAltmanMat (1,1) string = fullfile(fileparts(mfilename("fullpath")), "blandAltman_v001.mat")
        opts.NoiseModel (1,1) string = "shaped_xu"
        opts.Proportion (1,1) double = 0.90
        opts.SaveMat (1,1) logical = true
        opts.SaveTxt (1,1) logical = true
        opts.SavePng (1,1) logical = true
    end

    if ~isfile(opts.MetricsMat)
        error("analyzeTDICoverage:missingInput", "%s", "MetricsMat not found: " + opts.MetricsMat);
    end
    if ~isfile(opts.BlandAltmanMat)
        error("analyzeTDICoverage:missingInput", "%s", "BlandAltmanMat not found: " + opts.BlandAltmanMat);
    end
    S = load(opts.MetricsMat, "results");
    src = S.results;
    if ~isfield(src.byNoiseModel, opts.NoiseModel)
        error("analyzeTDICoverage:missingModel", "%s", "NoiseModel '" + opts.NoiseModel + "' not present in " + opts.MetricsMat);
    end
    pd = src.byNoiseModel.(opts.NoiseModel).perDataset;
    MDC = src.config.MDC;
    DEG_THRESH = src.config.DEG_THRESH;

    Sba = load(opts.BlandAltmanMat, "results");
    baPerDataset = Sba.results.perDataset;

    zScore = norminv((1 + opts.Proportion) / 2);

    nDatasets = numel(pd);
    perDataset = struct("dataset", {}, "n", {}, "bias", {}, "sdDiff", {}, ...
        "cpMDC", {}, "cp2xMDC", {}, "tdiParam", {}, "tdiEmp", {});

    for i = 1:nDatasets
        [diffVec, validRows] = pairedDiffs_local(pd(i).deltaBetaPred, pd(i).deltaBetaObs, MDC, DEG_THRESH);
        if numel(validRows) ~= pd(i).Nvalid
            error("analyzeTDICoverage:validityMismatch", "%s", ...
                pd(i).dataset + ": reconstructed valid-row count (" + numel(validRows) + ...
                ") does not match stored Nvalid (" + pd(i).Nvalid + "). Do not silently proceed.");
        end

        baRow = baPerDataset(string({baPerDataset.dataset}) == pd(i).dataset);
        if numel(baRow) ~= 1
            error("analyzeTDICoverage:noBlandAltmanMatch", "%s", ...
                pd(i).dataset + ": expected exactly one matching Bland-Altman row, found " + numel(baRow));
        end

        if isempty(diffVec)
            bias = NaN; sdDiff = NaN; cpMDC = NaN; cp2xMDC = NaN; tdiParam = NaN; tdiEmp = NaN;
        else
            bias = mean(diffVec);
            sdDiff = std(diffVec);
            cpMDC = 100 * mean(abs(diffVec) < MDC);
            cp2xMDC = 100 * mean(abs(diffVec) < 2 * MDC);
            tdiParam = zScore * sqrt(bias^2 + sdDiff^2);
            tdiEmp = prctile(abs(diffVec), 100 * opts.Proportion);

            % Cross-check against Finding #126's own numbers (independent
            % reconstruction, not a read of its saved output): bias, sdDiff
            % and CP(MDC) must agree to numerical precision.
            if abs(bias - baRow.bias) > 1e-9 || abs(sdDiff - baRow.sdDiff) > 1e-9 || abs(cpMDC - baRow.pctWithinMDC) > 1e-6
                error("analyzeTDICoverage:blandAltmanMismatch", "%s", ...
                    pd(i).dataset + ": independently reconstructed bias/sdDiff/CP(MDC) do not match " + ...
                    "blandAltman_v001.mat's stored values (Finding #126). Do not silently proceed.");
            end
        end

        perDataset(i).dataset = pd(i).dataset; %#ok<*AGROW>
        perDataset(i).n = numel(diffVec);
        perDataset(i).bias = bias;
        perDataset(i).sdDiff = sdDiff;
        perDataset(i).cpMDC = cpMDC;
        perDataset(i).cp2xMDC = cp2xMDC;
        perDataset(i).tdiParam = tdiParam;
        perDataset(i).tdiEmp = tdiEmp;
    end

    results = struct();
    results.perDataset = perDataset;
    results.config = struct("MetricsMat", opts.MetricsMat, "BlandAltmanMat", opts.BlandAltmanMat, ...
        "NoiseModel", opts.NoiseModel, "MDC", MDC, "DEG_THRESH", DEG_THRESH, "Proportion", opts.Proportion);
    results.runDate = string(datetime("now"));

    outDir = fileparts(mfilename("fullpath"));
    if opts.SaveTxt
        txtPath = fullfile(outDir, "tdiCoverage_v001_summary.txt");
        writeSummary_local(txtPath, perDataset, opts.Proportion);
    end
    if opts.SavePng
        figDir = fullfile(fileparts(outDir), "figures");
        if ~isfolder(figDir)
            mkdir(figDir);
        end
        pngPath = fullfile(figDir, "tdiCoverage_v001.png");
        plotTDICoverage_local(perDataset, MDC, opts.Proportion, pngPath);
    end
    if opts.SaveMat
        matPath = fullfile(outDir, "tdiCoverage_v001.mat");
        save(matPath, "results", "-v7.3");
    end
end

function [diffVec, validRows] = pairedDiffs_local(predMat, obsMat, MDC, DEG_THRESH)
    finiteRows = all(isfinite(predMat), 2) & all(isfinite(obsMat), 2);
    obsAbsMax = max(abs(obsMat), [], 2);
    degenRows = finiteRows & (obsAbsMax < DEG_THRESH);
    scenARows = finiteRows & (obsAbsMax < MDC);
    validRows = find(finiteRows & ~degenRows & ~scenARows);
    predValid = predMat(validRows, :);
    obsValid = obsMat(validRows, :);
    diffVec = predValid(:) - obsValid(:);
end

function writeSummary_local(txtPath, perDataset, prop)
    fid = fopen(txtPath, "w");
    if fid == -1
        error("analyzeTDICoverage:writeFail", "%s", "Could not open " + txtPath + " for write");
    end
    fprintf(fid, "Coverage Probability and Total Deviation Index, per dataset, shaped_xu, pooled (trial x pipeline-pair) points\n");
    fprintf(fid, "CP(MDC) cross-checked exact match to Finding #126 (blandAltman_v001.mat pctWithinMDC)\n");
    fprintf(fid, "TDI at proportion %.2f: param = Lin(2000) normal-theory approx; emp = empirical percentile of |diff|\n\n", prop);
    fprintf(fid, "%-14s %8s %10s %10s %12s %12s %12s %12s\n", ...
        "Dataset", "n", "bias", "sdDiff", "CP(MDC)%", "CP(2xMDC)%", "TDI_param", "TDI_emp");
    for i = 1:numel(perDataset)
        p = perDataset(i);
        fprintf(fid, "%-14s %8d %10.4f %10.4f %12.1f %12.1f %12.4f %12.4f\n", ...
            p.dataset, p.n, p.bias, p.sdDiff, p.cpMDC, p.cp2xMDC, p.tdiParam, p.tdiEmp);
    end
    fclose(fid);
end

function plotTDICoverage_local(perDataset, MDC, prop, pngPath)
    dsLabels = string({perDataset.dataset});
    cpMDC = [perDataset.cpMDC];
    tdiParam = [perDataset.tdiParam];
    tdiEmp = [perDataset.tdiEmp];

    figure("Visible", "off", "Position", [100 100 1300 550]);
    tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "compact");

    nexttile;
    bar(categorical(dsLabels, dsLabels), cpMDC, "FaceColor", [0.3 0.5 0.7]);
    yline(50, "k:");
    ylim([0 100]);
    ylabel("CP(MDC), %");
    title(sprintf("Coverage Probability: %% of pairwise \\Delta\\beta predictions within MDC (%.2f)", MDC));

    nexttile;
    bar(categorical(dsLabels, dsLabels), [tdiParam(:) tdiEmp(:)], "grouped");
    yline(MDC, "k--");
    ylabel(sprintf("TDI_{%.0f}", 100 * prop));
    legend(["parametric (Lin 2000)", "empirical percentile"], "Location", "northoutside");
    title(sprintf("TDI at p=%.2f (dashed line = MDC %.2f, for scale)", prop, MDC));

    sgtitle("Coverage Probability and TDI (shaped\\_xu, pooled trial x pair)");
    exportgraphics(gcf, pngPath, "Resolution", 200);
    close(gcf);
end
