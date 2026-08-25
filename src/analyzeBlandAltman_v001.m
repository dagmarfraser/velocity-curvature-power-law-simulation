function results = analyzeBlandAltman_v001(opts)
% analyzeBlandAltman_v001  Bland-Altman limits of agreement for the constellation
% predicted-vs-observed delta-beta contrasts.
%
% Prereg §7.5: "Bland-Altman limits of agreement plots will supplement whichever
% primary metric is selected, providing visual assessment of systematic bias and
% precision." Never built as a standalone analysis (docs/TODO_CCC_Toolbox_v001.md,
% docs/PaperDraft_BulletSkeleton_v022.md D6) until now.
%
% Pools all valid (trial x pipeline-pair) points per dataset -- the same
% deltaBetaPred/deltaBetaObs arrays that already feed cccMed/pearsonMed/maeMed
% (docs/CCC_MEASURES_REFERENCE.md) -- into one Bland-Altman analysis per dataset:
%   diff = pred - obs   (same sign convention as MAE)
%   mean = (pred + obs) / 2
%   bias = mean(diff); sdDiff = std(diff); LoA = bias +/- 1.96*sdDiff
%
% Validity mask reproduces constellationMetrics_v002.m's own Nvalid definition
% exactly (verified against the stored Nvalid field for all 7 datasets before
% this script was written): finite in both pred and obs, not degenerate
% (max|obs| < DEG_THRESH), not Scenario A (max|obs| < MDC).
%
% INPUTS (name-value)
%   MetricsMat   path to constellationMetrics_v002.mat (default: src/ colocated)
%   NoiseModel   default "shaped_xu" (canonical per Finding #75)
%   SaveMat      default true
%   SaveTxt      default true
%   SavePng      default true
%
% OUTPUT
%   results.perDataset(:).{dataset,n,bias,sdDiff,loaLo,loaHi,pctWithinMDC}
%
% Fraser, D.S. (2026)  v001
    arguments
        opts.MetricsMat (1,1) string = fullfile(fileparts(mfilename("fullpath")), "constellationMetrics_v002.mat")
        opts.NoiseModel (1,1) string = "shaped_xu"
        opts.SaveMat (1,1) logical = true
        opts.SaveTxt (1,1) logical = true
        opts.SavePng (1,1) logical = true
    end

    if ~isfile(opts.MetricsMat)
        error("analyzeBlandAltman:missingInput", "%s", "MetricsMat not found: " + opts.MetricsMat);
    end
    S = load(opts.MetricsMat, "results");
    src = S.results;

    if ~isfield(src.byNoiseModel, opts.NoiseModel)
        error("analyzeBlandAltman:missingModel", "%s", "NoiseModel '" + opts.NoiseModel + "' not present in " + opts.MetricsMat);
    end
    pd = src.byNoiseModel.(opts.NoiseModel).perDataset;
    MDC = src.config.MDC;
    DEG_THRESH = src.config.DEG_THRESH;

    nDatasets = numel(pd);
    perDataset = struct("dataset", {}, "n", {}, "bias", {}, "sdDiff", {}, ...
        "loaLo", {}, "loaHi", {}, "pctWithinMDC", {}, "meanVec", {}, "diffVec", {});

    for i = 1:nDatasets
        [meanVec, diffVec, validRows] = pairedPoints_local(pd(i).deltaBetaPred, pd(i).deltaBetaObs, MDC, DEG_THRESH);
        if numel(validRows) ~= pd(i).Nvalid
            error("analyzeBlandAltman:validityMismatch", "%s", ...
                pd(i).dataset + ": reconstructed valid-row count (" + numel(validRows) + ...
                ") does not match stored Nvalid (" + pd(i).Nvalid + "). Do not silently proceed.");
        end
        if isempty(diffVec)
            bias = NaN; sdDiff = NaN; loaLo = NaN; loaHi = NaN; pctWithin = NaN;
        else
            bias = mean(diffVec);
            sdDiff = std(diffVec);
            loaLo = bias - 1.96 * sdDiff;
            loaHi = bias + 1.96 * sdDiff;
            pctWithin = 100 * mean(abs(diffVec) < MDC);
        end
        perDataset(i).dataset = pd(i).dataset; %#ok<*AGROW>
        perDataset(i).n = numel(diffVec);
        perDataset(i).bias = bias;
        perDataset(i).sdDiff = sdDiff;
        perDataset(i).loaLo = loaLo;
        perDataset(i).loaHi = loaHi;
        perDataset(i).pctWithinMDC = pctWithin;
        perDataset(i).meanVec = meanVec;
        perDataset(i).diffVec = diffVec;
    end

    results = struct();
    results.perDataset = perDataset;
    results.config = struct("MetricsMat", opts.MetricsMat, "NoiseModel", opts.NoiseModel, ...
        "MDC", MDC, "DEG_THRESH", DEG_THRESH);
    results.runDate = string(datetime("now"));

    outDir = fileparts(mfilename("fullpath"));

    if opts.SaveTxt
        txtPath = fullfile(outDir, "blandAltman_v001_summary.txt");
        writeSummary_local(txtPath, perDataset);
    end

    if opts.SavePng
        figDir = fullfile(fileparts(outDir), "figures");
        if ~isfolder(figDir)
            mkdir(figDir);
        end
        pngPath = fullfile(figDir, "blandAltman_v001.png");
        plotBlandAltman_local(perDataset, pngPath);
    end

    if opts.SaveMat
        matPath = fullfile(outDir, "blandAltman_v001.mat");
        save(matPath, "results", "-v7.3");
    end
end

function [meanVec, diffVec, validRows] = pairedPoints_local(predMat, obsMat, MDC, DEG_THRESH)
    finiteRows = all(isfinite(predMat), 2) & all(isfinite(obsMat), 2);
    obsAbsMax = max(abs(obsMat), [], 2);
    degenRows = finiteRows & (obsAbsMax < DEG_THRESH);
    scenARows = finiteRows & (obsAbsMax < MDC);
    validRows = find(finiteRows & ~degenRows & ~scenARows);
    predValid = predMat(validRows, :);
    obsValid = obsMat(validRows, :);
    diffVec = predValid(:) - obsValid(:);
    meanVec = (predValid(:) + obsValid(:)) / 2;
end

function writeSummary_local(txtPath, perDataset)
    fid = fopen(txtPath, "w");
    if fid == -1
        error("analyzeBlandAltman:writeFail", "%s", "Could not open " + txtPath + " for write");
    end
    fprintf(fid, "Bland-Altman limits of agreement, per-dataset, shaped_xu, pooled (trial x pipeline-pair) points\n");
    fprintf(fid, "%-14s %8s %10s %10s %10s %10s %10s\n", "Dataset", "n", "bias", "sdDiff", "LoA_lo", "LoA_hi", "%|diff|<MDC");
    for i = 1:numel(perDataset)
        p = perDataset(i);
        fprintf(fid, "%-14s %8d %10.4f %10.4f %10.4f %10.4f %10.1f\n", ...
            p.dataset, p.n, p.bias, p.sdDiff, p.loaLo, p.loaHi, p.pctWithinMDC);
    end
    fclose(fid);
end

function plotBlandAltman_local(perDataset, pngPath)
    nDatasets = numel(perDataset);
    figure("Visible", "off", "Position", [100 100 1400 900]);
    tiledlayout(3, 3, "TileSpacing", "compact", "Padding", "compact");
    for i = 1:nDatasets
        p = perDataset(i);
        nexttile;
        if p.n == 0
            title(p.dataset + " (no valid trials)");
            continue
        end
        scatter(p.meanVec, p.diffVec, 4, [0.3 0.3 0.7], "filled", "MarkerFaceAlpha", 0.15);
        hold on;
        yline(0, "k:");
        yline(p.bias, "r-", "LineWidth", 1.2);
        yline(p.loaLo, "r--");
        yline(p.loaHi, "r--");
        hold off;
        title(sprintf("%s (n=%d)", p.dataset, p.n));
        xlabel("mean(pred,obs)"); ylabel("pred - obs");
        subtitle(sprintf("bias=%.3f, LoA=[%.3f, %.3f]", p.bias, p.loaLo, p.loaHi));
    end
    sgtitle("Bland-Altman: predicted vs observed pairwise \\Delta\\beta (shaped\\_xu)");
    exportgraphics(gcf, pngPath, "Resolution", 200);
    close(gcf);
end
