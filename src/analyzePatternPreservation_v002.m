function results = analyzePatternPreservation_v002(opts)
% analyzePatternPreservation_v002  Are pipeline rankings preserved?
%
% v002 fixes a real bug found running v001 against constellationMetrics_v004.mat
% (Fraser included, Pilot dropped -- Dagmar, 2026-08-07, "persona non
% grata"): v001's own tagMap dictionary never had a "Fraser" entry, so the
% moment the dataset loop reached Fraser's row, `tagMap(dsLabel)` threw
% "Key not found" -- this was missed earlier the same session when this
% script was checked and judged safe: that check confirmed there was no
% hardcoded NOISE CENTROID (correct, this script never needed one), but
% did not check for this SEPARATE hardcoded dataset-to-tag mapping, which
% v001 also has, for exactly the same structural reason analyzeSEMAgreement
% needed one (both scripts read raw loopClosureResults_<tag>_*.mat files
% directly, unlike Bland-Altman/TDI-Coverage which only ever read
% constellationMetrics's own pre-computed arrays).
%
% Changes from v001:
%   1. "Fraser" added to tagMap (Pilot simply absent now -- no other
%      change needed for Pilot's removal, since this script had no
%      Pilot-specific values beyond the tagMap entry itself).
%   2. Fraser's loop-closure file special-cased to the v008 suffix
%      (loopClosureResults_Fraser_all_shaped_xu_v008.mat -- a newer
%      runner version, not a different tier, same pattern as
%      constellationMetrics_v004.m and analyzeSEMAgreement_v002.m).
%   3. Default MetricsMat now points at constellationMetrics_v004.mat.
%
% Methodology unchanged from v001 -- see that file's header for the full
% description (per-trial Spearman rank correlation of the 6-pipeline
% predicted vs observed rankings; top-1/bottom-1 rank-agreement fractions).
%
% INPUTS (name-value)
%   MetricsMat   path to constellationMetrics_v004.mat (validity mask source)
%   SrcDir       directory holding loopClosureResults_<tag>_all_shaped_xu_v007.mat
%                (Fraser: _v008.mat, see above)
%   SaveMat      default true
%   SaveTxt      default true
%   SavePng      default true
%
% OUTPUT
%   results.perDataset(:).{dataset,n,rhoMed,rhoIQR,pctTop1Match,pctBottom1Match}
%
% Fraser, D.S. (2026)  v002
    arguments
        opts.MetricsMat (1,1) string = fullfile(fileparts(mfilename("fullpath")), "constellationMetrics_v004.mat")
        opts.SrcDir (1,1) string = fileparts(mfilename("fullpath"))
        opts.SaveMat (1,1) logical = true
        opts.SaveTxt (1,1) logical = true
        opts.SavePng (1,1) logical = true
    end

    if ~isfile(opts.MetricsMat)
        error("analyzePatternPreservation:missingInput", "%s", "MetricsMat not found: " + opts.MetricsMat);
    end
    Sc = load(opts.MetricsMat, "results");
    cm = Sc.results;
    sx = cm.byNoiseModel.shaped_xu.perDataset;
    MDC = cm.config.MDC;
    DEG_THRESH = cm.config.DEG_THRESH;

    tagMap = dictionary(["Fraser","Zarandi","Cook CTRL","Cook ASD","Hickman PLAC","Hickman HALO","Dhieb"], ...
                         ["Fraser","Zarandi","Cook_CTRL","Cook_ASD","Hickman_PLAC","Hickman_HALO","Dhieb"]);

    nDatasets = numel(sx);
    perDataset = struct("dataset", {}, "n", {}, "rhoMed", {}, "rhoQ25", {}, "rhoQ75", {}, ...
        "pctTop1Match", {}, "pctBottom1Match", {}, "rhoPerTrial", {});

    for i = 1:nDatasets
        dsLabel = sx(i).dataset;
        tag = tagMap(dsLabel);
        % Fraser: newer runner version (v008), not a different tier -- see header.
        if tag == "Fraser"
            lcPath = fullfile(opts.SrcDir, "loopClosureResults_" + tag + "_all_shaped_xu_v008.mat");
        else
            lcPath = fullfile(opts.SrcDir, "loopClosureResults_" + tag + "_all_shaped_xu_v007.mat");
        end
        if ~isfile(lcPath)
            error("analyzePatternPreservation:missingLoopClosure", "%s", ...
                dsLabel + ": expected file not found: " + lcPath);
        end
        Sl = load(lcPath, "results");
        lc = Sl.results;

        % Reconstruct the exact same validity mask constellationMetrics_v004.m uses,
        % verified against its stored Nvalid for every dataset before this script existed.
        obsFinite = all(isfinite(sx(i).deltaBetaObs), 2);
        predFinite = all(isfinite(sx(i).deltaBetaPred), 2);
        finiteRows = obsFinite & predFinite;
        obsAbsMax = max(abs(sx(i).deltaBetaObs), [], 2);
        degenRows = finiteRows & (obsAbsMax < DEG_THRESH);
        scenARows = finiteRows & (obsAbsMax < MDC);
        validRows = find(finiteRows & ~degenRows & ~scenARows);
        if numel(validRows) ~= sx(i).Nvalid
            error("analyzePatternPreservation:validityMismatch", "%s", ...
                dsLabel + ": reconstructed valid-row count (" + numel(validRows) + ...
                ") does not match stored Nvalid (" + sx(i).Nvalid + "). Do not silently proceed.");
        end
        if numel(lc) ~= size(sx(i).deltaBetaObs, 1)
            error("analyzePatternPreservation:rowCountMismatch", "%s", ...
                dsLabel + ": loopClosureResults trial count (" + numel(lc) + ...
                ") does not match constellationMetrics_v004 row count (" + size(sx(i).deltaBetaObs,1) + "). Row-order correspondence cannot be trusted.");
        end

        nValid = numel(validRows);
        rhoPerTrial = NaN(nValid, 1);
        top1Match = false(nValid, 1);
        bottom1Match = false(nValid, 1);
        for k = 1:nValid
            t = validRows(k);
            obsVec = lc(t).betaObs(:);
            predVec = median(lc(t).betaRecSlice, 1, "omitnan")';
            if numel(obsVec) ~= 6 || numel(predVec) ~= 6 || any(~isfinite(obsVec)) || any(~isfinite(predVec))
                error("analyzePatternPreservation:badTrialData", "%s", ...
                    dsLabel + " trial " + t + ": expected 6 finite pipeline values in both obsVec and predVec.");
            end
            rhoPerTrial(k) = corr(predVec, obsVec, "Type", "Spearman");
            [~, obsTopIdx] = max(obsVec); [~, predTopIdx] = max(predVec);
            [~, obsBotIdx] = min(obsVec); [~, predBotIdx] = min(predVec);
            top1Match(k) = (obsTopIdx == predTopIdx);
            bottom1Match(k) = (obsBotIdx == predBotIdx);
        end

        perDataset(i).dataset = dsLabel; %#ok<*AGROW>
        perDataset(i).n = nValid;
        perDataset(i).rhoMed = median(rhoPerTrial, "omitnan");
        perDataset(i).rhoQ25 = prctile(rhoPerTrial, 25);
        perDataset(i).rhoQ75 = prctile(rhoPerTrial, 75);
        perDataset(i).pctTop1Match = 100 * mean(top1Match);
        perDataset(i).pctBottom1Match = 100 * mean(bottom1Match);
        perDataset(i).rhoPerTrial = rhoPerTrial;
    end

    results = struct();
    results.perDataset = perDataset;
    results.config = struct("MetricsMat", opts.MetricsMat, "SrcDir", opts.SrcDir, "MDC", MDC, "DEG_THRESH", DEG_THRESH);
    results.runDate = string(datetime("now"));

    outDir = fileparts(mfilename("fullpath"));

    if opts.SaveTxt
        txtPath = fullfile(outDir, "patternPreservation_v002_summary.txt");
        writeSummary_local(txtPath, perDataset);
    end

    if opts.SavePng
        figDir = fullfile(fileparts(outDir), "figures");
        if ~isfolder(figDir)
            mkdir(figDir);
        end
        pngPath = fullfile(figDir, "patternPreservation_v002.png");
        plotPatternPreservation_local(perDataset, pngPath);
    end

    if opts.SaveMat
        matPath = fullfile(outDir, "patternPreservation_v002.mat");
        save(matPath, "results", "-v7.3");
    end
end

function writeSummary_local(txtPath, perDataset)
    fid = fopen(txtPath, "w");
    if fid == -1
        error("analyzePatternPreservation:writeFail", "%s", "Could not open " + txtPath + " for write");
    end
    fprintf(fid, "Pattern preservation (Spearman rank correlation of 6-pipeline rankings), per dataset, shaped_xu\n");
    fprintf(fid, "%-14s %8s %10s %18s %14s %14s\n", "Dataset", "n", "rhoMed", "[Q25,Q75]", "%top1match", "%bot1match");
    for i = 1:numel(perDataset)
        p = perDataset(i);
        fprintf(fid, "%-14s %8d %10.4f  [%6.4f,%6.4f] %14.1f %14.1f\n", ...
            p.dataset, p.n, p.rhoMed, p.rhoQ25, p.rhoQ75, p.pctTop1Match, p.pctBottom1Match);
    end
    fclose(fid);
end

function plotPatternPreservation_local(perDataset, pngPath)
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
        histogram(p.rhoPerTrial, "BinLimits", [-1 1], "NumBins", 21, "FaceColor", [0.3 0.3 0.7]);
        hold on;
        xline(p.rhoMed, "r-", "LineWidth", 1.5);
        xline(0, "k:");
        hold off;
        title(sprintf("%s (n=%d)", p.dataset, p.n));
        xlabel("Spearman \\rho (predicted vs observed pipeline ranking)"); ylabel("trial count");
        subtitle(sprintf("median=%.3f, top1match=%.0f%%", p.rhoMed, p.pctTop1Match));
    end
    sgtitle("Pattern preservation: Spearman rank correlation of 6-pipeline rankings (shaped\\_xu)");
    exportgraphics(gcf, pngPath, "Resolution", 200);
    close(gcf);
end
