function results = analyzeSEMAgreement_v002(opts)
% analyzeSEMAgreement_v002  Prereg SS7.5 "SEM agreement" metric, per pipeline.
%
% v002 changes from v001 (2026-08-07, Fraser D.S.):
%   1. Pilot's row REMOVED from centroidRows, not fixed. v001's row used
%      "3.91" for Pilot's sigma -- Finding #134 traced this to a wrong
%      pixel/mm conversion baked into that number itself (the true value,
%      recomputed directly from noiseCharacterisation_pilot.mat at the
%      confirmed-correct PX_PER_MM=9.73, is 2.594mm), and it was
%      deliberately left unfixed on Dagmar's explicit instruction ("leave
%      Pilot... moot once the swap is complete"). That swap is now
%      complete -- Pilot dropped outright, 2026-08-07 ("persona non
%      grata... we can iterate and leave it behind"), resolving
%      TODO_FraserIntegration_v001.md Part C4b. Removing the row, rather
%      than fixing its value, is the correct response to that decision:
%      fixing a number that no longer appears anywhere active would be
%      pointless work.
%   2. Fraser's row ADDED: (alpha=4.289, sigma=2.009mm, fs=240Hz), the
%      confirmed centroid from Findings #134/#135 -- not the stale-formula
%      figure, computed correctly from the start.
%   3. Loop-closure file lookup special-cased for Fraser: its all-trials
%      shaped_xu file is loopClosureResults_Fraser_all_shaped_xu_v008.mat
%      (a newer runner version, not a different tier -- same pattern as
%      constellationMetrics_v004.m and runConstellationRSA_HPC_v007.m).
%   4. Default MetricsMat now points at constellationMetrics_v004.mat
%      (Pilot-free by default, Fraser included -- see that file's header).
%
% Everything else -- Part A (grid SEMpredicted vs across-trial SEMobserved,
% pooled and beta-conditioned variants), Part B (grid SEMpredicted vs loop
% closure's own internal forward-simulation SEM), tolerance (0.0108, reused
% from the existing R5/V1 adequacy threshold, disclosed as a reuse not a
% prereg-specified value) -- unchanged from v001. See that file's header
% for the full methodology description; not repeated here.
%
% INPUTS (name-value)
%   SemMat        path to perCoordinateSEM_v2_001.mat
%   MetricsMat    path to constellationMetrics_v004.mat (pipeline label order)
%   SrcDir        directory holding loopClosureResults_<tag>_all_shaped_xu_v007.mat
%                 (Fraser: _v008.mat, see above)
%   Tolerance     default 0.0108
%   BetaBand      [lo hi] for the beta-conditioned variant, default [0.25 0.35]
%   SaveMat/SaveTxt   default true
%
% OUTPUT
%   results.partA(:).{dataset,pipeline,alpha,sigma,fs,semPredPooled,
%                      semPredBetaBand,semObs,absDiffPooled,absDiffBetaBand,
%                      passPooled,passBetaBand}
%   results.partB(:).{dataset,pipeline,semPredPooled,semLoopClosureMed,
%                      semLoopClosureIQR,absDiff}
%
% Fraser, D.S. (2026)  v002
    arguments
        opts.SemMat (1,1) string = fullfile(fileparts(mfilename("fullpath")), "perCoordinateSEM_v2_001.mat")
        opts.MetricsMat (1,1) string = fullfile(fileparts(mfilename("fullpath")), "constellationMetrics_v004.mat")
        opts.SrcDir (1,1) string = fileparts(mfilename("fullpath"))
        opts.Tolerance (1,1) double = 0.0108
        opts.BetaBand (1,2) double = [0.25 0.35]
        opts.SaveMat (1,1) logical = true
        opts.SaveTxt (1,1) logical = true
        opts.SavePng (1,1) logical = true
    end

    if ~isfile(opts.SemMat)
        error("analyzeSEMAgreement:missingInput", "%s", "SemMat not found: " + opts.SemMat);
    end
    if ~isfile(opts.MetricsMat)
        error("analyzeSEMAgreement:missingInput", "%s", "MetricsMat not found: " + opts.MetricsMat);
    end

    Ssem = load(opts.SemMat, "coordTable");
    T = Ssem.coordTable;
    Sc = load(opts.MetricsMat, "results");
    pipelineLabels = Sc.results.pipelineLabels;   % canonical betaObs/betaRecSlice column order
    nPipe = numel(pipelineLabels);

    allAlpha = sort(unique(T.alpha));
    allSigma = sort(unique(T.sigma));
    allFs    = sort(unique(T.fs));

    % v004 centroids (alpha, sigma_mm, fs). Pilot removed (see header).
    % Fraser added: confirmed centroid, Findings #134/#135.
    centroidRows = {
        "Fraser",        4.289, 2.009, 240;
        "Zarandi",       3.184, 4.77, 100;
        "Cook CTRL",     4.771, 8.15, 133;
        "Cook ASD",      5.062, 7.84, 133;
        "Hickman PLAC",  5.337, 7.17, 133;
        "Hickman HALO",  5.424, 7.42, 133;
        "Dhieb",         2.497, 7.50, 100
    };
    tagMap = dictionary(["Fraser","Zarandi","Cook CTRL","Cook ASD","Hickman PLAC","Hickman HALO","Dhieb"], ...
                         ["Fraser","Zarandi","Cook_CTRL","Cook_ASD","Hickman_PLAC","Hickman_HALO","Dhieb"]);

    nDatasets = size(centroidRows, 1);
    partA = struct("dataset", {}, "pipeline", {}, "alphaSnap", {}, "sigmaSnap", {}, "fsSnap", {}, ...
        "semPredPooled", {}, "semPredBetaBand", {}, "semObs", {}, "nObs", {}, ...
        "absDiffPooled", {}, "absDiffBetaBand", {}, "passPooled", {}, "passBetaBand", {});
    partB = struct("dataset", {}, "pipeline", {}, "semPredPooled", {}, ...
        "semLoopClosureMed", {}, "semLoopClosureQ25", {}, "semLoopClosureQ75", {}, ...
        "nTrialsUsed", {}, "absDiff", {});

    for di = 1:nDatasets
        dsLabel = centroidRows{di, 1};
        tA = centroidRows{di, 2}; tS = centroidRows{di, 3}; tF = centroidRows{di, 4};
        [~, ai] = min(abs(allAlpha - tA));
        [~, si] = min(abs(allSigma - tS));
        [~, fi] = min(abs(allFs - tF));
        sA = allAlpha(ai); sS = allSigma(si); sF = allFs(fi);

        tag = tagMap(dsLabel);
        % Fraser: newer runner version (v008), not a different tier -- see header.
        if tag == "Fraser"
            lcPath = fullfile(opts.SrcDir, "loopClosureResults_" + tag + "_all_shaped_xu_v008.mat");
        else
            lcPath = fullfile(opts.SrcDir, "loopClosureResults_" + tag + "_all_shaped_xu_v007.mat");
        end
        if ~isfile(lcPath)
            error("analyzeSEMAgreement:missingLoopClosure", "%s", ...
                dsLabel + ": expected file not found: " + lcPath);
        end
        Sl = load(lcPath, "results");
        lc = Sl.results;
        betaObsAll = vertcat(lc.betaObs);   % [Ntrial x 6], canonical pipeline order
        if size(betaObsAll, 2) ~= nPipe
            error("analyzeSEMAgreement:pipelineCountMismatch", "%s", ...
                dsLabel + ": betaObs has " + size(betaObsAll,2) + " columns, expected " + nPipe);
        end

        baseRows = (T.alpha == sA) & (T.sigma == sS) & (T.fs == sF);
        bandRows = baseRows & (T.betaGen >= opts.BetaBand(1)) & (T.betaGen <= opts.BetaBand(2));

        for p = 1:nPipe
            pipe = pipelineLabels(p);
            pooledRows = baseRows & (T.pipeline == pipe);
            bandRowsP = bandRows & (T.pipeline == pipe);
            if ~any(pooledRows)
                error("analyzeSEMAgreement:noGridRows", "%s", ...
                    dsLabel + " / " + pipe + ": no grid rows at snapped coordinate; cannot proceed silently.");
            end
            semPredPooled = mean(T.sem(pooledRows), "omitnan");
            if any(bandRowsP)
                semPredBand = mean(T.sem(bandRowsP), "omitnan");
            else
                semPredBand = NaN;
            end

            obsCol = betaObsAll(:, p);
            obsCol = obsCol(isfinite(obsCol));
            semObs = std(obsCol, "omitnan");
            nObs = numel(obsCol);

            absDiffPooled = abs(semPredPooled - semObs);
            absDiffBand = abs(semPredBand - semObs);

            row = struct("dataset", dsLabel, "pipeline", pipe, "alphaSnap", sA, "sigmaSnap", sS, "fsSnap", sF, ...
                "semPredPooled", semPredPooled, "semPredBetaBand", semPredBand, "semObs", semObs, "nObs", nObs, ...
                "absDiffPooled", absDiffPooled, "absDiffBetaBand", absDiffBand, ...
                "passPooled", absDiffPooled < opts.Tolerance, "passBetaBand", absDiffBand < opts.Tolerance);
            partA(end+1) = row; %#ok<AGROW>

            % Part B: loop closure's own internal forward-simulation SEM at
            % each trial's own betaGenStarMed anchor (N_REPS=20 replicates),
            % median across trials with a finite anchor and finite slice.
            semTrial = NaN(numel(lc), 1);
            for t = 1:numel(lc)
                slice = lc(t).betaRecSlice(:, p);
                if isfinite(lc(t).betaGenStarMed) && sum(isfinite(slice)) >= 2
                    semTrial(t) = std(slice, "omitnan");
                end
            end
            semTrial = semTrial(isfinite(semTrial));
            rowB = struct("dataset", dsLabel, "pipeline", pipe, "semPredPooled", semPredPooled, ...
                "semLoopClosureMed", median(semTrial, "omitnan"), ...
                "semLoopClosureQ25", prctile(semTrial, 25), "semLoopClosureQ75", prctile(semTrial, 75), ...
                "nTrialsUsed", numel(semTrial), "absDiff", abs(semPredPooled - median(semTrial, "omitnan")));
            partB(end+1) = rowB; %#ok<AGROW>
        end
    end

    results = struct();
    results.partA = partA;
    results.partB = partB;
    results.config = struct("SemMat", opts.SemMat, "MetricsMat", opts.MetricsMat, "SrcDir", opts.SrcDir, ...
        "Tolerance", opts.Tolerance, "BetaBand", opts.BetaBand);
    results.runDate = string(datetime("now"));

    outDir = fileparts(mfilename("fullpath"));
    if opts.SaveTxt
        txtPath = fullfile(outDir, "semAgreement_v002_summary.txt");
        writeSummary_local(txtPath, partA, partB, opts.Tolerance, opts.BetaBand);
    end
    if opts.SavePng
        figDir = fullfile(fileparts(outDir), "figures");
        if ~isfolder(figDir)
            mkdir(figDir);
        end
        pngPath = fullfile(figDir, "semAgreement_v002.png");
        plotSemAgreement_local(partA, partB, opts.Tolerance, pngPath);
    end
    if opts.SaveMat
        matPath = fullfile(outDir, "semAgreement_v002.mat");
        save(matPath, "results", "-v7.3");
    end
end

function plotSemAgreement_local(partA, partB, tol, pngPath)
    dsAll = string({partA.dataset});
    ds = unique(dsAll, "stable");
    nDatasets = numel(ds);
    figure("Visible", "off", "Position", [100 100 1500 950]);
    tiledlayout(3, 3, "TileSpacing", "compact", "Padding", "compact");
    for i = 1:nDatasets
        idxA = dsAll == ds(i);
        rowsA = partA(idxA);
        idxB = string({partB.dataset}) == ds(i);
        rowsB = partB(idxB);
        pipes = string({rowsA.pipeline});
        semPred = [rowsA.semPredPooled];
        semObs = [rowsA.semObs];
        semLC = [rowsB.semLoopClosureMed];

        nexttile;
        y = [semPred; semLC; semObs]';
        bar(categorical(pipes, pipes), y, "grouped");
        set(gca, "YScale", "log");
        yline(tol, "k--");
        title(sprintf("%s", ds(i)));
        ylabel("SEM (log scale)");
        if i == 1
            legend(["grid SEM_{pred}", "loop-closure internal SEM", "across-trial SEM_{obs}"], ...
                "Location", "northoutside", "FontSize", 7);
        end
    end
    sgtitle("SEM agreement: grid-predicted vs loop-closure-internal vs across-trial-observed (shaped\\_xu, log scale, dashed line = tolerance 0.0108)");
    exportgraphics(gcf, pngPath, "Resolution", 200);
    close(gcf);
end

function writeSummary_local(txtPath, partA, partB, tol, betaBand)
    fid = fopen(txtPath, "w");
    if fid == -1
        error("analyzeSEMAgreement:writeFail", "%s", "Could not open " + txtPath + " for write");
    end
    fprintf(fid, "Part A: SEM agreement |SEMpredicted(P) - SEMobserved(P)| < tolerance (%.4f)\n", tol);
    fprintf(fid, "  semPredPooled: grid sem averaged over all betaGen bins at snapped (alpha,sigma,fs)\n");
    fprintf(fid, "  semPredBetaBand: grid sem averaged over betaGen in [%.2f, %.2f] only\n", betaBand(1), betaBand(2));
    fprintf(fid, "  semObs: SD of all-trial betaObs(:,pipeline), NaN-omitted\n\n");
    fprintf(fid, "%-14s %-10s %8s %8s %10s %10s %8s %10s %10s %6s %6s\n", ...
        "Dataset", "Pipeline", "alpha", "sigma", "fs", "semObs", "nObs", "semPred", "semPredBB", "pass", "passBB");
    for i = 1:numel(partA)
        r = partA(i);
        fprintf(fid, "%-14s %-10s %8.3f %8.2f %10d %10.5f %8d %10.5f %10.5f %6d %6d\n", ...
            r.dataset, r.pipeline, r.alphaSnap, r.sigmaSnap, r.fsSnap, r.semObs, r.nObs, ...
            r.semPredPooled, r.semPredBetaBand, r.passPooled, r.passBetaBand);
    end

    fprintf(fid, "\nPart B: grid SEMpredicted(P) vs loop closure's own internal forward-sim SEM\n");
    fprintf(fid, "  semLoopClosureMed: median across trials of SD(betaRecSlice(:,pipeline)) at each trial's own betaGenStarMed anchor, N_REPS=20\n\n");
    fprintf(fid, "%-14s %-10s %10s %14s %20s %10s %8s\n", ...
        "Dataset", "Pipeline", "semPred", "semLC_med", "semLC_[Q25,Q75]", "absDiff", "nTrials");
    for i = 1:numel(partB)
        r = partB(i);
        fprintf(fid, "%-14s %-10s %10.5f %14.5f   [%7.5f,%7.5f] %10.5f %8d\n", ...
            r.dataset, r.pipeline, r.semPredPooled, r.semLoopClosureMed, ...
            r.semLoopClosureQ25, r.semLoopClosureQ75, r.absDiff, r.nTrialsUsed);
    end
    fclose(fid);
end
