function results = constellationMetrics_v002(opts)
% constellationMetrics_v002  Multi-model constellation metric family.
%
% Extends v001 to run all three noise models (xu, shaped_xu, bootstrap)
% simultaneously, producing a side-by-side comparison table. This supports
% two distinct uses of the constellation CCC, documented in detail in
% docs/CCC_MEASURES_REFERENCE.md:
%
%   NOISE MODEL SELECTION (loopCCC, stored in each mat, RELATIVE comparison):
%     Pool (betaObs_pp, median betaRecSlice_pp) across all trials and all six
%     pipelines; compute Lin's CCC. The between-pipeline main effect is constant
%     across noise models (it comes from the observed side), so it cancels in the
%     ranking shaped_xu > xu > bootstrap. This was used to SELECT shaped_xu as
%     the production surrogate. Because it is a RELATIVE comparison the pooling
%     inflation does not invalidate the selection. loopCCC should never be used
%     as an absolute validation claim.
%
%   VALIDATION (per-trial contrast CCC, PRIMARY here, ABSOLUTE claim):
%     Per trial, form the 15 inter-pipeline contrasts delta-beta; compare
%     predicted vs observed. Removes the per-trial common mode. Valid as an
%     absolute number. See docs/CCC_MEASURES_REFERENCE.md Section 4.
%
% PRIMARY metrics (per valid trial), predicted vs observed 15-contrast vectors:
%   CCC      Lin's concordance on signed delta-beta   (shape + bias penalty)
%   Pearson  r from the same call                     (shape, no bias penalty)
%   MAE      mean |delta-beta_pred - delta-beta_obs|  (magnitude error)
%   Mantel   label-permutation correlation of 6x6 |delta-beta| distance mats
%            (topological fidelity; prereg section 7.5/7.7)
%
% COHERENCE CHECK: loopCCC is independently recomputed from raw slices using
% the runner's own formula, so any plumbing error announces itself.
%
% INPUTS (name-value)
%   Datasets    string array (default: all in-scope ellipse datasets)
%   NoiseModels string array (default: ["xu","shaped_xu","bootstrap"])
%   NPerm       Mantel permutations (default 999)
%   RngSeed     default 1729
%   SaveMat     default true
%   SaveTxt     default true
%
% OUTPUT
%   results.byNoiseModel.<model>.perDataset  -- per-dataset structs (as v001)
%   results.comparisonTable                  -- wide table, all models side by side
%
% Reads (no re-simulation): loopClosureResults_<tag>_all_<model>_v007.mat,
% EXCEPT Fraser, which reads the _v008.mat runner suffix (Finding #137 --
% a newer runner version, not a different tier; matches constellationMetrics_v004.m's
% own special-case).
% Dhieb is auto-included once its full shaped_xu mat exists (skipped if absent).
%
% CORRECTION, 2026-08-26: default Datasets swapped Pilot -> Fraser. This
% script's own checked-in output (constellationMetrics_v002.mat) had drifted
% out of sync with its own already-published downstream dependents
% (blandAltman_v001.mat 07-Aug-2026 and the SEM-agreement/pattern-preservation/
% TDI-CP outputs built from it) -- those were generated from a Fraser-inclusive
% run of this file that was later overwritten by an archival Pilot-only
% snapshot (29-Jun-2026), reportedly during a decision to preserve v002 as
% "historical three-model rows" for a DIFFERENT citation (the xu/shaped_xu/
% bootstrap noise-model-selection narrative, which genuinely is Pilot-era and
% is not affected by this fix). Pilot itself is fully retired ("persona non
% grata", see claude.md) -- this correction resolves the drift by making
% Fraser the live default, matching v004. Pre-correction Pilot-based output
% preserved verbatim at constellationMetrics_v002_PILOT_ARCHIVE_20260629.mat.
%
% Fraser, D.S. (2026)  v002 (corrected 2026-08-26)

    arguments
        opts.Datasets (1,:) string = ["Fraser","Zarandi","Cook CTRL","Cook ASD", ...
                                      "Hickman PLAC","Hickman HALO","Dhieb"]
        opts.NoiseModels (1,:) string = ["xu","shaped_xu","bootstrap"]
        opts.NPerm (1,1) double {mustBeInteger, mustBePositive} = 999
        opts.RngSeed (1,1) double = 1729
        opts.SaveMat (1,1) logical = true
        opts.SaveTxt (1,1) logical = true
    end

    srcDir = fileparts(mfilename("fullpath"));
    addpath(genpath(fullfile(srcDir, "functions")));
    addpath(genpath(fullfile(srcDir, "req")));

    MDC        = 0.03;
    DEG_THRESH = 0.005;
    pipelineLabels = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];
    nPipe = numel(pipelineLabels);
    [pairP, pairQ] = pairIndices_local(nPipe);
    nPairs = numel(pairP);

    rng(opts.RngSeed, "twister");

    nModels = numel(opts.NoiseModels);
    byNoiseModel = struct();

    for mi = 1:nModels
        model  = opts.NoiseModels(mi);
        mField = matlab.lang.makeValidName(model);
        fprintf("\n=== Model: %s ===\n", model);
        fprintf("PRIMARY=per-trial contrast; loopCCC=SECONDARY(noise-model selection only).\n");

        perDataset = struct("dataset",{}, "noiseModel",{}, ...
            "Ntot",{}, "Nfinite",{}, "Nvalid",{}, "Ndegen",{}, "NscenA",{}, ...
            "cccMed",{}, "cccPool",{}, "cccPoolLo",{}, "cccPoolHi",{}, ...
            "pearsonMed",{}, "maeMed",{}, "maePool",{}, ...
            "mantelMed",{}, "mantelSigFrac",{}, "mantelPerTrial",{}, ...
            "loopCCCstored",{}, "loopCCCrecomp",{}, "loopCCCrecompLin",{}, ...
            "loopCCCmaxDelta",{}, "cccPerTrial",{}, ...
            "deltaBetaPred",{}, "deltaBetaObs",{});

        for d = 1:numel(opts.Datasets)
            ds  = opts.Datasets(d);
            tag = replace(ds, " ", "_");
            % Fraser uses the v008 runner suffix, not v007 (Finding #137) --
            % a newer runner version, not a different tier. Matches
            % constellationMetrics_v004.m's own special-case.
            if tag == "Fraser"
                f = fullfile(srcDir, sprintf( ...
                    "loopClosureResults_%s_all_%s_v008.mat", tag, model));
            else
                f = fullfile(srcDir, sprintf( ...
                    "loopClosureResults_%s_all_%s_v007.mat", tag, model));
            end
            if ~isfile(f)
                warning("constellationMetrics_v002:missingMat", "%s", ...
                    "Skipping " + ds + " [" + model + "] -- mat absent.");
                continue
            end

            S = load(f, "results", "loopCCC");
            r = S.results;  N = numel(r);
            loopCCCstored = NaN;
            if isfield(S,"loopCCC"); loopCCCstored = S.loopCCC; end

            [obsPool, predPool] = poolRawConstellation_local(r, nPipe);
            loopCCCrecomp    = linCCCshortcut_local(obsPool, predPool);
            ccLin            = linCCC_v001(obsPool, predPool);
            loopCCCrecompLin = ccLin.ccc;
            loopCCCmaxDelta  = max(abs([ ...
                loopCCCrecomp    - loopCCCstored, ...
                loopCCCrecompLin - loopCCCstored, ...
                loopCCCrecompLin - loopCCCrecomp]));

            dbObs  = NaN(N,nPairs); dbPred = NaN(N,nPairs);
            bObsM  = NaN(N,nPipe);  bPredM = NaN(N,nPipe);
            degen  = false(N,1); scenA = false(N,1); finite = false(N,1);

            for t = 1:N
                bo = r(t).betaObs(:)';
                if numel(bo)~=nPipe || any(~isfinite(bo)); continue; end
                slice = r(t).betaRecSlice;
                if isempty(slice)||size(slice,2)~=nPipe; continue; end
                bp = median(slice,1,"omitnan");
                if any(~isfinite(bp)); continue; end
                if ~isfinite(r(t).betaGenStarMed); continue; end
                finite(t)   = true;
                bObsM(t,:)  = bo;  bPredM(t,:) = bp;
                dbObs(t,:)  = bo(pairP) - bo(pairQ);
                dbPred(t,:) = bp(pairP) - bp(pairQ);
                degen(t)    = max(abs(dbObs(t,:))) < DEG_THRESH;
                scenA(t)    = all(abs(dbObs(t,:)) < MDC);
            end

            valid = finite & ~degen & ~scenA;
            vi = find(valid); nV = numel(vi);

            cccT=NaN(nV,1); rT=NaN(nV,1); maeT=NaN(nV,1);
            manT=NaN(nV,1); manP_=NaN(nV,1);
            for ii = 1:nV
                t  = vi(ii);
                cc = linCCC_v001(dbPred(t,:)', dbObs(t,:)');
                cccT(ii) = cc.ccc;  rT(ii) = cc.r;
                maeT(ii) = mean(abs(dbPred(t,:) - dbObs(t,:)));
                [manT(ii), manP_(ii)] = mantel_local( ...
                    bPredM(t,:), bObsM(t,:), pairP, pairQ, opts.NPerm);
            end

            cccPool=NaN; cccLo=NaN; cccHi=NaN; maePool=NaN;
            if nV >= 2
                allP = dbPred(vi,:); allO = dbObs(vi,:);
                cc   = linCCC_v001(allP(:), allO(:));
                cccPool = cc.ccc;
                if isfield(cc,"LowerCI"); cccLo=cc.LowerCI; cccHi=cc.UpperCI; end
                maePool = mean(abs(allP(:)-allO(:)));
            end

            rec = struct();
            rec.dataset          = ds;
            rec.noiseModel       = model;
            rec.Ntot             = N;
            rec.Nfinite          = sum(finite);
            rec.Nvalid           = nV;
            rec.Ndegen           = sum(degen & finite);
            rec.NscenA           = sum(scenA & finite);
            rec.cccMed           = median(cccT,"omitnan");
            rec.cccPool          = cccPool;
            rec.cccPoolLo        = cccLo;
            rec.cccPoolHi        = cccHi;
            rec.pearsonMed       = median(rT,"omitnan");
            rec.maeMed           = median(maeT,"omitnan");
            rec.maePool          = maePool;
            rec.mantelMed        = median(manT,"omitnan");
            rec.mantelSigFrac    = mean(manP_ < 0.05,"omitnan");
            rec.mantelPerTrial   = manT;
            rec.loopCCCstored    = loopCCCstored;
            rec.loopCCCrecomp    = loopCCCrecomp;
            rec.loopCCCrecompLin = loopCCCrecompLin;
            rec.loopCCCmaxDelta  = loopCCCmaxDelta;
            rec.cccPerTrial      = cccT;
            rec.deltaBetaPred    = dbPred;
            rec.deltaBetaObs     = dbObs;
            perDataset(end+1)    = rec; %#ok<AGROW>

            cohTag = "OK";
            if ~(loopCCCmaxDelta < 1e-3); cohTag = "FORMULA_DIFF"; end
            if nV >= 2
                fprintf("  %-13s Nv=%4d cccMed=%+.3f cccPool=%+.3f " + ...
                    "pear=%+.3f MAE=%.4f mant=%+.3f(%.0f%%) " + ...
                    "loopCCC=%.3f [coh:%s]\n", ...
                    ds, nV, rec.cccMed, rec.cccPool, rec.pearsonMed, ...
                    rec.maeMed, rec.mantelMed, 100*rec.mantelSigFrac, ...
                    loopCCCstored, cohTag);
            else
                fprintf("  %-13s Nv=%4d DEGENERATE (scenA=%d)  loopCCC=%.3f [coh:%s]\n", ...
                    ds, nV, rec.NscenA, loopCCCstored, cohTag);
            end
        end
        byNoiseModel.(mField).perDataset = perDataset;
    end

    comparisonTable = buildComparisonTable_local(byNoiseModel, opts.NoiseModels);

    fprintf("\n=== COMPARISON (loopCCC | cccMed | MAE | Mantel) ===\n");
    fprintf("%-13s %-10s %8s %8s %8s %8s\n", ...
        "dataset","model","loopCCC","cccMed","MAE","mantel");
    fprintf("%s\n", repmat('-',1,66));
    for mi = 1:nModels
        mf = matlab.lang.makeValidName(opts.NoiseModels(mi));
        if ~isfield(byNoiseModel,mf); continue; end
        pd = byNoiseModel.(mf).perDataset;
        for d = 1:numel(pd)
            if pd(d).Nvalid < 2
                fprintf("%-13s %-10s %8.3f %8s %8s %8s  [DEGENERATE Nvalid=%d]\n", ...
                    pd(d).dataset, pd(d).noiseModel, pd(d).loopCCCstored, ...
                    "---","---","---", pd(d).Nvalid);
            else
                fprintf("%-13s %-10s %8.3f %8.3f %8.4f %8.3f\n", ...
                    pd(d).dataset, pd(d).noiseModel, pd(d).loopCCCstored, ...
                    pd(d).cccMed, pd(d).maeMed, pd(d).mantelMed);
            end
        end
        fprintf("\n");
    end

    results = struct();
    results.byNoiseModel    = byNoiseModel;
    results.comparisonTable = comparisonTable;
    results.pipelineLabels  = pipelineLabels;
    results.pairP = pairP;  results.pairQ = pairQ;
    results.config = struct("NPerm",opts.NPerm,"RngSeed",opts.RngSeed, ...
        "MDC",MDC,"DEG_THRESH",DEG_THRESH,"NoiseModels",opts.NoiseModels);
    results.runDate = string(datetime("now"));

    if opts.SaveMat
        outMat = fullfile(srcDir,"constellationMetrics_v002.mat");
        save(outMat,"results","-v7.3");
        fprintf("Saved: %s\n", outMat);
    end
    if opts.SaveTxt
        outTxt = fullfile(srcDir,"constellationMetrics_v002_summary.txt");
        writeSummaryTxt_local(outTxt, byNoiseModel, opts, comparisonTable);
        fprintf("Saved: %s\n", outTxt);
    end
end

%% ============================ LOCAL HELPERS ================================

function [pairP, pairQ] = pairIndices_local(n)
% Lower-triangle ordered pairs P>Q (matches betaObs(:) column-major convention).
    pairP = zeros(1, n*(n-1)/2);
    pairQ = zeros(1, n*(n-1)/2);
    k = 0;
    for pp = 2:n
        for qq = 1:pp-1
            k = k + 1;
            pairP(k) = pp;
            pairQ(k) = qq;
        end
    end
end

function [mR, mP] = mantel_local(bPred, bObs, pairP, pairQ, nPerm)
% Mantel test between the two abs(delta-beta) distance structures of the six
% pipelines. Two-sided label-permutation p. Returns NaN when either distance
% vector is constant (degenerate constellation -- same failure mode as CCC).
    dPred = abs(bPred(pairP) - bPred(pairQ));
    dObs  = abs(bObs(pairP)  - bObs(pairQ));
    if std(dPred) == 0 || std(dObs) == 0
        mR = NaN; mP = NaN; return
    end
    oC  = dObs - mean(dObs);
    oN  = sqrt(sum(oC.^2));
    pC  = dPred - mean(dPred);
    mR  = sum(pC .* oC) / (sqrt(sum(pC.^2)) * oN);
    nGE = 0;
    nP  = numel(bPred);
    for p = 1:nPerm
        bp  = bPred(randperm(nP));
        dP  = abs(bp(pairP) - bp(pairQ));
        sP  = std(dP);
        if sP == 0; continue; end
        pCp = dP - mean(dP);
        rp  = sum(pCp .* oC) / (sqrt(sum(pCp.^2)) * oN);
        if abs(rp) >= abs(mR); nGE = nGE + 1; end
    end
    mP = (nGE + 1) / (nPerm + 1);
end

function [obsPool, predPool] = poolRawConstellation_local(r, nPipe)
% Pool (betaObs_pp, median betaRecSlice_pp) across all trials and all pipelines,
% replicating runLoopClosureFftnoise_v007/computeLoopCCC_local inclusion exactly
% (slice not all-NaN; both finite). This is the RAW pipeline-level pooling that
% underlies loopCCC -- used only for the coherence check, never as a primary.
    N = numel(r);
    obsPool  = NaN(N*nPipe, 1);
    predPool = NaN(N*nPipe, 1);
    c = 0;
    for t = 1:N
        slice = r(t).betaRecSlice;
        if isempty(slice) || all(isnan(slice(:))); continue; end
        bo = r(t).betaObs;
        for pp = 1:nPipe
            bObs  = bo(pp);
            bPred = median(slice(:, pp), "omitnan");
            if isfinite(bObs) && isfinite(bPred)
                c = c + 1;
                obsPool(c)  = bObs;
                predPool(c) = bPred;
            end
        end
    end
    obsPool  = obsPool(1:c);
    predPool = predPool(1:c);
end

function ccc = linCCCshortcut_local(x, y)
% Verbatim copy of runLoopClosureFftnoise_v007/linCCC_local so the recomputed
% loopCCC matches the stored value to numerical precision. Do not "improve" the
% formula here: its job is to reproduce the runner exactly.
    valid = isfinite(x) & isfinite(y);
    x = x(valid); y = y(valid);
    if numel(x) < 3; ccc = NaN; return; end
    mx  = mean(x);  my  = mean(y);
    sx2 = var(x, 0); sy2 = var(y, 0);
    sxy = mean((x - mx) .* (y - my));
    denom = sx2 + sy2 + (mx - my)^2;
    if denom < eps; ccc = NaN; return; end
    ccc = 2*sxy / denom;
end

function T = buildComparisonTable_local(byNoiseModel, models)
% Wide table: one row per (dataset x model).
    rows = struct("dataset",{}, "noiseModel",{}, "loopCCC",{}, ...
        "cccMed",{}, "cccPool",{}, "pearsonMed",{}, "maeMed",{}, ...
        "mantelMed",{}, "mantelSigFrac",{}, "Nvalid",{});
    for mi = 1:numel(models)
        mf = matlab.lang.makeValidName(models(mi));
        if ~isfield(byNoiseModel, mf); continue; end
        pd = byNoiseModel.(mf).perDataset;
        for d = 1:numel(pd)
            r = struct("dataset", pd(d).dataset, "noiseModel", pd(d).noiseModel, ...
                "loopCCC", pd(d).loopCCCstored, "cccMed", pd(d).cccMed, ...
                "cccPool", pd(d).cccPool, "pearsonMed", pd(d).pearsonMed, ...
                "maeMed", pd(d).maeMed, "mantelMed", pd(d).mantelMed, ...
                "mantelSigFrac", pd(d).mantelSigFrac, "Nvalid", pd(d).Nvalid);
            rows(end+1) = r; %#ok<AGROW>
        end
    end
    if isempty(rows); T = table(); return; end
    T = struct2table(rows);
end

function writeSummaryTxt_local(outTxt, byNoiseModel, opts, compTab)
    fid = fopen(outTxt, "w");
    if fid < 0
        error("constellationMetrics_v002:txtOpen", "%s", "Cannot open " + outTxt);
    end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, "constellationMetrics_v002 -- per-trial delta-beta, multi-model\n");
    fprintf(fid, "Generated: %s | NPerm=%d | Seed=%d\n", ...
        string(datetime("now")), opts.NPerm, opts.RngSeed);
    fprintf(fid, "Models: %s\n\n", strjoin(opts.NoiseModels, ", "));
    fprintf(fid, "PRIMARY = per-trial contrast CCC/Pearson/MAE/Mantel on delta-beta.\n");
    fprintf(fid, "loopCCC = pooled-across-pipelines SECONDARY (noise model selection,\n");
    fprintf(fid, "not absolute validation). cohDelta = max|recomp-stored| (~0 = OK).\n");
    fprintf(fid, "See docs/CCC_MEASURES_REFERENCE.md for full definitions.\n\n");
    fprintf(fid, "%-13s %-10s %6s %7s %8s %9s %8s %8s %8s | %8s %9s\n", ...
        "dataset","model","Nvald","cccMed","cccPool","pearsMed","maeMed", ...
        "mantMed","mantSig%","loopCCC","cohDelta");
    fprintf(fid, "%s\n", repmat('-', 1, 117));
    models = opts.NoiseModels;
    for mi = 1:numel(models)
        mf = matlab.lang.makeValidName(models(mi));
        if ~isfield(byNoiseModel, mf); continue; end
        pd = byNoiseModel.(mf).perDataset;
        for i = 1:numel(pd)
            if pd(i).Nvalid < 2
                fprintf(fid, "%-13s %-10s %6d %7s %8s %9s %8s %8s %8s | %8.3f %9s  [DEGENERATE]\n", ...
                    pd(i).dataset, pd(i).noiseModel, pd(i).Nvalid, ...
                    "---","---","---","---","---","---", ...
                    pd(i).loopCCCstored, "---");
            else
                fprintf(fid, "%-13s %-10s %6d %+7.3f %+8.3f %+9.3f %8.4f %+8.3f %8.0f | %8.3f %9.1e\n", ...
                    pd(i).dataset, pd(i).noiseModel, pd(i).Nvalid, ...
                    pd(i).cccMed, pd(i).cccPool, pd(i).pearsonMed, ...
                    pd(i).maeMed, pd(i).mantelMed, 100*pd(i).mantelSigFrac, ...
                    pd(i).loopCCCstored, pd(i).loopCCCmaxDelta);
            end
        end
        fprintf(fid, "\n");
    end
    if ~isempty(compTab) && istable(compTab)
        fprintf(fid, "\n--- Comparison table (struct2table) ---\n");
        fprintf(fid, "%s\n", formattable_local(compTab));
    end
end

function s = formattable_local(T)
    s = evalc("disp(T)");
end

