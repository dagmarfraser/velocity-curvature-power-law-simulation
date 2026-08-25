function results = constellationMetrics_v004(opts)
% constellationMetrics_v004  Multi-model constellation metric family.
%
% Only change from v003: Pilot dropped from the default dataset list.
% Decided (Dagmar, 2026-08-07): "persona non grata... we can iterate and
% leave it behind" -- resolves TODO_FraserIntegration_v001.md Part C4b,
% which had been explicitly open all session. Not a supplementary-check
% retention, not a fold-in -- dropped outright. v003 is retained on disk
% as the historical record of the moment Fraser was added while Pilot was
% still present (an unreviewed default choice at the time, not this
% decision); this version is the new canonical default.
%
% Datasets remains a name-value override, so a caller can still request
% Pilot explicitly if a historical comparison is ever genuinely needed --
% but the DEFAULT (what happens when no override is given) no longer
% includes it, which is the actual point: future callers should not have
% to remember to exclude Pilot manually.
%
% Everything else unchanged from v003: Fraser's v008 filename special-case,
% the three metric families (per-trial CCC/Pearson/MAE/Mantel as PRIMARY;
% loopCCC as SECONDARY), the coherence check, output structure.
%
% INPUTS (name-value)
%   Datasets    string array (default: Fraser, Zarandi, Cook CTRL, Cook ASD,
%               Hickman PLAC, Hickman HALO, Dhieb -- no Pilot)
%   NoiseModels string array (default: ["xu","shaped_xu","bootstrap"])
%   NPerm       Mantel permutations (default 999)
%   RngSeed     default 1729
%   SaveMat     default true
%   SaveTxt     default true
%
% OUTPUT
%   results.byNoiseModel.<model>.perDataset  -- per-dataset structs
%   results.comparisonTable                  -- wide table, all models side by side
%
% Reads: loopClosureResults_<tag>_all_<model>_v007.mat for every dataset
% except Fraser, which reads the _v008.mat file.
%
% Fraser, D.S. (2026)  v004

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
            % a newer runner version, not a different tier.
            if tag == "Fraser"
                f = fullfile(srcDir, sprintf( ...
                    "loopClosureResults_%s_all_%s_v008.mat", tag, model));
            else
                f = fullfile(srcDir, sprintf( ...
                    "loopClosureResults_%s_all_%s_v007.mat", tag, model));
            end
            if ~isfile(f)
                warning("constellationMetrics_v004:missingMat", "%s", ...
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
        outMat = fullfile(srcDir,"constellationMetrics_v004.mat");
        save(outMat,"results","-v7.3");
        fprintf("Saved: %s\n", outMat);
    end
    if opts.SaveTxt
        outTxt = fullfile(srcDir,"constellationMetrics_v004_summary.txt");
        writeSummaryTxt_local(outTxt, byNoiseModel, opts, comparisonTable);
        fprintf("Saved: %s\n", outTxt);
    end
end

%% ============================ LOCAL HELPERS ================================

function [pairP, pairQ] = pairIndices_local(n)
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
        error("constellationMetrics_v004:txtOpen", "%s", "Cannot open " + outTxt);
    end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, "constellationMetrics_v004 -- per-trial delta-beta, multi-model\n");
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
