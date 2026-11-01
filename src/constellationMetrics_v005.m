function [T, PT, results] = constellationMetrics_v005(opts)
%CONSTELLATIONMETRICS_V005  Constellation metric family on the v015 basis, each with a trial-agnostic NULL (replaces v004 / the v007-v008 main table; Z12, Addenda 52-54).
%   Per dataset (and pooled over in-domain / in-domain-without-Fraser / all seven):
%     PRIMARY (prereg S7.5 Step 3) per-trial Lin's CCC of the 15 predicted vs observed inter-pipeline contrasts: cccMed (+ bootstrap CI), cccPool, pearsonMed, maeMed, Mantel (RSA, exploratory)
%     SUPPLEMENTARY (prereg additional metrics / D6): pattern preservation (per-trial Spearman rho of the 6-pipeline beta rankings, top-1 / bottom-1 match, as analyzePatternPreservation_v002);
%       pooled (trial x pair) Bland-Altman bias and limits, Coverage Probability CP(MDC) and CP(2xMDC), TDI90 parametric and empirical (as analyzeTDICoverage_v001, same diff = pred - obs);
%       secondary loopCCC (stored and recomputed).
%     NULL (post hoc, not registered; Addendum 54): prediction replaced by the leave-one-out dataset-mean OBSERVED vector (contrasts for CCC / MAE / CP / TDI / bias; the 6 observed betas for rho and
%       top-1 / bottom-1). It carries no trial-specific information; metric - null shows what a stable pipeline signature alone achieves. No null is computed for Mantel.
%   Validity mask, extraction (median of betaRecSlice) and fail-loud checks as checkNoiseModelCCC_v001 / constellationMetrics_v004: N_REPS 200, N_BETA 25, seed 1729, identical trial IDs across models;
%   shaped_xu = loopClosureResults_<ds>_all_shaped_xu_v015.mat; xu / bootstrap = ..._v017.mat (D-H1). Dhieb and Zarandi are outside the validity domain (Findings #223, #242); Zarandi has 1 valid trial.
%   Writes ../results/constellationMetrics_v005.mat (results, T, PT) and logs/. Run from the folder holding the loopClosureResults files, or pass Dir.
%   USAGE  [T, PT] = constellationMetrics_v005;   constellationMetrics_v005(Dir="<folder>", NoiseModels=["shaped_xu" "xu" "bootstrap"])
%   Fraser, D.S. (2026)  v005
arguments
    opts.Dir (1,1) string = ""
    opts.Datasets (1,:) string = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Dhieb" "Zarandi"]
    opts.NoiseModels (1,:) string {mustBeMember(opts.NoiseModels, ["shaped_xu" "xu" "bootstrap"])} = "shaped_xu"
    opts.NPerm (1,1) double {mustBeInteger, mustBePositive} = 999
    opts.NBoot (1,1) double {mustBeInteger, mustBePositive} = 1000
    opts.OutTag (1,1) string = ""
end
srcDir = fileparts(mfilename('fullpath'));  addpath(genpath(fullfile(srcDir, 'functions')));
dataDir = opts.Dir;  if dataDir == "", dataDir = string(srcDir); end
resDir = fullfile(srcDir, '..', 'results');  logDir = fullfile(srcDir, 'logs');
for d = {resDir, logDir}, if ~isfolder(d{1}), mkdir(d{1}); end, end   % cell, not [char char]: the old form iterated over single characters and created stray one-letter folders
lg = fopen(fullfile(logDir, "constellationMetrics_v005_" + string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) + ".log"), 'w');  cleanup = onCleanup(@() fclose(lg)); %#ok<NASGU>
say = @(varargin) fprintf_both_local(lg, varargin{:});
MDC = 0.03;  DEG = 0.005;  nPipe = 6;  [pP, pQ] = pairs_local(nPipe);  z90 = norminv(0.95);  rng(1729, 'twister');
inDomain = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO"];
verOf = containers.Map(["shaped_xu" "xu" "bootstrap"], ["v015" "v017" "v017"]);
ids = containers.Map();  rows = {};  trials = struct('ds', {}, 'model', {}, 'cc', {}, 'nM', {}, 'rho', {}, 'nRho', {});
say("constellationMetrics_v005 | %s | data %s | NPerm %d\n", string(datetime('now')), dataDir, opts.NPerm);
for model = opts.NoiseModels
    for ds = opts.Datasets
        f = fullfile(dataDir, "loopClosureResults_" + ds + "_all_" + model + "_" + verOf(char(model)) + ".mat");
        if ~isfile(f), error('cm5:missing', 'No file: %s', f); end
        S0 = load(f, 'results', 'config', 'loopCCC');  r = S0.results;  N = numel(r);  c = S0.config;
        if c.N_REPS ~= 200 || c.N_BETA ~= 25 || c.RNG_SEED ~= 1729, error('cm5:config', '%s: N_REPS/N_BETA/seed differ from 200/25/1729', f); end
        tid = string({r.trialID})';
        if ~isKey(ids, char(ds)), ids(char(ds)) = tid; elseif ~isequal(ids(char(ds)), tid), error('cm5:trials', '%s: trial IDs differ from the first model file', f); end
        bO = NaN(N, nPipe);  bP = bO;  dO = NaN(N, numel(pP));  dP = dO;  ok = false(N, 1);  gStar = NaN(N, 1);
        for t = 1:N
            bo = r(t).betaObs(:)';  sl = r(t).betaRecSlice;  gStar(t) = r(t).betaGenStarMed;
            if numel(bo) ~= nPipe || any(~isfinite(bo)) || isempty(sl) || size(sl, 2) ~= nPipe || ~isfinite(gStar(t)), continue; end
            bp = median(sl, 1, 'omitnan');  if any(~isfinite(bp)), continue; end
            bO(t, :) = bo;  bP(t, :) = bp;  dO(t, :) = bo(pP) - bo(pQ);  dP(t, :) = bp(pP) - bp(pQ);  ok(t) = true;
        end
        degen = ok & max(abs(dO), [], 2) < DEG;  scenA = ok & all(abs(dO) < MDC, 2);
        vi = find(ok & ~degen & ~scenA);  n = numel(vi);  O = dO(vi, :);  Pm = dP(vi, :);  BO = bO(vi, :);  BP = bP(vi, :);
        cc = NaN(n, 1);  rr = cc;  mae = cc;  nM = cc;  nMae = cc;  rho = cc;  nRho = cc;  t1 = false(n, 1);  b1 = t1;  nt1 = t1;  nb1 = t1;  man = cc;  manP = cc;
        sumO = sum(O, 1);  sumB = sum(BO, 1);  Pn = NaN(n, size(O, 2));
        for i = 1:n
            q = linCCC_v001(Pm(i, :)', O(i, :)');  cc(i) = q.ccc;  rr(i) = q.r;  mae(i) = mean(abs(Pm(i, :) - O(i, :)));
            rho(i) = corr(BP(i, :)', BO(i, :)', 'Type', 'Spearman');  [~, a] = max(BO(i, :));  [~, b] = max(BP(i, :));  t1(i) = a == b;
            [~, a2] = min(BO(i, :));  [~, b2] = min(BP(i, :));  b1(i) = a2 == b2;
            [man(i), manP(i)] = mantel_local(BP(i, :), BO(i, :), pP, pQ, opts.NPerm);
            if n >= 2
                Pn(i, :) = (sumO - O(i, :)) / (n - 1);  bn = (sumB - BO(i, :)) / (n - 1);
                nM(i) = linCCC_v001(Pn(i, :)', O(i, :)').ccc;  nMae(i) = mean(abs(Pn(i, :) - O(i, :)));
                nRho(i) = corr(bn', BO(i, :)', 'Type', 'Spearman');  [~, nb] = max(bn);  nt1(i) = nb == a;  [~, nb2] = min(bn);  nb1(i) = nb2 == a2;
            end
        end
        lo = NaN;  hi = NaN;
        if nnz(isfinite(cc)) >= 10, bs = bootstrp(opts.NBoot, @(x) median(x, 'omitnan'), cc(isfinite(cc)));  lo = prctile(bs, 2.5);  hi = prctile(bs, 97.5); end
        cp = NaN;  cl = NaN;  ch = NaN;  sr = NaN;
        if n >= 2, q = linCCC_v001(Pm(:), O(:));  cp = q.ccc;  sr = q.r;  if isfield(q, 'LowerCI'), cl = q.LowerCI; ch = q.UpperCI; end, end
        dfM = Pm(:) - O(:);  dfN = Pn(:) - O(:);   % pooled differences, pred - obs (as Bland-Altman / TDI / CP scripts)
        stat = @(v) deal(mean(v), std(v), 100 * mean(abs(v) < MDC), 100 * mean(abs(v) < 2 * MDC), z90 * sqrt(mean(v)^2 + std(v)^2), prctile(abs(v), 90));
        [bi, sdd, cpm, cp2, tdp, tde] = stat(dfM);  [nbi, ~, ncp, ncp2, ~, ntde] = stat(dfN);
        fb = mean(cc > nM, 'omitnan');  nT1 = 100 * mean(nt1);  nB1 = 100 * mean(nb1);
        if n < 2, [nbi, ncp, ncp2, ntde, fb, nT1, nB1] = deal(NaN);  end   % no leave-one-out null from a single trial: NaN, not 0
        lcc = pooledObs_local(r, nPipe);
        rows{end+1, 1} = table(ds, model, N, nnz(ok), nnz(degen), nnz(scenA), n, median(cc, 'omitnan'), lo, hi, cp, cl, ch, sr, median(rr, 'omitnan'), median(mae, 'omitnan'), ... %#ok<AGROW>
            median(nM, 'omitnan'), median(cc - nM, 'omitnan'), fb, median(nMae, 'omitnan'), ...
            median(rho, 'omitnan'), prctile(rho, 25), prctile(rho, 75), 100 * mean(t1), 100 * mean(b1), median(nRho, 'omitnan'), nT1, nB1, ...
            median(man, 'omitnan'), mean(manP < 0.05, 'omitnan'), numel(dfM), bi, sdd, bi - 1.96 * sdd, bi + 1.96 * sdd, cpm, cp2, tdp, tde, nbi, ncp, ncp2, ntde, ...
            S0.loopCCC, lcc, median(gStar(vi), 'omitnan'), ismember(ds, inDomain), ...
            'VariableNames', ["dataset" "noiseModel" "Ntot" "Nfinite" "Ndegen" "NscenA" "Nvalid" "cccMed" "cccMedLo" "cccMedHi" "cccPool" "cccPoolLo" "cccPoolHi" "stackedR" "pearsonMed" "maeMed" ...
            "nullCccMed" "gainCccMed" "fracBeatNull" "nullMaeMed" "rhoMed" "rhoQ25" "rhoQ75" "pctTop1" "pctBottom1" "nullRhoMed" "nullPctTop1" "nullPctBottom1" ...
            "mantelMed" "mantelSigFrac" "nPairs" "baBias" "baSd" "loaLo" "loaHi" "cpMDC" "cp2xMDC" "tdiParam" "tdiEmp" "nullBias" "nullCpMDC" "nullCp2xMDC" "nullTdiEmp" ...
            "loopCCC" "loopCCCrecomp" "medBetaGenStar" "inDomain"]);
        trials(end+1) = struct('ds', ds, 'model', model, 'cc', cc, 'nM', nM, 'rho', rho, 'nRho', nRho); %#ok<AGROW>
        say("  %-13s %-9s Nv %4d  cccMed %.3f (null %.3f)  rho %.3f (null %.3f)  CP(MDC) %.1f%% (null %.1f%%)  TDI90 %.3f (null %.3f)\n", ds, model, n, median(cc, 'omitnan'), median(nM, 'omitnan'), ...
            median(rho, 'omitnan'), median(nRho, 'omitnan'), cpm, ncp, tde, ntde);
    end
end
T = vertcat(rows{:});
prow = {};
for model = opts.NoiseModels
    for sc = ["inDomain" "exFraser" "all7"]
        sel = arrayfun(@(q) q.model == model && scopeHas_local(q.ds, sc, inDomain), trials);
        cc = vertcat(trials(sel).cc);  nM = vertcat(trials(sel).nM);  rho = vertcat(trials(sel).rho);  nRho = vertcat(trials(sel).nRho);
        bs = bootstrp(opts.NBoot, @(x) median(x, 'omitnan'), cc(isfinite(cc)));
        prow{end+1, 1} = table(model, sc, nnz(isfinite(cc)), median(cc, 'omitnan'), prctile(bs, 2.5), prctile(bs, 97.5), median(nM, 'omitnan'), median(cc - nM, 'omitnan'), median(rho, 'omitnan'), median(nRho, 'omitnan'), ... %#ok<AGROW>
            'VariableNames', ["noiseModel" "scope" "Nvalid" "cccMed" "cccMedLo" "cccMedHi" "nullCccMed" "gainCccMed" "rhoMed" "nullRhoMed"]);
    end
end
PT = vertcat(prow{:});
say("\nper dataset (metric | null):\n%s\n", evalc("disp(T(:, [""dataset"" ""noiseModel"" ""Nvalid"" ""cccMed"" ""nullCccMed"" ""rhoMed"" ""nullRhoMed"" ""cpMDC"" ""nullCpMDC"" ""tdiEmp"" ""nullTdiEmp""]))"));
say("pooled medians of trials:\n%s\n", evalc("disp(PT)"));
results = struct('T', T, 'PT', PT, 'config', struct('MDC', MDC, 'DEG', DEG, 'NPerm', opts.NPerm, 'seed', 1729, 'models', opts.NoiseModels), 'runDate', string(datetime('now'))); %#ok<NASGU>
save(fullfile(resDir, "constellationMetrics_v005" + opts.OutTag + ".mat"), 'results', 'T', 'PT');
say("Saved results/constellationMetrics_v005%s.mat\n=== COMPLETED OK ===\n", opts.OutTag);
end

function tf = scopeHas_local(ds, sc, inDomain)
    switch sc
        case "inDomain", tf = ismember(ds, inDomain);
        case "exFraser", tf = ismember(ds, inDomain) && ds ~= "Fraser";
        otherwise,       tf = true;
    end
end

function lcc = pooledObs_local(r, nPipe)   % loopCCC recomputed from the raw pooled (trial x pipeline) betas, Lin's CCC
    o = NaN(numel(r) * nPipe, 1);  p = o;  k = 0;
    for t = 1:numel(r)
        sl = r(t).betaRecSlice;  if isempty(sl) || all(isnan(sl(:))), continue; end
        for j = 1:nPipe
            bo = r(t).betaObs(j);  bp = median(sl(:, j), 'omitnan');
            if isfinite(bo) && isfinite(bp), k = k + 1;  o(k) = bo;  p(k) = bp; end
        end
    end
    lcc = linCCC_v001(o(1:k), p(1:k)).ccc;
end

function [mR, mP] = mantel_local(bPred, bObs, pP, pQ, nPerm)
    dPred = abs(bPred(pP) - bPred(pQ));  dObs = abs(bObs(pP) - bObs(pQ));
    if std(dPred) == 0 || std(dObs) == 0, mR = NaN;  mP = NaN;  return, end
    oC = dObs - mean(dObs);  oN = sqrt(sum(oC .^ 2));  pC = dPred - mean(dPred);  mR = sum(pC .* oC) / (sqrt(sum(pC .^ 2)) * oN);  nGE = 0;
    for p = 1:nPerm
        bq = bPred(randperm(numel(bPred)));  dp = abs(bq(pP) - bq(pQ));   % permuted pipeline labels
        if std(dp) == 0, continue, end
        pc = dp - mean(dp);  rp = sum(pc .* oC) / (sqrt(sum(pc .^ 2)) * oN);  if abs(rp) >= abs(mR), nGE = nGE + 1; end
    end
    mP = (nGE + 1) / (nPerm + 1);
end

function fprintf_both_local(lg, varargin)
    fprintf(varargin{:});  fprintf(lg, varargin{:});
end

function [pP, pQ] = pairs_local(n)
    pP = [];  pQ = [];
    for p = 2:n, for q = 1:p-1, pP(end+1) = p; pQ(end+1) = q; end, end %#ok<AGROW>
end
