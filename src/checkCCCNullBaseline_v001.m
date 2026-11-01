function [P, S] = checkCCCNullBaseline_v001(opts)
%CHECKCCCNULLBASELINE_V001  Does cccMed reward trial-specific recovery, or a pipeline-pair signature shared by all trials?  (SESSION_LOG Addendum 53/54)
%   Per valid trial (same construction as checkNoiseModelCCC_v001: 15 inter-pipeline contrasts, Lin's CCC of model-implied vs observed) three scores are compared:
%     model   CCC(median betaRecSlice contrasts, observed contrasts)                          the registered cccMed
%     nullMean CCC(leave-one-out dataset-mean OBSERVED contrast vector, observed contrasts)   a trial-agnostic pipeline signature
%     nullRand CCC(observed contrasts of a random OTHER valid trial of the dataset, observed)  same idea, one random donor (seed 1729)
%   If the nulls are high, cccMed mostly measures the pipeline-pair signature; if near zero, it measures trial-specific recovery. 'gain' = median per-trial (model - nullMean).
%   Also the stacked (all trial x contrast pairs) CCC decomposed into Pearson r, SD ratio (pred/obs) and mean shift, per dataset, in-domain, in-domain without Fraser, and all seven.
%   Same data and fail-loud checks as checkNoiseModelCCC_v001 (N_REPS 200, N_BETA 25, seed 1729, identical trial IDs). Writes ../results/checkCCCNullBaseline_v001.mat and logs/.
%   Fraser, D.S. (2026)  v001
arguments
    opts.Dir (1,1) string = ""
    opts.Datasets (1,:) string = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Dhieb" "Zarandi"]
    opts.OutTag (1,1) string = ""
    opts.FromCache (1,1) logical = false   % reuse results/checkCCCNullBaseline_v001_POOL.mat (per-trial scores) and redo only the aggregation
end
srcDir = fileparts(mfilename('fullpath'));  addpath(genpath(fullfile(srcDir, 'functions')));
dataDir = opts.Dir;  if dataDir == "", dataDir = string(srcDir); end
resDir = fullfile(srcDir, '..', 'results');  logDir = fullfile(srcDir, 'logs');
for d = {resDir, logDir}, if ~isfolder(d{1}), mkdir(d{1}); end, end
lg = fopen(fullfile(logDir, "checkCCCNullBaseline_v001_" + string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) + ".log"), 'w');  cleanup = onCleanup(@() fclose(lg)); %#ok<NASGU>
say = @(varargin) fprintf_both_local(lg, varargin{:});
MDC = 0.03;  DEG = 0.005;  nPipe = 6;  [pP, pQ] = pairs_local(nPipe);  rng(1729);
models = ["shaped_xu" "xu" "bootstrap"];  vers = ["v015" "v017" "v017"];
inDomain = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO"];
POOL = struct('ds', {}, 'model', {}, 'cc', {}, 'nM', {}, 'nR', {}, 'dO', {}, 'dP', {});
say("checkCCCNullBaseline_v001 | %s | data %s\n", string(datetime('now')), dataDir);
cacheF = fullfile(resDir, 'checkCCCNullBaseline_v001_POOL.mat');
if opts.FromCache, C = load(cacheF, 'POOL');  POOL = C.POOL;  loadList = strings(1, 0);  else, loadList = opts.Datasets;  end
for d = loadList
    ids = [];
    for m = 1:numel(models)
        f = fullfile(dataDir, "loopClosureResults_" + d + "_all_" + models(m) + "_" + vers(m) + ".mat");
        if ~isfile(f), error('nullCCC:missing', 'No file: %s', f); end
        S0 = load(f, 'results', 'config');  r = S0.results;  N = numel(r);  c = S0.config;
        if c.N_REPS ~= 200 || c.N_BETA ~= 25 || c.RNG_SEED ~= 1729, error('nullCCC:config', '%s: N_REPS/N_BETA/seed differ from 200/25/1729', f); end
        tid = string({r.trialID})';  if m == 1, ids = tid; elseif ~isequal(tid, ids), error('nullCCC:trials', '%s: trial IDs differ from the shaped_xu file', f); end
        dbO = NaN(N, numel(pP));  dbP = dbO;  ok = false(N, 1);
        for t = 1:N
            bo = r(t).betaObs(:)';  sl = r(t).betaRecSlice;
            if numel(bo) ~= nPipe || any(~isfinite(bo)) || isempty(sl) || size(sl, 2) ~= nPipe || ~isfinite(r(t).betaGenStarMed), continue; end
            bp = median(sl, 1, 'omitnan');  if any(~isfinite(bp)), continue; end
            dbO(t, :) = bo(pP) - bo(pQ);  dbP(t, :) = bp(pP) - bp(pQ);  ok(t) = true;
        end
        valid = ok & max(abs(dbO), [], 2) >= DEG & ~all(abs(dbO) < MDC, 2);  vi = find(valid);  n = numel(vi);
        O = dbO(vi, :);  Pm = dbP(vi, :);  sumO = sum(O, 1);  cc = NaN(n, 1);  nM = cc;  nR = cc;
        for i = 1:n
            cc(i) = linCCC_v001(Pm(i, :)', O(i, :)').ccc;
            if n >= 2
                nM(i) = linCCC_v001(((sumO - O(i, :)) / (n - 1))', O(i, :)').ccc;
                j = randi(n - 1);  j = j + (j >= i);  nR(i) = linCCC_v001(O(j, :)', O(i, :)').ccc;
            end
        end
        POOL(end+1) = struct('ds', d, 'model', models(m), 'cc', cc, 'nM', nM, 'nR', nR, 'dO', O, 'dP', Pm); %#ok<AGROW>
        say("  %-13s %-9s valid %4d  model %.3f  nullMean %.3f  nullRand %.3f\n", d, models(m), n, median(cc, 'omitnan'), median(nM, 'omitnan'), median(nR, 'omitnan'));
    end
end
if ~opts.FromCache, save(cacheF, 'POOL', '-v7.3'); end
names = [opts.Datasets, "inDomain", "exFraser", "all7"];
pr = {};  sr = {};
for sc = names
    for m = 1:numel(models)
        sel = arrayfun(@(q) q.model == models(m) && inScope_local(q.ds, sc, inDomain), POOL);
        if ~any(sel), continue; end
        cc = vertcat(POOL(sel).cc);  nM = vertcat(POOL(sel).nM);  nR = vertcat(POOL(sel).nR);  O = vertcat(POOL(sel).dO);  Pm = vertcat(POOL(sel).dP);
        g = cc - nM;  lo = NaN;  hi = NaN;   % CI only with >= 10 valid trials (Zarandi has 1: reported without CI, not faked)
        if nnz(isfinite(cc)) >= 10, bs = bootstrp(1000, @(x) median(x, 'omitnan'), cc(isfinite(cc)));  lo = prctile(bs, 2.5);  hi = prctile(bs, 97.5); end
        pr{end+1, 1} = table(sc, models(m), numel(cc), median(cc, 'omitnan'), lo, hi, median(nM, 'omitnan'), median(nR, 'omitnan'), median(g, 'omitnan'), mean(g > 0, 'omitnan'), ... %#ok<AGROW>
            'VariableNames', ["scope" "noiseModel" "Nvalid" "cccMed" "ciLo" "ciHi" "nullMeanMed" "nullRandMed" "gainMed" "fracBeatNull"]);
        q = linCCC_v001(Pm(:), O(:));
        sr{end+1, 1} = table(sc, models(m), numel(O(:, 1)), q.ccc, q.r, std(Pm(:)) / std(O(:)), mean(Pm(:)) - mean(O(:)), ... %#ok<AGROW>
            'VariableNames', ["scope" "noiseModel" "Nvalid" "stackedCCC" "pearsonR" "sdRatio" "meanShift"]);
    end
end
P = vertcat(pr{:});  S = vertcat(sr{:});
say("\nper-trial CCC: model vs trial-agnostic nulls (gain = median of per-trial model - nullMean; fracBeatNull = share of trials where model > nullMean):\n%s\n", evalc("disp(P)"));
say("stacked CCC decomposition (sdRatio = SD pred / SD obs; meanShift = mean pred - mean obs):\n%s\n", evalc("disp(S)"));
runDate = string(datetime('now', 'Format', 'yyyy-MM-dd'));  seed = 1729;
save(fullfile(resDir, "checkCCCNullBaseline_v001" + opts.OutTag + ".mat"), 'P', 'S', 'runDate', 'seed');
say("Saved results/checkCCCNullBaseline_v001%s.mat\n=== COMPLETED OK ===\n", opts.OutTag);
end

function tf = inScope_local(ds, sc, inDomain)
    switch sc
        case "inDomain", tf = ismember(ds, inDomain);
        case "exFraser", tf = ismember(ds, inDomain) && ds ~= "Fraser";
        case "all7",     tf = true;
        otherwise,       tf = ds == sc;
    end
end

function fprintf_both_local(lg, varargin)
    fprintf(varargin{:});  fprintf(lg, varargin{:});
end

function [pP, pQ] = pairs_local(n)
    pP = [];  pQ = [];
    for p = 2:n, for q = 1:p-1, pP(end+1) = p; pQ(end+1) = q; end, end %#ok<AGROW>
end
