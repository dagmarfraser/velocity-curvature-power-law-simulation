function T = checkNoiseModelCCC_v001(opts)
%CHECKNOISEMODELCCC_V001  Three-noise-model primary-metric table on the v015 basis (replaces Supplement A's v007/v008 table; SESSION_LOG Addendum 52).
%   For each dataset and noise model the registered primary metric is recomputed exactly as checkPart3PrimaryMetric_v001 (per trial, the 15 inter-pipeline contrasts of
%   median(betaRecSlice) vs betaObs, Lin's CCC; cccMed = median over valid trials), plus the secondary loopCCC, the invertible share and the median model-implied exponent.
%   Files (all trials, N_REPS = 200, N_BETA = 25, seed 1729, same trials and pipelines; asserted):  shaped_xu: loopClosureResults_<ds>_all_shaped_xu_v015.mat;
%   xu and bootstrap: loopClosureResults_<ds>_all_{xu,bootstrap}_v017.mat (the D-H1 runs, runH1HPC_v002). Fails loud on a missing file or a config / trial-ID mismatch.
%   Run on the HPC (or any machine) from the folder holding the files; no parpool, a few minutes. Writes ../results/checkNoiseModelCCC_v001.mat and logs/checkNoiseModelCCC_v001_<stamp>.log.
%   USAGE:  checkNoiseModelCCC_v001            % src folder = this file's folder
%           checkNoiseModelCCC_v001(Dir="<folder with the loopClosureResults files>")
%   Fraser, D.S. (2026)  v001
arguments
    opts.Dir (1,1) string = ""
    opts.Datasets (1,:) string = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Dhieb" "Zarandi"]
    opts.OutTag (1,1) string = ""
end
srcDir = fileparts(mfilename('fullpath'));  addpath(genpath(fullfile(srcDir, 'functions')));
dataDir = opts.Dir;  if dataDir == "", dataDir = string(srcDir); end
resDir = fullfile(srcDir, '..', 'results');  logDir = fullfile(srcDir, 'logs');
for d = {resDir, logDir}, if ~isfolder(d{1}), mkdir(d{1}); end, end
lg = fopen(fullfile(logDir, "checkNoiseModelCCC_v001_" + string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) + ".log"), 'w');  cleanup = onCleanup(@() fclose(lg)); %#ok<NASGU>
say = @(varargin) fprintf_both_local(lg, varargin{:});
MDC = 0.03;  DEG = 0.005;  nPipe = 6;  [pP, pQ] = pairs_local(nPipe);
models = ["shaped_xu" "xu" "bootstrap"];  vers = ["v015" "v017" "v017"];
inDomain = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO"];   % Dhieb and Zarandi are outside the validity domain (Findings #223, #242)
rows = {};  POOL = struct('cc', {}, 'dO', {}, 'dP', {}, 'ds', {}, 'model', {});  say("checkNoiseModelCCC_v001 | %s | data %s\n", string(datetime('now')), dataDir);
for d = opts.Datasets
    ids = [];  cfg0 = [];
    for m = 1:numel(models)
        f = fullfile(dataDir, "loopClosureResults_" + d + "_all_" + models(m) + "_" + vers(m) + ".mat");
        if ~isfile(f), error('nmCCC:missing', 'No file: %s', f); end
        S = load(f, 'results', 'config', 'loopCCC');  r = S.results;  N = numel(r);
        c = S.config;
        if c.N_REPS ~= 200 || c.N_BETA ~= 25 || c.RNG_SEED ~= 1729, error('nmCCC:config', '%s: N_REPS %g, N_BETA %g, seed %g differ from 200 / 25 / 1729', f, c.N_REPS, c.N_BETA, c.RNG_SEED); end
        tid = string({r.trialID})';
        if m == 1, ids = tid; cfg0 = c; elseif ~isequal(tid, ids), error('nmCCC:trials', '%s: trial IDs differ from the shaped_xu file', f); end
        for fld = ["BW_CUTOFF" "SG_REF" "EDGE_CLIP" "FS_DS"]   % compared where both configs hold the field; a field absent from either is disclosed, not skipped silently
            if isfield(c, fld) && isfield(cfg0, fld)
                if ~isequal(c.(fld), cfg0.(fld)), error('nmCCC:config', '%s: %s differs from the shaped_xu file', f, fld); end
            elseif m > 1, say("  NOTE %s %s: config field %s not in both files (shaped_xu has it: %d, this file: %d); not compared\n", d, models(m), fld, isfield(cfg0, fld), isfield(c, fld));
            end
        end
        dbO = NaN(N, numel(pP));  dbP = dbO;  ok = false(N, 1);  inv = false(N, 1);  gMed = NaN(N, 1);
        for t = 1:N
            bo = r(t).betaObs(:)';  sl = r(t).betaRecSlice;  inv(t) = isfinite(r(t).betaGenStarMed);  gMed(t) = r(t).betaGenStarMed;
            if numel(bo) ~= nPipe || any(~isfinite(bo)) || isempty(sl) || size(sl, 2) ~= nPipe || ~inv(t), continue; end
            bp = median(sl, 1, 'omitnan');  if any(~isfinite(bp)), continue; end
            dbO(t, :) = bo(pP) - bo(pQ);  dbP(t, :) = bp(pP) - bp(pQ);  ok(t) = true;
        end
        valid = ok & max(abs(dbO), [], 2) >= DEG & ~all(abs(dbO) < MDC, 2);  vi = find(valid);  cc = NaN(numel(vi), 1);  pr = cc;  mae = cc;
        for i = 1:numel(vi), q = linCCC_v001(dbP(vi(i), :)', dbO(vi(i), :)');  cc(i) = q.ccc;  pr(i) = q.r;  mae(i) = mean(abs(dbP(vi(i), :) - dbO(vi(i), :))); end
        POOL(end+1) = struct('cc', {cc}, 'dO', {dbO(vi, :)}, 'dP', {dbP(vi, :)}, 'ds', {d}, 'model', {models(m)}); %#ok<AGROW>
        rows{end+1, 1} = table(d, models(m), vers(m), N, nnz(inv), numel(vi), median(cc, 'omitnan'), median(pr, 'omitnan'), median(mae, 'omitnan'), S.loopCCC, median(gMed, 'omitnan'), ismember(d, inDomain), ... %#ok<AGROW>
            'VariableNames', ["dataset" "noiseModel" "version" "Ntrials" "Ninvertible" "Nvalid" "cccMed" "pearsonMed" "maeMed" "loopCCC" "medBetaGenStar" "inDomain"]);
        say("  %-13s %-9s %s  N %4d  valid %4d  cccMed %.3f  loopCCC %.3f  median beta_gen* %.3f\n", d, models(m), vers(m), N, numel(vi), median(cc, 'omitnan'), S.loopCCC, median(gMed, 'omitnan'));
    end
end
T = vertcat(rows{:});
W = unstack(T(:, ["dataset" "noiseModel" "cccMed"]), "cccMed", "noiseModel");  W = W(:, ["dataset" models]);
say("\nprimary metric cccMed by noise model (v015 basis):\n%s\n", evalc("disp(W)"));
% Pooled across trials (not a random-effects estimate): median of per-trial CCC over all valid trials of the in-domain datasets (and over all seven, disclosed), with a trial-level
% bootstrap CI (seed 1729), plus the CCC of ALL stacked trial x contrast pairs. The stacked CCC is driven by between-trial (hidden) structure; the per-trial median is not.
rng(1729);  prow = {};
for m = 1:numel(models)
    for scope = ["inDomain" "all7"]
        sel = arrayfun(@(q) q.model == models(m) && (scope == "all7" || ismember(q.ds, inDomain)), POOL);
        cc = vertcat(POOL(sel).cc);  dO = vertcat(POOL(sel).dO);  dP = vertcat(POOL(sel).dP);  cc = cc(isfinite(cc));
        bs = bootstrp(2000, @(x) median(x, 'omitnan'), cc);  q = linCCC_v001(dP(:), dO(:));
        prow{end+1, 1} = table(models(m), scope, numel(cc), median(cc), prctile(bs, 2.5), prctile(bs, 97.5), q.ccc, ... %#ok<AGROW>
            'VariableNames', ["noiseModel" "scope" "Nvalid" "pooledCccMed" "ciLo" "ciHi" "stackedCCC"]);
    end
end
P = vertcat(prow{:});  say("\npooled across trials (median of per-trial CCC, bootstrap 95%% CI; stackedCCC = CCC of all stacked contrasts):\n%s\n", evalc("disp(P)"));
runDate = string(datetime('now', 'Format', 'yyyy-MM-dd'));
save(fullfile(resDir, "checkNoiseModelCCC_v001" + opts.OutTag + ".mat"), 'T', 'W', 'P', 'runDate');
say("Saved results/checkNoiseModelCCC_v001%s.mat\n=== COMPLETED OK ===\n", opts.OutTag);
end

function fprintf_both_local(lg, varargin)
    fprintf(varargin{:});  fprintf(lg, varargin{:});
end

function [pP, pQ] = pairs_local(n)
    pP = [];  pQ = [];
    for p = 2:n, for q = 1:p-1, pP(end+1) = p; pQ(end+1) = q; end, end %#ok<AGROW>
end
