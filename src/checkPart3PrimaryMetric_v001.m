function T = checkPart3PrimaryMetric_v001(ver)
% checkPart3PrimaryMetric_v001  Registered primary metric (per-trial contrast CCC, cccMed)
% recomputed from a chosen shaped_xu loop-closure version, no re-simulation.
% Same construction as constellationMetrics_v004 (CCC_MEASURES_REFERENCE, measure 4):
% per trial, 15 inter-pipeline contrasts of median(betaRecSlice) vs betaObs, Lin's CCC.
% Default ver = "v015" (checked, N_REPS = 200). v004's table used v007/v008 (N_REPS = 20).
% Fails loud on a missing file or a config mismatch. Writes ../results/.
arguments, ver (1,1) string = "v015"; end
srcDir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(srcDir, 'functions')));
MDC = 0.03; DEG = 0.005; nPipe = 6;
ds  = ["Fraser","Cook CTRL","Cook ASD","Hickman HALO","Hickman PLAC","Dhieb","Zarandi"];
[pP, pQ] = pairs_local(nPipe);
rows = cell(numel(ds), 1);
for d = 1:numel(ds)
    tag = replace(ds(d), " ", "_");
    f = fullfile(srcDir, sprintf('loopClosureResults_%s_all_shaped_xu_%s.mat', tag, ver));
    if ~isfile(f), error('checkPart3PrimaryMetric_v001:missing', 'No file: %s', f); end
    S = load(f, 'results', 'config', 'loopCCC');  r = S.results;  N = numel(r);
    dbO = NaN(N, numel(pP));  dbP = dbO;  ok = false(N, 1);  inv = false(N, 1);
    for t = 1:N
        bo = r(t).betaObs(:)';  sl = r(t).betaRecSlice;
        inv(t) = isfinite(r(t).betaGenStarMed);
        if numel(bo) ~= nPipe || any(~isfinite(bo)) || isempty(sl) || size(sl, 2) ~= nPipe || ~inv(t), continue; end
        bp = median(sl, 1, 'omitnan');  if any(~isfinite(bp)), continue; end
        dbO(t, :) = bo(pP) - bo(pQ);  dbP(t, :) = bp(pP) - bp(pQ);  ok(t) = true;
    end
    valid = ok & max(abs(dbO), [], 2) >= DEG & ~all(abs(dbO) < MDC, 2);
    vi = find(valid);  c = NaN(numel(vi), 1);  pr = c;  mae = c;
    for i = 1:numel(vi)
        cc = linCCC_v001(dbP(vi(i), :)', dbO(vi(i), :)');  c(i) = cc.ccc;  pr(i) = cc.r;
        mae(i) = mean(abs(dbP(vi(i), :) - dbO(vi(i), :)));
    end
    rows{d} = table(ds(d), N, nnz(inv), numel(vi), median(c, 'omitnan'), median(pr, 'omitnan'), ...
        median(mae, 'omitnan'), S.loopCCC, S.config.N_REPS, ...
        'VariableNames', ["dataset","Ntrials","Ninvertible","Nvalid","cccMed","pearsonMed","maeMed","loopCCC","NREPS"]);
end
T = vertcat(rows{:});
T.gradeRegistered = arrayfun(@gradeReg_local, T.cccMed);
T.gradeMcBride    = arrayfun(@gradeMcB_local, T.cccMed);
fprintf('Per-trial contrast CCC (primary), shaped_xu %s\n', ver);  disp(T);
save(fullfile(srcDir, '..', 'results', sprintf('checkPart3PrimaryMetric_v001_%s.mat', ver)), 'T');
end

function [pP, pQ] = pairs_local(n)
pP = [];  pQ = [];
for p = 2:n, for q = 1:p-1, pP(end+1) = p; pQ(end+1) = q; end, end %#ok<AGROW>
end

function g = gradeReg_local(x)   % prereg 7.5 convention
g = "n/a";  if isnan(x), return, end
g = "Poor";  if x >= 0.5, g = "Moderate"; end;  if x >= 0.7, g = "Good"; end;  if x >= 0.9, g = "Excellent"; end
end

function g = gradeMcB_local(x)   % McBride (2005) stricter scale
g = "n/a";  if isnan(x), return, end
g = "Poor";  if x >= 0.9, g = "Moderate"; end;  if x >= 0.95, g = "Substantial"; end;  if x >= 0.99, g = "Almost perfect"; end
end
