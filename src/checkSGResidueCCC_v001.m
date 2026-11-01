function T = checkSGResidueCCC_v001(ver)
% checkSGResidueCCC_v001  Sol's closure 1, second half: do the primary metric (per-trial contrast CCC, cccMed),
% the dataset constellation median of beta_gen* and the invertible share move when the SG curvature-residue
% trials (checkSGCurvatureResidue_v001, all 79 listed) are left out? Same construction as
% checkPart3PrimaryMetric_v001 (unchanged logic); the residue trials are removed from the whole trial set.
% Reads src/loopClosureResults_<dataset>_all_shaped_xu_<ver>.mat and results/checkSGCurvatureResidue_v001.mat.
% Writes results/checkSGResidueCCC_v001_<ver>.mat. Fails loud on a missing file.
arguments, ver (1,1) string = "v015"; end
srcDir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(srcDir, 'functions')));
MDC = 0.03; DEG = 0.005; nPipe = 6;
ds  = ["Fraser","Cook CTRL","Cook ASD","Hickman HALO","Hickman PLAC","Dhieb","Zarandi"];
[pP, pQ] = pairs_local(nPipe);
inRS = fullfile(srcDir, '..', 'results', 'checkSGCurvatureResidue_v001.mat');
if ~isfile(inRS), error('sgCCC:missing', 'No file: %s', inRS); end
RS = load(inRS, 'perTrial');  resid = strings(0, 1);
for fn = string(fieldnames(RS.perTrial))', resid = [resid; RS.perTrial.(fn).trialID]; end
rows = cell(numel(ds), 1);
for d = 1:numel(ds)
    tag = replace(ds(d), " ", "_");
    f = fullfile(srcDir, sprintf('loopClosureResults_%s_all_shaped_xu_%s.mat', tag, ver));
    if ~isfile(f), error('sgCCC:missing', 'No file: %s', f); end
    r = load(f, 'results').results;  N = numel(r);
    isR = ismember(string({r.trialID})', resid);
    dbO = NaN(N, numel(pP));  dbP = dbO;  ok = false(N, 1);  inv = false(N, 1);  gMed = NaN(N, 1);
    for t = 1:N
        bo = r(t).betaObs(:)';  sl = r(t).betaRecSlice;
        inv(t) = isfinite(r(t).betaGenStarMed);  gMed(t) = r(t).betaGenStarMed;
        if numel(bo) ~= nPipe || any(~isfinite(bo)) || isempty(sl) || size(sl, 2) ~= nPipe || ~inv(t), continue; end
        bp = median(sl, 1, 'omitnan');  if any(~isfinite(bp)), continue; end
        dbO(t, :) = bo(pP) - bo(pQ);  dbP(t, :) = bp(pP) - bp(pQ);  ok(t) = true;
    end
    valid = ok & max(abs(dbO), [], 2) >= DEG & ~all(abs(dbO) < MDC, 2);
    cA = ccc_local(dbP, dbO, valid);  cD = ccc_local(dbP, dbO, valid & ~isR);
    rows{d} = table(ds(d), N, nnz(isR), nnz(valid), nnz(valid & ~isR), cA, cD, cD - cA, ...
        median(gMed, 'omitnan'), median(gMed(~isR), 'omitnan'), 100 * mean(inv), 100 * mean(inv(~isR)), ...
        'VariableNames', ["dataset","N","nResid","NvalidAll","NvalidDrop","cccMedAll","cccMedDrop","dCCC", ...
        "medGenAll","medGenDrop","pctInvAll","pctInvDrop"]);
end
T = vertcat(rows{:});  T.dMedGen = T.medGenDrop - T.medGenAll;  T.dPctInv = T.pctInvDrop - T.pctInvAll;
fprintf('Residue-trial sensitivity, shaped_xu %s\n', ver);  disp(T);
fprintf('Largest |change|: cccMed %.4f, median beta_gen* %.4f, invertible share %.2f pp\n', ...
    max(abs(T.dCCC)), max(abs(T.dMedGen)), max(abs(T.dPctInv)));
save(fullfile(srcDir, '..', 'results', sprintf('checkSGResidueCCC_v001_%s.mat', ver)), 'T');
end

function m = ccc_local(dbP, dbO, valid)
vi = find(valid);  c = NaN(numel(vi), 1);
for i = 1:numel(vi), c(i) = linCCC_v001(dbP(vi(i), :)', dbO(vi(i), :)').ccc; end
m = median(c, 'omitnan');
end

function [pP, pQ] = pairs_local(n)
pP = [];  pQ = [];
for p = 2:n, for q = 1:p-1, pP(end+1) = p; pQ(end+1) = q; end, end %#ok<AGROW>
end
