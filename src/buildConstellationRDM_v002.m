function results = buildConstellationRDM_v002(deltaBetaPred, deltaBetaObs, validMask, options)
%BUILDCONSTELLATIONRDM_V002 Trial x trial RDMs from constellation vectors,
% Mantel-tested (prereg_v101.docx Section 7.5, "Second-order structural
% validation" / RSA).
%
% v002 change from v001 (2026-07-25/26, Fraser D.S.): v001's parfor
% permutation routine reconstructed the FULL n x n permuted distance
% matrix (`D2(permIdx,permIdx)`) on every single permutation, then
% extracted the lower triangle from that. At Pilot's full uncapped scale
% (n~4237), that is a fresh ~143MB, 17.9-million-element matrix allocated
% and filled per permutation -- wasteful, since only the ~9-million-
% element lower triangle is ever used. v002 precomputes the lower-
% triangle's (row, col) subscripts ONCE, then gathers each permutation's
% values directly via linear indexing into the original D2 -- no
% intermediate n x n matrix ever materialised inside the loop.
%
% This does NOT change what is computed, only how -- the consistency
% check against mantelTest's own serial statistic (unchanged, still
% present below) is exactly what would catch a mistake in this
% rewrite, not merely inspection.
%
% SEPARATE, LIKELY LARGER SUSPECT (not something this file can fix):
% every parfor worker needs its own copy of the full D2 matrix (a
% broadcast variable) regardless of this optimisation -- for Pilot,
% that is ~143MB x (number of workers). At 71 workers that is ~10GB for
% D2 alone, before any per-worker MATLAB runtime overhead or working
% memory. If the BlueBEAR job's memory request was interpreted as a
% TOTAL (e.g. a plain `--mem=4G` read as per-node rather than
% per-task), 71 workers could not even hold their own broadcast copies
% of D2, let alone compute anything -- this would present exactly as
% "hours, no output," via swapping, not as an error. Worth confirming
% the actual per-task memory available before assuming this rewrite
% alone fixes the incident that motivated it. See
% docs/TODO_ConstellationRSA_v001.md for the live incident notes.
%
% Everything else (Fail Loud checks, RDM construction, output fields,
% documentation of the class-level design) is unchanged from v001; see
% that file's header for the full original documentation.
%
% Fraser, D.S. (2026) v002

arguments
    deltaBetaPred (:, 15) double
    deltaBetaObs (:, 15) double
    validMask (:, 1) logical
    options.Method (1, 1) string {mustBeMember(options.Method, ["spearman", "pearson"])} = "spearman"
    options.NPerm (1, 1) double {mustBeInteger, mustBePositive} = 10000
    options.Exact (1, 1) string {mustBeMember(options.Exact, ["auto", "on", "off"])} = "off"
    options.UseParfor (1, 1) logical = false
    options.Label (1, 1) string = "unnamed"
end

if exist("mantelTest", "file") ~= 2
    error("buildConstellationRDM_v002:mantelTestNotFound", "%s", ...
        "mantelTest.m not found on path. addpath the mantel toolbox " + ...
        "(stagingForGithub/mantel/toolbox) before calling this function.");
end

idx = find(validMask);
n = numel(idx);
if n < 4
    error("buildConstellationRDM_v002:tooFewTrials", "%s", ...
        sprintf("Only %d valid trials for '%s'; need at least 4 to build an RDM.", ...
        n, options.Label));
end

P = deltaBetaPred(idx, :);
O = deltaBetaObs(idx, :);
if any(~isfinite(P(:))) || any(~isfinite(O(:)))
    error("buildConstellationRDM_v002:nonfiniteInput", "%s", ...
        sprintf("'%s': validMask selected trials contain non-finite " + ...
        "deltaBeta values -- pass a mask that excludes these.", options.Label));
end

RDM_pred = squareEuclideanRDM_local(P);
RDM_obs  = squareEuclideanRDM_local(O);

triIdx = triu(true(n), 1);
cvPred = std(RDM_pred(triIdx)) / mean(RDM_pred(triIdx));
cvObs  = std(RDM_obs(triIdx)) / mean(RDM_obs(triIdx));

if options.UseParfor
    % Cheap serial call (NPerm=1) purely to get mantelTest's own mantelR
    % for the consistency check below -- NPerm barely affects runtime,
    % only the O(n^2) statistic computation and one trivial permutation.
    [mantelR_toolbox, ~, ~] = mantelTest(RDM_pred, RDM_obs, ...
        Method=options.Method, NPerm=1, Exact="off");

    [mantelR, pValue, nPermUsed] = mantelPermuteParfor_local( ...
        RDM_pred, RDM_obs, options.Method, options.NPerm);

    if ~isequaln(mantelR, mantelR_toolbox) && ...
            abs(mantelR - mantelR_toolbox) > 1e-9
        error("buildConstellationRDM_v002:parforDrift", "%s", ...
            sprintf("'%s': parfor-path mantelR (%.10f) disagrees with " + ...
            "mantelTest's own serial mantelR (%.10f) by more than 1e-9 " + ...
            "-- the linear-indexing rewrite in this file has drifted " + ...
            "from mantelTest.m's math. Fix before trusting either result.", ...
            options.Label, mantelR, mantelR_toolbox));
    end
else
    [mantelR, pValue, info] = mantelTest(RDM_pred, RDM_obs, ...
        Method=options.Method, NPerm=options.NPerm, Exact=options.Exact);
    nPermUsed = info.nPermUsed;
end

results = struct();
results.label      = options.Label;
results.nTrials    = n;
results.mantelR    = mantelR;
results.pValue     = pValue;
results.nPermUsed  = nPermUsed;
results.usedParfor = options.UseParfor;
results.cvPred     = cvPred;
results.cvObs      = cvObs;
results.RDM_pred   = RDM_pred;
results.RDM_obs    = RDM_obs;

end

%% ============================ LOCAL HELPERS ================================

function D = squareEuclideanRDM_local(X)
% Hand-rolled pairwise Euclidean distance matrix, no Statistics and
% Machine Learning Toolbox pdist/squareform dependency (matches the
% mantel toolbox's own dependency-light ethos, even though this is
% project-internal code, not redistributed).
n = size(X, 1);
D = zeros(n, n);
for i = 1:n-1
    diffs = X(i+1:n, :) - X(i, :);
    d = sqrt(sum(diffs .^ 2, 2));
    D(i+1:n, i) = d;
    D(i, i+1:n) = d';
end
end

function [mantelR, pValue, nUsed] = mantelPermuteParfor_local(D1, D2, method, nPerm)
% Project-local parfor permutation routine, deliberately NOT part of the
% mantel toolbox (kept toolbox-serial per design decision). Replicates
% mantelTest.m's own math exactly (same lower-triangle extraction, same
% rank-once-for-the-fixed-side optimisation, same p-value formula
% (nGE+1)/(nUsed+1) with skipped degenerate permutations still counted
% in nUsed) -- verified against mantelTest's own serial mantelR by the
% caller's consistency check, not assumed correct on inspection alone.
%
% v002: no full n x n permuted matrix is ever built. (iLow, jLow) are
% the row/col subscripts of the lower triangle, computed once; each
% permutation gathers directly via linear indexing into the ORIGINAL D2
% -- D2(sub2ind([n n], permIdx(iLow), permIdx(jLow))) -- which is
% mathematically identical to D2(permIdx,permIdx) followed by the same
% lower-triangle extraction, just without the intermediate full matrix.
n = size(D1, 1);
lowerMask = tril(true(n), -1);
[iLow, jLow] = find(lowerMask);
v1 = D1(lowerMask);
v2 = D2(lowerMask);

if method == "spearman"
    v1rank = rankVector_local(v1);
else
    v1rank = v1;
end
C0 = corrcoef(v1rank, rankIfSpearman_local(v2, method));
mantelR = C0(1, 2);
absR = abs(mantelR);

hits = zeros(nPerm, 1);
parfor k = 1:nPerm
    permIdx = randperm(n);
    linIdx = sub2ind([n, n], permIdx(iLow), permIdx(jLow)); %#ok<PFBNS>
    v2perm = D2(linIdx);
    if std(v2perm) == 0
        hits(k) = 0;
        continue
    end
    v2p = rankIfSpearman_local(v2perm, method);
    Ck = corrcoef(v1rank, v2p);
    hits(k) = abs(Ck(1, 2)) >= absR;
end

nGE = sum(hits);
nUsed = nPerm;
pValue = (nGE + 1) / (nUsed + 1);
end

function y = rankIfSpearman_local(x, method)
if method == "spearman"
    y = rankVector_local(x);
else
    y = x;
end
end

function r = rankVector_local(x)
% Verbatim duplicate of mantel toolbox's own mantelRank (average-rank
% transform, ties handled), needed here since that function is private
% to mantelTest.m and not exported. Kept deliberately identical -- do
% not "improve" this independently of the toolbox's own copy, or the
% consistency check in buildConstellationRDM_v002 will start failing for
% the wrong reason.
x = x(:);
n = numel(x);
[sortedX, sortIdx] = sort(x);
% Vectorized average-rank tie correction (patched 2026-07-27, Fraser D.S.):
% the previous per-unique-value loop (tiedMask = ic==k, looped once per
% unique value) is O(n^2) -- fine at small n, but at RDM lower-triangle
% scale (n ~ 9e6 for Pilot's full n=4237 trials) it never completes in
% practical time. This groups by consecutive equal sorted values and
% averages the sorted-order rank position (1:n) within each group via
% accumarray -- mathematically identical output, O(n log n) via sort.
isNewGroup = [true; diff(sortedX) ~= 0];
groupID = cumsum(isNewGroup);
positions = (1:n)';
meanRankPerGroup = accumarray(groupID, positions) ./ accumarray(groupID, 1);
r = zeros(n, 1);
r(sortIdx) = meanRankPerGroup(groupID);
end
