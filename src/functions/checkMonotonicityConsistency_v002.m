function result = checkMonotonicityConsistency_v002(xGrid, yReplicates, betaObs, params)
%CHECKMONOTONICITYCONSISTENCY_V002 Same replicate-based empirical
% confidence check as v001, but each replicate is classified with
% findBothBranches_v008 (noise-tolerant trim, Finding #153's fix)
% instead of v004 (zero-tolerance trim).
%
% sigmaY (the noise-tolerance input v008 needs) is computed HERE, once,
% across all N_REPS replicates at each grid point -- shared context for
% every individual replicate's own trim step, the same role it played in
% the (unsuccessful) detection-side fixes v005/v006/v007, but applied to
% a different step this time.
%
% v001 is left completely unmodified.
%
% See also: checkMonotonicityConsistency_v001 (uses v004, zero-tolerance
% trim), findBothBranches_v008 (the per-replicate classifier this calls)
%
% Fraser, D.S. (2026)

nReps = size(yReplicates, 1);
perRepCovered = false(nReps, 1);

sigmaY = std(yReplicates, 0, 1, 'omitnan');  % 1 x N_BETA, across all replicates

for r = 1:nReps
    yRep = yReplicates(r, :);
    good = ~isnan(yRep) & ~isnan(sigmaY);
    if nnz(good) < 2, continue, end

    both = findBothBranches_v008(xGrid(good), yRep(good), sigmaY(good), params);
    allRuns = [both.rise, both.desc];
    for k = 1:numel(allRuns)
        if betaObs >= min(allRuns(k).y) - eps && betaObs <= max(allRuns(k).y) + eps
            perRepCovered(r) = true;
        end
    end
end

result.fracCovered   = mean(perRepCovered);
result.nCovered      = sum(perRepCovered);
result.nReps         = nReps;
result.perRepCovered = perRepCovered;

end
