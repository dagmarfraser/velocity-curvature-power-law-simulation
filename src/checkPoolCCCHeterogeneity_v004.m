function R = checkPoolCCCHeterogeneity_v004(opts)
% checkPoolCCCHeterogeneity_v004  Per-pipeline loopCCC (Finding #113) rebuilt on the v015 basis (N_REPS = 200).
% For each core dataset: Lin's CCC of betaObs vs median(betaRecSlice) over trials, separately per pipeline,
% plus the pooled all-six CCC (asserted equal to the stored loopCCC) and the OLS-excluded pool.
% Same pairing rule as computeLoopCCC_local in runLoopClosureFftnoise_v015 (trials with all-NaN betaRecSlice skipped,
% pairs finite). Saves results/checkPoolCCCHeterogeneity_v004.mat (R, T) and a CSV.
arguments
    opts.Datasets (1,:) string = ["Cook_CTRL","Cook_ASD","Hickman_PLAC","Hickman_HALO","Fraser"]
end
srcDir = fileparts(mfilename('fullpath'));  resDir = fullfile(fileparts(srcDir), 'results');
pl = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];  isOLS = contains(pl, "OLS");
rows = cell(0, 10);
for d = opts.Datasets
    f = fullfile(srcDir, sprintf('loopClosureResults_%s_all_shaped_xu_v015.mat', d));
    assert(isfile(f), 'poolCCC:missing', 'missing %s', f);
    S = load(f, 'results', 'config', 'pipelineLabels', 'loopCCC');
    assert(S.config.N_REPS == 200 && S.config.N_BETA == 25 && S.config.RNG_SEED == 1729, 'poolCCC:config', '%s not v015 N_REPS=200', d);
    assert(isequal(string(S.pipelineLabels(:))', pl), 'poolCCC:labels', 'pipeline order differs in %s', d);
    Rr = S.results;  n = numel(Rr);  O = nan(n, 6);  P = nan(n, 6);
    for s = 1:n
        if all(isnan(Rr(s).betaRecSlice(:))), continue; end
        O(s, :) = Rr(s).betaObs(:)';  P(s, :) = median(Rr(s).betaRecSlice, 1, 'omitnan');
    end
    ok = isfinite(O) & isfinite(P);
    pooled = lccc(O(ok), P(ok));
    assert(abs(pooled - S.loopCCC) < 1e-9, 'poolCCC:mismatch', '%s pooled %.6f vs stored loopCCC %.6f', d, pooled, S.loopCCC);
    okNoOLS = ok & ~isOLS;  pooledNoOLS = lccc(O(okNoOLS), P(okNoOLS));
    per = arrayfun(@(k) lccc(O(ok(:,k), k), P(ok(:,k), k)), 1:6);
    nPer = sum(ok, 1);
    rows(end+1, :) = {d, pooled, pooledNoOLS, per(1), per(2), per(3), per(4), per(5), per(6), min(nPer)}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'dataset','pooled','pooledNoOLS','BWFD_OLS','SG_OLS','BWFD_LMLS','SG_LMLS','BWFD_IRLS','SG_IRLS','minN'});
T.dataset = string(T.dataset);  R = T;
disp(T);
save(fullfile(resDir, 'checkPoolCCCHeterogeneity_v004.mat'), 'R', 'T');
writetable(T, fullfile(resDir, 'checkPoolCCCHeterogeneity_v004.csv'));
end

function c = lccc(x, y)
x = x(:); y = y(:);  assert(numel(x) >= 10, 'poolCCC:small', 'fewer than 10 pairs');
mx = mean(x); my = mean(y); sxy = mean((x-mx).*(y-my));
c = 2*sxy / (var(x,0) + var(y,0) + (mx-my)^2);
end
