function [T6a, TNC, info] = buildTable6aV015_v001(ver)
% Table 6a (beta_obs by pipeline: all-trial median; mean+-SD) and its companion (per-pipeline
% non-convergence = share of trials whose beta_gen* is NaN) from loopClosureResults_<ds>_all_shaped_xu_<ver>.mat.
% Default ver "v015" (N_REPS = 200). ver "v007" (Fraser "v008") reproduces the N_REPS = 20 table as a check.
% Output: results/buildTable6aV015_v001_<ver>.mat (T6a, TNC, info) and CSVs.
if nargin < 1, ver = "v015"; end
ver = string(ver);
srcDir = fileparts(mfilename('fullpath'));  resDir = fullfile(fileparts(srcDir), 'results');
ds = ["Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Fraser"];
pl = ["BWFD-OLS" "SG-OLS" "BWFD-IRLS" "SG-IRLS" "BWFD-LMLS" "SG-LMLS"];      % Table 6a row order
rows = {};  nc = nan(numel(pl), numel(ds));  N = zeros(1, numel(ds));
for j = 1:numel(ds)
    v = ver;  if ver == "v007" && ds(j) == "Fraser", v = "v008"; end
    f = fullfile(srcDir, sprintf('loopClosureResults_%s_all_shaped_xu_%s.mat', ds(j), v));
    assert(isfile(f), 'buildTable6a:missing', 'missing %s', f);
    S = load(f);  R = S.results;  N(j) = numel(R);
    if ver == "v015"
        assert(S.config.N_REPS == 200 && S.config.N_BETA == 25 && S.config.RNG_SEED == 1729, 'buildTable6a:config', '%s config not N_REPS=200/N_BETA=25/seed=1729', ds(j));
    end
    [tf, ix] = ismember(pl, string(S.pipelineLabels));  assert(all(tf), 'buildTable6a:labels', 'pipeline labels differ in %s', ds(j));
    BO = cell2mat(arrayfun(@(r) r.betaObs(:)', R(:), 'UniformOutput', false));       % trials x 6
    BG = cell2mat(arrayfun(@(r) r.betaGenStar(:)', R(:), 'UniformOutput', false));
    BO = BO(:, ix);  BG = BG(:, ix);
    for k = 1:numel(pl)
        x = BO(:, k);  x = x(isfinite(x));
        rows(end+1, :) = {ds(j), pl(k), numel(x), median(x), mean(x), std(x)};    %#ok<AGROW>
        nc(k, j) = 100 * mean(isnan(BG(:, k)));
    end
end
T6a = cell2table(rows, 'VariableNames', {'dataset','pipeline','nFinite','med','mean','sd'});
T6a.dataset = string(T6a.dataset);  T6a.pipeline = string(T6a.pipeline);
TNC = array2table(nc, 'VariableNames', cellstr(ds), 'RowNames', cellstr(pl));
info = struct('ver', ver, 'N', N, 'datasets', ds, 'pipelines', pl, 'runDate', datestr(now));
save(fullfile(resDir, sprintf('buildTable6aV015_v001_%s.mat', ver)), 'T6a', 'TNC', 'info');
writetable(T6a, fullfile(resDir, sprintf('buildTable6aV015_v001_%s_T6a.csv', ver)));
writetable(TNC, fullfile(resDir, sprintf('buildTable6aV015_v001_%s_TNC.csv', ver)), 'WriteRowNames', true);
fprintf('N = %s\n', mat2str(N));  disp(TNC);
end
