function plotSecondaryMetricsV015_v001(opts)
%PLOTSECONDARYMETRICSV015_V001  Figs 10 (Bland-Altman), 11 (pattern preservation) and 13 (CP / TDI), Supplement A, rebuilt on the v015 basis with the trial-agnostic comparator.
%   Replaces the v007/v008 versions (analyzeBlandAltman_v001, analyzePatternPreservation_v002, analyzeTDICoverage_v001). Reads loopClosureResults_<ds>_all_shaped_xu_v015.mat directly and
%   builds the validity mask, contrasts, predictions and the leave-one-out comparator exactly as constellationMetrics_v005; fails loud if N_REPS / N_BETA / seed differ from 200 / 25 / 1729, and if the
%   headline medians drawn here differ from results/constellationMetrics_v005.mat (cccMed not drawn; rho median, CP(MDC) and TDI90 are). Writes ../figures/fig10_blandAltman_v015.png, fig11_patternPreservation_v015.png,
%   fig13_cpTdi_v015.png. Fraser, D.S. (2026)  v001
arguments
    opts.Dir (1,1) string = ""
    opts.Datasets (1,:) string = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Dhieb" "Zarandi"]
    opts.MaxPoints (1,1) double = 6000
end
srcDir = fileparts(mfilename('fullpath'));  dataDir = opts.Dir;  if dataDir == "", dataDir = string(srcDir); end
figDir = fullfile(srcDir, '..', 'figures');  if ~isfolder(figDir), mkdir(figDir); end
V5 = load(fullfile(srcDir, '..', 'results', 'constellationMetrics_v005.mat'), 'T').T;
MDC = 0.03;  DEG = 0.005;  nPipe = 6;  rng(1729);  [pP, pQ] = pairs_local(nPipe);  nD = numel(opts.Datasets);
D = struct('ds', {}, 'O', {}, 'Pm', {}, 'Nu', {}, 'rho', {}, 'nRho', {});
for ds = opts.Datasets
    f = fullfile(dataDir, "loopClosureResults_" + ds + "_all_shaped_xu_v015.mat");  if ~isfile(f), error('plotSec:missing', 'No file: %s', f); end
    S0 = load(f, 'results', 'config');  c = S0.config;  r = S0.results;  N = numel(r);
    if c.N_REPS ~= 200 || c.N_BETA ~= 25 || c.RNG_SEED ~= 1729, error('plotSec:config', '%s: N_REPS/N_BETA/seed differ from 200/25/1729', f); end
    bO = NaN(N, nPipe);  bP = bO;  dO = NaN(N, numel(pP));  dP = dO;  ok = false(N, 1);
    for t = 1:N
        bo = r(t).betaObs(:)';  sl = r(t).betaRecSlice;
        if numel(bo) ~= nPipe || any(~isfinite(bo)) || isempty(sl) || size(sl, 2) ~= nPipe || ~isfinite(r(t).betaGenStarMed), continue; end
        bp = median(sl, 1, 'omitnan');  if any(~isfinite(bp)), continue; end
        bO(t, :) = bo;  bP(t, :) = bp;  dO(t, :) = bo(pP) - bo(pQ);  dP(t, :) = bp(pP) - bp(pQ);  ok(t) = true;
    end
    vi = find(ok & max(abs(dO), [], 2) >= DEG & ~all(abs(dO) < MDC, 2));  n = numel(vi);  O = dO(vi, :);  Pm = dP(vi, :);  BO = bO(vi, :);  BP = bP(vi, :);
    rho = NaN(n, 1);  nRho = rho;  Nu = NaN(n, size(O, 2));  sumO = sum(O, 1);  sumB = sum(BO, 1);
    for i = 1:n
        rho(i) = corr(BP(i, :)', BO(i, :)', 'Type', 'Spearman');
        if n >= 2, Nu(i, :) = (sumO - O(i, :)) / (n - 1);  nRho(i) = corr(((sumB - BO(i, :)) / (n - 1))', BO(i, :)', 'Type', 'Spearman'); end
    end
    D(end+1) = struct('ds', ds, 'O', O, 'Pm', Pm, 'Nu', Nu, 'rho', rho, 'nRho', nRho); %#ok<AGROW>
    v = V5(V5.dataset == ds & V5.noiseModel == "shaped_xu", :);  d = Pm(:) - O(:);
    if abs(median(rho, 'omitnan') - v.rhoMed) > 1e-9 || abs(100 * mean(abs(d) < MDC) - v.cpMDC) > 1e-6 || abs(prctile(abs(d), 90) - v.tdiEmp) > 1e-9
        error('plotSec:mismatch', '%s: rho median / CP / TDI drawn here differ from constellationMetrics_v005', ds);
    end
end
% proportional bias of the Bland-Altman differences: OLS slope of (pred - obs) on mean(pred, obs), and its Spearman correlation, over the pooled pairs (slope 2(a-1)/(a+1) for a prediction a times the observation; -2 is a constant prediction)
BA = table(strings(nD, 1), zeros(nD, 1), zeros(nD, 1), zeros(nD, 1), 'VariableNames', ["dataset" "nPairs" "slope" "spearman"]);
for k = 1:nD
    d = D(k).Pm(:) - D(k).O(:);  m = (D(k).Pm(:) + D(k).O(:)) / 2;  pf = polyfit(m, d, 1);
    BA(k, :) = {D(k).ds, numel(d), pf(1), corr(m, d, 'Type', 'Spearman')};
end
disp(BA);  save(fullfile(srcDir, '..', 'results', 'plotSecondaryMetricsV015_v001.mat'), 'BA');
lab = @(s) strrep(s, "_", " ");  cm = lines(2);
%% Fig 10  Bland-Altman
fh = figure('Visible', 'off', 'Position', [100 100 1500 950]);  tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:nD
    nexttile;  d = D(k).Pm(:) - D(k).O(:);  m = (D(k).Pm(:) + D(k).O(:)) / 2;  j = randperm(numel(d), min(numel(d), opts.MaxPoints));
    scatter(m(j), d(j), 4, cm(1, :), 'filled', 'MarkerFaceAlpha', 0.25);  hold on;  b = mean(d);  sd = std(d);
    yline(0, 'k-');  yline(b, 'b-');  yline([b - 1.96 * sd, b + 1.96 * sd], 'b--');  yline([-MDC MDC], 'r:');
    title(sprintf('%s (n = %d pairs)  bias %+.3f  CP(MDC) %.1f%%', lab(D(k).ds), numel(d), b, 100 * mean(abs(d) < MDC)), 'FontSize', 8);  xlabel('mean of predicted and observed \Delta\beta');  ylabel('predicted - observed');
end
sgtitle('Bland-Altman, pooled trial x contrast pairs (v015, shaped\_xu). Blue: bias and 95% limits; red dotted: \pm MDC');
exportgraphics(fh, fullfile(figDir, 'fig10_blandAltman_v015.png'), 'Resolution', 200);  close(fh);
%% Fig 11  pattern preservation
fh = figure('Visible', 'off', 'Position', [100 100 1500 950]);  tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:nD
    nexttile;  histogram(D(k).rho, -1:0.1:1, 'FaceColor', cm(1, :));  hold on;  xline(median(D(k).rho, 'omitnan'), 'b-', 'LineWidth', 1.5);
    if any(isfinite(D(k).nRho)), xline(median(D(k).nRho, 'omitnan'), 'r--', 'LineWidth', 1.5); end
    title(sprintf('%s (n = %d)  median \\rho %.3f, comparator %.3f', lab(D(k).ds), numel(D(k).rho), median(D(k).rho, 'omitnan'), median(D(k).nRho, 'omitnan')), 'FontSize', 8);  xlabel('Spearman \rho (predicted vs observed pipeline ranking)');  ylabel('trials');
end
sgtitle('Pattern preservation (v015, shaped\_xu). Blue: model median; red dashed: trial-agnostic comparator median');
exportgraphics(fh, fullfile(figDir, 'fig11_patternPreservation_v015.png'), 'Resolution', 200);  close(fh);
%% Fig 13  CP and TDI with comparator
cp = NaN(nD, 2);  td = cp;
for k = 1:nD
    dm = D(k).Pm(:) - D(k).O(:);  cp(k, 1) = 100 * mean(abs(dm) < MDC);  td(k, 1) = prctile(abs(dm), 90);
    if any(isfinite(D(k).Nu(:))), dn = D(k).Nu(:) - D(k).O(:);  cp(k, 2) = 100 * mean(abs(dn) < MDC);  td(k, 2) = prctile(abs(dn), 90); end
end
fh = figure('Visible', 'off', 'Position', [100 100 1400 520]);  tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;  bar(categorical(lab(opts.Datasets), lab(opts.Datasets)), cp);  ylabel('CP(MDC), %');  legend(["model" "comparator"], 'Location', 'southwest');  title('Coverage probability within MDC');
nexttile;  bar(categorical(lab(opts.Datasets), lab(opts.Datasets)), td);  yline(MDC, 'k--', 'MDC');  ylabel('TDI_{90} (empirical)');  title('Total deviation index, 90%');
sgtitle('CP and TDI (v015, shaped\_xu); Zarandi has one valid trial, so no comparator');
exportgraphics(fh, fullfile(figDir, 'fig13_cpTdi_v015.png'), 'Resolution', 200);  close(fh);
fprintf('Saved figures to %s\n', figDir);
end

function [pP, pQ] = pairs_local(n)
    pP = [];  pQ = [];
    for p = 2:n, for q = 1:p-1, pP(end+1) = p; pQ(end+1) = q; end, end %#ok<AGROW>
end
