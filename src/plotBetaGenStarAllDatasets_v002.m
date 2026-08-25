%% plotBetaGenStarAllDatasets_v002.m
%
% Seven-panel strip chart: beta_obs and beta_gen* (SG-IRLS) for all datasets,
% all trials, sorted by beta_obs. Mean +/- SD shown as horizontal bands.
%
% Reads shaped_xu v007 loop-closure results (post ceiling-fix, June 2026).
% Requires: loopClosureResults_<DATASET>_all_shaped_xu_v007.mat in src/
%
% v002: registry updated from _all_v003.mat (fftnoise) to
%       _all_shaped_xu_v007.mat (canonical post ceiling-fix).
%
% Fraser, D.S. (2026)  v002

cd(fileparts(mfilename('fullpath')));
addpath(genpath('functions')); addpath(genpath('req'));

PIPE = 6; % SG-IRLS

COL_V1 = [0.08 0.38 0.84];   % V1 — strong blue (consistent across all panels)
COL_V2 = [0.88 0.18 0.12];   % V2 — strong red  (consistent across all panels)

registry = {
    'Zarandi',       'loopClosureResults_Zarandi_all_shaped_xu_v007.mat',       [0.85 0.52 0.05], 'o';
    'Dhieb RETED',   'loopClosureResults_Dhieb_all_shaped_xu_v007.mat',         [0.55 0.18 0.65], 'd';
    'Cook CTRL',     'loopClosureResults_Cook_CTRL_all_shaped_xu_v007.mat',     [0.12 0.42 0.78], 's';
    'Cook ASD',      'loopClosureResults_Cook_ASD_all_shaped_xu_v007.mat',      [0.85 0.18 0.18], 's';
    'Hickman PLAC',  'loopClosureResults_Hickman_PLAC_all_shaped_xu_v007.mat',  [0.08 0.58 0.35], 'p';
    'Hickman HALO',  'loopClosureResults_Hickman_HALO_all_shaped_xu_v007.mat',  [0.75 0.35 0.10], 'p';
    'Pilot',         'loopClosureResults_Pilot_all_shaped_xu_v007.mat',         [0.10 0.65 0.60], '^';
};

nDS  = size(registry, 1);
nCol = 4;
nRow = ceil(nDS / nCol);

fig = figure('Name','beta_gen* SG-IRLS — all datasets, all trials', ...
    'Position', [30 30 420*nCol 360*nRow]);
tl  = tiledlayout(nRow, nCol, 'TileSpacing','compact', 'Padding','compact');
title(tl, {'\beta_{gen}^* inversion — SG-IRLS, all trials (mean \pm SD shown as horizontal bands)', ...
    'Grey dots = \beta_{obs}.  Colour dots = \beta_{gen}^*.  Sorted by \beta_{obs}.'}, ...
    'FontSize', 11);

for di = 1:nDS
    label = registry{di,1};
    fname = registry{di,2};
    col   = registry{di,3};

    S  = load(fname, 'results');
    R  = S.results;
    N  = numel(R);
    isPilot    = strcmp(label, 'Pilot');
    trialIDs   = string({R.trialID}');
    v1mask_raw = contains(trialIDs, '_V1_');
    v2mask_raw = ~v1mask_raw;
    bO = arrayfun(@(r) r.betaObs(PIPE),     R)';
    bS = arrayfun(@(r) r.betaGenStar(PIPE), R)';
    lo = arrayfun(@(r) r.ciLo(PIPE),        R)';
    hi = arrayfun(@(r) r.ciHi(PIPE),        R)';

    [bO_s, ord] = sort(bO, 'ascend');
    bS_s = bS(ord); lo_s = lo(ord); hi_s = hi(ord);
    x     = (1:N)';
    valid = isfinite(bS_s);

    % Statistics — guard empty valid
    mO  = mean(bO, 'omitnan');
    sdO = std(bO,  'omitnan');
    nC  = sum(valid);
    if nC > 1
        mS  = mean(bS_s(valid));
        sdS = std(bS_s(valid));
    else
        mS = NaN; sdS = NaN;
    end

    ax = nexttile; hold(ax, 'on');

    % 5-95% CI band removed — inverted lo/hi causes irregular polygon at high N

    % Mean +/- SD bands
    fill(ax, [0.5; N+0.5; N+0.5; 0.5], ...
        [mO-sdO; mO-sdO; mO+sdO; mO+sdO], ...
        [0.65 0.65 0.65], 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
    if isfinite(mS)
        fill(ax, [0.5; N+0.5; N+0.5; 0.5], ...
            [mS-sdS; mS-sdS; mS+sdS; mS+sdS], ...
            col, 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end

    % Scatter
    plot(ax, x, bO_s, '.', 'Color', [0.68 0.68 0.68], ...
        'MarkerSize', 3, 'DisplayName', '\beta_{obs}');
    if any(valid)
        if isPilot
            v1s = v1mask_raw(ord) & valid;
            v2s = v2mask_raw(ord) & valid;
            scatter(ax, x(v1s), bS_s(v1s), 6, COL_V1, 'filled', ...
                'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', 'none', ...
                'DisplayName', 'V1 \beta_{gen}^*');
            scatter(ax, x(v2s), bS_s(v2s), 6, COL_V2, 'filled', ...
                'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', 'none', ...
                'DisplayName', 'V2 \beta_{gen}^*');
        else
            scatter(ax, x(valid), bS_s(valid), 6, col, 'filled', ...
                'MarkerFaceAlpha', 0.50, 'MarkerEdgeColor', 'none', ...
                'DisplayName', '\beta_{gen}^*');
        end
    end

    % Mean lines
    yline(ax, mO, '-', 'Color', [0.50 0.50 0.50], 'LineWidth', 2.0, ...
        'HandleVisibility', 'off');
    if isfinite(mS)
        yline(ax, mS, '-', 'Color', col, 'LineWidth', 2.5, ...
            'HandleVisibility', 'off');
    end

    % Reference lines
    yline(ax, 1/3, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.9, ...
        'Label', '1/3', 'LabelHorizontalAlignment', 'left', ...
        'HandleVisibility', 'off');
    yline(ax, 0.5, ':',  'Color', [0.55 0.55 0.55], 'LineWidth', 0.9, ...
        'Label', '1/2', 'LabelHorizontalAlignment', 'left', ...
        'HandleVisibility', 'off');

    % Annotation
    if isfinite(mS)
        annStr = sprintf('obs  = %.3f \x00B1 %.3f\ngen* = %.3f \x00B1 %.3f\nconv = %d/%d (%.0f%%)', ...
            mO, sdO, mS, sdS, nC, N, 100*nC/N);
    else
        annStr = sprintf('obs  = %.3f \x00B1 %.3f\ngen* = N/A (conv=%d/%d)', ...
            mO, sdO, nC, N);
    end
    text(ax, N*0.97, 0.77, annStr, ...
        'FontSize', 8, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'EdgeColor', [0.80 0.80 0.80], 'FontName', 'Courier');

    title(ax, label, 'FontSize', 10);
    xlabel(ax, 'Trial (sorted by \beta_{obs})', 'FontSize', 8.5);
    ylabel(ax, '\beta', 'FontSize', 10);
    xlim(ax, [0.5 N+0.5]); ylim(ax, [-0.05 0.82]);
    legend(ax, 'Location', 'northwest', 'FontSize', 8);
    grid(ax, 'on'); box(ax, 'on'); hold(ax, 'off');
end

outFile = '../figures/betaGenStar_allDatasets_SGIRLS_allTrials.png';
print(fig, outFile, '-dpng', '-r150');
fprintf('Saved: %s\n', outFile);
