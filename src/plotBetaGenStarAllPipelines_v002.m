%% plotBetaGenStarAllPipelines_v002.m
%
% For each dataset: six-panel strip chart showing beta_obs and beta_gen*
% for all six pipelines, all trials, sorted by beta_obs.
% Mean +/- SD shown as horizontal bands.
%
% Reads shaped_xu v007 loop-closure results (post ceiling-fix, June 2026).
% Produces one figure per dataset, saved to ../figures/.
%
% v002: registry updated from _all_v003.mat (fftnoise) to
%       _all_shaped_xu_v007.mat (canonical post ceiling-fix).
%
% Requires: loopClosureResults_<DATASET>_all_shaped_xu_v007.mat in src/
%
% Fraser, D.S. (2026)  v002

cd(fileparts(mfilename('fullpath')));
addpath(genpath('functions')); addpath(genpath('req'));

PIPE_LABELS = ["BWFD-OLS","BWFD-LMLS","BWFD-IRLS","SG-OLS","SG-LMLS","SG-IRLS"];
PIPE_COLORS = [
    0.20 0.45 0.70;   % BWFD-OLS  — steel blue
    0.60 0.20 0.60;   % BWFD-LMLS — purple
    0.85 0.35 0.10;   % BWFD-IRLS — burnt orange
    0.10 0.60 0.40;   % SG-OLS    — teal green
    0.80 0.15 0.35;   % SG-LMLS   — crimson
    0.10 0.40 0.75;   % SG-IRLS   — blue (canonical)
];
nPipes = 6;

registry = {
    'Zarandi',      'loopClosureResults_Zarandi_all_shaped_xu_v007.mat',      [0.85 0.52 0.05];
    'Dhieb RETED',  'loopClosureResults_Dhieb_all_shaped_xu_v007.mat',        [0.55 0.18 0.65];
    'Cook CTRL',    'loopClosureResults_Cook_CTRL_all_shaped_xu_v007.mat',    [0.12 0.42 0.78];
    'Cook ASD',     'loopClosureResults_Cook_ASD_all_shaped_xu_v007.mat',     [0.85 0.18 0.18];
    'Hickman PLAC', 'loopClosureResults_Hickman_PLAC_all_shaped_xu_v007.mat', [0.08 0.58 0.35];
    'Hickman HALO', 'loopClosureResults_Hickman_HALO_all_shaped_xu_v007.mat', [0.75 0.35 0.10];
    'Pilot',        'loopClosureResults_Pilot_all_shaped_xu_v007.mat',        [0.10 0.65 0.60];
};

for di = 1:size(registry,1)
    label = registry{di,1};
    fname = registry{di,2};

    fprintf('Plotting %s...\n', label);
    S = load(fname, 'results');
    R = S.results;
    N = numel(R);

    % Visit mask for Pilot (V1/V2 coloured separately)
    isPilot    = strcmp(label, 'Pilot');
    trialIDs   = string({R.trialID}');
    v1mask_raw = contains(trialIDs, '_V1_');
    v2mask_raw = ~v1mask_raw;

    % Pre-extract all pipelines [N x 6]
    bO_all  = reshape([R.betaObs],     nPipes, N)';   % betaObs is [1x6]
    bS_all  = reshape([R.betaGenStar], nPipes, N)';   % betaGenStar is [6x1]
    lo_all  = reshape([R.ciLo],        nPipes, N)';
    hi_all  = reshape([R.ciHi],        nPipes, N)';

    fig = figure('Name', sprintf('beta_gen* all pipelines — %s', label), ...
        'Position', [30 30 420*3 360*2], 'Visible', 'off');
    tl  = tiledlayout(2, 3, 'TileSpacing','compact', 'Padding','compact');
    title(tl, {sprintf('\\beta_{gen}^* inversion — all pipelines, all trials: %s  (N=%d)', ...
        label, N), ...
        'Grey = \beta_{obs}.  Colour = \beta_{gen}^*.  Mean \pm SD bands.  Sorted by \beta_{obs}.'}, ...
        'FontSize', 10);

    for pp = 1:nPipes
        col  = PIPE_COLORS(pp,:);

COL_V1 = [0.08 0.38 0.84];   % V1 — strong blue (consistent across all panels)
COL_V2 = [0.88 0.18 0.12];   % V2 — strong red  (consistent across all panels)
        bO   = bO_all(:, pp);
        bS   = bS_all(:, pp);
        lo   = lo_all(:, pp);
        hi   = hi_all(:, pp);

        [bO_s, ord] = sort(bO, 'ascend');
        bS_s = bS(ord); lo_s = lo(ord); hi_s = hi(ord);
        x     = (1:N)';
        valid = isfinite(bS_s);

        mO  = mean(bO_s, 'omitnan');
        sdO = std(bO_s,  'omitnan');
        nC  = sum(valid);
        if nC > 1
            mS  = mean(bS_s(valid));
            sdS = std(bS_s(valid));
        else
            mS = NaN; sdS = NaN;
        end

        ax = nexttile; hold(ax, 'on');

        % CI band removed — inverted lo/hi causes irregular polygon at high N

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

        % Scatter — visit-split for Pilot, single colour otherwise
        plot(ax, x, bO_s, '.', 'Color', [0.68 0.68 0.68], ...
            'MarkerSize', 3, 'DisplayName', '\beta_{obs}');
        if any(valid)
            if isPilot
                v1s = v1mask_raw(ord) & valid;
                v2s = v2mask_raw(ord) & valid;
                scatter(ax, x(v1s), bS_s(v1s), 6, COL_V1, 'filled', ...
                    'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', 'none', ...
                    'DisplayName', 'V1 \beta_{gen}^*');
                COL_V2 = COL_V2;
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
            annStr = sprintf('obs  = %.3f \x00B1 %.3f\ngen* = %.3f \x00B1 %.3f\nconv=%d/%d (%.0f%%)', ...
                mO, sdO, mS, sdS, nC, N, 100*nC/N);
        else
            annStr = sprintf('obs  = %.3f \x00B1 %.3f\ngen* = N/A  conv=%d/%d', ...
                mO, sdO, nC, N);
        end
        text(ax, N*0.97, 0.77, annStr, ...
            'FontSize', 7.5, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'BackgroundColor', 'w', 'EdgeColor', [0.80 0.80 0.80], 'FontName', 'Courier');

        title(ax, PIPE_LABELS(pp), 'FontSize', 9.5, 'Color', col * 0.7);
        xlabel(ax, 'Trial (sorted by \beta_{obs})', 'FontSize', 8);
        ylabel(ax, '\beta', 'FontSize', 9);
        xlim(ax, [0.5 N+0.5]); ylim(ax, [-0.05 0.82]);
        legend(ax, 'Location', 'northwest', 'FontSize', 7.5);
        grid(ax, 'on'); box(ax, 'on'); hold(ax, 'off');
    end

    tag     = strrep(strrep(label,' ','_'), '/','_');
    outFile = sprintf('../figures/betaGenStar_allPipelines_%s_allTrials.png', tag);
    print(fig, outFile, '-dpng', '-r150');
    close(fig);
    fprintf('  Saved: %s\n', outFile);
end

fprintf('\nAll done.\n');
