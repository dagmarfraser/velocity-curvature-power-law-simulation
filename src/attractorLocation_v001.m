% attractorLocation_v001.m  —  Where does each pipeline's beta_rec actually converge?
%
% attractorBasin_v001 measured the width of the "basin" around beta = 1/3.
% That question assumes 1/3 IS the universal attractor. This script makes no
% such assumption. It asks: at maximum noise, does the beta_gen -> beta_rec
% mapping go flat (attractor regime), and if so, at what value?
%
% Three diagnostics:
%   Fig 1 — Recovery curves at 4 sigma levels (one figure per alpha).
%            2x3 tiled layout, one tile per pipeline. Coloured by sigma.
%            Solid horizontal line marks attractor at sigma = 10 mm.
%   Fig 2 — Attractor location vs alpha for all 6 pipelines at sigma = 10 mm.
%            Key question: do all pipelines share ~1/3, or diverge?
%   Fig 3 — Recovery slope (d beta_rec / d beta_gen, polyfit) vs sigma.
%            Slope = 1 is perfect recovery; slope = 0 is flat attractor.
%
% Console: attractor location table (pipeline x alpha) at sigma = 10 mm,
%          plus focused summary at alpha = 2 (red noise, empirically typical).
%
% CORRECTION 2026-08-26: was reading perCoordinateSEM_v001.mat, an older,
% smaller-grid artefact (4 alpha values 0-3, sigma capped at 10mm). The
% current v058 grid (perCoordinateSEM_v2_001.mat) has 31 alpha values (0-6)
% and sigma to 20mm -- the 10mm ceiling this script's own SIGMA_MAX assumed
% was correct for the OLD file, but v058 was specifically expanded because
% real datasets (Cook/Hickman) have empirical sigma ~14mm, above that old
% ceiling. This script was therefore blind to the exact noise regime that
% matters most for real data. VGF axis confirmed identical between the two
% files before swapping (14 values, max abs diff 0), so REF_VGF_I=8 still
% resolves to the same VGF=181.3.
%
% Requires: perCoordinateSEM_v2_001.mat (current v058-era grid)
% Fraser, D.S. (2026)

%% Configuration
REF_FS     = 120;    % Hz  — fs effect negligible; consistent with other scripts
REF_VGF_I  = 8;     % index into sorted VGF vector — exp(5.2) ≈ 181.3 (~1 Hz)
SIGMA_MAX  = 20.0;  % mm  — v058 grid maximum (extended from the pre-registered 10mm ceiling); used for attractor location estimate
ALPHA_SHOW = [0, 1, 2, 3];           % cardinal noise colours for Fig 1 (one fig each)
SIGMA_SHOW = [0.5, 2.0, 6.0, 15.0, 20.0]; % mm — noise levels overlaid per Fig 1 panel; 15mm is the nearest actual grid point to Cook/Hickman's own empirical sigma (~14mm; grid has 12 and 15, not 14), 20mm is the new grid max

%% Pipeline aesthetics (consistent across project)
pipeNames  = ["BWFD-OLS","BWFD-LMLS","BWFD-IRLS","SG-OLS","SG-LMLS","SG-IRLS"];
pipeColors = [0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.56;
              0.47 0.67 0.19; 0.30 0.75 0.93; 0.00 0.45 0.74];
nPipes     = numel(pipeNames);

%% Load
fprintf('Loading perCoordinateSEM_v2_001.mat ...\n');
load('perCoordinateSEM_v2_001.mat', 'coordTable');
T = coordTable;
fprintf('  %d coordinates loaded\n', height(T));

allVGF = sort(unique(T.VGF));
if REF_VGF_I > numel(allVGF)
    error('REF_VGF_I=%d exceeds available VGF levels (%d). Adjust.', ...
          REF_VGF_I, numel(allVGF));
end
refVGF = allVGF(REF_VGF_I);
fprintf('  Reference: VGF = %.1f (index %d), fs = %d Hz\n\n', ...
        refVGF, REF_VGF_I, REF_FS);

T = T(T.VGF == refVGF & T.fs == REF_FS, :);
fprintf('  After VGF/fs filter: %d coordinates\n\n', height(T));
if isempty(T)
    error('No data after VGF/fs filter. Check REF_FS and REF_VGF_I.');
end

allSigma = sort(unique(T.sigma));
allAlpha = sort(unique(T.alpha));
allBeta  = sort(unique(T.betaGen));

% Validate requested sigma/alpha are in grid
for s = SIGMA_SHOW
    if ~ismember(s, allSigma)
        error('SIGMA_SHOW value %.4f not in sigma grid.', s);
    end
end
for a = ALPHA_SHOW
    if ~ismember(a, allAlpha)
        error('ALPHA_SHOW value %.1f not in alpha grid.', a);
    end
end

%% Ensure figures directory exists
figDir = fullfile('..', 'figures');
if ~isfolder(figDir), mkdir(figDir); end

%% ================================================================
%  FIGURE 1: Recovery curves — 4 sigma levels, one figure per alpha
%  2x3 tile layout, one pipeline per tile.
%  ================================================================
alphaLabels = ["\alpha=0 (white)", "\alpha=1 (pink)", ...
               "\alpha=2 (red)", "\alpha=3 (black)"];

% Sigma colour ramp: light green (low noise) -> black (max noise)
sigmaColors = [0.40 0.75 0.40;   % sigma = 0.5 mm
               0.10 0.55 0.15;   % sigma = 2 mm
               0.00 0.25 0.05;   % sigma = 6 mm
               0.35 0.15 0.05;   % sigma = 15 mm (nearest grid point to Cook/Hickman empirical ~14mm)
               0.10 0.10 0.10];  % sigma = 20 mm (new grid max)

for ai = 1:numel(ALPHA_SHOW)
    alph = ALPHA_SHOW(ai);
    fig1 = figure('Name', sprintf('Fig1 alpha=%.0f', alph), ...
                  'Position', [50 50 1400 480]);
    tl = tiledlayout(2, 3, 'TileSpacing','compact', 'Padding','compact');

    for pi = 1:nPipes
        ax = nexttile; hold(ax, 'on');

        % Identity and 1/3 reference
        plot(ax, [0 2/3], [0 2/3], 'k--', 'LineWidth', 0.8, ...
             'HandleVisibility','off');
        yline(ax, 1/3, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7, ...
              'HandleVisibility','off');

        for si = 1:numel(SIGMA_SHOW)
            sig = SIGMA_SHOW(si);
            sel = sliceSelect(T, pipeNames(pi), sig, alph);
            if isempty(sel), continue; end
            plot(ax, sel.betaGen, sel.meanBetaRec, '-o', ...
                 'Color', sigmaColors(si,:), 'LineWidth', 1.4, ...
                 'MarkerSize', 3, 'MarkerFaceColor', sigmaColors(si,:), ...
                 'DisplayName', sprintf('\\sigma=%.4g mm', sig));
        end

        % Attractor location at max sigma — annotated horizontal line
        selMax = sliceSelect(T, pipeNames(pi), SIGMA_MAX, alph);
        if ~isempty(selMax)
            loc = mean(selMax.meanBetaRec);
            yline(ax, loc, '-', 'Color', pipeColors(pi,:), 'LineWidth', 2.0, ...
                  'HandleVisibility','off');
            text(ax, 0.01, loc + 0.018, sprintf('%.3f', loc), ...
                 'Color', pipeColors(pi,:), 'FontSize', 7.5, 'FontWeight','bold');
        end

        xlim(ax, [0 2/3]);
        ylim(ax, [-0.05 0.60]);
        xlabel(ax, '\beta_{gen}', 'FontSize', 8);
        if mod(pi-1,3) == 0, ylabel(ax, '\beta_{rec}', 'FontSize', 8); end
        title(ax, pipeNames(pi), 'FontSize', 9, ...
              'Color', pipeColors(pi,:), 'FontWeight','bold');
        grid(ax, 'on'); box(ax, 'on');
        if pi == nPipes
            lgd = legend(ax, 'Location','southeast', 'FontSize', 7);
            lgd.ItemTokenSize = [12 8];
        end
        hold(ax, 'off');
    end

    title(tl, sprintf('Recovery curves — %s (fs=%d Hz, VGF=%.0f)', ...
          alphaLabels(ai), REF_FS, refVGF), 'FontSize', 11);
    saveas(fig1, fullfile(figDir, ...
           sprintf('attractorLocation_fig1_alpha%d.png', round(alph))));
end

%% ================================================================
%  FIGURE 2: Attractor location vs alpha — all 6 pipelines, sigma = 10 mm
%  The key diagnostic: are all attractors at 1/3, or does location vary?
%  ================================================================
fig2 = figure('Name', 'Fig2: Attractor location vs alpha', ...
              'Position', [100 100 820 520]);
hold on;

yline(1/3,  'k:',  'LineWidth', 1.2, 'DisplayName', '\beta = 1/3 canonical');
yline(0.29, '--', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.9, ...
      'DisplayName', 'Prior note: BWFD-OLS ~0.29');
yline(0.21, '-.', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.9, ...
      'DisplayName', 'Prior note: SG-IRLS ~0.21');

lineStyles = {'-o','-s','-^','--o','--s','--^'};

for pi = 1:nPipes
    locs = NaN(numel(allAlpha), 1);
    for ai = 1:numel(allAlpha)
        sel = sliceSelect(T, pipeNames(pi), SIGMA_MAX, allAlpha(ai));
        if ~isempty(sel), locs(ai) = mean(sel.meanBetaRec); end
    end
    plot(allAlpha, locs, lineStyles{pi}, ...
         'Color', pipeColors(pi,:), 'LineWidth', 1.6, ...
         'MarkerSize', 4, 'MarkerFaceColor', pipeColors(pi,:), ...
         'DisplayName', pipeNames(pi));
end

xlabel('\alpha (noise colour)', 'FontSize', 11);
ylabel(sprintf('Attractor: mean(\\beta_{rec}) at \\sigma = %.0f mm', SIGMA_MAX), ...
       'FontSize', 11);
title(sprintf('Pipeline-specific attractor locations  (fs=%d Hz, VGF=%.0f, \\sigma=%.0f mm)', ...
    REF_FS, refVGF, SIGMA_MAX), 'FontSize', 11);
legend('Location','northeast', 'FontSize', 8);
ylim([-0.05 0.50]); xlim([0 3]);
grid on; box on; hold off;
saveas(fig2, fullfile(figDir, 'attractorLocation_fig2_vs_alpha.png'));

%% ================================================================
%  FIGURE 3: Recovery slope vs sigma — when does the curve go flat?
%  Slope = 1 is perfect recovery; slope = 0 is fully flat (attractor).
%  ================================================================
alphaForSlope = [0, 2];  % white and red noise side by side
slopeTileLabels = ["\alpha=0 (white noise)", "\alpha=2 (red noise)"];

fig3 = figure('Name', 'Fig3: Recovery slope vs sigma', ...
              'Position', [150 150 960 480]);
tl3 = tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','compact');

for ai = 1:numel(alphaForSlope)
    alph = alphaForSlope(ai);
    ax = nexttile; hold(ax, 'on');

    yline(ax, 1, 'k--', 'LineWidth', 0.9, 'HandleVisibility','off');  % perfect recovery
    yline(ax, 0, 'k:',  'LineWidth', 0.9, 'HandleVisibility','off');  % flat attractor
    text(ax, 0.15, 1.08, 'perfect recovery', 'FontSize', 7, 'Color', [0.4 0.4 0.4]);
    text(ax, 0.15, 0.04, 'flat attractor',   'FontSize', 7, 'Color', [0.4 0.4 0.4]);

    for pi = 1:nPipes
        slopes = NaN(numel(allSigma), 1);
        for si = 1:numel(allSigma)
            sel = sliceSelect(T, pipeNames(pi), allSigma(si), alph);
            if height(sel) >= 2
                p = polyfit(sel.betaGen, sel.meanBetaRec, 1);
                slopes(si) = p(1);
            end
        end
        plot(ax, allSigma, slopes, lineStyles{pi}, ...
             'Color', pipeColors(pi,:), 'LineWidth', 1.4, ...
             'MarkerSize', 4, 'MarkerFaceColor', pipeColors(pi,:), ...
             'DisplayName', pipeNames(pi));
    end

    % Mark empirical sigma ranges
    xline(ax, 4.8, '--', 'Color', [0.1 0.5 0.1], 'LineWidth', 0.8, ...
          'HandleVisibility','off', 'Label','Zarandi', ...
          'LabelOrientation','horizontal', 'FontSize', 7);
    xline(ax, 7.5, '--', 'Color', [0.6 0.1 0.1], 'LineWidth', 0.8, ...
          'HandleVisibility','off', 'Label','Cook/Hickman', ...
          'LabelOrientation','horizontal', 'FontSize', 7);

    xlabel(ax, '\sigma (mm)', 'FontSize', 10);
    ylabel(ax, 'd\beta_{rec}/d\beta_{gen}  (polyfit slope)', 'FontSize', 10);
    title(ax, slopeTileLabels(ai), 'FontSize', 10);
    xlim(ax, [0 SIGMA_MAX]); ylim(ax, [-0.3 1.3]);
    grid(ax, 'on'); box(ax, 'on');
    if ai == 2
        legend(ax, 'Location','northeast', 'FontSize', 8);
    end
    hold(ax, 'off');
end
title(tl3, 'Recovery slope vs noise magnitude  (slope=1: perfect; slope=0: attractor)', ...
      'FontSize', 11);
saveas(fig3, fullfile(figDir, 'attractorLocation_fig3_slope.png'));

%% ================================================================
%  CONSOLE: Attractor location table
%  ================================================================
fprintf('\n=== ATTRACTOR LOCATIONS: mean(beta_rec) at sigma=%.0f mm ===\n', SIGMA_MAX);
fprintf('    fs=%d Hz, VGF=%.1f\n', REF_FS, refVGF);
fprintf('    Attractor = mean(beta_rec) across all 21 beta_gen values.\n');
fprintf('    Low spread = flat curve = genuine attractor.\n\n');

% Header row
fprintf('%-14s', 'Pipeline');
% CORRECTION 2026-08-26: alphaColHeaders was a hardcoded 0:0.5:3, a leftover
% from the pre-v058 grid (alpha step 0.1). v058's own step is 0.2
% (defineParameterSpace.m, debug==5), so half of 0:0.5:3's values (0.5, 1.5,
% 2.5) no longer land on an actual grid point -- this crashed the now-removed
% colIdx computation (dead code, never used downstream) and would have
% printed silent NaN columns in the table below even without the crash.
% Derived directly from allAlpha (already loaded from T) rather than a
% hardcoded literal, so it can never float-point-mismatch against T.alpha's
% own stored values -- same reasoning plotBetaRecovery_v005.m's own
% snap-to-grid section already applies elsewhere in this project.
alphaLE3 = allAlpha(allAlpha <= 3 + 1e-9);
alphaColHeaders = alphaLE3(1:2:end);  % every other grid point, 0-3, for readability
for a = alphaColHeaders, fprintf(' a=%3.1f', a); end
fprintf('\n%s\n', repmat('-', 1, 14 + 7*numel(alphaColHeaders)));

for pi = 1:nPipes
    fprintf('%-14s', pipeNames(pi));
    for ci = 1:numel(alphaColHeaders)
        sel = sliceSelect(T, pipeNames(pi), SIGMA_MAX, alphaColHeaders(ci));
        if isempty(sel)
            fprintf(' %5s', '  NaN');
        else
            fprintf(' %5.3f', mean(sel.meanBetaRec));
        end
    end
    fprintf('\n');
end

% Focused summary at alpha = 2
fprintf('\n--- Summary at alpha=2 (red noise, most empirically representative) ---\n');
fprintf('%-14s  %8s  %8s  %8s  %8s\n', ...
        'Pipeline', 'AttrLoc', 'Spread', 'Slope', '|Loc-1/3|');
for pi = 1:nPipes
    sel = sliceSelect(T, pipeNames(pi), SIGMA_MAX, 2.0);
    if isempty(sel)
        fprintf('%-14s  (no data)\n', pipeNames(pi));
        continue
    end
    loc    = mean(sel.meanBetaRec);
    spread = max(sel.meanBetaRec) - min(sel.meanBetaRec);
    p      = polyfit(sel.betaGen, sel.meanBetaRec, 1);
    slope  = p(1);
    fprintf('%-14s  %8.4f  %8.4f  %8.4f  %8.4f\n', ...
            pipeNames(pi), loc, spread, slope, abs(loc - 1/3));
end

fprintf('\nNOTE: "spread" < 0.05 with slope near 0 indicates genuine attractor regime.\n');
fprintf('      Compare AttrLoc across pipelines: if they differ, there is no\n');
fprintf('      single shared attractor — each pipeline has its own compression floor.\n');
fprintf('\n=== DONE ===\n');

%% ================================================================
%  LOCAL FUNCTION: subset rows for a (pipeline, sigma, alpha) cell
%  ================================================================
function sel = sliceSelect(T, pipeName, sig, alph)
    sel = sortrows( ...
            T(T.pipeline == pipeName & T.sigma == sig & T.alpha == alph, :), ...
            'betaGen');
end
