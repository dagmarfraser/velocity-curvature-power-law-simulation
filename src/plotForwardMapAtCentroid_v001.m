function plotForwardMapAtCentroid_v001(tgtAlpha, tgtSigma, tgtFs, datasetLabel)
% PLOTFORWARDMAPCENTROID_V001  Forward map beta_gen->beta_rec at any (alpha,sigma,fs).
%
% Shows the bend-back non-monotonicity that causes invertibility failure.
% Reads aggregated meanBetaRec from perCoordinateSEM_v2_001.mat (no DB needed).
%
% USAGE:
%   plotForwardMapAtCentroid_v001()                         % Hickman PLAC defaults
%   plotForwardMapAtCentroid_v001(4.6, 5.677, 120, 'Hickman PLAC')
%   plotForwardMapAtCentroid_v001(3.649, 6.47, 120, 'Cook CTRL')
%   plotForwardMapAtCentroid_v001(2.531, 4.44,  120, 'Zarandi')
%
% Fraser, D.S. (2026)

    if nargin < 1; tgtAlpha      = 4.600; end
    if nargin < 2; tgtSigma      = 5.677; end
    if nargin < 3; tgtFs         = 120;   end
    if nargin < 4; datasetLabel  = 'Hickman PLAC'; end

    srcDir = fileparts(mfilename('fullpath'));
    cd(srcDir);

    %% Load
    fprintf('Loading perCoordinateSEM_v2_001.mat...\n');
    load('perCoordinateSEM_v2_001.mat', 'coordTable');
    T = coordTable;

    %% Snap to nearest grid node
    allAlpha = sort(unique(T.alpha));
    allSigma = sort(unique(T.sigma));
    allFs    = sort(unique(T.fs));
    [~,ai] = min(abs(allAlpha - tgtAlpha));
    [~,si] = min(abs(allSigma - tgtSigma));
    [~,fi] = min(abs(allFs    - tgtFs));
    snapAlpha = allAlpha(ai);
    snapSigma = allSigma(si);
    snapFs    = allFs(fi);

    fprintf('Target:  alpha=%.3f  sigma=%.3f mm  fs=%d Hz\n', tgtAlpha, tgtSigma, tgtFs);
    fprintf('Snapped: alpha=%.3f  sigma=%.3f mm  fs=%d Hz\n', snapAlpha, snapSigma, snapFs);

    %% Subset
    sub = T(T.alpha == snapAlpha & T.sigma == snapSigma & T.fs == snapFs, :);
    if isempty(sub)
        error('plotForwardMap:NoData', 'No data at snapped coordinate.');
    end

    %% Pipeline aesthetics
    PIPELINES  = ["BWFD-OLS","BWFD-LMLS","BWFD-IRLS","SG-OLS","SG-LMLS","SG-IRLS"];
    PIPE_COLS  = [0.85 0.33 0.10;
                  0.93 0.69 0.13;
                  0.49 0.18 0.56;
                  0.47 0.67 0.19;
                  0.30 0.75 0.93;
                  0.00 0.45 0.74];
    PIPE_LSTY  = {'-','--',':','-','--',':'};
    PIPE_LW    = [2.0 1.8 1.8 2.0 1.8 1.8];

    MDC          = 0.03;
    SEM_ADEQUATE = MDC / 2.77;
    PASS_THRESH  = 0.95;

    bg = sort(unique(sub.betaGen));
    nB = numel(bg);

    %% Build matrix: meanBetaRec [nBeta x nPipe], peak and monotonic segment
    brMat    = NaN(nB, numel(PIPELINES));
    peakIdx  = NaN(1, numel(PIPELINES));
    monoEnd  = NaN(1, numel(PIPELINES));   % last idx of rising segment

    for pi = 1:numel(PIPELINES)
        for bi = 1:nB
            % Average over VGF values at this coordinate (VGF marginalized)
            r = sub(sub.pipeline == PIPELINES(pi) & sub.betaGen == bg(bi), :);
            if ~isempty(r)
                brMat(bi, pi) = mean(r.meanBetaRec, 'omitnan');
            end
        end
        % Find peak (bend-back point)
        [~, pk] = max(brMat(:, pi));
        peakIdx(pi) = pk;
        % Find last monotonically increasing index before peak.
        % The beta_gen=0 region can show ill-conditioning artifacts (negative
        % curvature denominator instability). Skip the first two grid points
        % (beta_gen <= 0.05) when searching for the monotonic segment start.
        SKIP = 2;   % skip first SKIP grid points (beta_gen <= ~0.067)
        diffs = diff(brMat(1:pk, pi));
        significantDrop = find(diffs(SKIP:end) < -0.02, 1, 'last');
        if isempty(significantDrop)
            monoEnd(pi) = pk;
        else
            monoEnd(pi) = significantDrop + SKIP - 1;
        end
    end

    %% Figure
    fig = figure('Name', sprintf('Forward map: %s', datasetLabel), ...
                 'Position', [80 80 860 560]);
    ax = axes(fig); hold(ax, 'on');

    % Identity line
    plot(ax, [0 2/3], [0 2/3], 'k:', 'LineWidth', 1.2, 'HandleVisibility','off');
    text(ax, 0.68, 0.68, '\beta_{rec} = \beta_{gen}', 'FontSize', 8, ...
         'Color', [0.3 0.3 0.3], 'HandleVisibility','off');

    % MDC band around 1/3
    patch(ax, [0 max(bg) max(bg) 0], ...
          [1/3-MDC 1/3-MDC 1/3+MDC 1/3+MDC], ...
          [0.95 0.90 0.90], 'EdgeColor','none','FaceAlpha',0.5,...
          'HandleVisibility','off');
    text(ax, 0.02, 1/3+MDC+0.005, '\pm MDC around 1/3', ...
         'FontSize', 7.5, 'Color', [0.65 0.2 0.2], 'HandleVisibility','off');

    % Plot each pipeline curve
    for pi = 1:numel(PIPELINES)
        col = PIPE_COLS(pi,:);
        lsy = PIPE_LSTY{pi};
        lw  = PIPE_LW(pi);
        br  = brMat(:, pi);

        % Full curve (faded beyond peak)
        plot(ax, bg, br, lsy, 'Color', [col 0.35], 'LineWidth', lw*0.7, ...
             'HandleVisibility','off');

        % Monotonic rising segment (solid, full opacity)
        mEnd = monoEnd(pi);
        plot(ax, bg(1:mEnd), br(1:mEnd), lsy, 'Color', col, ...
             'LineWidth', lw, 'DisplayName', PIPELINES(pi));

        % Bend-back peak marker
        pk = peakIdx(pi);
        plot(ax, bg(pk), br(pk), 'v', 'Color', col, ...
             'MarkerFaceColor', col, 'MarkerSize', 7, ...
             'HandleVisibility','off');
    end

    % Vertical line at peak region (SG-IRLS)
    sgIdx = find(PIPELINES == "SG-IRLS");
    pkA   = bg(peakIdx(sgIdx));
    xline(ax, pkA, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, ...
          'HandleVisibility','off');
    text(ax, pkA+0.01, 0.05, sprintf('Bend-back\n\\beta_{gen}\\approx%.2f', pkA), ...
         'FontSize', 8, 'Color', [0.3 0.3 0.3], 'HandleVisibility','off');

    % Ambiguity zone shading
    pkMin = min(bg(peakIdx));
    pkMax = max(bg(peakIdx));
    patch(ax, [pkMin pkMax pkMax pkMin], [0 0 0.75 0.75], ...
          [0.85 0.85 0.85], 'EdgeColor','none','FaceAlpha',0.25,...
          'HandleVisibility','off');
    text(ax, (pkMin+pkMax)/2, 0.71, 'Ambiguous zone', ...
         'HorizontalAlignment','center','FontSize',8,'Color',[0.4 0.4 0.4],...
         'HandleVisibility','off');

    %% Annotations
    xlabel(ax, '\beta_{gen}  (ground truth)', 'FontSize', 12);
    ylabel(ax, '\beta_{rec}  (recovered)',    'FontSize', 12);
    title(ax, {sprintf('\\bf Forward map: %s', datasetLabel), ...
        sprintf('\\alpha=%.3f, \\sigma=%.2f mm, fs=%d Hz  (snapped: \\alpha=%.2f, \\sigma=%.2f mm)', ...
        tgtAlpha, tgtSigma, tgtFs, snapAlpha, snapSigma)}, 'FontSize', 10);

    lgd = legend(ax, 'Location', 'northwest', 'FontSize', 9);
    title(lgd, 'Solid = monotonic segment');

    % Console summary
    fprintf('\nBend-back peak per pipeline:\n');
    fprintf('%-12s  %8s  %8s  %8s\n','Pipeline','peak_bg','peak_br','monoEnd_bg');
    for pi = 1:numel(PIPELINES)
        pk = peakIdx(pi);
        fprintf('%-12s  %8.4f  %8.4f  %8.4f\n', ...
            PIPELINES(pi), bg(pk), brMat(pk,pi), bg(monoEnd(pi)));
    end

    xlim(ax, [0 max(bg)+0.02]); ylim(ax, [0 0.75]);
    grid(ax, 'on'); box(ax, 'on'); hold(ax, 'off');

    %% 4. Local monotonicity bounds — console + figure annotation
    fprintf('\nLocal monotonic LB (beta_rec>0.20, slope +ve for 3 pts):\n');
    fprintf('%-12s  %8s  %8s\n', 'Pipeline', 'LB_bg', 'peak_bg');
    localLB = NaN(1, numel(PIPELINES));
    for pi2 = 1:numel(PIPELINES)
        br2 = brMat(:, pi2);  pk2 = peakIdx(pi2);  lb2 = NaN;
        for k = 2:pk2-3
            if br2(k) > 0.20 && all(diff(br2(k:k+3)) > 0)
                lb2 = bg(k);  break;
            end
        end
        localLB(pi2) = lb2;
        fprintf('%-12s  %8.4f  %8.4f\n', PIPELINES(pi2), lb2, bg(pk2));
    end
    safeZoneLo = max(localLB, [], 'omitnan');
    safeZoneHi = min(bg(peakIdx));
    fprintf('Conservative monotonic zone: [%.3f, %.3f]\n', safeZoneLo, safeZoneHi);

    % Shade the safe local monotonic zone on the figure
    hold(ax, 'on');
    patch(ax, [safeZoneLo safeZoneHi safeZoneHi safeZoneLo], [0 0 0.75 0.75], ...
          [0.70 0.85 0.70], 'EdgeColor', 'none', 'FaceAlpha', 0.20, ...
          'HandleVisibility', 'off');
    text(ax, safeZoneLo + 0.01, 0.62, ...
         sprintf('Locally monotonic\n[%.2f, %.2f]', safeZoneLo, safeZoneHi), ...
         'FontSize', 8, 'Color', [0.12 0.48 0.12], 'HandleVisibility', 'off');
    hold(ax, 'off');

    %% 5. Save
    safeLabel = regexprep(datasetLabel, '[^a-zA-Z0-9]', '_');
    outFile = fullfile('..','figures', ...
        sprintf('forwardMapCentroid_%s_a%.2f_s%.2f_v001.png', safeLabel, snapAlpha, snapSigma));
    saveas(fig, outFile);
    fprintf('\nSaved: %s\n', outFile);

end
