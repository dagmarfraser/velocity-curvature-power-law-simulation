function plotLMMCoefficients_v002()
% plotLMMCoefficients_v002  Fig 1: LMM coefficient display -- three focused panels
%
% Panel A: All 9 main effects, horizontal forest plot ordered by |Estimate|.
% Panel B: Top 25 significant interactions by |Estimate|, coloured by
%          interaction order (2-way / 3-way / 4-way+). Labels truncated to 35 chars.
% Panel C: Significance summary bar (sig-main / sig-interaction / n.s.).
%
% Replaces v001 (192-row forest plot -- unreadable).
%
% Output: ../figures/plotLMMCoefficients_v002_mains.png
%         ../figures/plotLMMCoefficients_v002_interactions.png
%         ../figures/plotLMMCoefficients_v002_summary.png
%
% Run from src/. Requires L9_coefficients_v004.csv in src/.
% Fraser, Di Luca, Cook (2026) -- PreReg simulation study.

    csvPath = fullfile(fileparts(mfilename('fullpath')), 'L9_coefficients_v004.csv');
    if ~isfile(csvPath)
        error('plotLMMCoefficients_v002:missingFile', '%s', ...
            'L9_coefficients_v004.csv not found. Run from src/.');
    end

    T = readtable(csvPath, 'TextType', 'string');
    T.Label   = cleanLabels(T.Name);
    T.nColons = arrayfun(@(n) sum(char(n) == ':'), T.Name);
    T.Order   = T.nColons + 1;

    figDir = fullfile(fileparts(mfilename('fullpath')), '..', 'figures');

    colMain  = [0.18 0.45 0.69];
    col2way  = [0.80 0.33 0.00];
    col3way  = [0.47 0.67 0.19];
    col4plus = [0.55 0.27 0.07];
    colNS    = [0.70 0.70 0.70];

    %% Panel A: Main effects (9 terms including intercept)
    Tm = T(T.IsInteraction == 0, :);
    [~, ord] = sort(abs(Tm.Estimate), 'ascend');
    Tm = Tm(ord, :);
    nM = height(Tm);

    figA = figure('Units', 'centimeters', 'Position', [2 2 18 9], 'Color', 'w');
    axA  = axes(figA);
    hold(axA, 'on');

    for k = 1:nM
        col = colMain;
        if Tm.Significant(k) == 0, col = colNS; end
        line(axA, [Tm.Lower(k) Tm.Upper(k)], [k k], 'Color', col, 'LineWidth', 2.5);
        scatter(axA, Tm.Estimate(k), k, 50, col, 'filled');
    end
    xline(axA, 0, 'k--', 'LineWidth', 0.8, 'Alpha', 0.4);

    yticks(axA, 1:nM);
    yticklabels(axA, Tm.Label);
    axA.TickLabelInterpreter = 'none';
    axA.FontSize = 10;
    xlabel(axA, 'Standardised estimate (95% CI)', 'FontSize', 10);
    title(axA, 'Main effects  (N=17,458,535; R^2=0.9912)', ...
          'FontSize', 11, 'FontWeight', 'bold');
    box(axA, 'off');

    outA = fullfile(figDir, 'plotLMMCoefficients_v002_mains.png');
    exportgraphics(figA, outA, 'Resolution', 150);
    fprintf('Saved: %s\n', outA);
    close(figA);

    %% Panel B: Top 25 significant interactions, truncated labels
    Ti = T(T.IsInteraction == 1 & T.Significant == 1, :);
    [~, ord] = sort(abs(Ti.Estimate), 'descend');
    nTop = min(25, height(Ti));
    Ti = Ti(ord(1:nTop), :);
    % Re-sort ascending for forest plot (largest |estimate| at top)
    [~, ord2] = sort(abs(Ti.Estimate), 'ascend');
    Ti = Ti(ord2, :);
    nI = height(Ti);

    % Truncate labels to 35 chars, append order tag
    shortLabels = cell(nI, 1);
    for k = 1:nI
        lbl = char(Ti.Label(k));
        if numel(lbl) > 35, lbl = [lbl(1:32) '...']; end
        shortLabels{k} = sprintf('%s  [%d-way]', lbl, Ti.Order(k));
    end

    colours = zeros(nI, 3);
    for k = 1:nI
        switch Ti.Order(k)
            case 2,    colours(k,:) = col2way;
            case 3,    colours(k,:) = col3way;
            otherwise, colours(k,:) = col4plus;
        end
    end

    figB = figure('Units', 'centimeters', 'Position', [2 2 26 16], 'Color', 'w');
    axB  = axes(figB);
    hold(axB, 'on');

    for k = 1:nI
        line(axB, [Ti.Lower(k) Ti.Upper(k)], [k k], 'Color', colours(k,:), 'LineWidth', 2);
        scatter(axB, Ti.Estimate(k), k, 40, colours(k,:), 'filled');
    end
    xline(axB, 0, 'k--', 'LineWidth', 0.8, 'Alpha', 0.4);

    yticks(axB, 1:nI);
    yticklabels(axB, shortLabels);
    axB.TickLabelInterpreter = 'none';
    axB.FontSize = 7.5;
    axB.YAxis.FontSize = 7;
    xlabel(axB, 'Standardised estimate (95% CI)', 'FontSize', 9);
    title(axB, sprintf('Top %d significant interactions by |Estimate|', nTop), ...
          'FontSize', 10, 'FontWeight', 'bold');

    h2 = line(axB, nan, nan, 'Color', col2way,  'LineWidth', 2.5);
    h3 = line(axB, nan, nan, 'Color', col3way,  'LineWidth', 2.5);
    h4 = line(axB, nan, nan, 'Color', col4plus, 'LineWidth', 2.5);
    legend(axB, [h2 h3 h4], {'2-way', '3-way', '4-way+'}, ...
           'Location', 'southeast', 'FontSize', 8);
    box(axB, 'off');

    outB = fullfile(figDir, 'plotLMMCoefficients_v002_interactions.png');
    exportgraphics(figB, outB, 'Resolution', 150);
    fprintf('Saved: %s\n', outB);
    close(figB);

    %% Panel C: Significance summary bar
    nSigMain = sum(T.Significant == 1 & T.IsInteraction == 0);
    nSigIntr = sum(T.Significant == 1 & T.IsInteraction == 1);
    nNS      = sum(T.Significant == 0);
    vals     = [nSigMain, nSigIntr, nNS];
    lbls     = {'Sig. main', 'Sig. interaction', 'n.s.'};

    figC = figure('Units', 'centimeters', 'Position', [2 2 12 9], 'Color', 'w');
    axC  = axes(figC);
    bh   = bar(axC, 1:3, vals, 'FaceColor', 'flat', 'EdgeColor', 'none');
    bh.CData = [colMain; col2way; colNS];
    xticks(axC, 1:3);
    xticklabels(axC, lbls);
    axC.FontSize = 9;
    axC.XAxis.TickLabelInterpreter = 'none';
    ylabel(axC, 'Count', 'FontSize', 9);
    title(axC, sprintf('192 fixed effects: significance breakdown'), ...
          'FontSize', 9, 'FontWeight', 'bold');
    for k = 1:3
        text(axC, k, vals(k) + 0.8, num2str(vals(k)), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    end
    ylim(axC, [0 max(vals) * 1.15]);
    box(axC, 'off');

    outC = fullfile(figDir, 'plotLMMCoefficients_v002_summary.png');
    exportgraphics(figC, outC, 'Resolution', 150);
    fprintf('Saved: %s\n', outC);
    close(figC);
end

%% -- Helpers
function labels = cleanLabels(names)
    names  = cellstr(names);
    labels = strrep(names,  'betaGenerated',    'beta_gen');
    labels = strrep(labels, 'noiseMagnitude',   'sigma');
    labels = strrep(labels, 'noiseColor',       'alpha');
    labels = strrep(labels, 'samplingRate',     'fs');
    labels = strrep(labels, 'filterType_6',     'SG');
    labels = strrep(labels, 'regressionType_4', 'LMLS');
    labels = strrep(labels, 'regressionType_5', 'IRLS');
    labels = strrep(labels, '(Intercept)',      'Intercept');
end
