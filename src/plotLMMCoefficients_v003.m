function plotLMMCoefficients_v003()
% plotLMMCoefficients_v003  Fig 2: Stage 1 L9 LMM coefficients, three panels (as v002).
%   A: all 9 main effects (incl. intercept), forest plot by |Estimate|.
%   B: top 25 significant interactions by |Estimate|, coloured by order.
%   C: significance breakdown of the 192 fixed effects.
% v003: reads the fitted model's fixed effects from results/L9ModelInspect_v002.mat.
%   v002 read src/L9_coefficients_v004.csv, whose interactions do not match the fitted
%   model (Finding #237); its interaction panel plotted wrong values. CIs are Wald,
%   estimate +/- 1.96 SE, and p two-sided normal (t with df ~ 1.7e7 is identical).
%   Prints the counts behind Part 1 Results (a)'s "all nine main effects significant;
%   123 of 183 interaction terms significant".
% Run from src/.  Output: figures/plotLMMCoefficients_v003_{mains,interactions,summary}.png
% Fraser, Di Luca, Cook (2026)  v003

    SRC = fileparts(mfilename('fullpath'));  ROOT = fileparts(SRC);
    inMat  = fullfile(ROOT, 'results', 'L9ModelInspect_v002.mat');
    figDir = fullfile(ROOT, 'figures');
    R2_43  = 0.9912;                                   % Finding #43 (Adjusted R^2, L9 run)
    if ~isfile(inMat), error('plotLMM3:input', '%s', "Missing: " + inMat); end
    if ~isfolder(figDir), error('plotLMM3:outDir', '%s', "Missing folder: " + figDir); end
    I = load(inMat, 'out');  F = I.out.fixed;
    if height(F) ~= 192, error('plotLMM3:nCoef', '%s', sprintf('%d coefficients, expected 192', height(F))); end

    T = table(F.name, F.est, F.SE, 'VariableNames', {'Name', 'Estimate', 'SE'});
    T.Lower = T.Estimate - 1.96 * T.SE;  T.Upper = T.Estimate + 1.96 * T.SE;
    T.pValue = 2 * normcdf(-abs(T.Estimate ./ T.SE));
    T.Significant   = T.pValue < 0.05;
    T.IsInteraction = contains(T.Name, ":");
    T.Label = cleanLabels(T.Name);
    T.Order = count(T.Name, ":") + 1;
    fprintf('Main effects (incl. intercept): %d, significant %d | interactions: %d, significant %d\n', ...
        sum(~T.IsInteraction), sum(~T.IsInteraction & T.Significant), sum(T.IsInteraction), ...
        sum(T.IsInteraction & T.Significant));

    colMain = [0.18 0.45 0.69];  col2 = [0.80 0.33 0.00];  col3 = [0.47 0.67 0.19];
    col4 = [0.55 0.27 0.07];     colNS = [0.70 0.70 0.70];

    %% Panel A: main effects
    Tm = T(~T.IsInteraction, :);  [~, o] = sort(abs(Tm.Estimate), 'ascend');  Tm = Tm(o, :);
    figA = figure('Units', 'centimeters', 'Position', [2 2 18 9], 'Color', 'w');  axA = axes(figA);  hold(axA, 'on');
    for k = 1:height(Tm)
        c = colMain;  if ~Tm.Significant(k), c = colNS; end
        line(axA, [Tm.Lower(k) Tm.Upper(k)], [k k], 'Color', c, 'LineWidth', 2.5);
        scatter(axA, Tm.Estimate(k), k, 50, c, 'filled');
    end
    xline(axA, 0, 'k--', 'LineWidth', 0.8, 'Alpha', 0.4);
    yticks(axA, 1:height(Tm));  yticklabels(axA, Tm.Label);  axA.TickLabelInterpreter = 'none';  axA.FontSize = 10;
    xlabel(axA, 'Standardised estimate (95% CI)');
    title(axA, sprintf('Main effects  (N=%s; R^2=%.4f)', num2sepstr_local(I.out.nObs), R2_43), 'FontSize', 11);
    box(axA, 'off');  save_local(figA, fullfile(figDir, 'plotLMMCoefficients_v003_mains.png'));

    %% Panel B: top 25 significant interactions
    Ti = T(T.IsInteraction & T.Significant, :);  [~, o] = sort(abs(Ti.Estimate), 'descend');
    Ti = Ti(o(1:min(25, height(Ti))), :);  [~, o] = sort(abs(Ti.Estimate), 'ascend');  Ti = Ti(o, :);
    lbl = cell(height(Ti), 1);  cols = zeros(height(Ti), 3);
    for k = 1:height(Ti)
        s = char(Ti.Label(k));  if numel(s) > 35, s = [s(1:32) '...']; end
        lbl{k} = sprintf('%s  [%d-way]', s, Ti.Order(k));
        switch Ti.Order(k), case 2, cols(k, :) = col2; case 3, cols(k, :) = col3; otherwise, cols(k, :) = col4; end
    end
    figB = figure('Units', 'centimeters', 'Position', [2 2 26 16], 'Color', 'w');  axB = axes(figB);  hold(axB, 'on');
    for k = 1:height(Ti)
        line(axB, [Ti.Lower(k) Ti.Upper(k)], [k k], 'Color', cols(k, :), 'LineWidth', 2);
        scatter(axB, Ti.Estimate(k), k, 40, cols(k, :), 'filled');
    end
    xline(axB, 0, 'k--', 'LineWidth', 0.8, 'Alpha', 0.4);
    yticks(axB, 1:height(Ti));  yticklabels(axB, lbl);  axB.TickLabelInterpreter = 'none';  axB.YAxis.FontSize = 7;
    xlabel(axB, 'Standardised estimate (95% CI)');
    title(axB, sprintf('Top %d significant interactions by |Estimate|', height(Ti)), 'FontSize', 10);
    h = [line(axB, nan, nan, 'Color', col2, 'LineWidth', 2.5), line(axB, nan, nan, 'Color', col3, 'LineWidth', 2.5), ...
         line(axB, nan, nan, 'Color', col4, 'LineWidth', 2.5)];
    legend(axB, h, {'2-way', '3-way', '4-way+'}, 'Location', 'southeast', 'FontSize', 8);  box(axB, 'off');
    save_local(figB, fullfile(figDir, 'plotLMMCoefficients_v003_interactions.png'));

    %% Panel C: significance breakdown
    vals = [sum(T.Significant & ~T.IsInteraction), sum(T.Significant & T.IsInteraction), sum(~T.Significant)];
    figC = figure('Units', 'centimeters', 'Position', [2 2 12 9], 'Color', 'w');  axC = axes(figC);
    bh = bar(axC, 1:3, vals, 'FaceColor', 'flat', 'EdgeColor', 'none');  bh.CData = [colMain; col2; colNS];
    xticks(axC, 1:3);  xticklabels(axC, {'Sig. main', 'Sig. interaction', 'n.s.'});  ylabel(axC, 'Count');
    title(axC, sprintf('%d fixed effects: significance breakdown', height(T)), 'FontSize', 9);
    for k = 1:3, text(axC, k, vals(k) + 0.8, num2str(vals(k)), 'HorizontalAlignment', 'center', 'FontWeight', 'bold'); end
    ylim(axC, [0 max(vals) * 1.15]);  box(axC, 'off');
    save_local(figC, fullfile(figDir, 'plotLMMCoefficients_v003_summary.png'));
end

function labels = cleanLabels(names)
    labels = replace(string(names), ["betaGenerated" "noiseMagnitude" "noiseColor" "samplingRate" ...
        "filterType_6" "regressionType_4" "regressionType_5" "(Intercept)"], ...
        ["beta_gen" "sigma" "alpha" "fs" "SG" "LMLS" "IRLS" "Intercept"]);
end

function save_local(fig, f)
    set(findall(fig, 'Type', 'axes'), 'Toolbar', []);
    exportgraphics(fig, f, 'Resolution', 150);  fprintf('Saved: %s\n', f);  close(fig);
end

function s = num2sepstr_local(n)
    s = regexprep(sprintf('%d', n), '(\d)(?=(\d{3})+$)', '$1,');
end
