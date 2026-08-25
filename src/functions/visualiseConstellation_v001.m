function visualiseConstellation_v001(betaCanon, pairDiffs, bioResults, subIDs, opts)
% visualiseConstellation_v001  Generic 4-panel constellation consistency figure
%
% Panels:
%   A: Trial constellation profiles (6 pipelines), coloured by subject
%   B: Biological noise alpha vs sigma, coloured by subject
%   C: Within-subject CCC (horizontal bar, sorted)
%   D: Constellation spread vs composite noise impact index
%
% USAGE:
%   visualiseConstellation_v001(betaCanon, pairDiffs, bioResults, subIDs)
%   visualiseConstellation_v001(..., Title="Cook et al. (2026)", SigmaToMM=0.248)
%
% INPUTS:
%   betaCanon  - [nTrials x 6] recovered beta for 6 canonical pipelines
%   pairDiffs  - [nTrials x 15] pairwise pipeline differences
%   bioResults - table with .sigmaMean, .alphaMean, .subjectID columns
%   subIDs     - [nTrials x 1] string array of subject IDs
%
% NAME-VALUE:
%   Title       - sgtitle string (default: "Constellation consistency")
%   SigmaToMM   - multiplier to convert bioResults.sigmaMean to mm
%                  (default: 1.0; Zarandi cm -> mm = 10; Cook WACOM -> mm = 0.248)
%   SigmaUnit   - label for raw sigma axis (default: "data units")
%   PipeLabels  - [1x6] string of pipeline names
%                  (default: ["BW-OLS","BW-LMLS","BW-IRLS","SG-OLS","SG-LMLS","SG-IRLS"])
%
% Fraser, D.S. (2026)

arguments
    betaCanon   (:,6) double
    pairDiffs   (:,15) double
    bioResults  table
    subIDs      (:,1) string
    opts.Title       string = "Constellation consistency"
    opts.SigmaToMM   double = 1.0
    opts.SigmaUnit   string = "data units"
    opts.PipeLabels  (1,6) string = ["BW-OLS","BW-LMLS","BW-IRLS", ...
                                     "SG-OLS","SG-LMLS","SG-IRLS"]
end

uSubs  = unique(subIDs, 'stable');
nSubs  = numel(uSubs);
nPipes = 6;
cmap   = lines(nSubs);

%% Compute within-subject CCC
perSubCCCmean = NaN(nSubs, 1);
perSubCCCsd   = NaN(nSubs, 1);
perSubCCCmin  = NaN(nSubs, 1);
perSubCCCmax  = NaN(nSubs, 1);

for s = 1:nSubs
    mask  = subIDs == uSubs(s);
    subPD = pairDiffs(mask, :);
    nTri  = sum(mask);
    if nTri < 2
        continue;  % perSubCCC* already NaN
    end
    cccVals = NaN(nTri*(nTri-1)/2, 1);
    k = 0;
    for i = 1:nTri
        for j = (i+1):nTri
            k = k + 1;
            res = linCCC_v001(subPD(i,:)', subPD(j,:)');
            cccVals(k) = res.ccc;
        end
    end
    perSubCCCmean(s) = mean(cccVals);
    perSubCCCsd(s)   = std(cccVals);
    perSubCCCmin(s)  = min(cccVals);
    perSubCCCmax(s)  = max(cccVals);
end

trialSpread = max(betaCanon, [], 2) - min(betaCanon, [], 2);

%% Figure
figure('Position', [50 50 1400 900], 'Color', 'w');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel A: Trial constellation profiles ---
nexttile;
hold on;
hLeg = gobjects(nSubs, 1);
for s = 1:nSubs
    mask = subIDs == uSubs(s);
    idx  = find(mask);
    for k = 1:numel(idx)
        h = plot(1:nPipes, betaCanon(idx(k), :), '-o', ...
            'Color', [cmap(s,:) 0.4], 'MarkerSize', 3, ...
            'MarkerFaceColor', cmap(s,:), 'LineWidth', 0.8);
        if k == 1, hLeg(s) = h; end
    end
end
hold off;
set(gca, 'XTick', 1:nPipes, 'XTickLabel', opts.PipeLabels, 'XTickLabelRotation', 30);
ylabel('\beta_{rec}');
title('A: Trial constellation profiles (6 pipelines)');
xlim([0.5 6.5]);
yline(1/3, '--', '\beta = 1/3', 'Color', [0.5 0.5 0.5], 'Alpha', 0.6, ...
    'LabelHorizontalAlignment', 'left');
legend(hLeg, uSubs, 'Location', 'bestoutside', 'FontSize', 7, 'NumColumns', 2);

% --- Panel B: Noise parameter space ---
nexttile;
hold on;
hLeg2 = gobjects(nSubs, 1);
for s = 1:nSubs
    mask = subIDs == uSubs(s);
    hLeg2(s) = scatter(bioResults.sigmaMean(mask), bioResults.alphaMean(mask), ...
        30, cmap(s,:), 'filled', 'MarkerFaceAlpha', 0.7);
end
hold off;
xlabel(sprintf('\\sigma_{bio} (%s)', opts.SigmaUnit));
ylabel('\alpha_{bio} (spectral exponent)');
title('B: Biological noise per trial');
legend(hLeg2, uSubs, 'Location', 'bestoutside', 'FontSize', 7, 'NumColumns', 2);
for s = 1:nSubs
    mask = subIDs == uSubs(s);
    text(mean(bioResults.sigmaMean(mask)), mean(bioResults.alphaMean(mask)), ...
        "  " + uSubs(s), 'FontSize', 7, 'Color', cmap(s,:));
end

% --- Panel C: Within-subject CCC ---
nexttile;
hasCCC = any(~isnan(perSubCCCmean));
if hasCCC
    [~, sortIdx] = sort(perSubCCCmean, 'descend');
    barh(1:nSubs, perSubCCCmean(sortIdx));
    hold on;
    for s = 1:nSubs
        origIdx = sortIdx(s);
        if isnan(perSubCCCmean(origIdx)), continue; end
        sd = perSubCCCsd(origIdx);
        if isnan(sd), sd = 0; end
        % patch first (coloured background bar)
        patch([0 perSubCCCmean(origIdx) perSubCCCmean(origIdx) 0], ...
            [s-0.35 s-0.35 s+0.35 s+0.35], cmap(origIdx,:), ...
            'FaceAlpha', 0.7, 'EdgeColor', 'none');
        % errorbar on top, thicker so it reads through the patch
        errorbar(perSubCCCmean(origIdx), s, sd, sd, ...
            'horizontal', 'k', 'LineWidth', 1.5, 'CapSize', 5);
        % numeric label at the right end of the SD whisker
        text(min(perSubCCCmean(origIdx) + sd + 0.01, 0.99), s, ...
            sprintf('%.2f \\pm %.2f', perSubCCCmean(origIdx), sd), ...
            'FontSize', 7, 'VerticalAlignment', 'middle', ...
            'HorizontalAlignment', 'left', 'Color', [0.2 0.2 0.2]);
    end
    hold off;
    set(gca, 'YTick', 1:nSubs, 'YTickLabel', uSubs(sortIdx), 'YDir', 'reverse');
    xlabel('Within-subject CCC (mean \pm SD)');
    xlim([0 1.15]);  % widened from [0 1] so the numeric labels fitxline(0.70, ':', 'good', 'Color', [0.5 0.5 0.5], 'Alpha', 0.5, 'FontSize', 8);
    xline(0.90, ':', 'excellent', 'Color', [0.5 0.5 0.5], 'Alpha', 0.5, 'FontSize', 8);
else
    text(0.5, 0.5, {'Within-subject CCC unavailable:', ...
        'requires >= 2 trials per subject'}, ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 11, 'Color', [0.5 0.5 0.5]);
    set(gca, 'XTick', [], 'YTick', []);
end
title('C: Constellation consistency by subject');

% --- Panel D: Spread vs composite noise impact ---
nexttile;
boundaryMM = @(a) 0.375 .* exp(0.63 .* a);
sigmaMM    = bioResults.sigmaMean * opts.SigmaToMM;
noiseImpact = sigmaMM ./ boundaryMM(bioResults.alphaMean);

hold on;
for s = 1:nSubs
    mask = subIDs == uSubs(s);
    scatter(noiseImpact(mask), trialSpread(mask)*1000, ...
        30, cmap(s,:), 'filled', 'MarkerFaceAlpha', 0.7);
end
xline(1, ':', 'invertibility boundary', 'Color', [0.5 0.5 0.5], ...
    'Alpha', 0.5, 'FontSize', 8, 'LabelVerticalAlignment', 'bottom');
hold off;
xlabel('Noise impact (\sigma_{mm} / boundary(\alpha))');
ylabel('Constellation spread (\times10^{-3})');
title('D: Spread vs composite noise impact');

validIdx = ~isnan(noiseImpact) & ~isnan(trialSpread);
[rho, pval] = corr(noiseImpact(validIdx), trialSpread(validIdx));
text(0.05, 0.92, sprintf('r = %.3f, p = %.1e', rho, pval), ...
    'Units', 'normalized', 'FontSize', 9, 'Color', [0.3 0.3 0.3]);

sgtitle(opts.Title, 'FontSize', 14, 'FontWeight', 'bold');

%% Console summary
fprintf('\n=== CONSTELLATION SUMMARY: %s ===\n\n', opts.Title);
fprintf('%-8s %5s %8s %8s %8s %8s  |  %8s %8s %8s\n', ...
    'Sub', 'nTri', 'CCC_m', 'CCC_sd', 'CCC_min', 'CCC_max', 'alpha_m', 'sig_mm', 'spread');
fprintf('%s\n', repmat('-', 1, 84));
for s = 1:nSubs
    mask = subIDs == uSubs(s);
    fprintf('%-8s %5d %8.3f %8.3f %8.3f %8.3f  |  %8.2f %8.2f %8.4f\n', ...
        uSubs(s), sum(mask), perSubCCCmean(s), perSubCCCsd(s), ...
        perSubCCCmin(s), perSubCCCmax(s), ...
        mean(bioResults.alphaMean(mask)), ...
        mean(sigmaMM(mask)), ...
        mean(trialSpread(mask)));
end
fprintf('\nOverall within-subject CCC: mean=%.3f, median=%.3f (N=%d subjects)\n', ...
    mean(perSubCCCmean), median(perSubCCCmean), nSubs);
end
