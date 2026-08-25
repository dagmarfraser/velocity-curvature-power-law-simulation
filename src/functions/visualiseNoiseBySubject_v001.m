function visualiseNoiseBySubject_v001(bioResults, opts)
% visualiseNoiseBySubject_v001  Generic 4-panel per-subject noise variability
%
% Panels:
%   A: sigma_bio by subject (trial dots + mean/SD)
%   B: alpha_bio by subject (trial dots + mean/SD)
%   C: Within-subject CV for sigma and alpha
%   D: f0 (tracing frequency) by subject
%
% USAGE:
%   visualiseNoiseBySubject_v001(bioResults)
%   visualiseNoiseBySubject_v001(bioResults, Title="Cook", SigmaToMM=0.248)
%
% INPUTS:
%   bioResults - table with .sigmaMean, .alphaMean, .subjectID, .f0 columns
%
% NAME-VALUE:
%   Title      - sgtitle string (default: "Per-subject noise characterisation")
%   SigmaToMM  - multiplier to convert sigmaMean to mm (default: 1.0)
%   SigmaUnit  - label for raw sigma axis (default: "data units")
%
% Fraser, D.S. (2026)

arguments
    bioResults  table
    opts.Title      string = "Per-subject noise characterisation (pmtm)"
    opts.SigmaToMM  double = 1.0
    opts.SigmaUnit  string = "data units"
end

subIDs = bioResults.subjectID;
uSubs  = unique(subIDs, 'stable');
nSubs  = numel(uSubs);
cmap   = lines(nSubs);
sigmaMM = bioResults.sigmaMean * opts.SigmaToMM;

%% Figure
figure('Position', [50 50 1400 800], 'Color', 'w');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel A: sigma by subject ---
nexttile;
hold on;
for s = 1:nSubs
    mask = subIDs == uSubs(s);
    vals = sigmaMM(mask);
    nTri = sum(mask);
    jitter = 0.15 * (rand(nTri, 1) - 0.5);
    scatter(s + jitter, vals, 20, cmap(s,:), 'filled', 'MarkerFaceAlpha', 0.6);
    plot([s-0.3 s+0.3], [mean(vals) mean(vals)], '-', 'Color', cmap(s,:), 'LineWidth', 2);
    plot([s s], [mean(vals)-std(vals) mean(vals)+std(vals)], '-', ...
        'Color', cmap(s,:), 'LineWidth', 1.5);
end
hold off;
set(gca, 'XTick', 1:nSubs, 'XTickLabel', uSubs, 'XTickLabelRotation', 45);
ylabel('\sigma_{bio} (mm)');
title('A: Biological noise magnitude by subject');
xlim([0.3 nSubs+0.7]);

% --- Panel B: alpha by subject ---
nexttile;
hold on;
for s = 1:nSubs
    mask = subIDs == uSubs(s);
    vals = bioResults.alphaMean(mask);
    nTri = sum(mask);
    jitter = 0.15 * (rand(nTri, 1) - 0.5);
    scatter(s + jitter, vals, 20, cmap(s,:), 'filled', 'MarkerFaceAlpha', 0.6);
    plot([s-0.3 s+0.3], [mean(vals) mean(vals)], '-', 'Color', cmap(s,:), 'LineWidth', 2);
    plot([s s], [mean(vals)-std(vals) mean(vals)+std(vals)], '-', ...
        'Color', cmap(s,:), 'LineWidth', 1.5);
end
hold off;
set(gca, 'XTick', 1:nSubs, 'XTickLabel', uSubs, 'XTickLabelRotation', 45);
ylabel('\alpha_{bio} (spectral exponent)');
title('B: Biological noise colour by subject');
xlim([0.3 nSubs+0.7]);

% --- Panel C: Within-subject CV ---
nexttile;
cvSigma = NaN(nSubs, 1);
cvAlpha = NaN(nSubs, 1);
for s = 1:nSubs
    mask = subIDs == uSubs(s);
    sig = sigmaMM(mask);
    alp = bioResults.alphaMean(mask);
    cvSigma(s) = std(sig) / mean(sig) * 100;
    cvAlpha(s) = std(alp) / mean(alp) * 100;
end
b = bar([cvSigma, cvAlpha], 'grouped');
b(1).FaceColor = '#378ADD';
b(2).FaceColor = '#D85A30';
set(gca, 'XTick', 1:nSubs, 'XTickLabel', uSubs, 'XTickLabelRotation', 45);
ylabel('Within-subject CV (%)');
title('C: Measurement variability (CV)');
legend({'\sigma_{bio}', '\alpha_{bio}'}, 'Location', 'northeast');
xlim([0.3 nSubs+0.7]);

% --- Panel D: f0 by subject ---
nexttile;
hold on;
subMeanF0    = NaN(nSubs, 1);
subMeanAlpha = NaN(nSubs, 1);
subMeanSigma = NaN(nSubs, 1);
for s = 1:nSubs
    mask = subIDs == uSubs(s);
    vals = bioResults.f0(mask);
    nTri = sum(mask);
    jitter = 0.15 * (rand(nTri, 1) - 0.5);
    scatter(s + jitter, vals, 20, cmap(s,:), 'filled', 'MarkerFaceAlpha', 0.6);
    plot([s-0.3 s+0.3], [mean(vals) mean(vals)], '-', 'Color', cmap(s,:), 'LineWidth', 2);
    subMeanF0(s)    = mean(vals);
    subMeanAlpha(s) = mean(bioResults.alphaMean(mask));
    subMeanSigma(s) = mean(sigmaMM(mask));
end
hold off;
set(gca, 'XTick', 1:nSubs, 'XTickLabel', uSubs, 'XTickLabelRotation', 45);
ylabel('f_0 (Hz)');
title('D: Tracing frequency by subject');
xlim([0.3 nSubs+0.7]);

sgtitle(opts.Title, 'FontSize', 14, 'FontWeight', 'bold');

%% Console summary
fprintf('\n=== PER-SUBJECT NOISE SUMMARY: %s ===\n\n', opts.Title);
fprintf('%-8s %5s %8s %7s %8s %7s %8s %7s\n', ...
    'Sub', 'nTri', 'sig_mm', 'CV_sig', 'alpha', 'CV_alp', 'f0_Hz', 'CV_f0');
fprintf('%s\n', repmat('-', 1, 65));
for s = 1:nSubs
    mask = subIDs == uSubs(s);
    sig = sigmaMM(mask);
    alp = bioResults.alphaMean(mask);
    f0  = bioResults.f0(mask);
    fprintf('%-8s %5d %8.2f %6.1f%% %8.2f %6.1f%% %8.3f %6.1f%%\n', ...
        uSubs(s), sum(mask), mean(sig), std(sig)/mean(sig)*100, ...
        mean(alp), std(alp)/mean(alp)*100, ...
        mean(f0), std(f0)/mean(f0)*100);
end

[r1, p1] = corr(subMeanF0, subMeanAlpha);
[r2, p2] = corr(subMeanF0, subMeanSigma);
fprintf('\nCross-subject: f0 vs alpha r=%+.3f p=%.3f | f0 vs sigma r=%+.3f p=%.3f\n', ...
    r1, p1, r2, p2);
end
