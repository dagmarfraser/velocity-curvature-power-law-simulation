% checkEmpiricalTempoVsGridVGF_v001.m  Where each dataset's tempo sits relative
% to the registered grid's VGF range, per trial and by summary statistic.
%
% PURPOSE: resolves coherence-pass item A16 and the Session 124 units flag.
%   (1) Units: every published empirical VGF (#207, #208, D17, Part 1 (c)) is
%       f0 x K_conv, with K_conv from the grid's own reference ellipse in px
%       (baselineShp6_120Hz.mat). f0 is in Hz, so these values are already on
%       the grid's VGF axis; no mm-to-px conversion applies. Not vgfObsMM.
%   (2) A16: v002 Part 1 (c) ("Cook CTRL 28% below the floor", Cook ASD not
%       listed) used MEAN f0 (checkVGFSEMCentroid_v002 L125); D17/#207/#208
%       ("Cook CTRL 46.4, Cook ASD 66.0") used MEDIAN f0. Both are reported.
%   (3) Per-trial share of each dataset outside the grid's tempo range.
%
% CAVEAT: K_conv is evaluated at BETA_REF = 1/3, as in every consumer. At fixed
%   VGF the grid's tempo rises with beta_gen (Finding #224), so "inside the
%   range" means inside the grid's tempo range at beta_gen = 1/3.
%
% SELF-CHECKS (Fail Loud): K_conv reproduces 175.2636; median-f0 VGFs reproduce
%   results/datasetOperatingPointCoverage_v001.mat (#207) to 1e-9; Cook CTRL's
%   mean-f0 VGF reproduces v002's "28% below the floor".
%
% INPUTS:  functions/baselineShp6_120Hz.mat; noiseCharacterisation_<name>.mat;
%          results/datasetOperatingPointCoverage_v001.mat (check only).
% OUTPUT:  results/checkEmpiricalTempoVsGridVGF_v001.mat; console table.
% USAGE:   run from src/: checkEmpiricalTempoVsGridVGF_v001

%% CONFIG
BETA_REF   = 1/3;
VGF_LO     = exp(4.5);      % registered grid VGF floor   (90.017)
VGF_HI     = exp(5.8);      % registered grid VGF ceiling (330.30)
KCONV_PUB  = 175.2636;      % published K_conv (#207, #208, Part 1 (c))
DATASETS   = ["Fraser","Cook_CTRL","Cook_ASD","Hickman_PLAC","Hickman_HALO","Zarandi","Dhieb"];
NOISEFILES = ["fraser","cook","cookASD","hickmanPLAC","hickmanHALO","zarandi","dhieb"];

srcDir = fileparts(mfilename('fullpath'));

%% K_conv from the grid's own reference ellipse (as datasetOperatingPointCoverage_v001 L45-55)
shapeFile = fullfile(srcDir, 'functions', 'baselineShp6_120Hz.mat');
if ~isfile(shapeFile)
    error('tempoVsGrid:NoShapeFile', '%s', sprintf('FAILED PATH: %s', shapeFile));
end
S = load(shapeFile, 'pathXYresample', 'K');
perimeterPx = sum(sqrt(sum(diff(S.pathXYresample, 1, 1).^2, 2)));
kappa = S.K(:); kappa = kappa(kappa > 0 & isfinite(kappa));
K_conv = perimeterPx / mean(kappa .^ (-BETA_REF));
if abs(K_conv - KCONV_PUB) > 1e-3
    error('tempoVsGrid:Kconv', '%s', sprintf('K_conv %.6f does not reproduce published %.4f', K_conv, KCONV_PUB));
end
f0Lo = VGF_LO / K_conv; f0Hi = VGF_HI / K_conv;

%% Per-dataset summary and per-trial shares
n = numel(DATASETS);
results = table('Size', [n 10], ...
    'VariableTypes', ["string","double","double","double","double","double","double","double","double","double"], ...
    'VariableNames', ["dataset","nTrials","meanF0","medianF0","vgfMeanF0","vgfMedianF0", ...
                      "pctBelowMean","pctBelowMedian","pctTrialsBelow","pctTrialsAbove"]);
for d = 1:n
    f = fullfile(srcDir, "noiseCharacterisation_" + NOISEFILES(d) + ".mat");
    if ~isfile(f)
        error('tempoVsGrid:NoNoiseMat', '%s', sprintf('FAILED PATH: %s', f));
    end
    b = load(f, 'bioResults').bioResults;
    f0 = double(b.f0(:));
    nBad = sum(~isfinite(f0));
    if nBad > 0
        warning('tempoVsGrid:NonFiniteF0', '%s', sprintf('%s: %d non-finite f0 excluded', DATASETS(d), nBad));
    end
    f0 = f0(isfinite(f0));
    vMean = mean(f0) * K_conv; vMed = median(f0) * K_conv;
    results(d,:) = {DATASETS(d), numel(f0), mean(f0), median(f0), vMean, vMed, ...
        100 * (VGF_LO - vMean) / VGF_LO, 100 * (VGF_LO - vMed) / VGF_LO, ...
        100 * mean(f0 < f0Lo), 100 * mean(f0 > f0Hi)};
end

%% Self-check: median-f0 VGFs reproduce #207's saved table
refFile = fullfile(srcDir, 'results', 'datasetOperatingPointCoverage_v001.mat');
if ~isfile(refFile)
    error('tempoVsGrid:NoRef', '%s', sprintf('FAILED PATH: %s', refFile));
end
R = load(refFile, 'results').results;
for d = 1:n
    j = find(R.dataset == DATASETS(d));
    if numel(j) ~= 1
        error('tempoVsGrid:RefRow', '%s', sprintf('%s: %d rows in #207 table', DATASETS(d), numel(j)));
    end
    if abs(R.VGF(j) - results.vgfMedianF0(d)) > 1e-9
        error('tempoVsGrid:RefMismatch', '%s', sprintf('%s: median VGF %.6f vs #207 %.6f', ...
            DATASETS(d), results.vgfMedianF0(d), R.VGF(j)));
    end
end

%% Self-check: v002's "Cook CTRL 28% below the floor" is the mean-f0 figure
iC = results.dataset == "Cook_CTRL";
if round(results.pctBelowMean(iC)) ~= 28
    error('tempoVsGrid:V002Check', '%s', sprintf('Cook CTRL mean-f0 figure is %.2f%% below, not 28%%', ...
        results.pctBelowMean(iC)));
end

%% Report
fprintf('K_conv = %.4f (BETA_REF = 1/3); grid VGF %.2f-%.2f = f0 %.3f-%.3f Hz\n', ...
    K_conv, VGF_LO, VGF_HI, f0Lo, f0Hi);
fprintf('Self-checks passed: K_conv; median VGFs = #207 table; Cook CTRL mean = 28%% below.\n\n');
fprintf('%-13s %5s %7s %7s %8s %8s %8s %8s %7s %7s\n', 'dataset', 'n', 'meanF0', 'medF0', ...
    'VGF(mn)', 'VGF(md)', '%bel(mn)', '%bel(md)', '%trBel', '%trAbv');
for d = 1:n
    fprintf('%-13s %5d %7.4f %7.4f %8.1f %8.1f %8.1f %8.1f %7.1f %7.1f\n', results.dataset(d), ...
        results.nTrials(d), results.meanF0(d), results.medianF0(d), results.vgfMeanF0(d), ...
        results.vgfMedianF0(d), results.pctBelowMean(d), results.pctBelowMedian(d), ...
        results.pctTrialsBelow(d), results.pctTrialsAbove(d));
end
fprintf('(%%bel: summary VGF''s distance below the floor, negative = above the floor.)\n');

matOut = fullfile(srcDir, 'results', 'checkEmpiricalTempoVsGridVGF_v001.mat');
save(matOut, 'results', 'K_conv', 'VGF_LO', 'VGF_HI', 'BETA_REF');
fprintf('Saved: %s\n', matOut);
