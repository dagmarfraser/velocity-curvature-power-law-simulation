% checkVGFSEMCentroid_v003.m  Adds the %-of-own-f0 bias/SEM figures and
% cross-dataset summary that the manuscript (Part 1 Results (c)) and
% Finding #194 actually cite -- v002 printed VGF bias/SEM in mm/s and Hz
% per cell, but never divided by each dataset's own f0 or aggregated
% across pipelines/datasets. Those percentages and the six-pipeline-mean
% summary were computed by hand in conversation, not by any script --
% this closes that gap so every number in Part 1 Results (c) traces to a
% printed line, not to arithmetic done off-script.
%
% Adds two things v002 did not have:
%   1. Per-cell columns: bias as %|own f0|, SEM as % own f0.
%   2. Per-dataset six-pipeline-mean summary (INCLUDED cells only), plus
%      the cross-dataset min-max range of those per-dataset means -- this
%      is what "25-74%" and "1.7-12.1%" in the manuscript/Finding #194
%      actually refer to.
%
% Everything else (centroid derivation, K_conv, snapping, beta gate,
% extrapolation flag) is unchanged from v002 -- see that file's header
% for the full method note; not repeated here.
%
% OUTPUT: console table only. No files written.
%
% USAGE: checkVGFSEMCentroid_v003
%
% Fraser, D.S. (2026)  v003

clearvars
srcDir = fileparts(mfilename('fullpath'));
if isempty(srcDir), srcDir = pwd; end
cd(srcDir);
addpath(genpath(fullfile(srcDir, 'functions')));

fprintf('=== VGF SEM AT EACH DATASET''S OWN TEMPO, WITH %% OF OWN F0 (v003) ===\n\n');

%% ---------------------------------------------------------------------
%% CONFIG: identical to v002
%% ---------------------------------------------------------------------
CONFIG = {
    'Fraser',       'noiseCharacterisation_fraser.mat',       240, 0.095988;
    'Cook CTRL',    'noiseCharacterisation_cook.mat',         133, 0.248;
    'Cook ASD',     'noiseCharacterisation_cookASD.mat',      133, 0.248;
    'Hickman PLAC', 'noiseCharacterisation_hickmanPLAC.mat',  133, 0.248;
    'Hickman HALO', 'noiseCharacterisation_hickmanHALO.mat',  133, 0.248;
    'Zarandi',      'noiseCharacterisation_zarandi.mat',      100, 10.0;
    'Dhieb',        'noiseCharacterisation_dhieb.mat',        100, 0.1478;
};
nDatasets = size(CONFIG, 1);

%% ---------------------------------------------------------------------
%% K_conv, derived fresh (identical method to v002/plotVGFRecovery_v001.m)
%% ---------------------------------------------------------------------
BETA_REF = 1/3;
shapeFile = fullfile(srcDir, 'functions', 'baselineShp6_120Hz.mat');
if ~isfile(shapeFile)
    shapeFile = fullfile(srcDir, 'functions', 'baselineShp6_60Hz.mat');
end
if ~isfile(shapeFile)
    error('checkVGFSEMCentroid:NoShapeFile', '%s', ...
        'baselineShp6_*Hz.mat not found in src/functions/.');
end
S = load(shapeFile, 'pathXYresample', 'K');
xy = S.pathXYresample;
diffs = diff(xy, 1, 1);
perimeter = sum(sqrt(sum(diffs.^2, 2)));
kappa = S.K(:);
kappa = kappa(kappa > 0);
mean_kappa_term = mean(kappa .^ (-BETA_REF));
K_conv = perimeter / mean_kappa_term;

VGF_LO = exp(4.5);
VGF_HI = exp(5.8);

fprintf('K_conv (derived fresh from %s) = %.6f\n\n', shapeFile, K_conv);

%% ---------------------------------------------------------------------
%% Per-dataset centroids (alpha, sigma, f0, VGF-equivalent)
%% ---------------------------------------------------------------------
centroids = table('Size', [nDatasets 6], ...
    'VariableTypes', {'string','double','double','double','double','logical'}, ...
    'VariableNames', {'label','alpha','sigmaMM','fs','meanF0','extrapolated'});

for i = 1:nDatasets
    ncFile = fullfile(srcDir, CONFIG{i,2});
    if ~isfile(ncFile)
        error('checkVGFSEMCentroid:NotFound', '%s', sprintf('FAILED PATH: %s not found.', ncFile));
    end
    D = load(ncFile, 'bioResults');
    T = D.bioResults;

    meanAlpha = mean(T.ira_alphaMean, 'omitnan');
    meanSigmaNative = mean(T.sigmaMean, 'omitnan');
    meanF0 = mean(T.f0, 'omitnan');
    vgfEquiv = meanF0 * K_conv;

    centroids.label(i) = CONFIG{i,1};
    centroids.alpha(i) = meanAlpha;
    centroids.sigmaMM(i) = meanSigmaNative * CONFIG{i,4};
    centroids.fs(i) = CONFIG{i,3};
    centroids.meanF0(i) = meanF0;
    centroids.extrapolated(i) = vgfEquiv < VGF_LO || vgfEquiv > VGF_HI;
end

%% ---------------------------------------------------------------------
%% Load beta and VGF per-coordinate tables
%% ---------------------------------------------------------------------
betaFile = fullfile(srcDir, 'perCoordinateSEM_v2_001.mat');
if ~isfile(betaFile)
    error('checkVGFSEMCentroid:notFound', '%s', sprintf('FAILED PATH: %s not found.', betaFile));
end
vgfFile = fullfile(srcDir, 'perCoordinateSEM_VGF_v001.mat');
if ~isfile(vgfFile)
    error('checkVGFSEMCentroid:notFound', '%s', sprintf(...
        'FAILED PATH: %s not found. Run computePerCoordinateSEM_VGF_v001() first.', vgfFile));
end

Tbeta = load(betaFile, 'coordTable'); Tbeta = Tbeta.coordTable;
Tvgf  = load(vgfFile,  'coordTable'); Tvgf  = Tvgf.coordTable;

SEM_ADEQUATE = semAdequacyThreshold_v001();

allAlphaB = sort(unique(Tbeta.alpha));
allSigmaB = sort(unique(Tbeta.sigma));
allFsB    = sort(unique(Tbeta.fs));
allVGF    = sort(unique(Tvgf.VGF));

%% ---------------------------------------------------------------------
%% Per-dataset, per-pipeline lookup, now with %-of-own-f0 columns
%% ---------------------------------------------------------------------
fprintf('%-14s %-10s %8s | %8s %8s | %8s %8s | %s\n', ...
    'Dataset', 'Pipeline', 'adequate', 'bias(Hz)', 'bias(%f0)', 'sem(Hz)', 'sem(%f0)', 'Included?');
fprintf('%s\n', repmat('-', 1, 90));

% Per-dataset accumulators for the six-pipeline-mean summary
biasPctByDataset = cell(nDatasets, 1);
semPctByDataset  = cell(nDatasets, 1);

for di = 1:nDatasets
    label = centroids.label(di);
    tA = centroids.alpha(di);
    tS = centroids.sigmaMM(di);
    tF = centroids.fs(di);
    f0 = centroids.meanF0(di);
    tV = f0 * K_conv;
    extrapFlag = centroids.extrapolated(di);

    [~, ai] = min(abs(allAlphaB - tA));
    [~, si] = min(abs(allSigmaB - tS));
    [~, fi] = min(abs(allFsB - tF));
    [~, vi] = min(abs(allVGF - tV));
    sA = allAlphaB(ai); sS = allSigmaB(si); sF = allFsB(fi); sV = allVGF(vi);

    subBeta = Tbeta(Tbeta.alpha == sA & Tbeta.sigma == sS & Tbeta.fs == sF, :);
    subVGF  = Tvgf(Tvgf.alpha == sA & Tvgf.sigma == sS & Tvgf.fs == sF & Tvgf.VGF == sV, :);

    if isempty(subBeta)
        fprintf('%-14s  WARNING: no beta rows at snapped coordinate -- skipped.\n', label);
        continue
    end

    [Gb, pipNamesB] = findgroups(subBeta.pipeline);
    betaSEMperPipe = splitapply(@(x) mean(x, 'omitnan'), subBeta.sem, Gb);

    biasPctThisDataset = [];
    semPctThisDataset  = [];

    for pIdx = 1:numel(pipNamesB)
        pipe = pipNamesB(pIdx);
        betaSEM = betaSEMperPipe(pIdx);
        adequate = betaSEM < SEM_ADEQUATE;

        vgfRow = subVGF(subVGF.pipeline == pipe, :);
        if isempty(vgfRow)
            biasHz = NaN; semHz = NaN;
        else
            biasHz = mean(vgfRow.meanBias, 'omitnan') / K_conv;
            semHz  = mean(vgfRow.sem, 'omitnan') / K_conv;
        end
        biasPct = 100 * abs(biasHz) / f0;
        semPct  = 100 * semHz / f0;

        if adequate
            includedStr = 'INCLUDED';
            biasPctThisDataset(end+1) = biasPct; %#ok<AGROW>
            semPctThisDataset(end+1)  = semPct;  %#ok<AGROW>
        else
            includedStr = 'EXCLUDED (beta gate)';
        end
        if extrapFlag
            includedStr = [includedStr, ' [VGF EXTRAPOLATED]']; %#ok<AGROW>
        end

        fprintf('%-14s %-10s %8s | %8.3f %8.1f | %8.3f %8.1f | %s\n', ...
            label, pipe, string(adequate), biasHz, biasPct, semHz, semPct, includedStr);
    end

    biasPctByDataset{di} = biasPctThisDataset;
    semPctByDataset{di}  = semPctThisDataset;
end

%% ---------------------------------------------------------------------
%% Per-dataset six-pipeline-mean summary, and cross-dataset range
%% ---------------------------------------------------------------------
fprintf('\n--- Six-pipeline-mean summary, INCLUDED cells only ---\n');
fprintf('%-14s %12s %12s\n', 'Dataset', 'mean bias%f0', 'mean sem%f0');
datasetMeanBias = nan(nDatasets, 1);
datasetMeanSem  = nan(nDatasets, 1);
for di = 1:nDatasets
    if isempty(biasPctByDataset{di})
        fprintf('%-14s %12s %12s  (no included cells)\n', centroids.label(di), '--', '--');
        continue
    end
    datasetMeanBias(di) = mean(biasPctByDataset{di});
    datasetMeanSem(di)  = mean(semPctByDataset{di});
    fprintf('%-14s %12.1f %12.1f\n', centroids.label(di), datasetMeanBias(di), datasetMeanSem(di));
end

validBias = datasetMeanBias(~isnan(datasetMeanBias));
validSem  = datasetMeanSem(~isnan(datasetMeanSem));
fprintf('\nCross-dataset range (six-pipeline means, INCLUDED datasets only):\n');
fprintf('  bias: %.1f%%-%.1f%%\n', min(validBias), max(validBias));
fprintf('  sem:  %.1f%%-%.1f%%\n', min(validSem), max(validSem));

fprintf('\n=== DONE ===\n');
