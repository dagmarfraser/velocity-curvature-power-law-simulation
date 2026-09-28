% checkVGFSEMCentroid_v002.m  VGF SEM/bias at each dataset's OWN tempo,
% not averaged across the whole VGF grid -- gated by the pre-reg's beta
% SEM adequacy condition, same as v001.
%
% SUPERSEDES v001's VGF lookup method. v001 snapped on (alpha, sigma, fs)
% only and averaged VGF bias/SEM across every VGF_gen node in the grid
% (14 values, 90-330 mm/s) -- appropriate for beta (which the paper
% already evaluates this way, see checkFraserSEMCentroid_v001.m /
% knockDownFlags_v001.m), but not for VGF itself: VGF is a monotone
% cipher for orbital frequency f0 AT THE SIMULATION'S OWN FIXED ELLIPSE
% GEOMETRY (D1 desiderata, docs/LMM_VGF_Desiderata_v003.md), which does
% not match any real dataset's own ellipse size. Averaging across the
% full VGF range therefore mixed in tempos the dataset never produced.
%
% v002 instead: (1) derives each dataset's own real f0 straight from its
% noiseCharacterisation_<name>.mat (the SAME source and SAME mean()
% aggregation already used for the canonical alpha/sigma figures --
% confirmed to reproduce docs/EMPIRICAL_DATASETS.md's published values
% exactly, 2026-09-13); (2) converts that f0 to the simulation's own VGF
% axis via K_conv, computed fresh from baselineShp6_120Hz.mat (the exact
% method in plotVGFRecovery_v001.m, not a hardcoded constant); (3) snaps
% the VGF lookup on THAT value, in addition to the existing alpha/sigma/fs
% snap, so the reported VGF bias/SEM is read at the dataset's own tempo,
% not smeared across tempos it never exhibited.
%
% REPRODUCIBILITY: every derived number below (K_conv, each dataset's
% alpha/sigma/VGF-equivalent centroid) is computed from saved project
% files in this run, not copied from a table or a prior session's
% printed numbers. Only fs (sampling rate) and each dataset's
% SigmaToMM conversion factor are literal CONFIG constants below -- both
% are fixed recording-hardware properties, not derived quantities, and
% each SigmaToMM is cited against its own importer file/line (matching
% docs/EMPIRICAL_DATASETS.md's own audit convention) rather than
% invented.
%
% EXTRAPOLATION: the simulation's VGF grid spans exp(4.5)-exp(5.8) =
% 90.0-330.3 mm/s. A dataset whose own tempo maps outside that range has
% no real neighbouring grid node -- the nearest-neighbour snap still
% returns a number (the grid's own floor or ceiling), but it is
% extrapolated, not interpolated, and is flagged EXTRAPOLATED below
% rather than presented as an ordinary lookup.
%
% OUTPUT: console table only. No files written (same convention as
% checkFraserSEMCentroid_v001.m / v001 of this script).
%
% USAGE: checkVGFSEMCentroid_v002
%
% Fraser, D.S. (2026)  v002

clearvars
srcDir = fileparts(mfilename('fullpath'));
if isempty(srcDir), srcDir = pwd; end
cd(srcDir);
addpath(genpath(fullfile(srcDir, 'functions')));

fprintf('=== VGF SEM AT EACH DATASET''S OWN TEMPO, GATED BY BETA SEM ADEQUACY (v002) ===\n\n');

%% ---------------------------------------------------------------------
%% CONFIG: fixed recording-hardware properties only (not derived).
%% SigmaToMM factors cited against their own importer source (see
%% docs/EMPIRICAL_DATASETS.md "Unit conversion reference" for the same
%% audit, done independently there); fs is each study's own sampling rate.
%% ---------------------------------------------------------------------
% {label, noiseCharFile, fs, SigmaToMM, SigmaToMM source}
CONFIG = {
    'Fraser',       'noiseCharacterisation_fraser.mat',       240, 0.095988, 'importDB_fraser_v001.m PX_PER_MM=10.41793';
    'Cook CTRL',    'noiseCharacterisation_cook.mat',         133, 0.248,    'importDB_cook_v002.m line 452';
    'Cook ASD',     'noiseCharacterisation_cookASD.mat',      133, 0.248,    'importDB_cook_v002.m line 452';
    'Hickman PLAC', 'noiseCharacterisation_hickmanPLAC.mat',  133, 0.248,    'importDB_hickman_v003.m line 238';
    'Hickman HALO', 'noiseCharacterisation_hickmanHALO.mat',  133, 0.248,    'importDB_hickman_v003.m line 238';
    'Zarandi',      'noiseCharacterisation_zarandi.mat',      100, 10.0,     'importDB_zarandi_v001.m line 79';
    'Dhieb',        'noiseCharacterisation_dhieb.mat',        100, 0.1478,   'importDB_dhieb_v001.m';
};

%% ---------------------------------------------------------------------
%% Derive K_conv fresh from the simulation's own reference ellipse
%% geometry (identical method to plotVGFRecovery_v001.m -- not copied
%% from that script's own printed value, recomputed from its source mat).
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

fprintf('K_conv (derived fresh from %s) = %.6f\n', shapeFile, K_conv);
fprintf('VGF grid range: [%.3f, %.3f] mm/s\n\n', VGF_LO, VGF_HI);

%% ---------------------------------------------------------------------
%% Derive each dataset's own centroid (alpha, sigma, VGF-equivalent) from
%% its own noiseCharacterisation mat -- NOT hardcoded from a doc.
%% ---------------------------------------------------------------------
nDatasets = size(CONFIG, 1);
centroids = table('Size', [nDatasets 6], ...
    'VariableTypes', {'string','double','double','double','double','logical'}, ...
    'VariableNames', {'label','alpha','sigmaMM','fs','vgfEquiv','extrapolated'});

fprintf('--- Per-dataset centroids, derived from noiseCharacterisation mats ---\n');
fprintf('%-14s %8s %8s %6s %10s %10s %s\n', 'Dataset', 'alpha', 'sigmaMM', 'fs', 'meanF0', 'VGFequiv', 'in-grid?');
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

    sigmaMM = meanSigmaNative * CONFIG{i,4};
    vgfEquiv = meanF0 * K_conv;
    extrapolated = vgfEquiv < VGF_LO || vgfEquiv > VGF_HI;

    centroids.label(i) = CONFIG{i,1};
    centroids.alpha(i) = meanAlpha;
    centroids.sigmaMM(i) = sigmaMM;
    centroids.fs(i) = CONFIG{i,3};
    centroids.vgfEquiv(i) = vgfEquiv;
    centroids.extrapolated(i) = extrapolated;

    fprintf('%-14s %8.4f %8.3f %6d %10.4f %10.3f %s\n', ...
        CONFIG{i,1}, meanAlpha, sigmaMM, CONFIG{i,3}, meanF0, vgfEquiv, string(~extrapolated));
end
fprintf('\n');

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

SEM_ADEQUATE = semAdequacyThreshold_v001();  % MDC/2.77 = 0.0108303, exact

allAlphaB = sort(unique(Tbeta.alpha));
allSigmaB = sort(unique(Tbeta.sigma));
allFsB    = sort(unique(Tbeta.fs));
allVGF    = sort(unique(Tvgf.VGF));

%% ---------------------------------------------------------------------
%% Per-dataset, per-pipeline lookup
%% ---------------------------------------------------------------------
fprintf('--- VGF bias/SEM at each dataset''s OWN snapped (alpha,sigma,fs,VGF) coordinate ---\n');
fprintf('(beta SEM/adequacy still uses the canonical mean-over-ALL-VGF-and-betaGen method, unchanged)\n\n');
fprintf('%-14s %-10s %7s %8s | %10s %8s %8s %8s | %s\n', ...
    'Dataset', 'Pipeline', 'betaSEM', 'adequate', 'VGFbias', 'VGFbias', 'VGFsem', 'VGFsem', 'Included?');
fprintf('%-14s %-10s %7s %8s | %10s %8s %8s %8s | %s\n', ...
    '', '', '', '', '(mm/s)', '(Hz)', '(mm/s)', '(Hz)', '');
fprintf('%s\n', repmat('-', 1, 110));

nIncluded = 0;
nExcluded = 0;

for di = 1:nDatasets
    label = centroids.label(di);
    tA = centroids.alpha(di);
    tS = centroids.sigmaMM(di);
    tF = centroids.fs(di);
    tV = centroids.vgfEquiv(di);
    extrapFlag = centroids.extrapolated(di);

    [~, ai] = min(abs(allAlphaB - tA));
    [~, si] = min(abs(allSigmaB - tS));
    [~, fi] = min(abs(allFsB - tF));
    [~, vi] = min(abs(allVGF - tV));
    sA = allAlphaB(ai); sS = allSigmaB(si); sF = allFsB(fi); sV = allVGF(vi);

    % beta: canonical method, unchanged -- alpha/sigma/fs only, averaged
    % over ALL betaGen and ALL VGF at that noise coordinate.
    subBeta = Tbeta(Tbeta.alpha == sA & Tbeta.sigma == sS & Tbeta.fs == sF, :);
    % VGF: NEW -- alpha/sigma/fs/VGF all four snapped, averaged over
    % betaGen only.
    subVGF = Tvgf(Tvgf.alpha == sA & Tvgf.sigma == sS & Tvgf.fs == sF & Tvgf.VGF == sV, :);

    if isempty(subBeta)
        fprintf('%-14s  WARNING: no beta coordTable rows at snapped (%.3f, %.2f, %d) -- skipped.\n', ...
            label, sA, sS, sF);
        continue
    end

    [Gb, pipNamesB] = findgroups(subBeta.pipeline);
    betaSEMperPipe = splitapply(@(x) mean(x, 'omitnan'), subBeta.sem, Gb);

    for pIdx = 1:numel(pipNamesB)
        pipe = pipNamesB(pIdx);
        betaSEM = betaSEMperPipe(pIdx);
        adequate = betaSEM < SEM_ADEQUATE;

        vgfRow = subVGF(subVGF.pipeline == pipe, :);
        if isempty(vgfRow)
            vgfBiasMM = NaN; vgfSemMM = NaN;
        else
            vgfBiasMM = mean(vgfRow.meanBias, 'omitnan');
            vgfSemMM  = mean(vgfRow.sem, 'omitnan');
        end
        vgfBiasHz = vgfBiasMM / K_conv;
        vgfSemHz  = vgfSemMM / K_conv;

        if adequate
            includedStr = 'INCLUDED';
            nIncluded = nIncluded + 1;
        else
            includedStr = 'EXCLUDED (beta gate)';
            nExcluded = nExcluded + 1;
        end
        if extrapFlag
            includedStr = [includedStr, ' [VGF EXTRAPOLATED]']; %#ok<AGROW>
        end

        fprintf('%-14s %-10s %7.4f %8s | %10.3f %8.3f %8.3f %8.3f | %s\n', ...
            label, pipe, betaSEM, string(adequate), vgfBiasMM, vgfBiasHz, vgfSemMM, vgfSemHz, includedStr);
    end
end

fprintf('\n%d / %d dataset-pipeline cells pass the pre-reg''s beta-SEM gate (Sections 3.2, 4.2).\n', ...
    nIncluded, nIncluded + nExcluded);
extrapDatasets = centroids.label(centroids.extrapolated);
if ~isempty(extrapDatasets)
    fprintf('VGF EXTRAPOLATED for: %s (own tempo maps outside the simulated VGF grid [%.1f, %.1f] mm/s).\n', ...
        strjoin(extrapDatasets, ', '), VGF_LO, VGF_HI);
else
    fprintf('No dataset''s own tempo falls outside the simulated VGF grid.\n');
end
fprintf('Excluded/extrapolated cells are shown, not hidden (Fail Loud).\n');

fprintf('\n=== DONE ===\n');
