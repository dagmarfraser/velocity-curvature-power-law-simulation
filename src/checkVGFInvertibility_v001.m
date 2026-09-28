% checkVGFInvertibility_v001.m  Is the VGF_gen -> VGF_rec forward map
% monotonic (hence invertible) at each dataset's own noise centroid?
%
% Mirrors the logic of this project's existing beta invertibility gate
% (checkInvertibilityForEmpirical_v2_002.m), applied to VGF instead of
% beta, using perCoordinateSEM_VGF_v001.m's own meanVGFrec column.
%
% At each dataset's snapped (alpha, sigma, fs) noise coordinate, and for
% each of the 6 pipelines, this sweeps the 14 VGF_gen grid nodes
% (averaging over all betaGen at that coordinate, same convention as the
% rest of this project's centroid lookups), and checks whether
% mean(VGF_rec) rises monotonically with VGF_gen. A non-monotonic map
% means the same recovered VGF_rec could plausibly have come from more
% than one true VGF_gen -- i.e. VGF is not safely invertible there, on
% top of (and independent from) the systematic bias already found.
%
% REPRODUCIBILITY: dataset centroids (alpha, sigma) are derived fresh
% from each noiseCharacterisation_<name>.mat, same method as
% checkVGFSEMCentroid_v002.m -- not copied from that script's own
% printed output.
%
% OUTPUT: console table only. No files written.
%
% USAGE: checkVGFInvertibility_v001
%
% Fraser, D.S. (2026)  v001

clearvars
srcDir = fileparts(mfilename('fullpath'));
if isempty(srcDir), srcDir = pwd; end
cd(srcDir);
addpath(genpath(fullfile(srcDir, 'functions')));

fprintf('=== VGF FORWARD-MAP MONOTONICITY AT EACH DATASET''S OWN NOISE CENTROID ===\n\n');

%% ---------------------------------------------------------------------
%% CONFIG: same fixed recording-hardware properties as checkVGFSEMCentroid_v002.m
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
%% Derive each dataset's own (alpha, sigma) centroid, fresh from source
%% ---------------------------------------------------------------------
alphaC = nan(nDatasets,1); sigmaC = nan(nDatasets,1); fsC = nan(nDatasets,1);
for i = 1:nDatasets
    ncFile = fullfile(srcDir, CONFIG{i,2});
    if ~isfile(ncFile)
        error('checkVGFInvertibility:NotFound', '%s', sprintf('FAILED PATH: %s not found.', ncFile));
    end
    D = load(ncFile, 'bioResults');
    T = D.bioResults;
    alphaC(i) = mean(T.ira_alphaMean, 'omitnan');
    sigmaC(i) = mean(T.sigmaMean, 'omitnan') * CONFIG{i,4};
    fsC(i)    = CONFIG{i,3};
end

%% ---------------------------------------------------------------------
%% Load VGF per-coordinate table
%% ---------------------------------------------------------------------
vgfFile = fullfile(srcDir, 'perCoordinateSEM_VGF_v001.mat');
if ~isfile(vgfFile)
    error('checkVGFInvertibility:NotFound', '%s', sprintf('FAILED PATH: %s not found.', vgfFile));
end
Tvgf = load(vgfFile, 'coordTable'); Tvgf = Tvgf.coordTable;

allAlpha = sort(unique(Tvgf.alpha));
allSigma = sort(unique(Tvgf.sigma));
allFs    = sort(unique(Tvgf.fs));
allVGF   = sort(unique(Tvgf.VGF));
pipes    = categories(Tvgf.pipeline);

fprintf('%-14s %-10s %10s %8s %s\n', 'Dataset', 'Pipeline', 'nNodes', 'nMono', 'Monotonic?');
fprintf('%s\n', repmat('-', 1, 60));

for di = 1:nDatasets
    label = CONFIG{di,1};
    [~, ai] = min(abs(allAlpha - alphaC(di)));
    [~, si] = min(abs(allSigma - sigmaC(di)));
    [~, fi] = min(abs(allFs - fsC(di)));
    sA = allAlpha(ai); sS = allSigma(si); sF = allFs(fi);

    sub = Tvgf(Tvgf.alpha == sA & Tvgf.sigma == sS & Tvgf.fs == sF, :);
    if isempty(sub)
        fprintf('%-14s  WARNING: no rows at snapped (%.3f, %.2f, %d) -- skipped.\n', label, sA, sS, sF);
        continue
    end

    for pIdx = 1:numel(pipes)
        pipe = pipes{pIdx};
        subP = sub(sub.pipeline == pipe, :);
        % Average over betaGen at each VGF node, then sort by VGF ascending
        [G, vgfVals] = findgroups(subP.VGF);
        meanRec = splitapply(@(x) mean(x, 'omitnan'), subP.meanVGFrec, G);
        [vgfSorted, order] = sort(vgfVals);
        recSorted = meanRec(order);

        steps = diff(recSorted);
        nSteps = numel(steps);
        nMono = sum(steps >= 0);
        isMono = all(steps >= 0);

        fprintf('%-14s %-10s %10d %8d %s\n', label, pipe, numel(vgfSorted), nMono, string(isMono));
        if ~isMono
            badIdx = find(steps < 0);
            fprintf('    non-monotonic step(s) at VGF=%.1f->%.1f (VGF_rec %.2f->%.2f), etc. (%d/%d steps bad)\n', ...
                vgfSorted(badIdx(1)), vgfSorted(badIdx(1)+1), recSorted(badIdx(1)), recSorted(badIdx(1)+1), ...
                numel(badIdx), nSteps);
        end
    end
end

fprintf('\n=== DONE ===\n');
