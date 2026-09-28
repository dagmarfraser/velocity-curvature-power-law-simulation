function results = vgfMatchedSEMSensitivity_v001(opts)
% vgfMatchedSEMSensitivity_v001  Does the published 36-of-42 SEM-adequacy
% classification (Fig 0) survive evaluating each dataset at its own real
% VGF/tempo instead of averaged across the grid's full VGF range?
%
% Prompted by GPT SOL's revised critique (2026-09-17): "the corrected
% Figure 0 lookup averages across VGF levels... If SEM varies materially
% with VGF, the 36-of-42 result is not literally evaluated at each full
% empirical operating point." Sources Finding #208.
%
% METHOD: canonical (alpha, sigma, fs) centroids per dataset
% (docs/EMPIRICAL_DATASETS.md v004 table); "existing method" SEM = mean
% of perCoordinateSEM_v2_001.mat's sem column over ALL betaGen and ALL
% VGF grid nodes at the snapped (alpha,sigma,fs,pipeline) coordinate,
% replicating knockDownFlags_v001.m's own lookup exactly. "VGF-matched"
% SEM additionally snaps to the dataset's own real VGF (f0 x K_conv,
% computed fresh via the same registry as
% datasetOperatingPointCoverage_v001.m) before averaging over betaGen
% only. Adequacy threshold 0.0108 (MDC/2.77), matching Finding #128/R5.
%
% USAGE:
%   R = vgfMatchedSEMSensitivity_v001()
%
% Sanity check: asserts the existing-method column reproduces the
% published 36-of-42 split exactly before trusting the VGF-matched
% column at all.
%
% Fraser, D.S. (2026)
% See also: Finding #206 (Zarandi/LMLS fold proximity -- same cells flip
% here), Finding #207/#208, knockDownFlags_v001.m (the lookup this
% extends)

    arguments
        opts.Save (1,1) logical = true
    end

    srcDir = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(srcDir, 'functions')));

    semFile = fullfile(srcDir, 'perCoordinateSEM_v2_001.mat');
    if ~isfile(semFile)
        error('vgfSEMSens:NoSemFile', '%s', 'perCoordinateSEM_v2_001.mat not found.');
    end
    S = load(semFile, 'coordTable');
    T = S.coordTable;
    allAlpha = sort(unique(T.alpha)); allSigma = sort(unique(T.sigma));
    allFs = sort(unique(T.fs)); allVGF = sort(unique(T.VGF));
    pipes = unique(T.pipeline);
    THRESH = 0.0108;

    % Canonical (alpha, sigma, fs) centroids -- docs/EMPIRICAL_DATASETS.md
    % v004 table, the exact coordinates behind the published 36-of-42.
    names  = ["Fraser","Cook_CTRL","Cook_ASD","Hickman_PLAC","Hickman_HALO","Zarandi","Dhieb"];
    alphaC = [4.289, 4.771, 5.062, 5.337, 5.424, 3.184, 2.497];
    sigmaC = [2.009, 8.15,  7.84,  7.17,  7.42,  4.77,  7.50];
    fsC    = [240,   133,   133,   133,   133,   100,   100];

    % Each dataset's own real VGF (f0 x K_conv) -- same registry/importers
    % as datasetOperatingPointCoverage_v001.m, recomputed fresh here.
    BETA_REF = 1/3;
    shapeFile = fullfile(srcDir, 'functions', 'baselineShp6_120Hz.mat');
    if ~isfile(shapeFile), shapeFile = fullfile(srcDir, 'functions', 'baselineShp6_60Hz.mat'); end
    Sshp = load(shapeFile, 'pathXYresample', 'K');
    perimeter_px = sum(sqrt(sum(diff(Sshp.pathXYresample, 1, 1).^2, 2)));
    kappa = Sshp.K(:); kappa = kappa(kappa > 0 & isfinite(kappa));
    K_conv = perimeter_px / mean(kappa .^ (-BETA_REF));

    registry = { ...
        struct('name', "Fraser", 'noiseMat', fullfile(srcDir, "noiseCharacterisation_fraser.mat")); ...
        struct('name', "Cook_CTRL", 'noiseMat', fullfile(srcDir, "noiseCharacterisation_cook.mat")); ...
        struct('name', "Cook_ASD", 'noiseMat', fullfile(srcDir, "noiseCharacterisation_cookASD.mat")); ...
        struct('name', "Hickman_PLAC", 'noiseMat', fullfile(srcDir, "noiseCharacterisation_hickmanPLAC.mat")); ...
        struct('name', "Hickman_HALO", 'noiseMat', fullfile(srcDir, "noiseCharacterisation_hickmanHALO.mat")); ...
        struct('name', "Zarandi", 'noiseMat', fullfile(srcDir, "noiseCharacterisation_zarandi.mat")); ...
        struct('name', "Dhieb", 'noiseMat', fullfile(srcDir, "noiseCharacterisation_dhieb.mat")); ...
        };
    vgfOwn = nan(numel(names),1);
    for k = 1:numel(registry)
        nz = load(registry{k}.noiseMat, 'bioResults');
        vgfOwn(k) = median(double(nz.bioResults.f0), 'omitnan') * K_conv;
    end

    rowsOut = {};
    for d = 1:numel(names)
        [~,ai] = min(abs(allAlpha-alphaC(d)));
        [~,si] = min(abs(allSigma-sigmaC(d)));
        [~,fi] = min(abs(allFs-fsC(d)));
        [~,vi] = min(abs(allVGF-vgfOwn(d)));
        sA=allAlpha(ai); sS=allSigma(si); sF=allFs(fi); sV=allVGF(vi);
        for p = 1:numel(pipes)
            subAll = T(T.alpha==sA & T.sigma==sS & T.fs==sF & T.pipeline==pipes(p), :);
            semAvg = mean(subAll.sem, 'omitnan');
            subVGF = subAll(subAll.VGF==sV, :);
            semVGFmatch = mean(subVGF.sem, 'omitnan');
            adeqA = semAvg < THRESH;
            adeqB = semVGFmatch < THRESH;
            rowsOut(end+1,:) = {names(d), pipes(p), semAvg, semVGFmatch, adeqA, adeqB, sV, adeqA~=adeqB}; %#ok<AGROW>
        end
    end

    results = cell2table(rowsOut, 'VariableNames', ...
        {'dataset','pipeline','semAvgVGF','semMatchedVGF','adequateAvg','adequateMatched','snappedVGF','flip'});

    fprintf('=== vgfMatchedSEMSensitivity_v001 (Finding #208) ===\n');
    fprintf('%-14s %-10s %10s %10s %6s %6s\n','dataset','pipeline','SEM(avg)','SEM(VGFm)','adeqA','adeqB');
    for i = 1:height(results)
        marker = ''; if results.flip(i), marker = '  <== FLIP'; end
        fprintf('%-14s %-10s %10.5f %10.5f %6d %6d%s\n', results.dataset(i), results.pipeline(i), ...
            results.semAvgVGF(i), results.semMatchedVGF(i), results.adequateAvg(i), results.adequateMatched(i), marker);
    end
    fprintf('\nTotal cells: %d, adequate (existing method): %d, flips under VGF-matching: %d\n', ...
        height(results), sum(results.adequateAvg), sum(results.flip));

    %% Sanity check against the published figure
    nAdeqExisting = sum(results.adequateAvg);
    if nAdeqExisting ~= 36
        warning('vgfSEMSens:SanityCheckFailed', ...
            'Existing-method adequate count %d does not match the published 36-of-42. Do not trust the VGF-matched column until this is resolved.', nAdeqExisting);
    else
        fprintf('\nSanity check PASSED: existing-method column reproduces the published 36-of-42 exactly.\n');
    end

    %% Save
    if opts.Save
        resDir = fullfile(srcDir, 'results');
        if ~exist(resDir, 'dir'), mkdir(resDir); end
        matOut = fullfile(resDir, 'vgfMatchedSEMSensitivity_v001.mat');
        save(matOut, 'results', 'K_conv', '-v7.3');
        fprintf('\nSaved: %s\n', matOut);
    end
end
