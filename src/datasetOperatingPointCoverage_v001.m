function results = datasetOperatingPointCoverage_v001(opts)
% datasetOperatingPointCoverage_v001  Nearest-neighbour coverage of the
% seven empirical datasets in the simulation grid's own (alpha, sigma, fs,
% VGF) operating-point space.
%
% Sources Finding #207 (docs/FINDINGS_REFERENCE.md): explains Finding
% #128's Fraser Part-B SEM outlier as a validation-architecture coverage
% gap rather than a dataset-specific defect. Four datasets (Cook CTRL/ASD,
% Hickman PLAC/HALO) cluster tightly enough to mutually corroborate one
% another in this space; the other three (Dhieb, Fraser, Zarandi) each sit
% relatively isolated, and each is independently the one dataset in the
% corpus showing anomalous validation behaviour elsewhere (Dhieb: EPP
% exclusion D2/D13; Zarandi: Scenario-A saturation; Fraser: Finding #128).
%
% METHOD: for each dataset, take the median per-trial (alpha, sigma, fs,
% VGF) -- alpha via ira_alphaMean (bioResults), sigma via template sigma
% times each dataset's own sigmaToMM, fs from the raw trial import, VGF
% via f0 x K_conv (same conversion as queryMonotonicSlopeAtEmpirical_v003.m).
% Normalise each dimension by monotonicSegments_v2_002.mat's own grid
% range, then take pairwise Euclidean distance and each dataset's nearest
% neighbour.
%
% USAGE:
%   R = datasetOperatingPointCoverage_v001()
%   R = datasetOperatingPointCoverage_v001(Save=false)
%
% Sanity check: asserts (as a warning, not a hard failure -- the numbers
% are still returned) that Fraser's NN distance reproduces Finding #207's
% documented 0.6640 (tolerance 1e-3). If the underlying mats have changed
% since 2026-09-17, this fails loud rather than silently reporting a
% different number under the same finding's name.
%
% Fraser, D.S. (2026)
% See also: queryMonotonicSlopeAtEmpirical_v003 (source of the alpha/sigma
% methodology this reuses), zarandiLMLSFoldDiagnostic_v001 (Finding #206,
% the companion diagnostic from the same session)

    arguments
        opts.Save (1,1) logical = true
    end

    srcDir = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(srcDir, 'functions')));

    %% K_conv (identical derivation to queryMonotonicSlopeAtEmpirical_v003)
    BETA_REF = 1/3;
    shapeFile = fullfile(srcDir, 'functions', 'baselineShp6_120Hz.mat');
    if ~isfile(shapeFile), shapeFile = fullfile(srcDir, 'functions', 'baselineShp6_60Hz.mat'); end
    if ~isfile(shapeFile)
        error('datasetOpCoverage:NoShapeFile', '%s', 'baselineShp6_*Hz.mat not found in src/functions/.');
    end
    Sshp = load(shapeFile, 'pathXYresample', 'K');
    perimeter_px = sum(sqrt(sum(diff(Sshp.pathXYresample, 1, 1).^2, 2)));
    kappa = Sshp.K(:); kappa = kappa(kappa > 0 & isfinite(kappa));
    K_conv = perimeter_px / mean(kappa .^ (-BETA_REF));

    monoFile = fullfile(srcDir, 'monotonicSegments_v2_002.mat');
    if ~isfile(monoFile)
        error('datasetOpCoverage:NoMonoBundle', '%s', 'monotonicSegments_v2_002.mat not found.');
    end
    mono = load(monoFile, 'alphaGrid', 'sigmaGrid', 'fsGrid', 'VGFGrid');

    %% Registry -- identical to queryMonotonicSlopeAtEmpirical_v003's own
    AlphaSource = "ira_alphaMean";
    registry = { ...
        struct('name', "Fraser", 'constellationMat', fullfile(srcDir, "constellationFraser_v001.mat"), 'noiseMat', fullfile(srcDir, "noiseCharacterisation_fraser.mat"), 'importerHandle', @importDB_fraser_v001, 'importerArgs', {{}}, 'sigmaToMM', 1/10.41793); ...
        struct('name', "Zarandi", 'constellationMat', fullfile(srcDir, "constellationZarandi_v001.mat"), 'noiseMat', fullfile(srcDir, "noiseCharacterisation_zarandi.mat"), 'importerHandle', @importDB_zarandi_v001, 'importerArgs', {{ "Verbose", false }}, 'sigmaToMM', 10.0); ...
        struct('name', "Cook_CTRL", 'constellationMat', fullfile(srcDir, "constellationCook_v001.mat"), 'noiseMat', fullfile(srcDir, "noiseCharacterisation_cook.mat"), 'importerHandle', @importDB_cook_v002, 'importerArgs', {{ "Group", "CTRL", "Tasks", 7, "Verbose", false }}, 'sigmaToMM', 0.248); ...
        struct('name', "Cook_ASD", 'constellationMat', fullfile(srcDir, "constellationCookASD_v001.mat"), 'noiseMat', fullfile(srcDir, "noiseCharacterisation_cookASD.mat"), 'importerHandle', @importDB_cook_v002, 'importerArgs', {{ "Group", "ASD", "Tasks", 7, "Verbose", false }}, 'sigmaToMM', 0.248); ...
        struct('name', "Dhieb", 'constellationMat', fullfile(srcDir, "constellationDhieb_v001.mat"), 'noiseMat', fullfile(srcDir, "noiseCharacterisation_dhieb.mat"), 'importerHandle', @importDB_dhieb_v001, 'importerArgs', {{}}, 'sigmaToMM', 0.1478); ...
        struct('name', "Hickman_PLAC", 'constellationMat', fullfile(srcDir, "constellationHickmanPLAC_v001.mat"), 'noiseMat', fullfile(srcDir, "noiseCharacterisation_hickmanPLAC.mat"), 'importerHandle', @importDB_hickman_v002, 'importerArgs', {{ "Study", 2, "Group", "PLAC", "Shapes", 3, "Verbose", false }}, 'sigmaToMM', 0.248); ...
        struct('name', "Hickman_HALO", 'constellationMat', fullfile(srcDir, "constellationHickmanHALO_v001.mat"), 'noiseMat', fullfile(srcDir, "noiseCharacterisation_hickmanHALO.mat"), 'importerHandle', @importDB_hickman_v002, 'importerArgs', {{ "Study", 2, "Group", "HALO", "Shapes", 3, "Verbose", false }}, 'sigmaToMM', 0.248); ...
        };
    nSets = numel(registry);

    %% Per-dataset median (alpha, sigma, fs, VGF)
    names = strings(nSets,1);
    alphaMed = nan(nSets,1); sigmaMed = nan(nSets,1); fsMed = nan(nSets,1); vgfMed = nan(nSets,1);

    for k = 1:nSets
        spec = registry{k};
        if ~isfile(spec.constellationMat) || ~isfile(spec.noiseMat)
            error('datasetOpCoverage:MissingFile', '%s', ...
                sprintf('Constellation or noise mat missing for %s.', spec.name));
        end
        nz = load(spec.noiseMat, 'bioResults');
        bio = nz.bioResults;
        trials = spec.importerHandle(spec.importerArgs{:});
        if numel(trials) ~= height(bio)
            error('datasetOpCoverage:LengthMismatch', '%s', ...
                sprintf('%s: importer returned %d trials but bioResults has %d rows.', ...
                spec.name, numel(trials), height(bio)));
        end
        fsPerTri = arrayfun(@(t) double(t.fs), trials).';
        alphaEmp = double(bio.(AlphaSource));
        sigmaEmpMM = double(bio.sigmaMean) * spec.sigmaToMM;
        vgfEmp = double(bio.f0) * K_conv;

        names(k) = spec.name;
        alphaMed(k) = median(alphaEmp, 'omitnan');
        sigmaMed(k) = median(sigmaEmpMM, 'omitnan');
        fsMed(k) = median(fsPerTri, 'omitnan');
        vgfMed(k) = median(vgfEmp, 'omitnan');
    end

    %% Normalise by grid range, pairwise distance, nearest neighbour
    rngA = range(mono.alphaGrid); rngS = range(mono.sigmaGrid);
    rngF = range(mono.fsGrid);    rngV = range(mono.VGFGrid);
    Xn = [ (alphaMed - min(mono.alphaGrid)) / rngA, ...
           (sigmaMed - min(mono.sigmaGrid)) / rngS, ...
           (fsMed    - min(mono.fsGrid))    / rngF, ...
           (vgfMed   - min(mono.VGFGrid))   / rngV ];

    Dmat = squareform(pdist(Xn));
    Dmat(Dmat == 0) = NaN;
    [nnDist, nnIdx] = min(Dmat, [], 2, 'omitnan');

    results = table(names, alphaMed, sigmaMed, fsMed, vgfMed, nnDist, names(nnIdx), ...
        'VariableNames', {'dataset','alpha','sigmaMM','fs','VGF','nnDist','nearestNeighbour'});
    results = sortrows(results, 'nnDist', 'descend');

    fprintf('=== datasetOperatingPointCoverage_v001 (Finding #207) ===\n');
    fprintf('%-14s %8s %8s %6s %8s %10s %14s\n', 'dataset','alpha','sigmaMM','fs','VGF','nnDist','nearest');
    for r = 1:height(results)
        fprintf('%-14s %8.3f %8.3f %6.0f %8.1f %10.4f %14s\n', results.dataset(r), results.alpha(r), ...
            results.sigmaMM(r), results.fs(r), results.VGF(r), results.nnDist(r), results.nearestNeighbour(r));
    end

    %% Sanity check against Finding #207's documented figure
    fraserRow = results(results.dataset == "Fraser", :);
    fraserNN = fraserRow.nnDist;
    if abs(fraserNN - 0.6640) > 1e-3
        warning('datasetOpCoverage:SanityCheckFailed', ...
            'Fraser NN distance %.4f does not match Finding #207''s documented 0.6640 (tol 1e-3). Source mats may have changed since 2026-09-17 -- do not cite this run under Finding #207 without checking why.', fraserNN);
    else
        fprintf('\nSanity check PASSED: Fraser NN distance %.4f matches Finding #207 (0.6640) within tolerance.\n', fraserNN);
    end

    %% Save
    if opts.Save
        resDir = fullfile(srcDir, 'results');
        if ~exist(resDir, 'dir'), mkdir(resDir); end
        matOut = fullfile(resDir, 'datasetOperatingPointCoverage_v001.mat');
        save(matOut, 'results', 'Xn', 'Dmat', 'K_conv', 'monoFile', '-v7.3');
        fprintf('\nSaved: %s\n', matOut);
    end
end
