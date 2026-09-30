function results = datasetOperatingPointCoverage_v003(opts)
% datasetOperatingPointCoverage_v003  Nearest-neighbour coverage of the seven empirical
% datasets in the grid's own (alpha, sigma, fs, VGF) operating-point space (Finding #207).
%
% v003 (2026-09-30) supersedes v002 only in the out-of-domain rule: Zarandi and Dhieb now use their own
% median betaGenStarMed too (flagged inDomain = false), not the 1/3 reference; in-domain results are identical.
% v001 placed each dataset on the VGF axis as median f0 x K_conv, K_conv = 175.2636 at
% beta = 1/3. That constant is 5.28% below the grid's true tempo constant K(1/3) = 185.02
% (testGridKConv_v001), and K falls steeply with beta. v002 runs three variants on
% identical inputs so each change is attributable:
%   V0 RETIRED   median f0 x 175.2636         regression anchor: must reproduce Finding #207's
%                                             Fraser NN distance 0.6640 (tolerance 1e-3).
%   V1 CONSTANT  median f0 x K(1/3)           the constant fixed only.
%   V2 PRIMARY   median f0 (all trials) x K(beta_ds), beta_ds from datasetOwnTempoExponent_v002:
%                the dataset's own median betaGenStarMed for all seven (out-of-domain flagged, not replaced).
% Method otherwise unchanged from v001: per-dataset median alpha (ira_alphaMean), sigma (template
% sigma x sigmaToMM), fs (importer), normalised by monotonicSegments_v2_002's own grid ranges,
% pairwise Euclidean distance, nearest neighbour. The output reports whether any dataset's
% nearest neighbour, or the ordering of the isolation distances, changes between variants.
%
% Writes results/datasetOperatingPointCoverage_v003.mat (the project's results/, NOT src/results/
% where v001 wrote). USAGE: from the project root: R = datasetOperatingPointCoverage_v003()
%                            R = datasetOperatingPointCoverage_v003(Save=false)
% Fraser, D.S. (2026)  v003
    arguments
        opts.Save (1,1) logical = true
    end
    ROOT   = fileparts(fileparts(mfilename('fullpath')));
    srcDir = fullfile(ROOT, 'src');
    addpath(ROOT); addpath(srcDir); addpath(genpath(fullfile(srcDir, 'functions')));

    KCONV_RETIRED = 175.2636;  K13_PUB = 185.0242;  FRASER_NN_PUB = 0.6640;
    K13 = gridKConv_v001(1/3);
    if abs(K13 - K13_PUB) > 1e-3
        error('datasetOpCoverage2:K13', '%s', sprintf('K(1/3) %.4f does not reproduce %.4f', K13, K13_PUB));
    end

    monoFile = fullfile(srcDir, 'monotonicSegments_v2_002.mat');
    if ~isfile(monoFile), error('datasetOpCoverage2:NoMonoBundle', '%s', 'monotonicSegments_v2_002.mat not found.'); end
    mono = load(monoFile, 'alphaGrid', 'sigmaGrid', 'fsGrid', 'VGFGrid');

    %% Registry -- identical to v001 (and queryMonotonicSlopeAtEmpirical_v003)
    AlphaSource = "ira_alphaMean";
    registry = { ...
        struct('name', "Fraser", 'noiseMat', fullfile(srcDir, "noiseCharacterisation_fraser.mat"), 'importerHandle', @importDB_fraser_v001, 'importerArgs', {{}}, 'sigmaToMM', 1/10.41793); ...
        struct('name', "Zarandi", 'noiseMat', fullfile(srcDir, "noiseCharacterisation_zarandi.mat"), 'importerHandle', @importDB_zarandi_v001, 'importerArgs', {{ "Verbose", false }}, 'sigmaToMM', 10.0); ...
        struct('name', "Cook_CTRL", 'noiseMat', fullfile(srcDir, "noiseCharacterisation_cook.mat"), 'importerHandle', @importDB_cook_v002, 'importerArgs', {{ "Group", "CTRL", "Tasks", 7, "Verbose", false }}, 'sigmaToMM', 0.248); ...
        struct('name', "Cook_ASD", 'noiseMat', fullfile(srcDir, "noiseCharacterisation_cookASD.mat"), 'importerHandle', @importDB_cook_v002, 'importerArgs', {{ "Group", "ASD", "Tasks", 7, "Verbose", false }}, 'sigmaToMM', 0.248); ...
        struct('name', "Dhieb", 'noiseMat', fullfile(srcDir, "noiseCharacterisation_dhieb.mat"), 'importerHandle', @importDB_dhieb_v001, 'importerArgs', {{}}, 'sigmaToMM', 0.1478); ...
        struct('name', "Hickman_PLAC", 'noiseMat', fullfile(srcDir, "noiseCharacterisation_hickmanPLAC.mat"), 'importerHandle', @importDB_hickman_v002, 'importerArgs', {{ "Study", 2, "Group", "PLAC", "Shapes", 3, "Verbose", false }}, 'sigmaToMM', 0.248); ...
        struct('name', "Hickman_HALO", 'noiseMat', fullfile(srcDir, "noiseCharacterisation_hickmanHALO.mat"), 'importerHandle', @importDB_hickman_v002, 'importerArgs', {{ "Study", 2, "Group", "HALO", "Shapes", 3, "Verbose", false }}, 'sigmaToMM', 0.248); ...
        };
    nSets = numel(registry);

    %% Per-dataset medians, and each dataset's own tempo statistics
    names = strings(nSets,1);
    alphaMed = nan(nSets,1); sigmaMed = nan(nSets,1); fsMed = nan(nSets,1); f0MedBio = nan(nSets,1);
    betaDs = nan(nSets,1); f0Med15 = nan(nSets,1); inDom = false(nSets,1); rule = strings(nSets,1);
    for k = 1:nSets
        spec = registry{k};
        if ~isfile(spec.noiseMat), error('datasetOpCoverage2:MissingFile', '%s', 'FAILED PATH: ' + spec.noiseMat); end
        bio = load(spec.noiseMat, 'bioResults').bioResults;
        trials = spec.importerHandle(spec.importerArgs{:});
        if numel(trials) ~= height(bio)
            error('datasetOpCoverage2:LengthMismatch', '%s', sprintf('%s: importer returned %d trials but bioResults has %d rows.', ...
                spec.name, numel(trials), height(bio)));
        end
        names(k)    = spec.name;
        alphaMed(k) = median(double(bio.(AlphaSource)), 'omitnan');
        sigmaMed(k) = median(double(bio.sigmaMean) * spec.sigmaToMM, 'omitnan');
        fsMed(k)    = median(arrayfun(@(t) double(t.fs), trials), 'omitnan');
        f0MedBio(k) = median(double(bio.f0), 'omitnan');
        D = datasetOwnTempoExponent_v002(spec.name);
        betaDs(k) = D.betaDs;  f0Med15(k) = D.f0Median;  inDom(k) = D.inDomain;  rule(k) = D.betaRule;
        if abs(f0Med15(k) / f0MedBio(k) - 1) > 1e-6
            error('datasetOpCoverage2:F0Source', '%s', sprintf('%s: median f0 %.6f (v015) vs %.6f (noiseCharacterisation)', ...
                spec.name, f0Med15(k), f0MedBio(k)));
        end
    end

    %% Variants
    rngA = range(mono.alphaGrid); rngS = range(mono.sigmaGrid); rngF = range(mono.fsGrid); rngV = range(mono.VGFGrid);
    vgfV = {f0MedBio * KCONV_RETIRED, f0MedBio * K13, f0Med15 .* gridKConv_v001(betaDs)};
    vname = ["V0 RETIRED (175.2636)", "V1 CONSTANT (K(1/3))", "V2 PRIMARY (K(beta_ds))"];
    allRes = table();  nnByVariant = strings(nSets, numel(vgfV));  R = cell(1, numel(vgfV));
    for v = 1:numel(vgfV)
        Xn = [(alphaMed - min(mono.alphaGrid)) / rngA, (sigmaMed - min(mono.sigmaGrid)) / rngS, ...
              (fsMed - min(mono.fsGrid)) / rngF, (vgfV{v} - min(mono.VGFGrid)) / rngV];
        Dmat = squareform(pdist(Xn));  Dmat(Dmat == 0) = NaN;
        [nnDist, nnIdx] = min(Dmat, [], 2, 'omitnan');
        T = table(names, inDom, alphaMed, sigmaMed, fsMed, betaDs, vgfV{v}, nnDist, names(nnIdx), ...
            'VariableNames', {'dataset','inDomain','alpha','sigmaMM','fs','betaDs','VGF','nnDist','nearestNeighbour'});
        nnByVariant(:, v) = names(nnIdx);
        T.variant = repmat(v - 1, height(T), 1);  allRes = [allRes; T]; %#ok<AGROW>
        R{v} = T;
    end

    %% Regression anchor
    fr0 = R{1}.nnDist(R{1}.dataset == "Fraser");
    if abs(fr0 - FRASER_NN_PUB) > 1e-3
        error('datasetOpCoverage2:Anchor', '%s', sprintf('V0 Fraser NN distance %.4f does not reproduce Finding #207 (%.4f)', fr0, FRASER_NN_PUB));
    end
    fprintf('REGRESSION ANCHOR passed: V0 Fraser NN distance %.4f reproduces Finding #207 (%.4f).\n', fr0, FRASER_NN_PUB);

    %% Report
    for v = 1:numel(vgfV)
        T = sortrows(R{v}, 'nnDist', 'descend');
        fprintf('\n=== %s ===\n', vname(v));
        fprintf('%-13s %3s %7s %7s %5s %7s %8s %9s %s\n', 'dataset','dom','alpha','sigmaMM','fs','betaDs','VGF','nnDist','nearest');
        for r = 1:height(T)
            fprintf('%-13s %3d %7.3f %7.3f %5.0f %7.3f %8.1f %9.4f %s\n', T.dataset(r), T.inDomain(r), T.alpha(r), T.sigmaMM(r), ...
                T.fs(r), T.betaDs(r), T.VGF(r), T.nnDist(r), T.nearestNeighbour(r));
        end
    end
    fprintf('\nNEAREST-NEIGHBOUR STRUCTURE AGAINST V0:\n');
    for v = 2:numel(vgfV)
        chg = nnByVariant(:, v) ~= nnByVariant(:, 1);
        o0 = R{1}.dataset(argsort(R{1}.nnDist));  o1 = R{v}.dataset(argsort(R{v}.nnDist));
        if any(chg), fprintf('  %s: nearest neighbour CHANGES for %s\n', vname(v), strjoin(names(chg)', ', '));
        else, fprintf('  %s: all seven nearest neighbours unchanged\n', vname(v)); end
        fprintf('     isolation order (most to least) %s: %s\n', 'V0', strjoin(flip(o0)', ' > '));
        fprintf('     isolation order (most to least) %s: %s\n', 'V' + string(v - 1), strjoin(flip(o1)', ' > '));
    end

    results = R{end};   % the primary variant
    if opts.Save
        resDir = fullfile(ROOT, 'results');
        if ~exist(resDir, 'dir'), error('datasetOpCoverage2:NoResultsDir', '%s', 'results/ not found at the project root'); end
        matOut = fullfile(resDir, 'datasetOperatingPointCoverage_v003.mat');
        constants = struct('K13', K13, 'KCONV_RETIRED', KCONV_RETIRED, 'runDate', string(datetime('now', 'Format', 'yyyy-MM-dd')));
        save(matOut, 'results', 'allRes', 'nnByVariant', 'vname', 'constants', 'monoFile', '-v7.3');
        fprintf('\nSaved: %s\n', matOut);
    end
end

function i = argsort(x)
    [~, i] = sort(x, 'ascend');
end
