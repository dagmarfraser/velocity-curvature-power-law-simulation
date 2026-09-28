% knockDownFlags_v001.m  Re-run R5/R6 analyses at v004 noise centroids.
%
% Knocks down the three **FLAG** markers written into the skeleton during the
% characteriseBiologicalNoise_v004 update:
%
%   FLAG 1 : SEM values at empirical centroids (R5/R6):
%     Cook CTRL  SEM was 0.0075 at old α=3.65; now α=4.77
%     Cook ASD   SEM was 0.0056 at old α=4.12; now α=5.06
%     Hickman PLAC  SEM was 0.0013 at old α=4.60; now α=5.34
%     Hickman HALO  SEM was 0.0013 at old α=4.61; now α=5.42
%
%   FLAG 2 : Fig 6 invertibility coverage heatmap:
%     checkInvertibilityForEmpirical_v2_001 reads ira_alphaMean from
%     noiseCharacterisation mats (now v004). Re-run produces updated
%     empiricalCoverage_v2_001.png and invertibilityCoverage_v2_001.mat.
%
%   FLAG 3 : Fig 6 centroid forward-map figures:
%     Regenerate forwardMapCentroid_*.png at v004 (α, σ_mm) centroids:
%       Zarandi      α=3.18  σ=4.77 mm  fs=100 Hz
%       Cook CTRL    α=4.77  σ=8.15 mm  fs=133 Hz
%       Hickman PLAC α=5.34  σ=7.17 mm  fs=133 Hz
%
% v004 centroids from rerunNoiseCharacterisation_v004.m:
%   Zarandi:       ira_alphaMean=3.184  sigmaMean=0.477 cm  → 4.77 mm
%   Cook CTRL:     ira_alphaMean=4.771  sigmaMean=32.866 px → 8.15 mm
%   Cook ASD:      ira_alphaMean=5.062  sigmaMean=31.630 px → 7.84 mm
%   Hickman PLAC:  ira_alphaMean=5.337  sigmaMean=28.892 px → 7.17 mm
%   Hickman HALO:  ira_alphaMean=5.424  sigmaMean=29.895 px → 7.42 mm
%
% USAGE (run from src/ on iMac):
%   knockDownFlags_v001
%
% Fraser, D.S. (2026)

clearvars
srcDir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(srcDir, 'functions')));
cd(srcDir);

fprintf('=== KNOCK DOWN FLAGS: v004 centroid re-runs ===\n\n');

%% =========================================================
%% FLAG 1: SEM values at new centroids
%% =========================================================
fprintf('--- FLAG 1: SEM lookup at v004 centroids ---\n');

semFile = fullfile(srcDir, 'perCoordinateSEM_v2_001.mat');
if ~isfile(semFile)
    fprintf('  WARNING: %s not found : skip SEM lookup.\n', semFile);
else
    load(semFile, 'coordTable');
    T = coordTable;

    allAlpha = sort(unique(T.alpha));
    allSigma = sort(unique(T.sigma));
    allFs    = sort(unique(T.fs));

    % v004 centroids: [alpha_v004, sigma_mm, fs, label]
    centroids = {
        4.771,  8.15, 133, 'Cook CTRL   (was α=3.65 SEM=0.0075)';
        5.062,  7.84, 133, 'Cook ASD    (was α=4.12 SEM=0.0056)';
        5.337,  7.17, 133, 'Hickman PLAC (was α=4.60 SEM=0.0013)';
        5.424,  7.42, 133, 'Hickman HALO (was α=4.61 SEM=0.0013)';
        3.184,  4.77, 100, 'Zarandi      (was α=2.53)';
        4.289,  3.91, 240, 'Pilot        (was α=2.56)';
    };

    fprintf('\n  %-28s  %6s  %6s  %6s  %8s\n', ...
        'Dataset', 'alpha', 'sigma', 'fs', 'SEM(SG-IRLS)');
    fprintf('  %s\n', repmat('-',1,65));

    for ci = 1:size(centroids,1)
        tA = centroids{ci,1};
        tS = centroids{ci,2};
        tF = centroids{ci,3};
        lb = centroids{ci,4};

        [~,ai] = min(abs(allAlpha - tA));
        [~,si] = min(abs(allSigma - tS));
        [~,fi] = min(abs(allFs    - tF));
        sA = allAlpha(ai); sS = allSigma(si); sF = allFs(fi);

        % SEM for SG-IRLS (pipeline 6 in canonical order)
        sub = T(T.alpha==sA & T.sigma==sS & T.fs==sF & T.pipeline=="SG-IRLS", :);
        if isempty(sub)
            semVal = NaN;
        else
            semVal = mean(sub.sem, 'omitnan');
        end

        % Also get SEM across all pipelines
        subAll = T(T.alpha==sA & T.sigma==sS & T.fs==sF, :);
        if ~isempty(subAll)
            [G, pipNames] = findgroups(subAll.pipeline);
            semPerPipe = splitapply(@(x) mean(x,'omitnan'), subAll.sem, G);
            semStr = strjoin(arrayfun(@(i) sprintf('%s=%.4f', ...
                pipNames(i), semPerPipe(i)), 1:numel(pipNames), ...
                'UniformOutput',false), '  ');
        else
            semStr = 'no data';
        end

        fprintf('  %-28s  %6.3f  %6.2f  %6d  %8.4f\n', lb, sA, sS, sF, semVal);
        if ~isempty(subAll)
            fprintf('    Snap: α=%.2f σ=%.2f fs=%d | all pipes: %s\n', ...
                sA, sS, sF, semStr);
        end
    end
    fprintf('\nCopy these SEM values to replace the FLAG lines in skeleton R5.\n');
end

%% =========================================================
%% FLAG 2: Invertibility coverage at v004 centroids
%% =========================================================
fprintf('\n--- FLAG 2: checkInvertibilityForEmpirical_v2_001 ---\n');
fprintf('Running (reads updated noiseCharacterisation mats automatically)...\n');
try
    run(fullfile(srcDir, 'checkInvertibilityForEmpirical_v2_001.m'));
    fprintf('FLAG 2 DONE : empiricalCoverage_v2_001.png regenerated.\n');
catch ME
    fprintf('ERROR: %s\n', ME.message);
    fprintf('Run checkInvertibilityForEmpirical_v2_001.m manually.\n');
end

%% =========================================================
%% FLAG 3: Centroid forward-map figures at v004 centroids
%% =========================================================
fprintf('\n--- FLAG 3: plotForwardMapAtCentroid_v001 ---\n');

centroidFigs = {
    3.184,  4.77, 100, 'Zarandi';
    4.771,  8.15, 133, 'Cook_CTRL';
    5.337,  7.17, 133, 'Hickman_PLAC';
};

for ci = 1:size(centroidFigs,1)
    tA = centroidFigs{ci,1};
    tS = centroidFigs{ci,2};
    tF = centroidFigs{ci,3};
    lb = centroidFigs{ci,4};
    fprintf('  Generating: %s (α=%.3f, σ=%.2f mm, fs=%d Hz)...\n', lb, tA, tS, tF);
    try
        plotForwardMapAtCentroid_v001(tA, tS, tF, lb);
        fprintf('  OK\n');
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
    end
end

fprintf('\n=== ALL FLAGS KNOCKED DOWN ===\n');
fprintf('Next: update skeleton R5/R6 SEM values from FLAG 1 output above.\n');
fprintf('      Update Fig 6 filenames to new centroid coordinates.\n');
fprintf('      Run validateIrasaRecovery_v001(Part=2) if not already done.\n');
