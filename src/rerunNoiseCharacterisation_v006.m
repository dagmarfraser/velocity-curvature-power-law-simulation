function rerunNoiseCharacterisation_v006(opts)
% rerunNoiseCharacterisation_v006  v005 + the FraserN48 registry row (2026-11-03).
%
% WHAT v006 CHANGES IN v005, AND ONLY THIS
%   1. Registry row 'FraserN48': importDB_fraser_v002 (N=48 corpus: gaps split,
%      direct touches dropped, ID101/102/214 out, ID401-403 in), saved to
%      noiseCharacterisation_fraserN48.mat. The 'Fraser' row and its file
%      (the v015/v019 input) are untouched.
%   2. Fail loud: an importer error now stops the run (v005 printed it and
%      skipped the dataset).
%   3. No silent overwrite: an existing noise file is never replaced unless
%      Overwrite=true is passed.
%   Method unchanged: characteriseBiologicalNoise_v004(Method="irasa").
%   FraserN48 sigma is in UIKit points (5.209 per mm on the M4 13" iPad);
%   runLoopClosureFftnoise_v020 converts with sigToMM = 1/5.209.
%   USAGE: rerunNoiseCharacterisation_v006(DS='FraserN48')
%
% v005's own header follows unchanged.
% rerunNoiseCharacterisation_v005  Adds Fraser to the v004 dataset registry.
%
% v004 fixed the raw-data IRASA artefact across the original 7 in-scope
% datasets (closed out, TODO_characteriseBiologicalNoise_v004.md archived).
% v005 is a distinct purpose: register Fraser (Analysis_2026 replacement
% cohort) using the same characteriseBiologicalNoise_v004 method, unchanged.
% Pilot's registry entry is KEPT, not removed -- whether to drop Pilot,
% retain it as a supplementary check, or fold it in some other way is an
% open decision (TODO_GitHubReleasePrep_v001.md A1), not something to
% resolve by deleting its registry row here.
%
% NOTE ON UNITS: characteriseBiologicalNoise_v004's sigmaMean output is in
% whatever units trial.x/trial.y are in -- raw pixels for both Pilot and
% Fraser (importDB_pilot_v001 / importDB_fraser_v001 store raw px, do not
% pre-convert to mm). The px-to-mm conversion (Pilot: 9.73, Fraser:
% 10.41793 -- see data/fraser/README.md point 3) is applied downstream in
% runLoopClosureFftnoise_v007.m's per-dataset sigToMM, NOT here. Any sigma
% figure quoted from this script's own printed summary or saved mat is
% therefore in PIXELS, not mm, for both datasets -- convert before quoting
% alongside other datasets' mm-native sigma values.
%
% USAGE:
%   rerunNoiseCharacterisation_v006(DS='Fraser')   % Fraser only (typical use)
%   rerunNoiseCharacterisation_v006()               % all 8 registered datasets
%
% Fraser, D.S. (2026)
% See also: rerunNoiseCharacterisation_v004, characteriseBiologicalNoise_v004

    arguments
        opts.DS      (1,1) string = "all"
        opts.DryRun  (1,1) logical = false
        opts.Overwrite (1,1) logical = false   % v006
        opts.Verbose (1,1) logical = true
    end

    srcDir  = fileparts(mfilename('fullpath'));
    dataDir = fullfile(srcDir, '..', 'data');
    addpath(genpath(fullfile(srcDir, 'functions')));

    % Dataset registry: name, importer call, mat filename, v003 ira_alphaMean
    % (NaN where no v003 predecessor exists -- Fraser is new, not a v003->v004
    % migration, so the delta-based sanity check below is skipped for it.)
    DS = {
        'Pilot',        @pilotImport_local,    'noiseCharacterisation_pilot.mat',       2.571;
        'Cook CTRL',    @cookCTRLImport_local,  'noiseCharacterisation_cook.mat',        3.649;
        'Cook ASD',     @cookASDImport_local,   'noiseCharacterisation_cookASD.mat',     4.115;
        'Hickman PLAC', @hickPLACImport_local,  'noiseCharacterisation_hickmanPLAC.mat', 4.562;
        'Hickman HALO', @hickHALOImport_local,  'noiseCharacterisation_hickmanHALO.mat', 4.536;
        'Zarandi',      @zarandiImport_local,   'noiseCharacterisation_zarandi.mat',     2.531;
        'Dhieb',        @dhiebImport_local,     'noiseCharacterisation_dhieb.mat',       3.437;
        'Fraser',       @fraserImport_local,    'noiseCharacterisation_fraser.mat',      NaN;
        'FraserN48',    @fraserN48Import_local, 'noiseCharacterisation_fraserN48.mat',   NaN;   % v006
    };

    if ~strcmp(opts.DS, "all")
        mask = strcmpi(DS(:,1), opts.DS);
        if ~any(mask)
            error('rerunNoiseCharacterisation_v006:unknownDS', '%s', ...
                ['Unknown dataset: ' char(opts.DS) '. Options: ' strjoin(DS(:,1)', ', ')]);
        end
        DS = DS(mask, :);
    end

    nDS   = size(DS, 1);
    tTotal = tic;

    fprintf('=== NOISE CHARACTERISATION (v006 registry): characteriseBiologicalNoise_v004 ===\n');
    fprintf('Datasets: %d   DryRun: %d\n\n', nDS, opts.DryRun);

    results_summary = struct();

    for di = 1:nDS
        dsName   = DS{di,1};
        importFn = DS{di,2};
        matName  = DS{di,3};
        alphaV3  = DS{di,4};
        matPath  = fullfile(srcDir, matName);

        fprintf('--- [%d/%d] %s ---\n', di, nDS, dsName);
        tDS = tic;

        %% Load trials
        fprintf('  Loading trials...\n');
        if isfile(matPath) && ~opts.Overwrite && ~opts.DryRun   % v006: never replace silently
            error('rerunNoiseCharacterisation_v006:exists', '%s', ...
                ['Noise file exists: ' matPath '. Pass Overwrite=true to replace it.']);
        end
        trials = importFn(dataDir);   % v006: an importer error stops the run (v005 skipped the dataset)
        fprintf('  N trials: %d\n', numel(trials));

        %% Run characteriseBiologicalNoise_v004 (method itself unchanged)
        fprintf('  Running characteriseBiologicalNoise_v004(Method="irasa")...\n');
        bioResults = characteriseBiologicalNoise_v004(trials, ...
            Method="irasa", Verbose=opts.Verbose);

        %% Report
        alphaV4 = mean(bioResults.ira_alphaMean, 'omitnan');
        sigmaM  = mean(bioResults.sigmaMean, 'omitnan');   % PIXELS, see header note

        if isnan(alphaV3)
            fprintf('\n  ira_alphaMean: %.3f  (no v003 predecessor -- new dataset)\n', alphaV4);
        else
            delta = alphaV4 - alphaV3;
            fprintf('\n  ira_alphaMean: %.3f  (v003 was %.3f  delta=%.3f)\n', ...
                alphaV4, alphaV3, delta);
            if delta < 0.5
                fprintf('  WARNING: delta < 0.5 : v004 alpha not substantially higher than v003.\n');
            else
                fprintf('  OK: delta=%.3f confirms harmonic contamination removed.\n', delta);
            end
        end
        fprintf('  sigmaMean:     %.3f px (NOT mm -- see header note)\n', sigmaM);

        results_summary(di).name     = dsName;
        results_summary(di).alphaV4  = alphaV4;
        results_summary(di).sigmaPx  = sigmaM;
        results_summary(di).nTrials  = numel(trials);

        %% Save
        if ~opts.DryRun
            save(matPath, 'bioResults', '-v7.3');
            fprintf('  Saved: %s\n', matPath);
        else
            fprintf('  DryRun: NOT saved.\n');
        end

        fprintf('  Time: %.0f s\n\n', toc(tDS));
    end

    %% Final summary
    fprintf('=== SUMMARY ===\n');
    fprintf('%-14s  %8s  %10s  %8s\n', 'Dataset', 'alpha', 'sigma(px)', 'nTrials');
    fprintf('%s\n', repmat('-', 1, 46));
    for di = 1:numel(results_summary)
        r = results_summary(di);
        if isempty(r.name), continue; end
        fprintf('%-14s  %8.3f  %10.3f  %8d\n', r.name, r.alphaV4, r.sigmaPx, r.nTrials);
    end
    fprintf('\nTotal time: %.0f s\n', toc(tTotal));
end

%% ========== DATASET IMPORTERS ==========

function trials = pilotImport_local(dataDir)
    trials = importDB_pilot_v001(fullfile(dataDir, 'pilot'), ...
        'Shapes', 3, 'Visits', [1 2], 'Verbose', false);
end

function trials = cookCTRLImport_local(~)
    trials = importDB_cook_v002('Group', 'CTRL', 'Tasks', 7, 'Verbose', false);
end

function trials = cookASDImport_local(~)
    trials = importDB_cook_v002('Group', 'ASD', 'Tasks', 7, 'Verbose', false);
end

function trials = hickPLACImport_local(~)
    trials = importDB_hickman_v003('Study', 2, 'Group', 'PLAC', 'Verbose', false);
end

function trials = hickHALOImport_local(~)
    trials = importDB_hickman_v003('Study', 2, 'Group', 'HALO', 'Verbose', false);
end

function trials = zarandiImport_local(~)
    trials = importDB_zarandi_v001('Verbose', false);
end

function trials = dhiebImport_local(~)
    trials = importDB_dhieb_v001('Verbose', false);
end

function trials = fraserImport_local(dataDir)
    trials = importDB_fraser_v001(fullfile(dataDir, 'fraser'), ...
        'Shapes', 3, 'Visits', [1 2], 'Verbose', false);
end

function trials = fraserN48Import_local(~)   % v006
    trials = importDB_fraser_v002('Verbose', false);
end
