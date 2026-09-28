%% runSimexBetaCorpus_v002.m
%
% Batch runner for simexBetaRecovery_v001 across empirical trials.
% Generalises v001 (Hickman Study 2 only) to any of Fraser / Cook CTRL /
% Cook ASD / Hickman PLAC / Hickman HALO / Zarandi / Dhieb, via a
% per-dataset spec table (CONFIG.DatasetSpecs). Each spec supplies its own
% noise-characterisation mat and a zero-arg closure over that dataset's own
% importer, so each importer's own name-value contract (Study/Group for
% Hickman, Group for Cook, none for Fraser/Zarandi/Dhieb) stays local to
% its own closure rather than forcing a common interface across importers
% that don't share one. Parallel-over-trials structure, pipeline
% definitions, NaN-noise skip logic, and the OLS-intercept fix are
% unchanged from v001 -- only dataset selection is new.
%
% Pilot and Peters CONG are deliberately NOT specced here: Pilot is
% dropped project-wide ("persona non grata", 2026-08-07), and Peters CONG
% is closed (EPP-OFF retest never available). Neither is a live dataset
% for this project; adding them would need the same importer-contract
% check the others below already got before being added.
%
% SCOPE, this version: CONFIG.ActiveDatasets selects which specs actually
% run. Default is the two fast Cook arms (~9 min combined at Hickman's
% measured 0.445 s/cell) -- a smoke test of the generalisation itself
% before committing to Fraser's much larger corpus (est. ~2 hours at
% N_total=2,774 valid trials). Widen ActiveDatasets once this default
% run is confirmed clean.
%
% CONFIG.MaxTrialsPerGroup defaults to a small number so a first run
% times a real subset on whichever machine runs this before committing
% to the full corpus. Set to Inf for the real run.
%
% Created 2026-09-21. Dagmar Scott Fraser, d.s.fraser@bham.ac.uk

%% ===================== CONFIG ==============================
CONFIG.ProjectRoot      = fileparts(fileparts(mfilename('fullpath'))); % assumes this file lives in src/
CONFIG.MaxTrialsPerGroup = Inf;                 % Inf for the real run; small for a timing dry-run first
CONFIG.NReplicates      = 50;
CONFIG.NBootstrap       = 200;
CONFIG.LambdaGrid       = [0 0.5 1 1.5 2];
CONFIG.NWorkers         = max(1, feature('numcores') - 1);   % leave one core free for the machine's own use

% --- Which specs actually run this call. Edit to widen. ---
% Phase 1 (default): Cook CTRL + Cook ASD, ~656 cells, ~9 min -- fast
% smoke test of the v001->v002 generalisation itself.
% Phase 2: add "fraser" (2,774 valid trials, est. ~2h -- decide N_total
% vs N_fin=2,744 before running; this file uses N_total, i.e. every valid
% trial the importer + noise mat agree on, same convention as v001 used
% for Hickman).
% Phase 3 (lower priority -- see scoping note): "zarandi", "dhieb". Both
% already excluded from the paper's headline identifiability table for
% reasons SIMEX won't change (Zarandi validation degeneracy; Dhieb's
% Windows-mouse-precision prefiltering masks the true noise floor below
% ~10 Hz, undercutting SIMEX's own noise-reinjection premise). Include
% only with that caveat carried into any write-up.
CONFIG.ActiveDatasets = ["hickmanPLAC", "hickmanHALO"];

CONFIG.Pipelines = {  % {diffFilterType, diffFilterParams, regressType, label}
    2, [2 10 1], 3, 'BWFD-OLS'
    2, [2 10 1], 4, 'BWFD-LMLS'
    2, [2 10 1], 5, 'BWFD-IRLS'
    6, [4 17],   3, 'SG-OLS'
    6, [4 17],   4, 'SG-LMLS'
    6, [4 17],   5, 'SG-IRLS'
};
CONFIG.LMseeds = [1 1/3];

CONFIG.OutputMatFile = fullfile(CONFIG.ProjectRoot, 'results', ...
    sprintf('simexBetaCorpus_%s_%s.mat', ...
        strjoin(CONFIG.ActiveDatasets, '-'), datestr(now, 'yyyymmdd_HHMMSS')));
%% ============================================================

tRunStart = tic;

addpath(fullfile(CONFIG.ProjectRoot, 'src', 'functions'));
addpath(fullfile(CONFIG.ProjectRoot, 'src', 'req', ...
    'dagmarfraser-fractal-noise-6cfae90', 'toolbox'));

if ~isfolder(fullfile(CONFIG.ProjectRoot, 'results'))
    mkdir(fullfile(CONFIG.ProjectRoot, 'results'));
end

%% ----- Dataset spec table -----
% Each spec: .label (dataset name for output rows), .group (arm name for
% output rows; '' for single-arm datasets), .noiseFile (full path,
% resolved below), .importerFn (zero-arg closure returning the trials
% struct array). Verified field-compatible against simexBetaRecovery_v001's
% contract (.x .y .fs .trialID .SigmaToMM) and against each
% noiseCharacterisation_*.mat's bioResults schema (trialID, alphaMean,
% sigmaMean) directly, 2026-09-21 -- not assumed.
specsAll = struct('label', {}, 'group', {}, 'noiseFile', {}, 'importerFn', {});

specsAll(end+1) = struct('label', 'hickman', 'group', 'PLAC', ...
    'noiseFile', 'noiseCharacterisation_hickmanPLAC.mat', ...
    'importerFn', @() importDB_hickman_v003(Study=2, Group="PLAC", Verbose=false));
specsAll(end+1) = struct('label', 'hickman', 'group', 'HALO', ...
    'noiseFile', 'noiseCharacterisation_hickmanHALO.mat', ...
    'importerFn', @() importDB_hickman_v003(Study=2, Group="HALO", Verbose=false));
specsAll(end+1) = struct('label', 'cook', 'group', 'CTRL', ...
    'noiseFile', 'noiseCharacterisation_cook.mat', ...
    'importerFn', @() importDB_cook_v002(Group="CTRL", Verbose=false));
specsAll(end+1) = struct('label', 'cook', 'group', 'ASD', ...
    'noiseFile', 'noiseCharacterisation_cookASD.mat', ...
    'importerFn', @() importDB_cook_v002(Group="ASD", Verbose=false));
specsAll(end+1) = struct('label', 'fraser', 'group', '', ...
    'noiseFile', 'noiseCharacterisation_fraser.mat', ...
    'importerFn', @() importDB_fraser_v001(Verbose=false));
specsAll(end+1) = struct('label', 'zarandi', 'group', '', ...
    'noiseFile', 'noiseCharacterisation_zarandi.mat', ...
    'importerFn', @() importDB_zarandi_v001(Verbose=false));
specsAll(end+1) = struct('label', 'dhieb', 'group', '', ...
    'noiseFile', 'noiseCharacterisation_dhieb.mat', ...
    'importerFn', @() importDB_dhieb_v001(Verbose=false));

% Map CONFIG.ActiveDatasets strings onto specsAll rows
activeKeyOf = @(s) lower(strrep([s.label s.group], ' ', ''));
allKeys     = arrayfun(activeKeyOf, specsAll, 'UniformOutput', false);

specs = struct('label', {}, 'group', {}, 'noiseFile', {}, 'importerFn', {});
for iAct = 1:numel(CONFIG.ActiveDatasets)
    wantKey = lower(CONFIG.ActiveDatasets(iAct));
    matchIdx = find(strcmp(allKeys, wantKey));
    if isempty(matchIdx)
        error('runSimexBetaCorpus_v002:UnknownDataset', '%s', ...
            sprintf(['ActiveDatasets entry "%s" does not match any spec key. ' ...
                     'Known keys: %s'], wantKey, strjoin(allKeys, ', ')));
    end
    specs(end+1) = specsAll(matchIdx); %#ok<AGROW>
end

%% ----- Build the job list: one entry per trial, each carrying its own noise -----
allJobs = struct('trialID', {}, 'dataset', {}, 'group', {}, 'x', {}, 'y', {}, ...
    'fs', {}, 'alpha', {}, 'sigma', {});
skippedNaNTrials = {};

for iSpec = 1:numel(specs)
    spec = specs(iSpec);

    noiseFile = fullfile(CONFIG.ProjectRoot, 'src', spec.noiseFile);
    if ~isfile(noiseFile)
        error('runSimexBetaCorpus_v002:NoiseFileMissing', '%s', ...
            sprintf('Expected noise file not found: %s', noiseFile));
    end
    S = load(noiseFile);
    bioResults = S.bioResults;

    trials = spec.importerFn();
    if isempty(trials)
        error('runSimexBetaCorpus_v002:NoTrials', '%s', ...
            sprintf('Importer for %s/%s returned zero trials.', spec.label, spec.group));
    end

    nTrials = min(numel(trials), CONFIG.MaxTrialsPerGroup);
    for iTrial = 1:nTrials
        tr  = trials(iTrial);
        row = bioResults(strcmp(bioResults.trialID, tr.trialID), :);
        if height(row) ~= 1
            warning('runSimexBetaCorpus_v002:NoiseRowMismatch', ...
                '%s/%s trial %s: expected exactly 1 matching noise row, found %d. Skipping.', ...
                spec.label, spec.group, tr.trialID, height(row));
            continue
        end
        if isnan(row.alphaMean) || isnan(row.sigmaMean)
            % Not a negative-sigma problem (mustBeNonnegative in
            % simexBetaRecovery_v001 rejects NaN too, with a message
            % that suggests the wrong cause). Skipped here, once, with
            % the real reason, instead of failing 6x downstream per
            % pipeline with a misleading "must be nonnegative" error.
            skippedNaNTrials{end+1} = sprintf('%s/%s:%s', spec.label, spec.group, tr.trialID); %#ok<AGROW>
            continue
        end

        job          = struct();
        job.trialID  = tr.trialID;
        job.dataset  = spec.label;
        job.group    = spec.group;
        job.x        = tr.x * tr.SigmaToMM;
        job.y        = tr.y * tr.SigmaToMM;
        job.fs       = tr.fs;
        job.alpha    = row.alphaMean;
        job.sigma    = row.sigmaMean * tr.SigmaToMM;
        allJobs(end+1) = job; %#ok<AGROW>
    end
end

nJobs = numel(allJobs);
if ~isempty(skippedNaNTrials)
    warning('runSimexBetaCorpus_v002:SkippedNaNNoiseTrials', ...
        ['%d trial(s) skipped for NaN alpha/sigma in the noise ' ...
         'characterisation table (not a data-quality issue this ' ...
         'runner introduces -- worth checking why these specific ' ...
         'trials have no valid noise characterisation): %s'], ...
        numel(skippedNaNTrials), strjoin(skippedNaNTrials, ', '));
end
fprintf('%d trials queued across %d spec(s), %d pipelines each = %d cells.\n', ...
    nJobs, numel(specs), size(CONFIG.Pipelines, 1), nJobs * size(CONFIG.Pipelines, 1));

%% ----- Run, parallel over trials -----
if isempty(gcp('nocreate'))
    parpool('local', CONFIG.NWorkers);
end

pipelines  = CONFIG.Pipelines;
nPipelines = size(pipelines, 1);
lambdaGrid = CONFIG.LambdaGrid;
nRep       = CONFIG.NReplicates;
nBoot      = CONFIG.NBootstrap;
lmSeeds    = CONFIG.LMseeds;

progressQueue = parallel.pool.DataQueue;
afterEach(progressQueue, @(~) fprintf('.'));

trialResults = cell(nJobs, 1);

parfor iJob = 1:nJobs
    job = allJobs(iJob); %#ok<PFBNS>
    pipelineOut = cell(nPipelines, 1);

    for iPipe = 1:nPipelines
        diffType   = pipelines{iPipe, 1};
        diffParams = pipelines{iPipe, 2};
        regType    = pipelines{iPipe, 3};
        label      = pipelines{iPipe, 4};

        try
            r = simexBetaRecovery_v001(job.x, job.y, job.fs, job.alpha, job.sigma, ...
                diffType, diffParams, regType, lmSeeds, ...
                LambdaGrid=lambdaGrid, NReplicates=nRep, NBootstrap=nBoot);
            pipelineOut{iPipe} = struct('pipeline', label, 'success', true, ...
                'errorMessage', '', 'results', r);
        catch ME
            pipelineOut{iPipe} = struct('pipeline', label, 'success', false, ...
                'errorMessage', ME.message, 'results', []);
        end
    end

    trialResults{iJob} = struct('trialID', job.trialID, 'dataset', job.dataset, ...
        'group', job.group, 'pipelineResults', {pipelineOut});

    send(progressQueue, iJob);
end
fprintf('\n');

%% ----- Flatten to a summary table; failures kept visible, never dropped -----
rows = {};
for iJob = 1:nJobs
    tRes = trialResults{iJob};
    for iPipe = 1:nPipelines
        p = tRes.pipelineResults{iPipe};
        if p.success
            rows(end+1, :) = {tRes.trialID, tRes.dataset, tRes.group, p.pipeline, true, '', ...
                p.results.betaMean(1), p.results.betaSimex, ...
                p.results.betaSimexCI(1), p.results.betaSimexCI(2), ...
                min(p.results.nValid)}; %#ok<AGROW>
        else
            rows(end+1, :) = {tRes.trialID, tRes.dataset, tRes.group, p.pipeline, false, ...
                p.errorMessage, NaN, NaN, NaN, NaN, NaN}; %#ok<AGROW>
        end
    end
end
summaryTable = cell2table(rows, 'VariableNames', ...
    {'trialID', 'dataset', 'group', 'pipeline', 'success', 'errorMessage', ...
     'betaObserved', 'betaSimex', 'betaSimexCI_lo', 'betaSimexCI_hi', 'minNValid'});

nFailed = sum(~summaryTable.success);
if nFailed > 0
    warning('runSimexBetaCorpus_v002:SomeCellsFailed', ...
        '%d of %d trial-pipeline cells failed; see summaryTable.errorMessage. Not silently dropped.', ...
        nFailed, height(summaryTable));
end

elapsedSeconds = toc(tRunStart);
save(CONFIG.OutputMatFile, 'summaryTable', 'trialResults', 'CONFIG', ...
    'elapsedSeconds', 'skippedNaNTrials', '-v7.3');

fprintf('Saved: %s\n', CONFIG.OutputMatFile);
fprintf('%d cells run (%d failed) in %.1f min (%.2f s/cell).\n', ...
    height(summaryTable), nFailed, elapsedSeconds / 60, ...
    elapsedSeconds / height(summaryTable));
