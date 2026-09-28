%% runSimexBetaCorpus_v001.m
%
% Batch runner for simexBetaRecovery_v001 across empirical trials.
% Parallel over TRIALS (parfor); each worker runs all configured
% pipelines for its own trial sequentially. This is the natural unit of
% independent work -- MATLAB does not support nested parfor, so
% parallelising simexBetaRecovery_v001's internal replicate loop would
% conflict with parallelising across trials, and would also mean
% duplicating each trial's x/y data six-fold across workers instead of
% once. The diagnostic function itself is left untouched.
%
% SCOPE, this version: Hickman Study 2 only (PLAC and/or HALO), the
% dataset already verified end-to-end in the single-trial smoke test
% (trialID matched against noiseCharacterisation_hickman*.mat, x/y
% converted to mm via each trial's own SigmaToMM, not a guessed
% constant). Extending to Zarandi/Cook/Dhieb/Fraser/Pilot/Dagenais needs
% each one's own importer and SigmaToMM convention checked the same way
% before being added here -- do not assume they match this one.
%
% CONFIG.MaxTrialsPerGroup defaults to a small number so a first run
% times a real subset on whichever machine runs this before committing
% to the full corpus. Set to Inf for the real run.
%
% CORRECTED 2026-09-20: the OLS pipeline rows were wired to
% regressDataEBR case 2 (fitlm, intercept forced to zero), not case 3
% (fitlm, intercept). The prereg's own OLS definition is log-linear
% regression WITH intercept -- case 2 assumes VGF=1, which real trial
% data does not satisfy, and produced betaSimex ~1.3-1.4 for
% BWFD-OLS/SG-OLS in the first run on 20 real trials. Fraser et al.
% (2025) already document this exact failure mode for forced-zero-
% intercept regression. LMLS/IRLS rows were unaffected (case 3 was
% never involved). Re-verified against the same 20-trial subset after
% the fix.
%
% Created 2026-09-21. Dagmar Scott Fraser, d.s.fraser@bham.ac.uk

%% ===================== CONFIG ==============================
CONFIG.ProjectRoot      = fileparts(fileparts(mfilename('fullpath'))); % assumes this file lives in src/
CONFIG.Groups           = {'PLAC', 'HALO'};   % Hickman Study 2 groups to include
CONFIG.MaxTrialsPerGroup = Inf;                 % Inf for the real run; small for a timing dry-run first
CONFIG.NReplicates      = 50;
CONFIG.NBootstrap       = 200;
CONFIG.LambdaGrid       = [0 0.5 1 1.5 2];
CONFIG.NWorkers         = min(10, feature('numcores'));
CONFIG.OutputMatFile    = fullfile(CONFIG.ProjectRoot, 'results', ...
    sprintf('simexBetaCorpus_hickman_%s.mat', datestr(now, 'yyyymmdd_HHMMSS')));

CONFIG.Pipelines = {  % {diffFilterType, diffFilterParams, regressType, label}
    2, [2 10 1], 3, 'BWFD-OLS'
    2, [2 10 1], 4, 'BWFD-LMLS'
    2, [2 10 1], 5, 'BWFD-IRLS'
    6, [4 17],   3, 'SG-OLS'
    6, [4 17],   4, 'SG-LMLS'
    6, [4 17],   5, 'SG-IRLS'
};
CONFIG.LMseeds = [1 1/3];
%% ============================================================

tRunStart = tic;

addpath(fullfile(CONFIG.ProjectRoot, 'src', 'functions'));
addpath(fullfile(CONFIG.ProjectRoot, 'src', 'req', ...
    'dagmarfraser-fractal-noise-6cfae90', 'toolbox'));

if ~isfolder(fullfile(CONFIG.ProjectRoot, 'results'))
    mkdir(fullfile(CONFIG.ProjectRoot, 'results'));
end

%% ----- Build the job list: one entry per trial, each carrying its own noise -----
allJobs = struct('trialID', {}, 'group', {}, 'x', {}, 'y', {}, 'fs', {}, ...
    'alpha', {}, 'sigma', {});
skippedNaNTrials = {};

for iGroup = 1:numel(CONFIG.Groups)
    grp = CONFIG.Groups{iGroup};

    noiseFile = fullfile(CONFIG.ProjectRoot, 'src', ...
        sprintf('noiseCharacterisation_hickman%s.mat', grp));
    if ~isfile(noiseFile)
        error('runSimexBetaCorpus_v001:NoiseFileMissing', '%s', ...
            sprintf('Expected noise file not found: %s', noiseFile));
    end
    S = load(noiseFile);
    bioResults = S.bioResults;

    trials = importDB_hickman_v003(Study=2, Group=string(grp), Verbose=false);
    if isempty(trials)
        error('runSimexBetaCorpus_v001:NoTrials', '%s', ...
            sprintf('importDB_hickman_v003 returned zero trials for group %s.', grp));
    end

    nTrials = min(numel(trials), CONFIG.MaxTrialsPerGroup);
    for iTrial = 1:nTrials
        tr  = trials(iTrial);
        row = bioResults(strcmp(bioResults.trialID, tr.trialID), :);
        if height(row) ~= 1
            warning('runSimexBetaCorpus_v001:NoiseRowMismatch', ...
                'Trial %s: expected exactly 1 matching noise row, found %d. Skipping.', ...
                tr.trialID, height(row));
            continue
        end
        if isnan(row.alphaMean) || isnan(row.sigmaMean)
            % Not a negative-sigma problem (mustBeNonnegative in
            % simexBetaRecovery_v001 rejects NaN too, with a message
            % that suggests the wrong cause). Skipped here, once, with
            % the real reason, instead of failing 6x downstream per
            % pipeline with a misleading "must be nonnegative" error.
            skippedNaNTrials{end+1} = tr.trialID; %#ok<AGROW>
            continue
        end

        job          = struct();
        job.trialID  = tr.trialID;
        job.group    = grp;
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
    warning('runSimexBetaCorpus_v001:SkippedNaNNoiseTrials', ...
        ['%d trial(s) skipped for NaN alpha/sigma in the noise ' ...
         'characterisation table (not a data-quality issue this ' ...
         'runner introduces -- worth checking why these specific ' ...
         'trials have no valid noise characterisation): %s'], ...
        numel(skippedNaNTrials), strjoin(skippedNaNTrials, ', '));
end
fprintf('%d trials queued across %d group(s), %d pipelines each = %d cells.\n', ...
    nJobs, numel(CONFIG.Groups), size(CONFIG.Pipelines, 1), nJobs * size(CONFIG.Pipelines, 1));

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

    trialResults{iJob} = struct('trialID', job.trialID, 'group', job.group, ...
        'pipelineResults', {pipelineOut});

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
            rows(end+1, :) = {tRes.trialID, tRes.group, p.pipeline, true, '', ...
                p.results.betaMean(1), p.results.betaSimex, ...
                p.results.betaSimexCI(1), p.results.betaSimexCI(2), ...
                min(p.results.nValid)}; %#ok<AGROW>
        else
            rows(end+1, :) = {tRes.trialID, tRes.group, p.pipeline, false, ...
                p.errorMessage, NaN, NaN, NaN, NaN, NaN}; %#ok<AGROW>
        end
    end
end
summaryTable = cell2table(rows, 'VariableNames', ...
    {'trialID', 'group', 'pipeline', 'success', 'errorMessage', ...
     'betaObserved', 'betaSimex', 'betaSimexCI_lo', 'betaSimexCI_hi', 'minNValid'});

nFailed = sum(~summaryTable.success);
if nFailed > 0
    warning('runSimexBetaCorpus_v001:SomeCellsFailed', ...
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
