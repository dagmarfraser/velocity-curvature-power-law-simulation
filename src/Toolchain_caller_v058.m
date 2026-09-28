% TOOLCHAIN_CALLER_V058  Pre-registration deviation: expanded noise grid
%
% Changes from v057 (frozen prereg production run):
%   1. debug = 5  ->  defineParameterSpace returns expanded noise/beta axes:
%        generatedBetas  0:(2/3)/20:0.7  (22 values; step=1/30; 0.7=21/30 exact;
%                        extends past 2/3 to clip descending beta_rec~1/3 branch)
%        noiseTypes      0:0.2:6   (31 values; coarser 0.2 step covers alpha 0-6;
%                        was 0:0.1:3, 31 values in v057)
%        noiseMagnitudes [..., 12, 15, 20] mm  (21 values; was 18)
%   2. Output DB: powerlaw_multiverse_v058.db  (~18.0M configs, 1.2x v057)
%   3. versionTc = 'Tc_v0058'
%
% Motivation (prereg deviation documented in defineParameterSpace.m debug=5 block):
%   Empirical alpha Cook/Hickman 4.1-4.6 exceeds current simulation ceiling 3.0.
%   Template-bias-corrected sigma ~14 mm exceeds ceiling 10 mm.
%   Invertibility analysis (Finding #18) shows 100% NoRise at empirical coords
%   when snapped to alpha=3; expanded grid resolves the snapping artefact,
%   enabling constellationCCC_v002 and Stage 5 CCC validation.
%   All other parameters (beta_gen, VGF, fs, shapes, filters, regressors,
%   trials) are identical to the frozen prereg (v057/debug=0).
%
% All v057 infrastructure (parallel pool, checkpoint/resume, JDBC workarounds,
% batch inserts) is unchanged.  Run from BlueBEAR interactive MATLAB session.
%
% Created April 2026
% Correspondence Dagmar Scott Fraser  d.s.fraser@bham.ac.uk
%
% Requires:
%   Database Toolbox, Parallel Computing Toolbox,
%   Curve Fitting Toolbox, Statistics and Machine Learning Toolbox

try
    rng('default');

    % ===== CONFIGURATION TOGGLES =====
    % debug=5 : expanded noise grid (prereg deviation), full parallel run
    debug        = 5;
    useFastBatch = true;

    callerDir = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(callerDir, 'functions')));
    addpath(genpath(fullfile(callerDir, 'req')));
    addpath(genpath(fullfile(callerDir, 'utils')));
    addpath(callerDir);
    commandwindow;

    [conn, dbFile, masterDir, jobIdentifier, isHPC] = setupEnvironment(debug);
    paramSpace = defineParameterSpace(debug);
    [parallelSettings, debugCores] = determineParallelSettings(debug, isHPC);
    cfg = createConfigSettings();
    cfg.dbTestMode = false;   % debug=5 is a real run, not a DB logic test

    [configIDs, startIdx, totalConfigs, generationTime] = ...
        initializeJobWithMethodChoice(conn, jobIdentifier, paramSpace, useFastBatch);

    displayParameterSpaceInfo(paramSpace, totalConfigs, parallelSettings, debugCores, isHPC, generationTime);

    [parallelMetrics, parRun, resourceTimer, poolObj] = ...
        createParallelPoolJustInTime(parallelSettings, debugCores, isHPC);

    processConfigurations(conn, configIDs, startIdx, totalConfigs, debugCores, ...
        parRun, parallelMetrics, cfg, jobIdentifier, isHPC, ...
        dbFile, resourceTimer, poolObj);

    markJobComplete(conn, jobIdentifier);

    if ~debugCores && isfield(parallelMetrics, 'taskTimes') && ~isempty(parallelMetrics.taskTimes)
        displayPerformanceReport(parallelMetrics, parRun, masterDir, jobIdentifier, isHPC);
    end

    try
        close(conn);
    catch ME
        warning(ME.identifier, '%s', ME.message);
    end

    fprintf('Processing complete! Results stored in: %s\n', dbFile);
    disp('Run extractL9Results or equivalent to analyse results.');

catch ME
    handleFatalError(ME);
end

%% ========================================================================
%% ALL INTERNAL FUNCTIONS (identical to v057 except setupEnvironment DB name
%% and createConfigSettings versionTc string)
%% ========================================================================

function [conn, dbFile, masterDir, jobIdentifier, isHPC] = setupEnvironment(debug)
thisFile  = mfilename('fullpath');
srcDir    = fileparts(thisFile);
masterDir = fileparts(srcDir);

% v058: points at powerlaw_multiverse_v058.db
if debug
    dbFile = fullfile(masterDir, 'results', 'powerlaw_debug_v058.db');
else
    dbFile = fullfile(masterDir, 'results', 'powerlaw_multiverse_v058.db');
end

setupPowerLawDB(dbFile);

try
    conn = sqlite(dbFile);
catch ME
    error('Failed to connect to SQLite database: %s', ME.message);
end

addpath(genpath(fullfile(masterDir, 'src', 'functions')));

slurm_job_id = getenv('SLURM_JOB_ID');
isHPC        = ~isempty(slurm_job_id);

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
if isHPC
    hostname = getenv('HOSTNAME');
    if isempty(hostname), hostname = 'hpc'; end
    currentJobIdentifier = sprintf('hpc_%s_%s', hostname, timestamp);
    fprintf('Running on HPC host: %s\n', hostname);
else
    currentJobIdentifier = sprintf('local_%s', timestamp);
    fprintf('Running locally\n');
end

try
    sqlquery = ['SELECT job_id, last_completed_idx, total_configs FROM job_checkpoints ' ...
                'WHERE status = "running" ORDER BY last_update_time DESC'];
    result = fetch(conn, sqlquery);

    if ~isempty(result)
        if istable(result)
            incompleteJobs = result;
        else
            incompleteJobs = cell2table(result, ...
                'VariableNames', {'job_id', 'last_completed_idx', 'total_configs'});
        end

        fprintf('\n===== INCOMPLETE JOBS FOUND =====\n');
        fprintf('There are %d incomplete jobs in the database.\n', height(incompleteJobs));

        rowsToShow = min(5, height(incompleteJobs));
        for i = 1:rowsToShow
            jobID   = incompleteJobs.job_id{i};
            lastIdx = double(incompleteJobs.last_completed_idx(i));
            totIdx  = double(incompleteJobs.total_configs(i));
            fprintf('[%d] %s - Progress: %d/%d (%.1f%%)\n', ...
                i, jobID, lastIdx, totIdx, 100*lastIdx/totIdx);
        end

        answer = input('\nDo you want to resume an incomplete job? (y/n): ', 's');
        if lower(answer) == 'y'
            jobToResume = 1;
            if height(incompleteJobs) > 1
                jobToResume = input('Enter the number of the job to resume [1]: ');
                if isempty(jobToResume), jobToResume = 1; end
            end
            if jobToResume > 0 && jobToResume <= height(incompleteJobs)
                jobIdentifier = incompleteJobs.job_id{jobToResume};
                fprintf('Resuming job: %s\n', jobIdentifier);
            else
                fprintf('Invalid selection. Using: %s\n', currentJobIdentifier);
                jobIdentifier = currentJobIdentifier;
            end
        else
            jobIdentifier = currentJobIdentifier;
        end
    else
        fprintf('No incomplete jobs found.\n');
        jobIdentifier = currentJobIdentifier;
    end
catch ME
    warning(ME.identifier, '%s', ME.message);
    jobIdentifier = currentJobIdentifier;
end

fprintf('Job identifier: %s\n', jobIdentifier);
end

function [parallelSettings, debugCores] = determineParallelSettings(debug, isHPC)
% debug=5 is a full production run — enable parallel
debugCores = false;

if isHPC
    poolType   = 'Processes';
    maxWorkers = feature('numcores');
    fprintf('HPC environment detected - will use ProcessPool\n');
else
    poolType   = 'Threads';
    maxWorkers = min(max(1, feature('numcores') - 1), 16);
    fprintf('Local environment detected - will use ThreadPool\n');
end

parallelSettings            = struct();
parallelSettings.poolType   = poolType;
parallelSettings.maxWorkers = maxWorkers;
parallelSettings.isHPC      = isHPC;
fprintf('Parallel pool will be created just-in-time with %d workers\n', maxWorkers);
end

function cfg = createConfigSettings()
cfg             = struct();
cfg.versionTc   = 'Tc_v0058';   % v058 marker
cfg.orbitCount  = 10;
cfg.saveAll     = 0;
cfg.display     = [0 0 0];
cfg.canvas      = [1920; 1080];
cfg.edgeClip    = 50;
cfg.curvatureChoice = 1;
cfg.limitBreak  = 0;
cfg.displayGraphs = 0;
cfg.MaticSpline = 0;
cfg.resample    = 20;
cfg.pixelScale  = 480 / 100;
cfg.variantPowerLaw = 1;
cfg.rethrowErrors   = false;
cfg.dbTestMode      = false;
end

function [configIDs, startIdx, totalConfigs, generationTime] = ...
    initializeJobWithMethodChoice(conn, jobIdentifier, paramSpace, useFastBatch)
setupCheckpointTable(conn);
jobExists = checkJobExists(conn, jobIdentifier);

if jobExists
    fprintf('Resuming job %s from checkpoint...\n', jobIdentifier);
    [configIDs, lastCompletedIdx] = getResumeInfo(conn, jobIdentifier);
    totalConfigs    = length(configIDs);
    generationTime  = 0;
    if lastCompletedIdx >= totalConfigs
        fprintf('Job %s already completed.\n', jobIdentifier);
        close(conn);
        error('Job already completed');
    end
    startIdx = lastCompletedIdx + 1;
    fprintf('Resuming from configuration %d of %d\n', startIdx, totalConfigs);
else
    fprintf('Starting new job %s...\n', jobIdentifier);
    generationTimer = tic;
    try
        if useFastBatch
            fprintf('Using fast streaming batch method...\n');
            configIDs = generateParameterConfigsDB_batch(conn, paramSpace);
        else
            fprintf('Using original individual INSERT method...\n');
            configIDs = generateParameterConfigsDB_fixed(conn, paramSpace);
        end
        generationTime = toc(generationTimer);
        totalConfigs   = length(configIDs);
        fprintf('Generated %d configurations in %s (%.0f/s)\n', ...
            totalConfigs, formatTime(generationTime), totalConfigs/generationTime);
    catch ME
        isMissing = strcmp(ME.identifier, 'MATLAB:UndefinedFunction') || ...
                    contains(ME.message, 'Undefined function');
        if useFastBatch && isMissing
            warning('MATLAB:PowerLaw:BatchMethodFailed', ...
                'Batch function not found, falling back to original method...');
            generationTimer = tic;
            configIDs      = generateParameterConfigsDB_fixed(conn, paramSpace);
            generationTime  = toc(generationTimer);
            totalConfigs    = length(configIDs);
        else
            rethrow(ME);
        end
    end
    createJobEntry(conn, jobIdentifier, configIDs);
    startIdx = 1;
end
end

function displayParameterSpaceInfo(paramSpace, totalConfigs, parallelSettings, debugCores, isHPC, generationTime)
fprintf('------------------------------------------------------\n');
fprintf('Multiverse analysis v058 (EXPANDED NOISE GRID): %d total configurations\n', totalConfigs);
fprintf('Parameter space: %d fs x %d beta x %d VGF x %d sigma x %d alpha x %d filters x %d regressors\n', ...
    length(paramSpace.samplingRates), length(paramSpace.generatedBetas), ...
    length(paramSpace.vgfValues),     length(paramSpace.noiseMagnitudes), ...
    length(paramSpace.noiseTypes),    length(paramSpace.filterTypes), ...
    length(paramSpace.regressTypes));
fprintf('Noise extension:  alpha 0->6 (%d values)  |  sigma up to 20 mm (%d values)\n', ...
    length(paramSpace.noiseTypes), length(paramSpace.noiseMagnitudes));
if generationTime > 0
    fprintf('Config generation: %s (%.0f configs/s)\n', ...
        formatTime(generationTime), totalConfigs/generationTime);
end
if ~debugCores
    fprintf('Pool: %s  %d workers\n', parallelSettings.poolType, parallelSettings.maxWorkers);
end
fprintf('Checkpoint interval: 200 configurations\n');
fprintf('------------------------------------------------------\n');
end

function [parallelMetrics, parRun, resourceTimer, poolObj] = ...
    createParallelPoolJustInTime(parallelSettings, debugCores, isHPC)
parallelMetrics = struct('poolType', '', 'startupTime', 0, 'taskTimes', []);
resourceTimer   = [];
poolObj         = [];

if ~debugCores
    fprintf('\n===== CREATING PARALLEL POOL JUST-IN-TIME =====\n');
    setupTimer    = tic;
    existingPool  = gcp('nocreate');
    needNewPool   = false;

    if isempty(existingPool)
        needNewPool = true;
    else
        existingType = class(existingPool);
        if isHPC && ~contains(existingType, 'ProcessPool')
            delete(existingPool); needNewPool = true;
        elseif ~isHPC && ~contains(existingType, 'ThreadPool')
            delete(existingPool); needNewPool = true;
        else
            fprintf('Using existing %s pool (%d workers)\n', existingType, existingPool.NumWorkers);
            poolObj = existingPool;
        end
    end

    if needNewPool
        fprintf('Creating %s pool with %d workers...\n', ...
            parallelSettings.poolType, parallelSettings.maxWorkers);
        poolObj = parpool(parallelSettings.poolType, parallelSettings.maxWorkers);
    end

    parallelMetrics.poolType    = parallelSettings.poolType;
    parallelMetrics.startupTime = toc(setupTimer);
    parRun                      = poolObj.NumWorkers;

    resourceTimer = timer('ExecutionMode', 'fixedRate', 'Period', 60, ...
                          'TimerFcn', @checkSystemResources);
    start(resourceTimer);
    fprintf('Parallel pool ready: %d workers\n', parRun);
    fprintf('================================================\n\n');
else
    parRun = 0;
end
end

function processConfigurations(conn, configIDs, startIdx, totalConfigs, debugCores, ...
    parRun, parallelMetrics, cfg, jobIdentifier, isHPC, dbFile, resourceTimer, poolObj)
remainingConfigs = totalConfigs - startIdx + 1;
fprintf('Processing %d remaining / %d total configurations...\n', remainingConfigs, totalConfigs);
if totalConfigs == 0 || remainingConfigs <= 0
    warning('No configurations to process!');  return;
end
if ~debugCores
    if isempty(gcp('nocreate'))
        error('Parallel pool is no longer available.');
    end
end

chunkSize     = min(1000, remainingConfigs);
startChunkIdx = ceil(startIdx / chunkSize);
numChunks     = ceil(totalConfigs / chunkSize);

for chunkIdx = startChunkIdx:numChunks
    processChunk(conn, configIDs, startIdx, totalConfigs, debugCores, ...
        parRun, parallelMetrics, cfg, jobIdentifier, isHPC, ...
        dbFile, chunkIdx, chunkSize, numChunks, poolObj);
end

if ~debugCores && ~isempty(resourceTimer) && isvalid(resourceTimer)
    stop(resourceTimer);
    delete(resourceTimer);
end
end

function processChunk(conn, configIDs, startIdx, totalConfigs, debugCores, ...
    parRun, parallelMetrics, cfg, jobIdentifier, isHPC, ...
    dbFile, chunkIdx, chunkSize, numChunks, poolObj)
startInChunk  = max(startIdx - (chunkIdx-1)*chunkSize, 1);
endInChunk    = min(chunkIdx*chunkSize, totalConfigs);
chunkStartIdx = (chunkIdx-1)*chunkSize + startInChunk;
chunkEndIdx   = min(endInChunk, length(configIDs));

if chunkStartIdx > chunkEndIdx || chunkStartIdx > length(configIDs)
    warning('Skipping empty chunk %d', chunkIdx);  return;
end

chunkIDs           = configIDs(chunkStartIdx:chunkEndIdx);
checkpointInterval = min(200, length(chunkIDs));

fprintf('Processing chunk %d/%d (configs %d-%d)...\n', ...
    chunkIdx, numChunks, chunkStartIdx, chunkEndIdx);

if debugCores
    processChunkSerial(conn, chunkIDs, chunkStartIdx, totalConfigs, cfg, ...
        checkpointInterval, jobIdentifier);
else
    processChunkParallel(conn, chunkIDs, chunkStartIdx, totalConfigs, ...
        parallelMetrics, cfg, jobIdentifier, isHPC, ...
        dbFile, checkpointInterval, poolObj, startIdx);
end
fprintf('Completed chunk %d/%d\n', chunkIdx, numChunks);
end

function processChunkSerial(conn, chunkIDs, chunkStartIdx, totalConfigs, cfg, ...
    checkpointInterval, jobIdentifier)
tic;
for i = 1:length(chunkIDs)
    localIdx = chunkStartIdx + i - 1;
    configID = chunkIDs(i);
    params   = getConfigParamsDB_minimal(conn, configID);

    localCfg = createWorkerConfig(params, cfg);

    try
        [DATA, beta, VGF, duration, errMadirolas, errCurvature] = ...
            Toolchain_func_v032(localIdx, totalConfigs, localCfg);
        results = struct('beta', beta, 'vgf', VGF, 'duration', duration, ...
            'err_madirolas', errMadirolas, 'err_curvature', errCurvature, ...
            'success', DATA);
        storeResultDB(conn, configID, 0, results);
    catch ME
        warning('MATLAB:PowerLaw:ConfigError', 'Config %d: %s', configID, ME.message);
        storeResultDB(conn, configID, 0, struct('error_message', ME.message, 'success', false));
    end

    if mod(i, checkpointInterval) == 0 || i == length(chunkIDs)
        updateCheckpoint(conn, jobIdentifier, localIdx);
        fprintf('Checkpoint %d/%d  elapsed %s\n', localIdx, totalConfigs, formatTime(toc));
    end
end
end

function processChunkParallel(conn, chunkIDs, chunkStartIdx, totalConfigs, ...
    parallelMetrics, cfg, jobIdentifier, isHPC, ...
    dbFile, checkpointInterval, poolObj, startIdxGlobal)
if ~isfield(parallelMetrics, 'taskTimes'), parallelMetrics.taskTimes = []; end
numSubChunks = ceil(length(chunkIDs) / checkpointInterval);

for subChunkIdx = 1:numSubChunks
    subStart   = (subChunkIdx-1)*checkpointInterval + 1;
    subEnd     = min(subChunkIdx*checkpointInterval, length(chunkIDs));
    subChunkIDs = chunkIDs(subStart:subEnd);
    tic;

    if ~isHPC && strcmp(parallelMetrics.poolType, 'Threads')
        [resultBatch, taskTimings] = processSubChunkThreadPoolBatch(conn, subChunkIDs, ...
            chunkStartIdx, subStart, totalConfigs, dbFile, cfg);
    else
        [resultBatch, taskTimings] = processSubChunkProcessPoolBatch(subChunkIDs, ...
            chunkStartIdx, subStart, totalConfigs, dbFile, cfg);
    end

    storeResultsBatch(conn, resultBatch);
    parallelMetrics.taskTimes = [parallelMetrics.taskTimes, taskTimings];
    batchTime = toc;

    globalIdx = chunkStartIdx + subChunkIdx*checkpointInterval - 1;
    globalIdx = min(globalIdx, totalConfigs);
    updateCheckpoint(conn, jobIdentifier, globalIdx);
    reportSubChunkPerformance(taskTimings, batchTime, globalIdx, totalConfigs, ...
        subChunkIDs, startIdxGlobal);
end
end

function [resultBatch, taskTimings] = processSubChunkThreadPoolBatch(conn, subChunkIDs, ...
    chunkStartIdx, subStartIdx, totalConfigs, dbFile, cfg)
paramsBatch = cell(length(subChunkIDs), 1);
for i = 1:length(subChunkIDs)
    try
        paramsBatch{i} = getConfigParamsDB_minimal(conn, subChunkIDs(i));
    catch ME
        warning('MATLAB:PowerLaw:DBParamError', 'Pre-fetch config %d: %s', subChunkIDs(i), ME.message);
    end
end

resultBatch  = cell(length(subChunkIDs), 1);
taskTimings  = zeros(1, length(subChunkIDs));

parfor i = 1:length(subChunkIDs)
    configID = subChunkIDs(i);
    workerID = labindex;
    localIdx = chunkStartIdx + subStartIdx + i - 2;
    result   = struct('configID', configID, 'workerID', workerID, 'success', false);
    try
        params    = paramsBatch{i};
        if isempty(params), error('Parameters not available'); end
        t0        = tic;
        localCfg  = createWorkerConfig(params, cfg);
        [DATA, beta, VGF, duration, errMadirolas, errCurvature] = ...
            Toolchain_func_v032(localIdx, totalConfigs, localCfg);
        result.success       = DATA;
        result.beta          = beta;
        result.vgf           = VGF;
        result.duration      = duration;
        result.err_madirolas = errMadirolas;
        result.err_curvature = errCurvature;
        result.processing_time = toc(t0);
        taskTimings(i)       = result.processing_time;
    catch ME
        warning('MATLAB:PowerLaw:WorkerError', 'Worker %d config %d: %s', workerID, configID, ME.message);
        result.error_message = ME.message;
    end
    resultBatch{i} = result;
end
resultBatch = resultBatch(~cellfun(@isempty, resultBatch));
end

function [resultBatch, taskTimings] = processSubChunkProcessPoolBatch(subChunkIDs, ...
    chunkStartIdx, subStartIdx, totalConfigs, dbFile, cfg)
resultBatch = cell(length(subChunkIDs), 1);
taskTimings = zeros(1, length(subChunkIDs));

parfor i = 1:length(subChunkIDs)
    configID = subChunkIDs(i);
    workerID = labindex;
    localIdx = chunkStartIdx + subStartIdx + i - 2;
    result   = struct('configID', configID, 'workerID', workerID, 'success', false);
    try
        localConn = sqlite(dbFile);
        params    = getConfigParamsDB_minimal(localConn, configID);
        close(localConn);
        t0        = tic;
        localCfg  = createWorkerConfig(params, cfg);
        [DATA, beta, VGF, duration, errMadirolas, errCurvature] = ...
            Toolchain_func_v032(localIdx, totalConfigs, localCfg);
        result.success       = DATA;
        result.beta          = beta;
        result.vgf           = VGF;
        result.duration      = duration;
        result.err_madirolas = errMadirolas;
        result.err_curvature = errCurvature;
        result.processing_time = toc(t0);
        taskTimings(i)       = result.processing_time;
    catch ME
        warning('MATLAB:PowerLaw:WorkerError', 'Worker %d config %d: %s', workerID, configID, ME.message);
        result.error_message = ME.message;
    end
    resultBatch{i} = result;
end
resultBatch = resultBatch(~cellfun(@isempty, resultBatch));
end

function localCfg = createWorkerConfig(params, cfg)
localCfg              = cfg;
localCfg.ShapeChoice  = double(params.shape_type);
localCfg.fs           = double(params.sampling_rate);
localCfg.powerLaw     = double(params.generated_beta);
localCfg.yGain        = double(params.vgf_value);
localCfg.noiseType    = double(params.noise_type);
localCfg.filterType   = double(params.filter_type);
localCfg.filterParams = params.filter_params;
localCfg.regressType  = double(params.regress_type);
localCfg.TrialNum     = double(params.trial_num);
localCfg.noiseStdDev  = double(params.noise_magnitude) * cfg.pixelScale;
end

function storeResultsBatch(conn, resultBatch)
fprintf('Storing %d results...\n', length(resultBatch));
batchSize  = 50;
numBatches = ceil(length(resultBatch) / batchSize);
hasTransaction = false;
try
    execute(conn, 'BEGIN IMMEDIATE TRANSACTION');
    hasTransaction = true;
catch ME
    if ~contains(ME.message, 'within a transaction')
        warning(ME.identifier, '%s', ME.message);
    end
end
for batchIdx = 1:numBatches
    s = (batchIdx-1)*batchSize + 1;
    e = min(batchIdx*batchSize, length(resultBatch));
    for i = s:e
        r = resultBatch{i};
        try
            storeResultDB(conn, r.configID, r.workerID, r);
        catch ME
            warning('MATLAB:PowerLaw:DBStoreError', 'Config %d: %s', r.configID, ME.message);
        end
    end
end
if hasTransaction
    try
        execute(conn, 'COMMIT');
    catch ME
        warning(ME.identifier, '%s', ME.message);
        try, execute(conn, 'ROLLBACK'); catch, end
    end
end
end

function reportSubChunkPerformance(taskTimings, batchTime, globalIdx, totalConfigs, subChunkIDs, startIdxGlobal)
validT = taskTimings(taskTimings > 0);
if ~isempty(validT)
    completed = globalIdx - startIdxGlobal + 1;
    elapsed   = toc;
    rate      = completed / elapsed;
    remaining = (totalConfigs - globalIdx) / rate;
    fprintf('Progress: %d/%d (%.1f%%)  elapsed %s  remaining %s\n', ...
        globalIdx, totalConfigs, 100*globalIdx/totalConfigs, ...
        formatTime(elapsed), formatTime(remaining));
end
end

function displayPerformanceReport(parallelMetrics, parRun, masterDir, jobIdentifier, isHPC)
fprintf('------------------------------------------------------\n');
fprintf('PARALLEL PERFORMANCE SUMMARY  (%s, %d workers)\n', parallelMetrics.poolType, parRun);
taskTimes = parallelMetrics.taskTimes(parallelMetrics.taskTimes > 0);
if ~isempty(taskTimes)
    fprintf('Tasks: %d  |  mean %.2f s  median %.2f s  min %.2f s  max %.2f s\n', ...
        length(taskTimes), mean(taskTimes), median(taskTimes), min(taskTimes), max(taskTimes));
    efficiency = (sum(taskTimes)/parRun) / toc * 100;
    fprintf('Parallel efficiency: %.1f%%\n', efficiency);
end
fprintf('------------------------------------------------------\n');
end

function markJobComplete(conn, jobIdentifier)
try
    execute(conn, sprintf( ...
        "UPDATE job_checkpoints SET status='completed', last_update_time=CURRENT_TIMESTAMP WHERE job_id='%s'", ...
        jobIdentifier));
    fprintf('Job %s marked complete\n', jobIdentifier);
catch ME
    warning('MATLAB:PowerLaw:CheckpointError', '%s', ME.message);
end
end

function handleFatalError(ME)
fprintf('Fatal error: %s\n', ME.message);
end

function result = ternary(condition, trueValue, falseValue)
if condition, result = trueValue; else, result = falseValue; end
end

function timeStr = formatTime(seconds)
h = floor(seconds/3600);
m = floor(mod(seconds,3600)/60);
s = mod(seconds,60);
if h > 0,      timeStr = sprintf('%dh%dm%.0fs', h, m, s);
elseif m > 0,  timeStr = sprintf('%dm%.0fs', m, s);
else,          timeStr = sprintf('%.1fs', s);
end
end

function checkSystemResources(~, ~)
[u, ~] = memory;
pct = u.MemUsedMATLAB / u.MaxPossibleArrayBytes * 100;
if pct > 75
    fprintf('[Resources] Memory: %.1f%%\n', pct);
    if pct > 85
        warning('MATLAB:PowerLaw:HighMemory', 'Memory high (%.1f%%).', pct);
    end
end
end
