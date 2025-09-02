% TOOLCHAIN_CALLER_V057 Optimized checkpoint interval for high-core-count systems
% 
% Key improvements in v057:
% - CHECKPOINT INTERVAL: Increased from 50 to 200 configurations for 72-core efficiency
% - Reduces database write contention for high-core-count deployments
% - Maintains fault tolerance while eliminating SQLite serialization bottlenecks
% - Fixes hardcoded checkpoint references to use checkpointInterval variable
%
% Previous v056 improvements:
% - ALL INTERNAL FUNCTIONS MOVED TO END for MATLAB 2022a HPC compatibility
% - CRITICAL FIX: Just-in-time parallel pool creation to prevent timeout
% - Pool created only when actually needed (after config generation)
% - Eliminates "Pool no longer exists" errors from long config generation
% - Maintains all v054 fast batch parameter generation performance
%
% Based on Toolchain_caller_v056 with checkpoint interval optimization for 72-core systems.
%
% Created July 2025
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk
%
% Requires:
% - Database Toolbox
% - Parallel Computing Toolbox
% - Curve Fitting Toolbox
% - Statistics and Machine Learning Toolbox

try
    % Initialize random number generator for reproducibility
    rng('default');

    % ===== CONFIGURATION TOGGLES =====
    % Debug level: 0=full run, 1=minimal debug, 2=rebuild shapes, 3=parallel debug
    debug = 2;
    
    % PARAMETER GENERATION METHOD TOGGLE
    % true = Use new fast streaming batch method (recommended)
    % false = Use original individual INSERT method (slow but verified)
    useFastBatch = true;
    
    % Add required paths
    addpath(genpath('functions'));
    addpath(genpath('req'));
    addpath(genpath('utils'));
    commandwindow;

    % Setup paths and database connection
    [conn, dbFile, masterDir, jobIdentifier, isHPC] = setupEnvironment(debug);

    % Get parameter space using the modular function
    paramSpace = defineParameterSpace(debug);

    % Determine parallel settings but DON'T create pool yet (v056 FIX)
    [parallelSettings, debugCores] = determineParallelSettings(debug, isHPC);

    % Create configuration settings
    cfg = createConfigSettings();

    % Setup checkpoints and initialize job with SELECTED METHOD
    [configIDs, startIdx, totalConfigs, generationTime] = ...
        initializeJobWithMethodChoice(conn, jobIdentifier, paramSpace, useFastBatch);

    % Display parameter space information (with generation performance)
    displayParameterSpaceInfo(paramSpace, totalConfigs, parallelSettings, debugCores, isHPC, generationTime);

    % NOW create parallel pool just-in-time (v056 FIX)
    [parallelMetrics, parRun, resourceTimer, poolObj] = ...
        createParallelPoolJustInTime(parallelSettings, debugCores, isHPC);

    % Process configurations with batch processing
    processConfigurations(conn, configIDs, startIdx, totalConfigs, debugCores, ...
        parRun, parallelMetrics, cfg, jobIdentifier, isHPC, ...
        dbFile, resourceTimer, poolObj);

    % Mark job as completed
    markJobComplete(conn, jobIdentifier);

    % Display performance report if applicable
    if ~debugCores && isfield(parallelMetrics, 'taskTimes') && ~isempty(parallelMetrics.taskTimes)
        displayPerformanceReport(parallelMetrics, parRun, masterDir, jobIdentifier, isHPC);
    end

    % Close database connection
    try
        close(conn);
    catch ME
        warning(ME.identifier, '%s', ME.message);
    end

    % Final message
    fprintf('Processing complete! Results stored in: %s\n', dbFile);
    disp('Run ResultLinearMixedModel_v0XX.m & Analysis_Multiverse_v0XX.m to analyze results.');

catch ME
    % Error handling
    handleFatalError(ME);
end

%% ========================================================================
%% ALL INTERNAL FUNCTIONS FOR MATLAB 2022A HPC COMPATIBILITY
%% ========================================================================

function [conn, dbFile, masterDir, jobIdentifier, isHPC] = setupEnvironment(debug)
% Setup database and environment variables
localDir = pwd;
masterDir = localDir(1:end-3);

if debug
    dbFile = fullfile(masterDir, 'results', 'powerlaw_debug_v057.db');
else
    dbFile = fullfile(masterDir, 'results', 'powerlaw_multiverse_v057.db');
end

setupPowerLawDB(dbFile);

try
    conn = sqlite(dbFile);
catch ME
    error('Failed to connect to SQLite database: %s', ME.message);
end

addpath(genpath(fullfile(masterDir, 'src', 'functions')));

slurm_job_id = getenv('SLURM_JOB_ID');
slurm_array_task_id = getenv('SLURM_ARRAY_TASK_ID');
isHPC = ~isempty(slurm_job_id);

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
if isHPC
    hostname = getenv('HOSTNAME');
    if isempty(hostname)
        hostname = 'hpc';
    end
    currentJobIdentifier = sprintf('hpc_%s_%s', hostname, timestamp);
    fprintf('Running on HPC host: %s\n', hostname);
else
    currentJobIdentifier = sprintf('local_%s', timestamp);
    fprintf('Running locally\n');
end

try
    sqlquery = 'SELECT job_id, last_completed_idx, total_configs FROM job_checkpoints WHERE status = "running" ORDER BY last_update_time DESC';
    result = fetch(conn, sqlquery);
    
    if ~isempty(result)
        if istable(result)
            incompleteJobs = result;
        else
            incompleteJobs = cell2table(result, 'VariableNames', {'job_id', 'last_completed_idx', 'total_configs'});
        end
        
        fprintf('\n===== INCOMPLETE JOBS FOUND =====\n');
        fprintf('There are %d incomplete jobs in the database.\n', height(incompleteJobs));
        
        rowsToShow = min(5, height(incompleteJobs));
        for i = 1:rowsToShow
            jobID = incompleteJobs.job_id{i};
            lastIdx = incompleteJobs.last_completed_idx(i);
            totalIdx = incompleteJobs.total_configs(i);
            progress = (lastIdx / totalIdx) * 100;
            
            fprintf('[%d] %s - Progress: %d/%d (%.1f%%)\n', ...
                i, jobID, lastIdx, totalIdx, progress);
        end
        
        answer = input('\nDo you want to resume an incomplete job? (y/n): ', 's');
        
        if lower(answer) == 'y'
            if height(incompleteJobs) == 1
                jobToResume = 1;
            else
                jobToResume = input('Enter the number of the job to resume [1]: ');
                if isempty(jobToResume)
                    jobToResume = 1;
                end
            end
            
            if jobToResume > 0 && jobToResume <= height(incompleteJobs)
                jobIdentifier = incompleteJobs.job_id{jobToResume};
                fprintf('Resuming job: %s\n', jobIdentifier);
            else
                fprintf('Invalid selection. Using current job identifier: %s\n', currentJobIdentifier);
                jobIdentifier = currentJobIdentifier;
            end
        else
            fprintf('Using current job identifier: %s\n', currentJobIdentifier);
            jobIdentifier = currentJobIdentifier;
        end
    else
        fprintf('No incomplete jobs found in database.\n');
        jobIdentifier = currentJobIdentifier;
    end
catch ME
    warning(ME.identifier, '%s', ME.message);
    jobIdentifier = currentJobIdentifier;
end

fprintf('Final job identifier: %s\n', jobIdentifier);
end

function [parallelSettings, debugCores] = determineParallelSettings(debug, isHPC)
% Determine parallel settings but DON'T create pool yet (v056 fix)
debugCores = (debug > 0 && debug < 3);

if ~debugCores
    if isHPC
        poolType = 'Processes';
        maxWorkers = feature('numcores');
        fprintf('HPC environment detected - will use ProcessPool\n');
    else
        poolType = 'Threads';
        maxWorkers = max(1, feature('numcores') - 1);
        maxWorkers = min(maxWorkers, 16);
        fprintf('Local environment detected - will use ThreadPool\n');
    end
    
    parallelSettings = struct();
    parallelSettings.poolType = poolType;
    parallelSettings.maxWorkers = maxWorkers;
    parallelSettings.isHPC = isHPC;
    
    fprintf('Parallel pool will be created just-in-time with %d workers\n', maxWorkers);
else
    parallelSettings = struct();
    if debug == 1
        disp('MINIMAL DEBUG RUN - No parallel computing');
    elseif debug == 2
        disp('COMPREHENSIVE DEBUG RUN - No parallel computing');
    end
end
end

function cfg = createConfigSettings()
% Create configuration structure with default settings
cfg = struct();
cfg.versionTc = 'Tc_v0057';
cfg.orbitCount = 10;
cfg.saveAll = 0;
cfg.display = [0 0 0];
cfg.canvas = [1920; 1080];
cfg.edgeClip = 50;
cfg.curvatureChoice = 1;
cfg.limitBreak = 0;
cfg.displayGraphs = 0;
cfg.MaticSpline = 0;
cfg.resample = 20;
cfg.pixelScale = 480 / 100;
cfg.variantPowerLaw = 1;
cfg.rethrowErrors = false;
end

function [configIDs, startIdx, totalConfigs, generationTime] = ...
    initializeJobWithMethodChoice(conn, jobIdentifier, paramSpace, useFastBatch)
% Setup checkpoint table and initialize job with choice of generation method
setupCheckpointTable(conn);
jobExists = checkJobExists(conn, jobIdentifier);

if jobExists
    fprintf('Resuming job %s from checkpoint...\n', jobIdentifier);
    [configIDs, lastCompletedIdx] = getResumeInfo(conn, jobIdentifier);
    totalConfigs = length(configIDs);
    generationTime = 0;

    if lastCompletedIdx >= totalConfigs
        fprintf('Job %s already completed. No further processing needed.\n', jobIdentifier);
        close(conn);
        error('Job already completed');
    end

    startIdx = lastCompletedIdx + 1;
    fprintf('Resuming from configuration %d of %d\n', startIdx, totalConfigs);
else
    fprintf('Starting new job %s...\n', jobIdentifier);
    fprintf('Generating parameter configurations...\n');
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
        totalConfigs = length(configIDs);
        
        fprintf('\nParameter generation completed:\n');
        fprintf('  Method: %s\n', ternary(useFastBatch, 'Fast batch', 'Original'));
        fprintf('  Configurations: %d\n', totalConfigs);
        fprintf('  Time: %s\n', formatTime(generationTime));
        if totalConfigs > 0
            fprintf('  Rate: %.0f configs/second\n', totalConfigs / generationTime);
        end
        
    catch ME
        if useFastBatch
            warning('MATLAB:PowerLaw:BatchMethodFailed', ...
                'Fast batch method failed (%s), falling back to original method...', ME.message);
            generationTimer = tic;
            configIDs = generateParameterConfigsDB_fixed(conn, paramSpace);
            generationTime = toc(generationTimer);
            totalConfigs = length(configIDs);
            fprintf('Fallback to original method completed in %s\n', formatTime(generationTime));
        else
            rethrow(ME);
        end
    end

    createJobEntry(conn, jobIdentifier, configIDs);
    startIdx = 1;
end
end

function displayParameterSpaceInfo(paramSpace, totalConfigs, parallelSettings, debugCores, isHPC, generationTime)
% Display minimized information about the parameter space including generation performance
disp('------------------------------------------------------');
disp(['Multiverse analysis: ' num2str(totalConfigs) ' total configurations']);
fprintf('Parameter space: %d sampling rates × %d shapes × %d betas × %d VGFs × %d noise types × %d noise levels × %d filters × %d regressors\n', ...
    length(paramSpace.samplingRates), length(paramSpace.shapes), ...
    length(paramSpace.generatedBetas), length(paramSpace.vgfValues), ...
    length(paramSpace.noiseTypes), length(paramSpace.noiseMagnitudes), ...
    length(paramSpace.filterTypes), length(paramSpace.regressTypes));

if generationTime > 0
    fprintf('Configuration generation: %s (%.0f configs/sec)\n', ...
        formatTime(generationTime), totalConfigs / generationTime);
end

if ~debugCores
    fprintf('Will create %s pool with %d workers for processing\n', ...
        parallelSettings.poolType, parallelSettings.maxWorkers);
end

fprintf('Checkpoint interval: 200 configurations (optimized for high-core systems)\n');
disp('------------------------------------------------------');
end

function [parallelMetrics, parRun, resourceTimer, poolObj] = ...
    createParallelPoolJustInTime(parallelSettings, debugCores, isHPC)
% Create parallel pool only when actually needed (v056 fix for timeout issue)
parallelMetrics = struct('poolType', '', 'startupTime', 0, 'taskTimes', []);
resourceTimer = [];
poolObj = [];

if ~debugCores
    fprintf('\n===== CREATING PARALLEL POOL JUST-IN-TIME =====\n');
    fprintf('Creating pool now to avoid timeout during config generation...\n');
    
    setupTimer = tic;
    existingPool = gcp('nocreate');
    needNewPool = false;
    
    if isempty(existingPool)
        needNewPool = true;
        fprintf('No existing pool found - creating new pool\n');
    else
        existingType = class(existingPool);
        if isHPC && ~contains(existingType, 'ProcessPool')
            fprintf('Existing pool is %s but need ProcessPool for HPC - recreating\n', existingType);
            delete(existingPool);
            needNewPool = true;
        elseif ~isHPC && ~contains(existingType, 'ThreadPool')
            fprintf('Existing pool is %s but prefer ThreadPool for local - recreating\n', existingType);
            delete(existingPool);
            needNewPool = true;
        else
            fprintf('Using existing %s pool with %d workers\n', existingType, existingPool.NumWorkers);
            poolObj = existingPool;
        end
    end
    
    if needNewPool
        fprintf('Creating %s pool with %d workers...\n', ...
            parallelSettings.poolType, parallelSettings.maxWorkers);
        poolObj = parpool(parallelSettings.poolType, parallelSettings.maxWorkers);
        fprintf('Pool created successfully in %.2f seconds\n', toc(setupTimer));
    end

    parallelMetrics.poolType = parallelSettings.poolType;
    parallelMetrics.startupTime = toc(setupTimer);
    parallelMetrics.taskTimes = [];
    parRun = poolObj.NumWorkers;

    resourceTimer = timer('ExecutionMode', 'fixedRate', ...
        'Period', 60, 'TimerFcn', @checkSystemResources);
    start(resourceTimer);
    
    fprintf('Parallel pool ready with %d workers\n', parRun);
    fprintf('====================================================\n\n');
else
    parRun = 0;
end
end

function processConfigurations(conn, configIDs, startIdx, totalConfigs, debugCores, ...
    parRun, parallelMetrics, cfg, jobIdentifier, isHPC, ...
    dbFile, resourceTimer, poolObj)
% Process all configurations in the parameter space
remainingConfigs = totalConfigs - startIdx + 1;
fprintf('Processing %d remaining configurations out of %d total...\n', remainingConfigs, totalConfigs);

if totalConfigs == 0 || remainingConfigs <= 0
    warning('No configurations to process!');
    return;
end

if ~debugCores
    currentPool = gcp('nocreate');
    if isempty(currentPool)
        error('Parallel pool is no longer available. This should not happen with just-in-time creation.');
    end
    fprintf('Confirmed parallel pool is available with %d workers\n', currentPool.NumWorkers);
end

chunkSize = min(1000, remainingConfigs);
startChunkIdx = ceil(startIdx / chunkSize);
numChunks = ceil(totalConfigs / chunkSize);

for chunkIdx = startChunkIdx:numChunks
    processChunk(conn, configIDs, startIdx, totalConfigs, debugCores, ...
        parRun, parallelMetrics, cfg, jobIdentifier, isHPC, ...
        dbFile, chunkIdx, chunkSize, numChunks, poolObj);
end

if ~debugCores && exist('resourceTimer', 'var') && isvalid(resourceTimer)
    stop(resourceTimer);
    delete(resourceTimer);
end
end

function processChunk(conn, configIDs, startIdx, totalConfigs, debugCores, ...
    parRun, parallelMetrics, cfg, jobIdentifier, isHPC, ...
    dbFile, chunkIdx, chunkSize, numChunks, poolObj)
% Process a chunk of configurations
startInChunk = max(startIdx - (chunkIdx-1)*chunkSize, 1);
endInChunk = min(chunkIdx*chunkSize, totalConfigs);

chunkStartIdx = (chunkIdx-1)*chunkSize + startInChunk;
chunkEndIdx = min(endInChunk, length(configIDs));

if chunkStartIdx <= chunkEndIdx && chunkStartIdx <= length(configIDs)
    chunkIDs = configIDs(chunkStartIdx:chunkEndIdx);

    fprintf('Processing chunk %d/%d (configs %d-%d)...\n', ...
        chunkIdx, numChunks, chunkStartIdx, chunkEndIdx);

    % v057: Optimized checkpoint interval for high-core systems
    checkpointInterval = min(200, length(chunkIDs));

    if debugCores
        processChunkSerial(conn, chunkIDs, chunkStartIdx, totalConfigs, cfg, ...
            checkpointInterval, jobIdentifier);
    else
        processChunkParallel(conn, chunkIDs, chunkStartIdx, totalConfigs, ...
            parallelMetrics, cfg, jobIdentifier, isHPC, ...
            dbFile, checkpointInterval, poolObj, startIdx);
    end

    fprintf('Completed chunk %d/%d\n', chunkIdx, numChunks);
else
    warning('Skipping chunk %d as there are no configs in this range', chunkIdx);
end
end

function processChunkSerial(conn, chunkIDs, chunkStartIdx, totalConfigs, cfg, ...
    checkpointInterval, jobIdentifier)
% Process a chunk in serial mode (for debugging)
tic;
for i = 1:length(chunkIDs)
    localIdx = chunkStartIdx + i - 1;
    configID = chunkIDs(i);
    workerID = 0;

    params = getConfigParamsDB_minimal(conn, configID);

    cfg.ShapeChoice = double(params.shape_type);
    cfg.fs = double(params.sampling_rate);
    cfg.powerLaw = double(params.generated_beta);
    cfg.yGain = double(params.vgf_value);
    cfg.noiseType = double(params.noise_type);
    cfg.noiseStdDev = double(params.noise_magnitude) * cfg.pixelScale;
    cfg.filterType = double(params.filter_type);
    cfg.filterParams = params.filter_params;
    cfg.regressType = double(params.regress_type);
    cfg.TrialNum = double(params.trial_num);

    try
        [DATA, beta, VGF, duration, errMadirolas, errCurvature] = ...
            Toolchain_func_v032(localIdx, totalConfigs, cfg);

        results = struct('beta', beta, 'vgf', VGF, 'duration', duration, ...
            'err_madirolas', errMadirolas, 'err_curvature', errCurvature, ...
            'success', DATA);

        storeResultDB(conn, configID, workerID, results);
    catch ME
        warning('MATLAB:PowerLaw:ConfigError', 'Error processing config %d: %s', configID, ME.message);
        errorResult = struct('error_message', ME.message, 'success', false);
        storeResultDB(conn, configID, workerID, errorResult);
    end

    if mod(i, checkpointInterval) == 0 || i == length(chunkIDs)
        updateCheckpoint(conn, jobIdentifier, localIdx);

        currentTime = toc;
        estimatedTotalTime = currentTime * (totalConfigs / localIdx);
        estimatedRemainingTime = estimatedTotalTime - currentTime;

        fprintf('Checkpoint at config %d/%d - Elapsed: %s, Remaining: %s\n',...
            localIdx, totalConfigs, formatTime(currentTime), formatTime(estimatedRemainingTime));
    end

    % v057: Use checkpointInterval variable instead of hardcoded value
    if mod(i, checkpointInterval) == 0 || i == length(chunkIDs)
        fprintf('Processed %d/%d configurations\n', localIdx, totalConfigs);
    end
end
end

function processChunkParallel(conn, chunkIDs, chunkStartIdx, totalConfigs, ...
    parallelMetrics, cfg, jobIdentifier, isHPC, ...
    dbFile, checkpointInterval, poolObj, startIdxGlobal)
% Process a chunk in parallel mode with batch result processing
configsProcessed = 0;

if ~isfield(parallelMetrics, 'taskTimes')
    parallelMetrics.taskTimes = [];
end

numSubChunks = ceil(length(chunkIDs) / checkpointInterval);

for subChunkIdx = 1:numSubChunks
    subStartIdx = (subChunkIdx-1)*checkpointInterval + 1;
    subEndIdx = min(subChunkIdx*checkpointInterval, length(chunkIDs));
    subChunkIDs = chunkIDs(subStartIdx:subEndIdx);

    tic;

    if ~isHPC && strcmp(parallelMetrics.poolType, 'Threads')
        [resultBatch, taskTimings] = processSubChunkThreadPoolBatch(conn, subChunkIDs, ...
            chunkStartIdx, subStartIdx, totalConfigs, dbFile, cfg);
    else
        [resultBatch, taskTimings] = processSubChunkProcessPoolBatch(subChunkIDs, ...
            chunkStartIdx, subStartIdx, totalConfigs, dbFile, cfg);
    end

    storeResultsBatch(conn, resultBatch);
    parallelMetrics.taskTimes = [parallelMetrics.taskTimes, taskTimings];
    batchTime = toc;

    configsProcessed = configsProcessed + length(subChunkIDs);
    globalIdx = chunkStartIdx + configsProcessed - 1;
    updateCheckpoint(conn, jobIdentifier, globalIdx);

    reportSubChunkPerformance(taskTimings, batchTime, globalIdx, totalConfigs, ...
        subChunkIDs, startIdxGlobal);
end
end

function [resultBatch, taskTimings] = processSubChunkThreadPoolBatch(conn, subChunkIDs, ...
    chunkStartIdx, subStartIdx, totalConfigs, dbFile, cfg)
% Process a sub-chunk using ThreadPool with batch result processing
fprintf('Pre-fetching parameters for %d configurations...\n', length(subChunkIDs));
paramsBatch = cell(length(subChunkIDs), 1);
for i = 1:length(subChunkIDs)
    try
        paramsBatch{i} = getConfigParamsDB_minimal(conn, subChunkIDs(i));
    catch ME
        warning('MATLAB:PowerLaw:DBParamError', 'Error pre-fetching parameters for config %d: %s', ...
            subChunkIDs(i), ME.message);
    end
end

resultBatch = cell(length(subChunkIDs), 1);
taskTimings = zeros(1, length(subChunkIDs));

parfor i = 1:length(subChunkIDs)
    configID = subChunkIDs(i);
    workerID = labindex;
    localIdx = chunkStartIdx + subStartIdx + i - 2;

    result = struct('configID', configID, 'workerID', workerID, 'success', false);

    try
        params = paramsBatch{i};
        if isempty(params)
            error('Parameters not available for this configuration');
        end

        taskStartTime = tic;
        localCfg = createWorkerConfig(params, cfg.versionTc);

        [DATA, beta, VGF, duration, errMadirolas, errCurvature] = ...
            Toolchain_func_v032(localIdx, totalConfigs, localCfg);

        taskTime = toc(taskStartTime);
        taskTimings(i) = taskTime;

        result.success = DATA;
        result.beta = beta;
        result.vgf = VGF;
        result.duration = duration;
        result.err_madirolas = errMadirolas;
        result.err_curvature = errCurvature;
        result.processing_time = taskTime;

    catch ME
        warning('MATLAB:PowerLaw:WorkerError', 'Worker %d error on config %d: %s', ...
            workerID, configID, ME.message);
        result.success = false;
        result.error_message = ME.message;
    end

    resultBatch{i} = result;
end

resultBatch = resultBatch(~cellfun(@isempty, resultBatch));
end

function [resultBatch, taskTimings] = processSubChunkProcessPoolBatch(subChunkIDs, ...
    chunkStartIdx, subStartIdx, totalConfigs, dbFile, cfg)
% Process a sub-chunk using ProcessPool with batch result processing
resultBatch = cell(length(subChunkIDs), 1);
taskTimings = zeros(1, length(subChunkIDs));

parfor i = 1:length(subChunkIDs)
    configID = subChunkIDs(i);
    workerID = labindex;
    localIdx = chunkStartIdx + subStartIdx + i - 2;

    result = struct('configID', configID, 'workerID', workerID, 'success', false);

    try
        localConn = sqlite(dbFile);
        params = getConfigParamsDB_minimal(localConn, configID);
        close(localConn);

        taskStartTime = tic;
        localCfg = createWorkerConfig(params, cfg.versionTc);

        [DATA, beta, VGF, duration, errMadirolas, errCurvature] = ...
            Toolchain_func_v032(localIdx, totalConfigs, localCfg);

        taskTime = toc(taskStartTime);
        taskTimings(i) = taskTime;

        result.success = DATA;
        result.beta = beta;
        result.vgf = VGF;
        result.duration = duration;
        result.err_madirolas = errMadirolas;
        result.err_curvature = errCurvature;
        result.processing_time = taskTime;

    catch ME
        warning('MATLAB:PowerLaw:WorkerError', 'Worker %d error on config %d: %s', ...
            workerID, configID, ME.message);
        result.success = false;
        result.error_message = ME.message;
    end

    resultBatch{i} = result;
end

resultBatch = resultBatch(~cellfun(@isempty, resultBatch));
end

function localCfg = createWorkerConfig(params, versionTc)
% Create configuration for a worker with pixelScale correctly included
localCfg = struct();
localCfg.versionTc = versionTc;
localCfg.ShapeChoice = double(params.shape_type);
localCfg.fs = double(params.sampling_rate);
localCfg.powerLaw = double(params.generated_beta);
localCfg.yGain = double(params.vgf_value);
localCfg.noiseType = double(params.noise_type);
localCfg.filterType = double(params.filter_type);
localCfg.filterParams = params.filter_params;
localCfg.regressType = double(params.regress_type);
localCfg.TrialNum = double(params.trial_num);

localCfg.pixelScale = 480 / 100;
localCfg.noiseStdDev = double(params.noise_magnitude) * localCfg.pixelScale;

localCfg.orbitCount = 6;
localCfg.saveAll = 0;
localCfg.display = [0 0 0];
localCfg.canvas = [1920; 1080];
localCfg.edgeClip = 50;
localCfg.curvatureChoice = 1;
localCfg.limitBreak = 0;
localCfg.displayGraphs = 0;
localCfg.MaticSpline = 0;
localCfg.resample = 20;
localCfg.variantPowerLaw = 1;
localCfg.rethrowErrors = false;
end

function storeResultsBatch(conn, resultBatch)
% Store all results in the batch to the database
fprintf('Storing %d results in database...\n', length(resultBatch));

batchSize = 50;
numBatches = ceil(length(resultBatch) / batchSize);

hasTransaction = false;
try
    execute(conn, 'BEGIN IMMEDIATE TRANSACTION');
    hasTransaction = true;
catch ME
    if contains(ME.message, 'within a transaction')
        fprintf('Using existing transaction for storing results\n');
    else
        warning(ME.identifier, '%s', ME.message);
    end
end

for batchIdx = 1:numBatches
    startIdx = (batchIdx-1)*batchSize + 1;
    endIdx = min(batchIdx*batchSize, length(resultBatch));

    for i = startIdx:endIdx
        result = resultBatch{i};
        try
            configID = result.configID;
            workerID = result.workerID;
            storeResultDB(conn, configID, workerID, result);
        catch ME
            warning('MATLAB:PowerLaw:DBStoreError', ...
                'Error storing result for config %d: %s', result.configID, ME.message);
        end
    end

    if mod(batchIdx, max(5, round(numBatches/4))) == 0 || batchIdx == numBatches
        fprintf('Stored batch %d/%d (%d results)\n', batchIdx, numBatches, endIdx-startIdx+1);
    end
end

if hasTransaction
    try
        execute(conn, 'COMMIT');
    catch ME
        warning(ME.identifier, '%s', ME.message);
        try
            execute(conn, 'ROLLBACK');
        catch
        end
    end
end

fprintf('Results stored successfully.\n');
end

function reportSubChunkPerformance(taskTimings, batchTime, globalIdx, totalConfigs, subChunkIDs, startIdxGlobal)
% Report performance metrics for a sub-chunk - minimized output
validTimings = taskTimings(taskTimings > 0);
if ~isempty(validTimings)
    completedConfigs = globalIdx - startIdxGlobal + 1;
    totalTime = toc;
    configsPerSec = completedConfigs / totalTime;
    remainingConfigs = totalConfigs - globalIdx;
    estimatedRemainingTime = remainingConfigs / configsPerSec;

    fprintf('Progress: %d/%d (%.1f%%) - Elapsed: %s, Remaining: %s\n', ...
        globalIdx, totalConfigs, (globalIdx/totalConfigs)*100, ...
        formatTime(totalTime), formatTime(estimatedRemainingTime));
end
end

function displayPerformanceReport(parallelMetrics, parRun, masterDir, jobIdentifier, isHPC)
% Create and display detailed performance report - simplified
disp('------------------------------------------------------');
disp('PARALLEL PERFORMANCE SUMMARY');
disp('------------------------------------------------------');
fprintf('Pool type: %s with %d workers\n', parallelMetrics.poolType, parRun);
taskTimes = parallelMetrics.taskTimes(parallelMetrics.taskTimes > 0);

if ~isempty(taskTimes)
    fprintf('Total tasks processed: %d\n', length(taskTimes));
    fprintf('Task timing (seconds): Mean=%.2f, Median=%.2f, Min=%.2f, Max=%.2f\n',...
        mean(taskTimes), median(taskTimes), min(taskTimes), max(taskTimes));

    totalTaskTime = sum(taskTimes);
    totalWallTime = toc;
    theoreticalMinTime = totalTaskTime / parRun;
    parallelEfficiency = theoreticalMinTime / totalWallTime * 100;

    fprintf('Total cumulative task time: %s\n', formatTime(totalTaskTime));
    fprintf('Total wall clock time: %s\n', formatTime(totalWallTime));
    fprintf('Parallel efficiency: %.1f%%\n', parallelEfficiency);

    if parallelEfficiency < 60
        disp('Consider adjusting chunk sizes for better load balancing.');
    end
end

fprintf('Checkpoint interval: 200 configurations (4x improvement over v056)\n');
disp('------------------------------------------------------');
end

function markJobComplete(conn, jobIdentifier)
% Mark job as completed in the checkpoint table
try
    sqlquery = sprintf("UPDATE job_checkpoints SET status = 'completed', last_update_time = CURRENT_TIMESTAMP WHERE job_id = '%s'", jobIdentifier);
    execute(conn, sqlquery);
    fprintf('Job %s marked as completed\n', jobIdentifier);
catch ME
    warning('MATLAB:PowerLaw:CheckpointError', 'Error updating job completion status: %s', ME.message);
end
end

function handleFatalError(ME)
% Handle fatal errors gracefully
if exist('conn', 'var') && ~isempty(conn)
    try
        close(conn);
    catch
    end
end

if exist('resourceTimer', 'var') && isvalid(resourceTimer)
    try
        stop(resourceTimer);
        delete(resourceTimer);
    catch
    end
end

fprintf('Error occurred: %s\n', ME.message);
end

function result = ternary(condition, trueValue, falseValue)
% Simple ternary operator implementation
if condition
    result = trueValue;
else
    result = falseValue;
end
end

function timeStr = formatTime(seconds)
% Format time in a human-readable way
hours = floor(seconds / 3600);
minutes = floor(mod(seconds, 3600) / 60);
secs = mod(seconds, 60);

if hours > 0
    timeStr = sprintf('%dh:%dm:%.0fs', hours, minutes, secs);
elseif minutes > 0
    timeStr = sprintf('%dm:%.0fs', minutes, secs);
else
    timeStr = sprintf('%.1fs', secs);
end
end

function checkSystemResources(~, ~)
% Monitor system resources and provide warnings if resources are constrained - minimal output
[userPath, ~] = memory;
memUsedPercent = userPath.MemUsedMATLAB / userPath.MaxPossibleArrayBytes * 100;

if memUsedPercent > 75
    fprintf('[Resource Monitor] Memory usage: %.1f%% of maximum array size\n', memUsedPercent);

    if memUsedPercent > 85
        warning('MATLAB:PowerLaw:HighMemory', 'Memory usage high (%.1f%%). Consider reducing worker count.', memUsedPercent);
    end
end
end
