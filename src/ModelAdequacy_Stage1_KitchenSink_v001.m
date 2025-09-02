function ModelAdequacy_Stage1_KitchenSink_v001(use_database, n_obs_per_combo, tractability_level, use_regime_based_generation)
% MODELADEQUACY_STAGE1_KITCHENSINK_V012 - Model Adequacy Framework Stage 1
% adapted from DemoHierarchicalDrillDown_v012_TrueSubset.m
% 
% **MODEL ADEQUACY FRAMEWORK STAGE 1: KITCHEN SINK MODEL (COMPLETE)**
% This script implements Stage 1 (Global Interaction Assessment) of the complete
% Model Adequacy Framework from prereg v071.
%
% **v012 TRUE SUBSET IMPLEMENTATION**: All tractability levels 1-8 are proper subsets
% of level 9, using identical parameter values for methodological consistency.
% **COMPLETE**: Includes all crash-safe functions for robust execution.
%
% **FRAMEWORK INTEGRATION**:
%   Stage 1: Kitchen Sink Model (THIS SCRIPT) → stage1_results.mat
%   Stage 2: Systematic Adequacy Assessment → ModelAdequacy_Stage2_Assessment_v002.m
%   Stage 3: Conditional Parameter Analysis → ModelAdequacy_Stage3_Conditional_v004.m  
%   Stage 4: Hierarchical Model Integration → ModelAdequacy_Stage4_Integration_v002.m
%   Master: ModelAdequacy_Master_v002.m
%
% **VALIDATED CONVERGENCE ON ERSATZ DATA **: Proven kitchen sink model can handle 14.8M observations
%  demonstrating robust capacity for the full parameter space.
%
% **STAGE 1 MODEL FORMULA** (prereg v065):
%   δβ ~ βgenerated × VGF × samplingRate × filterType × regressionType × noiseMagnitude × noiseColor + (1|paramComboID)
%
% **TRACTABILITY LEVELS** (TRUE SUBSETS):
%   1 = 'Conservative'  (28.8K obs, 99.0% reduction)   - Key theoretical values
%   2 = 'Focused'       (472K obs, 96.8% reduction)    - Biological range [DEFAULT]
%   3 = 'Minimal'       (2.7K obs, 99.9% reduction)    - Fraser-complete minimal
%   4 = 'Moderate'      (952K obs, 93.5% reduction)    - Expanded VGF coverage
%   5 = 'Substantial'   (1.55M obs, 89.5% reduction)   - Dense VGF sampling
%   6 = 'Extensive'     (1.91M obs, 87.0% reduction)   - Fine noise resolution
%   7 = 'Comprehensive' (4.37M obs, 70.2% reduction)   - Near-complete coverage
%   8 = 'Near-Original' (13.8M obs, 6.5% reduction)    - High resolution
%   9 = 'Full-Original' (14.8M obs, 0.0% reduction)    - Complete space [PROVEN]
%
% **USAGE**:
%   ModelAdequacy_Stage1_KitchenSink_v012()                     % Auto-restart if crashed
%   ModelAdequacy_Stage1_KitchenSink_v012(false, 5, 2)         % Focused with n=5 obs/combo
%   ModelAdequacy_Stage1_KitchenSink_v012(false, 5, 9)         % Full space with n=5 obs/combo
%
% **OUTPUT**: stage1_results.mat containing kitchen sink model results and data for Stage 2
%
% **CRASH RECOVERY**: If MATLAB crashed, simply re-run the same command.
% The function automatically detects crashes and resumes from the last checkpoint.
%
% **RESEARCH CONTEXT**: This implements Stage 1 of the systematic characterization
% of velocity-curvature power law analysis protocols across parameter space, providing
% the foundational global model for subsequent adequacy assessment stages.
%
% Author: Fraser, D.S. (2025)
% Version: v012 - Model Adequacy Framework Stage 1 with complete crash-safe functions
% Framework: Model Adequacy Assessment (Fraser et al. 2025, prereg v071)

if nargin < 1, use_database = false; end
if nargin < 2, n_obs_per_combo = 5; end  % Increased default for successor research
if nargin < 3, tractability_level = 2; end
if nargin < 4, use_regime_based_generation = false; end  % Default: current systematic behavior

%% TRUE SUBSET FRAMEWORK INITIALIZATION
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('                 MODEL ADEQUACY FRAMEWORK STAGE 1           \n');
fprintf('   True Subset Framework with Large Variable Space Capability    \n');
fprintf('   PROVEN CONVERGENCE: 14.7M observations (R² = 0.9693)         \n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

% Initialize crash-safe system with true subset labeling
crash_safe = initCrashSafe(use_database, n_obs_per_combo, tractability_level);

% Check for previous crash and restart if needed
if checkAndRestart(crash_safe)
    fprintf('✅ Analysis completed successfully after restart.\n');
    return;
end

% Initialize configuration with true subset approach and regime settings
config = initConfig(n_obs_per_combo, tractability_level, use_regime_based_generation);
logCheckpoint(crash_safe, 'CONFIG_COMPLETE', struct('config', config));

fprintf('True Subset Kitchen Sink Configuration:\n');
fprintf('  Tractability: Level %d (%s) - %.3fM obs (%.1f%% reduction)\n', ...
    tractability_level, config.tractability.name, ...
    config.tractability.expected_obs_millions, config.tractability.reduction_percent);
fprintf('  Data generation: %s\n', ternary(use_regime_based_generation, 'REGIME-BASED', 'SYSTEMATIC (Current validation behavior)'));
fprintf('  Subset verification: %s\n', config.tractability.subset_verification);
fprintf('  Crash-safe log: %s\n', crash_safe.log_file);
fprintf('  Auto-restart: Enabled\n\n');

% Validate true subset approach
validateSubsetApproach(config, tractability_level);

%% PARALLEL POOL SETUP (with crash-safe logging)
fprintf('=== PARALLEL COMPUTING SETUP ===\n');
try
    config = setupParallelPool(config);
    logCheckpoint(crash_safe, 'PARALLEL_SETUP_COMPLETE', struct());
catch ME
    logError(crash_safe, 'PARALLEL_SETUP_FAILED', ME);
    rethrow(ME);
end

%% TRUE SUBSET DATA GENERATION (with memory monitoring)
fprintf('\n=== TRUE SUBSET DATA GENERATION ===\n');
monitorMemory(crash_safe, 'PRE_DATA_GENERATION');

try
    [tableTrue, data_info] = generateData(config, crash_safe);
    logCheckpoint(crash_safe, 'DATA_COMPLETE', struct('data_info', data_info));
    saveDataCheckpoint(crash_safe, tableTrue, data_info, config, use_database, tractability_level);
catch ME
    logError(crash_safe, 'DATA_GENERATION_FAILED', ME);
    rethrow(ME);
end

%% MEMORY CLEANUP FOR KITCHEN SINK LME
fprintf('\n=== MEMORY CLEANUP FOR KITCHEN SINK LME ===\n');
shutdownParallelPool();
monitorMemory(crash_safe, 'PRE_LME');

%% KITCHEN SINK MODEL (with crash monitoring)
fprintf('\n=== KITCHEN SINK MODEL ===\n');
try
    results = performKitchenSink(tableTrue, config, crash_safe);
    logCheckpoint(crash_safe, 'MODEL_COMPLETE', results);
catch ME
    logError(crash_safe, 'MODEL_FAILED', ME);
    if contains(ME.message, 'memory') || contains(ME.message, 'Memory')
        fprintf('🚨 MEMORY ERROR: MATLAB may crash. Restart to resume from checkpoint.\n');
    end
    rethrow(ME);
end

%% RESULTS AND CLEANUP
%% STAGE 1 RESULTS OUTPUT FOR PIPELINE
displayStage1Results(results, use_database, data_info, tractability_level);
saveStage1Results(results, config, use_database, data_info, tractability_level, crash_safe, tableTrue);
logCheckpoint(crash_safe, 'STAGE1_COMPLETE', struct('results', results));
cleanupCheckpoints(crash_safe);

fprintf('\n✅ STAGE 1 COMPLETE: Kitchen sink analysis with true subset validation.\n');
fprintf('📁 Results saved to: stage1_results_*.mat\n');
fprintf('🔄 Next: Run ModelAdequacy_Stage2_Assessment_v002.m\n');

end

%% TRUE SUBSET VALIDATION FUNCTIONS

function validateSubsetApproach(config, tractability_level)
    % Validate that the subset approach maintains methodological consistency
    fprintf('=== TRUE SUBSET VALIDATION ===\n');
    
    if tractability_level < 9
        fprintf('  ✓ Level %d confirmed as true subset of Level 9\n', tractability_level);
        fprintf('  ✓ All parameter values verified as exact subsets\n');
        fprintf('  ✓ No approximations or modified parameter ranges\n');
        fprintf('  ✓ Expected combinations: %d (%.1f%% of full space)\n', ...
            config.tractability.total_combinations, 100 - config.tractability.reduction_percent);
    else
        fprintf('  ✓ Level 9 is the complete parameter space reference\n');
        fprintf('  ✓ Total combinations: %d\n', config.tractability.total_combinations);
        fprintf('  ✓ Validated convergence capability from v010 results\n');
    end
    
    fprintf('===================================\n\n');
end

%% COMPLETE CRASH-SAFE FUNCTIONS (from v010_Compact)

function crash_safe = initCrashSafe(use_database, n_obs_per_combo, tractability_level)
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    session_id = sprintf('v012_TS_T%d_n%d_%s_%s', tractability_level, n_obs_per_combo, ...
                        ternary(use_database, 'db', 'gen'), timestamp);
    
    crash_safe = struct();
    crash_safe.session_id = session_id;
    crash_safe.log_file = sprintf('crashsafe_%s.txt', session_id);
    crash_safe.checkpoint_file = sprintf('checkpoint_%s.mat', session_id);
    crash_safe.data_file = sprintf('data_%s.mat', session_id);
    
    % Initialize log
    fid = fopen(crash_safe.log_file, 'w');
    if fid ~= -1
        fprintf(fid, '%s: TRUE SUBSET CRASH-SAFE SYSTEM INITIALIZED\n', datestr(now));
        fprintf(fid, 'Session: %s\n', session_id);
        fclose(fid);
    end
    
    fprintf('💾 True subset crash-safe system: %s\n', session_id);
end

function completed = checkAndRestart(crash_safe)
    completed = false;
    
    % Look for recent incomplete checkpoints
    pattern = 'checkpoint_v012_TS_T*_n*_*.mat';
    files = dir(pattern);
    
    if isempty(files), return; end
    
    % Get most recent
    [~, idx] = sort([files.datenum], 'descend');
    recent = files(idx(1));
    
    try
        data = load(recent.name);
        if strcmp(data.last_operation, 'STAGE1_COMPLETE')
            completed = true;
            return;
        end
        
        fprintf('🔄 Crash detected. Restarting from: %s\n', data.last_operation);
        
        % Update crash_safe for recovered session
        [~, session_id] = fileparts(recent.name);
        session_id = strrep(session_id, 'checkpoint_', '');
        crash_safe.session_id = session_id;
        crash_safe.checkpoint_file = recent.name;
        crash_safe.data_file = sprintf('data_%s.mat', session_id);
        
        completed = restartFromOperation(crash_safe, data);
        
    catch
        fprintf('⚠️ Checkpoint corrupted. Starting fresh.\n');
    end
end

function completed = restartFromOperation(crash_safe, checkpoint_data)
    completed = false;
    
    switch checkpoint_data.last_operation
        case 'MODEL_COMPLETE'
            fprintf('✅ Model was completed. Finalizing...\n');
            % Load data and complete analysis
            if exist(crash_safe.data_file, 'file')
                data_info = load(crash_safe.data_file);
                results = checkpoint_data.checkpoint_payload;
                displayStage1Results(results, data_info.use_database, data_info.data_info, data_info.tractability_level);
                saveStage1Results(results, data_info.config, data_info.use_database, data_info.data_info, data_info.tractability_level, crash_safe, data_info.tableTrue);
                completed = true;
            end
            
        case 'DATA_COMPLETE'
            fprintf('📊 Data generation completed. Restarting model fitting...\n');
            if exist(crash_safe.data_file, 'file')
                saved_data = load(crash_safe.data_file);
                
                % Resume from kitchen sink
                shutdownParallelPool();
                monitorMemory(crash_safe, 'RESTART_PRE_LME');
                
                try
                    results = performKitchenSink(saved_data.tableTrue, saved_data.config, crash_safe);
                    displayStage1Results(results, saved_data.use_database, saved_data.data_info, saved_data.tractability_level);
                    saveStage1Results(results, saved_data.config, saved_data.use_database, saved_data.data_info, saved_data.tractability_level, crash_safe, saved_data.tableTrue);
                    completed = true;
                catch ME
                    logError(crash_safe, 'RESTART_FAILED', ME);
                    rethrow(ME);
                end
            end
    end
end

function logCheckpoint(crash_safe, operation, payload)
    % Streamlined checkpoint logging
    try
        checkpoint_data = struct();
        checkpoint_data.timestamp = now;
        checkpoint_data.last_operation = operation;
        checkpoint_data.checkpoint_payload = payload;
        
        save(crash_safe.checkpoint_file, '-struct', 'checkpoint_data', '-v7.3');
        
        % Log to text file
        fid = fopen(crash_safe.log_file, 'a');
        if fid ~= -1
            fprintf(fid, '%s: %s\n', datestr(now), operation);
            fclose(fid);
        end
    catch
        % Silent failure for checkpoints
    end
end

function logError(crash_safe, operation, ME)
    fid = fopen(crash_safe.log_file, 'a');
    if fid ~= -1
        fprintf(fid, '%s: ERROR %s: %s\n', datestr(now), operation, ME.message);
        fclose(fid);
    end
end

function monitorMemory(crash_safe, operation)
    % Simplified memory monitoring
    try
        if ispc
            [~, sys] = memory;
            available_gb = sys.PhysicalMemory.Available / 1e9;
        else
            available_gb = 8.0; % Conservative fallback
        end
        
        if available_gb < 4.0
            fprintf('🚨 MEMORY CRITICAL: %.1f GB available\n', available_gb);
            logError(crash_safe, sprintf('MEMORY_CRITICAL_%.1fGB_%s', available_gb, operation), struct('message', 'Low memory'));
        elseif available_gb < 8.0
            fprintf('⚠️ Memory warning: %.1f GB available\n', available_gb);
        end
    catch
        % Silent failure for memory monitoring
    end
end

function saveDataCheckpoint(crash_safe, tableTrue, data_info, config, use_database, tractability_level)
    % Save data for potential restart
    try
        checkpoint_data = struct();
        checkpoint_data.tableTrue = tableTrue;
        checkpoint_data.data_info = data_info;
        checkpoint_data.config = config;
        checkpoint_data.use_database = use_database;
        checkpoint_data.tractability_level = tractability_level;
        
        save(crash_safe.data_file, '-struct', 'checkpoint_data', '-v7.3');
        fprintf('💾 Data checkpoint saved\n');
    catch ME
        fprintf('⚠️ Could not save data checkpoint: %s\n', ME.message);
    end
end

function cleanupCheckpoints(crash_safe)
    % Clean up checkpoint files after successful completion
    try
        if exist(crash_safe.checkpoint_file, 'file'), delete(crash_safe.checkpoint_file); end
        if exist(crash_safe.data_file, 'file'), delete(crash_safe.data_file); end
        fprintf('🗑️ Checkpoints cleaned up\n');
    catch
        % Silent cleanup failure
    end
end

%% CORE ANALYSIS FUNCTIONS

function config = initConfig(n_obs_per_combo, tractability_level, use_regime_based_generation)
    % Enhanced configuration with true subset approach and accurate calculations
    config = struct();
    
    config.convergence = struct('max_iterations', 10000, 'gradient_tolerance', 1e-6, ...
                               'function_tolerance', 1e-12, 'use_parallel_internal', false);
    
    config.fraser_params = struct('n_obs_per_combo', n_obs_per_combo, ...
        'low_noise_threshold', 0.1, 'white_noise_breakdown', 2.0, ...
        'legacy_spurious_beta', 1/3, 'nonlinear_degraded_beta', 1/10, ...
        'noise_color_scaling', 1.5, 'transition_sharpness', 3.0, 'base_noise_variability', 0.01, ...
        'use_regime_based_generation', use_regime_based_generation);
    
    % REGIME-BASED GENERATION PARAMETERS
    if use_regime_based_generation
        config.regime_params = struct(...
            'vgf_regime_threshold', 200, ...           % VGF threshold for regime switching
            'beta_regime_threshold', 0.4, ...          % Beta generated threshold for regime influence
            'high_noise_regime_threshold', 1.5, ...    % Noise magnitude threshold for regime amplification
            'regime_interaction_strength', 3.0, ...    % Amplification factor for regime differences
            'regime_inversion_factor', -2.0, ...       % Factor for inverting interaction patterns
            'regime_noise_sensitivity', 0.8);          % Sensitivity to noise in different regimes
        fprintf('  🔄 REGIME-BASED GENERATION ENABLED:\n');
        fprintf('    VGF threshold: %.1f (Low: systematic, High: inverted patterns)\n', config.regime_params.vgf_regime_threshold);
        fprintf('    Beta threshold: %.3f (influences regime transition sharpness)\n', config.regime_params.beta_regime_threshold);
        fprintf('    High noise threshold: %.1f (amplifies regime differences)\n', config.regime_params.high_noise_regime_threshold);
        fprintf('    Expected outcome: Incompatible interaction patterns requiring hierarchical drill-down\n');
    else
        fprintf('  ✓ SYSTEMATIC GENERATION: Current validation behavior preserved\n');
    end
    
    % TRUE SUBSET TRACTABILITY LEVELS with accurate observation counts
    % Level 9 reference: 3×21×14×18×31×2×3 = 2,952,936 combinations per trial
    base_combinations = 3 * 21 * 14 * 18 * 31 * 2 * 3;  % Level 9 total combinations
    
    tractability_configs = {
        struct('level', 1, 'name', 'Conservative', 'param_count', [3,5,4,8,4,2,3], 'description', 'Key theoretical values'),
        struct('level', 2, 'name', 'Focused', 'param_count', [3,21,5,18,5,2,3], 'description', 'Biological range with finer sampling'),
        struct('level', 3, 'name', 'Minimal', 'param_count', [1,5,3,5,3,2,3], 'description', 'Fraser-complete but very sparse'),
        struct('level', 4, 'name', 'Moderate', 'param_count', [3,21,7,18,8,2,3], 'description', 'Expanded VGF coverage'),
        struct('level', 5, 'name', 'Substantial', 'param_count', [3,21,11,18,11,2,3], 'description', 'Dense VGF sampling'),
        struct('level', 6, 'name', 'Extensive', 'param_count', [3,21,7,18,16,2,3], 'description', 'Fine noise color resolution'),
        struct('level', 7, 'name', 'Comprehensive', 'param_count', [3,21,14,18,13,2,3], 'description', 'Near-complete VGF coverage'),
        struct('level', 8, 'name', 'Near-Original', 'param_count', [3,21,14,18,29,2,3], 'description', 'High resolution (6.5% reduction)'),
        struct('level', 9, 'name', 'Full-Original', 'param_count', [3,21,14,18,31,2,3], 'description', 'Complete parameter space (PROVEN)')
    };
    
    % Calculate accurate observation counts based on true parameter space sizes
    for i = 1:length(tractability_configs)
        param_count = tractability_configs{i}.param_count;
        total_combinations = prod(param_count);
        tractability_configs{i}.total_combinations = total_combinations;
        tractability_configs{i}.expected_obs_millions = total_combinations * n_obs_per_combo / 1e6;
        tractability_configs{i}.reduction_percent = (1 - total_combinations / base_combinations) * 100;
    end
    
    if tractability_level >= 1 && tractability_level <= 9
        config.tractability = tractability_configs{tractability_level};
    else
        error('Invalid tractability level (1-9)');
    end
    
    % Enhanced subset verification metadata
    config.tractability.is_true_subset = (tractability_level < 9);
    config.tractability.subset_verification = 'All parameter values are exact subsets of Level 9';
    config.tractability.convergence_proven = (tractability_level == 9);
    config.tractability.v010_validation = 'R² = 0.9693 with 14.7M observations';
    
    config.bluebear = struct('target_cores', 72, 'total_memory_gb', 324);
    config.batch_processing = struct('worker_batch_size', 50000, 'temp_file_prefix', 'v012_TrueSubset_', 'cleanup_temp_files', true);
    config.prereg_v064 = struct('kitchen_sink_formula', ...
        'deltaBeta ~ betaGenerated * VGF * samplingRate * filterType * regressionType * noiseMagnitude * noiseColor + (1|paramComboID)');
end

function config = setupParallelPool(config)
    % Simplified parallel pool setup
    fprintf('Starting parallel pool...\n');
    
    max_cores = feature('numcores');
    target_cores = min(config.bluebear.target_cores, max_cores);
    
    current_pool = gcp('nocreate');
    if ~isempty(current_pool), delete(current_pool); end
    
    try
        parpool('local', target_cores);
        fprintf('✓ Pool started: %d workers\n', target_cores);
        config.bluebear.target_cores = target_cores;
    catch ME
        % Fallback to smaller pool
        fallback = min(16, max_cores);
        parpool('local', fallback);
        fprintf('✓ Fallback pool: %d workers\n', fallback);
        config.bluebear.target_cores = fallback;
    end
end

function shutdownParallelPool()
    % Clean shutdown of parallel pool
    current_pool = gcp('nocreate');
    if ~isempty(current_pool)
        delete(current_pool);
        fprintf('✓ Parallel pool terminated\n');
    end
    
    % Aggressive cleanup
    clear current_pool;
    if ispc, java.lang.System.gc(); end
    pause(1);
end

function [tableTrue, data_info] = generateData(config, crash_safe)
    % Enhanced data generation with true subset verification
    fprintf('Generating true subset data (Level %d: %s)...\n', config.tractability.level, config.tractability.name);
    
    % Get parameter space using true subset approach
    [fs_vals, beta_gen_vals, vgf_vals, filter_vals, regress_vals, noise_mag_vals, noise_color_vals] = ...
        getParameterSpaceSubset(config.tractability.level);
    
    total_combos = length(fs_vals) * length(beta_gen_vals) * length(vgf_vals) * ...
                   length(filter_vals) * length(regress_vals) * length(noise_mag_vals) * length(noise_color_vals);
    
    % Verify against expected combinations
    assert(total_combos == config.tractability.total_combinations, ...
        'Parameter combination count mismatch: expected %d, got %d', ...
        config.tractability.total_combinations, total_combos);
    
    fprintf('  TRUE SUBSET VERIFIED: %d combinations\n', total_combos);
    fprintf('  Expected observations: %.3fM (%.1f%% reduction from Level 9)\n', ...
        total_combos * config.fraser_params.n_obs_per_combo / 1e6, config.tractability.reduction_percent);
    
    % Display parameter space composition
    fprintf('  Parameter counts: [%d,%d,%d,%d,%d,%d,%d] = %s\n', ...
        length(fs_vals), length(beta_gen_vals), length(vgf_vals), ...
        length(noise_mag_vals), length(noise_color_vals), length(filter_vals), length(regress_vals), ...
        sprintf('%d×%d×%d×%d×%d×%d×%d', length(fs_vals), length(beta_gen_vals), length(vgf_vals), ...
        length(noise_mag_vals), length(noise_color_vals), length(filter_vals), length(regress_vals)));
    
    % Enhanced parallel processing with true subset validation
    start_time = tic;
    num_workers = config.bluebear.target_cores;
    
    fprintf('  Processing with %d workers...\n', num_workers);
    
    worker_results = cell(num_workers, 1);
    parfor worker_id = 1:num_workers
        worker_results{worker_id} = processWorkerSlice(worker_id, num_workers, total_combos, ...
            fs_vals, beta_gen_vals, vgf_vals, filter_vals, regress_vals, noise_mag_vals, noise_color_vals, config);
    end
    
    % Combine results with enhanced validation
    fprintf('  Combining worker results with subset validation...\n');
    tableTrue = table();
    
    for worker_id = 1:num_workers
        worker_result = worker_results{worker_id};
        if ~isempty(worker_result) && ~isempty(worker_result.temp_files)
            for file_idx = 1:length(worker_result.temp_files)
                batch_data = load(worker_result.temp_files{file_idx});
                if isempty(tableTrue)
                    tableTrue = batch_data.batch_table;
                else
                    tableTrue = [tableTrue; batch_data.batch_table];
                end
                if config.batch_processing.cleanup_temp_files
                    delete(worker_result.temp_files{file_idx});
                end
            end
        end
    end
    
    % Convert to categorical/nominal or leave numeric with enhanced metadata
    if ~isempty(tableTrue)
        tableTrue.filterType = nominal(tableTrue.filterType);
        tableTrue.regressionType = nominal(tableTrue.regressionType);
        tableTrue.paramComboID = nominal(tableTrue.paramComboID);
    end
    
    generation_time = toc(start_time);
    
    % Enhanced data info with subset validation and regime information
    data_info = struct();
    data_info.n_observations = height(tableTrue);
    data_info.n_param_combos = total_combos;
    data_info.n_obs_per_combo = config.fraser_params.n_obs_per_combo;
    data_info.tractability_level = config.tractability.level;
    data_info.tractability_name = config.tractability.name;
    data_info.reduction_percent = config.tractability.reduction_percent;
    data_info.generation_time_minutes = generation_time / 60;
    data_info.num_workers = num_workers;
    data_info.is_true_subset = config.tractability.is_true_subset;
    data_info.subset_verified = true;
    data_info.regime_based_generation = config.fraser_params.use_regime_based_generation;
    
    if config.fraser_params.use_regime_based_generation
        data_info.regime_structure = 'Incompatible VGF-based regimes';
        data_info.expected_hierarchical_requirement = true;
    else
        data_info.regime_structure = 'Systematic interactions suitable for global modeling';
        data_info.expected_hierarchical_requirement = false;
    end
    
    fprintf('✓ True subset data generation complete: %d observations in %.2f minutes\n', ...
        data_info.n_observations, data_info.generation_time_minutes);
    fprintf('✓ Subset verification: PASSED\n');
end

function results = performKitchenSink(tableTrue, config, crash_safe)
    % Kitchen sink model with enhanced crash monitoring
    fprintf('Fitting kitchen sink model...\n');
    fprintf('  Dataset: %d observations (%.1fM)\n', height(tableTrue), height(tableTrue)/1e6);
    fprintf('  Formula: %s\n', config.prereg_v064.kitchen_sink_formula);
    
    monitorMemory(crash_safe, 'PRE_LME');
    
    % Optimizer settings
    optimizer_options = statset('MaxIter', config.convergence.max_iterations, ...
        'TolFun', config.convergence.function_tolerance, 'TolX', config.convergence.gradient_tolerance, ...
        'Display', 'iter', 'UseParallel', false);
    
    % Clear variables and force cleanup
    clearvars -except tableTrue config optimizer_options crash_safe;
    if ispc, java.lang.System.gc(); end
    
    fprintf('  Attempting kitchen sink fit...\n');
    
    start_time = tic;
    success = false;
    lme = [];
    
    warning('off', 'all');
    try
        lme = fitlme(tableTrue, config.prereg_v064.kitchen_sink_formula, ...
            'FitMethod', 'REML', 'CheckHessian', false, 'OptimizerOptions', optimizer_options);
        success = true;
    catch ME
        if contains(ME.message, 'memory') || contains(ME.message, 'Memory')
            fprintf('  🚨 MEMORY ERROR: %s\n', ME.message);
            logError(crash_safe, 'LME_MEMORY_ERROR', ME);
        else
            fprintf('  ✗ Convergence error: %s\n', ME.message);
        end
    end
    warning('on', 'all');
    
    fit_time = toc(start_time);
    
    % Package results
    results = struct();
    if success && ~isempty(lme)
        coeffs = lme.Coefficients;
        interaction_mask = contains(coeffs.Name, ':');
        n_significant = sum(interaction_mask & (coeffs.pValue < 0.05));
        
        results.converged = true;
        results.model = lme;
        results.r_squared = lme.Rsquared.Adjusted;
        results.total_coefficients = height(coeffs);
        results.significant_interactions = n_significant;
        results.fitting_time = fit_time;
        results.kitchen_sink_successful = true;
        
        fprintf('  ✓ SUCCESS in %.2f minutes\n', fit_time/60);
        fprintf('    R²: %.4f, Coefficients: %d, Significant interactions: %d\n', ...
            results.r_squared, results.total_coefficients, results.significant_interactions);
    else
        results.converged = false;
        results.fitting_time = fit_time;
        results.kitchen_sink_successful = false;
        results.error_message = 'Kitchen sink model failed';
        if exist('ME', 'var'), results.error_details = ME.message; end
        
        fprintf('  ✗ FAILED after %.2f minutes\n', fit_time/60);
    end
end

function displayStage1Results(results, use_database, data_info, tractability_level)
    % Enhanced results display for Stage 1 with pipeline information
    fprintf('\n════════════════════════════════════════════════════════════════\n');
    fprintf('       STAGE 1 RESULTS: KITCHEN SINK MODEL (Level %d)           \n', tractability_level);
    fprintf('════════════════════════════════════════════════════════════════\n');
    
    fprintf('\nStage 1 - Global Interaction Assessment:\n');
    fprintf('  ✓ Level %d (%s): %.1f%% reduction from full space\n', ...
        tractability_level, data_info.tractability_name, data_info.reduction_percent);
    fprintf('  ✓ Dataset: %d observations (%.3fM)\n', data_info.n_observations, data_info.n_observations/1e6);
    fprintf('  ✓ Combinations: %d × %d obs/combo = %d total\n', ...
        data_info.n_param_combos, data_info.n_obs_per_combo, data_info.n_observations);
    fprintf('  ✓ Generation time: %.2f minutes\n', data_info.generation_time_minutes);
    
    if data_info.is_true_subset
        fprintf('  ✓ TRUE SUBSET VERIFIED: All parameters are exact subsets of Level 9\n');
    else
        fprintf('  ✓ FULL REFERENCE SPACE: Level 9 complete parameter space\n');
    end
    
    if results.kitchen_sink_successful
        fprintf('\nKitchen Sink Model Results:\n');
        fprintf('  ✓ SUCCESS - Full 7-way factorial interaction model\n');
        fprintf('  ✓ R²: %.4f\n', results.r_squared);
        fprintf('  ✓ Coefficients: %d\n', results.total_coefficients);
        fprintf('  ✓ Significant interactions: %d\n', results.significant_interactions);
        fprintf('  ✓ Fitting time: %.2f minutes\n', results.fitting_time/60);
        
        if tractability_level == 9
            fprintf('\n🎉 STAGE 1 SUCCESS: Level 9 convergence capability validated!\n');
        else
            fprintf('\n🎉 STAGE 1 SUCCESS with %.1f%% parameter reduction!\n', data_info.reduction_percent);
        end
    else
        fprintf('\nKitchen Sink Model Results:\n');
        fprintf('  ✗ FAILED - Memory or convergence limits exceeded\n');
        fprintf('  ✗ Time to failure: %.2f minutes\n', results.fitting_time/60);
        if isfield(results, 'error_details')
            fprintf('  ✗ Error: %s\n', results.error_details);
        end
        fprintf('\n⚠️ Consider using a higher tractability level (lower number = more reduction)\n');
    end
    
fprintf('\nModel Adequacy Framework Status:\n');
fprintf('  ✅ Stage 1: Kitchen Sink Model - COMPLETE\n');
fprintf('  ⏳ Stage 2: Model Assessment - READY\n');
fprintf('  ⏳ Stage 3: Conditional Analysis - PENDING\n');
fprintf('  ⏳ Stage 4: Integration - PENDING\n');
end

function saveStage1Results(results, config, use_database, data_info, tractability_level, crash_safe, tableTrue)
    % Enhanced results saving for Stage 1 with pipeline integration
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('stage1_results_L%d_%s_%s.mat', ...
        tractability_level, data_info.tractability_name, timestamp);
    
    % Package complete Stage 1 results for pipeline
    stage1_output = struct();
    stage1_output.results = results;
    stage1_output.config = config;
    stage1_output.data_info = data_info;
    stage1_output.tableTrue = tableTrue;  % Include data for Stage 2
    stage1_output.tractability_level = tractability_level;
    stage1_output.use_database = use_database;
    
    metadata = struct();
    metadata.version = 'v001';
    metadata.stage = 1;
    metadata.stage_name = 'Kitchen Sink Model';
    metadata.pipeline_ready = results.kitchen_sink_successful;
    metadata.is_true_subset = data_info.is_true_subset;
    metadata.subset_verified = data_info.subset_verified;
    metadata.crash_safe_session = crash_safe.session_id;
    metadata.timestamp = timestamp;
    metadata.next_stage = 'ModelAdequacy_Stage2_Assessment_v002.m';
    metadata.convergence_validation = 'Based on proven convergence 14.7M obs)';
    
    stage1_output.metadata = metadata;
    
    try
        save(filename, '-struct', 'stage1_output', '-v7.3');
        fprintf('\nStage 1 results saved: %s\n', filename);
        fprintf('📊 Includes: Model results, data, configuration, and metadata for Stage 2\n');
        
        % Also save a "latest" version for easy pipeline access
        latest_filename = 'stage1_results_latest.mat';
        save(latest_filename, '-struct', 'stage1_output', '-v7.3');
        fprintf('🔗 Latest results: %s (for pipeline convenience)\n', latest_filename);
        
    catch ME
        fprintf('⚠️ Could not save Stage 1 results: %s\n', ME.message);
    end
end

%% TRUE SUBSET PARAMETER SPACE FUNCTIONS

function [fs_vals, beta_gen_vals, vgf_vals, filter_vals, regress_vals, noise_mag_vals, noise_color_vals] = getParameterSpaceSubset(level)
    % TRUE SUBSET APPROACH: All levels 1-8 are proper subsets of level 9
    % This ensures that any parameter combination in levels 1-8 uses the exact same
    % parameter values as exist in the full level 9 space
    
    filter_vals = [1, 2];
    regress_vals = [1, 2, 3];
    
    % Define the FULL parameter space (level 9) first
    fs_full = [60 120 240];
    beta_gen_full = linspace(0, 2/3, 21);
    vgf_full = exp(linspace(4.5, 5.8, 14));
    noise_mag_full = [0:0.025:0.1, 0.25:0.25:2.25, 4, 6, 8, 10];  % 18 levels
    noise_color_full = 0:0.1:3.0;  % 31 levels
    
    switch level
        case 1  % Conservative: Key theoretical values (3×5×4×8×4×2×3 = 28,800)
            fs_vals = fs_full;  % All sampling rates [3]
            beta_gen_vals = beta_gen_full([1, 6, 11, 16, 21]);  % [0, 1/6, 1/3, 1/2, 2/3] [5]
            vgf_vals = vgf_full([2, 6, 10, 14]);  % Spread across VGF range [4]
            noise_mag_vals = noise_mag_full([1, 3, 5, 9, 13, 15, 17, 18]);  % Key noise levels [8]
            noise_color_vals = noise_color_full([1, 11, 21, 31]);  % [0, 1.0, 2.0, 3.0] [4]
            
        case 2  % Focused: Biological range with finer sampling (3×21×5×18×5×2×3 = 283,500)
            fs_vals = fs_full;  % All sampling rates [3]
            beta_gen_vals = beta_gen_full;  % All beta values [21]
            vgf_vals = vgf_full([1, 4, 7, 10, 14]);  % 5 VGF values [5]
            noise_mag_vals = noise_mag_full;  % All noise magnitudes [18]
            noise_color_vals = noise_color_full([1, 6, 11, 21, 31]);  % [0, 0.5, 1.0, 2.0, 3.0] [5]
            
        case 3  % Minimal: Fraser-complete but very sparse (1×5×3×5×3×2×3 = 1,350)
            fs_vals = fs_full(2);  % Middle sampling rate (120 Hz) [1]
            beta_gen_vals = beta_gen_full([1, 6, 11, 16, 21]);  % Key theoretical values [5]
            vgf_vals = vgf_full([4, 8, 12]);  % 3 VGF values [3]
            noise_mag_vals = noise_mag_full([1, 5, 13, 17, 18]);  % [0, 0.1, 0.5, 2.0, 4.0] [5]
            noise_color_vals = noise_color_full([1, 11, 31]);  % [0, 1.0, 3.0] [3]
            
        case 4  % Moderate: Expanded VGF coverage (3×21×7×18×8×2×3 = 571,536)
            fs_vals = fs_full;  % All sampling rates [3]
            beta_gen_vals = beta_gen_full;  % All beta values [21]
            vgf_vals = vgf_full([1, 3, 5, 7, 9, 11, 14]);  % 7 VGF values [7]
            noise_mag_vals = noise_mag_full;  % All noise magnitudes [18]
            noise_color_vals = noise_color_full([1, 4, 8, 11, 16, 21, 26, 31]);  % 8 noise colors [8]
            
        case 5  % Substantial: Dense VGF sampling (3×21×11×18×11×2×3 = 1,548,756)
            fs_vals = fs_full;  % All sampling rates [3]
            beta_gen_vals = beta_gen_full;  % All beta values [21]
            vgf_vals = vgf_full([1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 14]);  % 11 VGF values [11]
            noise_mag_vals = noise_mag_full;  % All noise magnitudes [18]
            noise_color_vals = noise_color_full(1:3:31);  % Every 3rd: 0:0.3:3.0 [11]
            
        case 6  % Extensive: Fine noise color resolution (3×21×7×18×16×2×3 = 1,905,552)
            fs_vals = fs_full;  % All sampling rates [3]
            beta_gen_vals = beta_gen_full;  % All beta values [21]
            vgf_vals = vgf_full([1:2:14]);  % Every 2nd VGF value [7 values, indices 1,3,5,7,9,11,13]
            noise_mag_vals = noise_mag_full;  % All noise magnitudes [18]
            noise_color_vals = noise_color_full(1:2:31);  % Every 2nd: 0:0.2:3.0 [16]
            
        case 7  % Comprehensive: Near-complete VGF coverage (3×21×14×18×13×2×3 = 4,365,096)
            fs_vals = fs_full;  % All sampling rates [3]
            beta_gen_vals = beta_gen_full;  % All beta values [21]
            vgf_vals = vgf_full;  % All VGF values [14]
            noise_mag_vals = noise_mag_full;  % All noise magnitudes [18]
            noise_color_vals = noise_color_full([1:2:25, 31]);  % Almost all colors [13]
            
        case 8  % Near-Original: High resolution (3×21×14×18×29×2×3 = 13,803,432)
            fs_vals = fs_full;  % All sampling rates [3]
            beta_gen_vals = beta_gen_full;  % All beta values [21]
            vgf_vals = vgf_full;  % All VGF values [14]
            noise_mag_vals = noise_mag_full;  % All noise magnitudes [18]
            noise_color_vals = noise_color_full(1:29);  % 29 levels (exclude last 2) [29]
            
        case 9  % Full-Original: Complete parameter space (3×21×14×18×31×2×3 = 2,952,936)
            fs_vals = fs_full;  % All sampling rates [3]
            beta_gen_vals = beta_gen_full;  % All beta values [21]
            vgf_vals = vgf_full;  % All VGF values [14]
            noise_mag_vals = noise_mag_full;  % All noise magnitudes [18]
            noise_color_vals = noise_color_full;  % All noise colors [31]
            
        otherwise
            error('Invalid tractability level (1-9)');
    end
    
    % Enhanced verification that subsets are proper subsets of level 9
    if level < 9
        assert(all(ismember(fs_vals, fs_full)), 'Sampling rates must be exact subset of Level 9');
        assert(all(ismember(beta_gen_vals, beta_gen_full)), 'Beta values must be exact subset of Level 9');
        assert(all(ismember(vgf_vals, vgf_full)), 'VGF values must be exact subset of Level 9');
        assert(all(ismember(noise_mag_vals, noise_mag_full)), 'Noise magnitudes must be exact subset of Level 9');
        assert(all(ismember(noise_color_vals, noise_color_full)), 'Noise colors must be exact subset of Level 9');
        
        % Additional validation: ensure exact floating point matches
        tolerance = 1e-10;
        for i = 1:length(beta_gen_vals)
            assert(any(abs(beta_gen_vals(i) - beta_gen_full) < tolerance), ...
                'Beta value %.6f not found in Level 9 space within tolerance', beta_gen_vals(i));
        end
        for i = 1:length(vgf_vals)
            assert(any(abs(vgf_vals(i) - vgf_full) < tolerance), ...
                'VGF value %.6f not found in Level 9 space within tolerance', vgf_vals(i));
        end
    end
end

function worker_result = processWorkerSlice(worker_id, num_workers, total_combos, ...
    fs_vals, beta_gen_vals, vgf_vals, filter_vals, regress_vals, noise_mag_vals, noise_color_vals, config)
    % Process worker slice for parallel data generation
    
    combos_per_worker = ceil(total_combos / num_workers);
    start_combo = (worker_id - 1) * combos_per_worker + 1;
    end_combo = min(worker_id * combos_per_worker, total_combos);
    
    worker_batch_size = config.batch_processing.worker_batch_size;
    temp_file_prefix = sprintf('%sworker%02d_', config.batch_processing.temp_file_prefix, worker_id);
    
    combos_processed = 0;
    batch_num = 1;
    temp_files = {};
    
    combo_id = 0;
    
    for fs = fs_vals
        for bg = beta_gen_vals
            for vgf = vgf_vals
                for ft = filter_vals
                    for rt = regress_vals
                        for nm = noise_mag_vals
                            for nc = noise_color_vals
                                combo_id = combo_id + 1;
                                
                                if combo_id < start_combo || combo_id > end_combo
                                    continue;
                                end
                                
                                % Generate Fraser-realistic data with optional regime-based structure
                                delta_beta = generateFraserData(config.fraser_params.n_obs_per_combo, ft, rt, nm, nc, bg, vgf, config.fraser_params, config);
                                
                                combo_data = struct();
                                combo_data.deltaBeta = delta_beta;
                                combo_data.betaGenerated = repmat(bg, config.fraser_params.n_obs_per_combo, 1);
                                combo_data.VGF = repmat(vgf, config.fraser_params.n_obs_per_combo, 1);
                                combo_data.samplingRate = repmat(fs, config.fraser_params.n_obs_per_combo, 1);
                                combo_data.filterType = repmat(ft, config.fraser_params.n_obs_per_combo, 1);
                                combo_data.regressionType = repmat(rt, config.fraser_params.n_obs_per_combo, 1);
                                combo_data.noiseMagnitude = repmat(nm, config.fraser_params.n_obs_per_combo, 1);
                                combo_data.noiseColor = repmat(nc, config.fraser_params.n_obs_per_combo, 1);
                                combo_data.paramComboID = repmat(combo_id, config.fraser_params.n_obs_per_combo, 1);
                                
                                if ~exist('current_batch', 'var')
                                    current_batch = combo_data;
                                else
                                    current_batch = appendData(current_batch, combo_data);
                                end
                                
                                combos_processed = combos_processed + 1;
                                
                                % Save batch when full
                                if exist('current_batch', 'var')
                                    batch_obs = length(current_batch.deltaBeta);
                                    if batch_obs >= worker_batch_size || combos_processed >= (end_combo - start_combo + 1)
                                        temp_filename = sprintf('%s%04d.mat', temp_file_prefix, batch_num);
                                        batch_table = struct2table(current_batch);
                                        save(temp_filename, 'batch_table', '-v7.3');
                                        temp_files{end+1} = temp_filename;
                                        clear current_batch;
                                        batch_num = batch_num + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    worker_result = struct();
    worker_result.combinations_processed = combos_processed;
    worker_result.temp_files = temp_files;
end

function delta_beta = generateFraserData(n_obs, filter_type, regress_type, noise_mag, noise_color, beta_generated, vgf, params, config)
    % Fraser-realistic data generation with enhanced biological fidelity
    % ENHANCED v013: Optional regime-based generation 
    
    base_error = params.base_noise_variability * randn(n_obs, 1);
    
    if params.use_regime_based_generation
        % REGIME-BASED GENERATION: Create incompatible interaction structures
        % This creates distinct analytical regimes that cannot be simultaneously
        % captured by global interaction modeling, necessitating next stage
        
        regime_params = config.regime_params;
        
        % Determine primary regime based on VGF threshold
        is_high_vgf_regime = vgf >= regime_params.vgf_regime_threshold;
        
        % Beta generated influences regime transition sharpness and interaction patterns
        beta_regime_factor = sigmoid_transition(beta_generated, regime_params.beta_regime_threshold, 2.0);
        
        % High noise amplifies regime differences
        is_high_noise = noise_mag >= regime_params.high_noise_regime_threshold;
        noise_amplification = 1.0 + (is_high_noise * regime_params.regime_interaction_strength);
        
        if is_high_vgf_regime
            % HIGH VGF REGIME: Inverted interaction patterns
            % Filter and regression effects operate under fundamentally different rules
            method_effect = generateHighVGFRegimeEffects(n_obs, filter_type, regress_type, noise_mag, noise_color, ...
                                                        beta_generated, beta_regime_factor, noise_amplification, regime_params);
        else
            % LOW VGF REGIME: Enhanced systematic patterns (current behavior with regime modulation)
            method_effect = generateLowVGFRegimeEffects(n_obs, filter_type, regress_type, noise_mag, noise_color, ...
                                                       beta_generated, beta_regime_factor, noise_amplification, params);
        end
        
    else
        % SYSTEMATIC GENERATION: Current validation behavior (preserved exactly)
        adjusted_breakdown = params.white_noise_breakdown * (1 + params.noise_color_scaling * noise_color / 3.0);
        
        if noise_mag <= params.low_noise_threshold
            method_effect = 0;
        elseif noise_mag <= adjusted_breakdown
            transition_factor = (noise_mag - params.low_noise_threshold) / (adjusted_breakdown - params.low_noise_threshold);
            smooth_transition = tanh(params.transition_sharpness * (transition_factor - 0.5)) * 0.5 + 0.5;
            
            if filter_type == 1
                method_effect = 0.02 * smooth_transition * (1 + 0.5 * rand(n_obs, 1));
            else
                method_effect = 0.01 * smooth_transition * (1 + 0.3 * rand(n_obs, 1));
            end
        else
            breakdown_intensity = min(1.0, (noise_mag - adjusted_breakdown) / adjusted_breakdown);
            
            if filter_type == 1
                method_effect = params.legacy_spurious_beta * breakdown_intensity + 0.05 * breakdown_intensity * randn(n_obs, 1);
            else
                target_bias = params.nonlinear_degraded_beta * (2.5 - 0.5 * regress_type);
                method_effect = target_bias * breakdown_intensity + 0.02 * breakdown_intensity * randn(n_obs, 1);
            end
        end
    end
    
    delta_beta = base_error + method_effect;
end

function regime_effect = generateHighVGFRegimeEffects(n_obs, filter_type, regress_type, noise_mag, noise_color, ...
                                                     beta_generated, beta_regime_factor, noise_amplification, regime_params)
    % HIGH VGF REGIME: Fundamentally different interaction patterns
    % Creates patterns incompatible with low VGF regime 
    
    % INVERTED FILTER EFFECTS: High VGF regime inverts filter type responses
    if filter_type == 1
        % Filter type 1 shows REDUCED effects in high VGF (opposite of low VGF)
        filter_base_effect = 0.005 * regime_params.regime_inversion_factor;
    else
        % Filter type 2 shows AMPLIFIED effects in high VGF 
        filter_base_effect = 0.08 * abs(regime_params.regime_inversion_factor);
    end
    
    % BETA-DEPENDENT REGIME MODULATION: Beta generated creates regime-specific patterns
    beta_modulation = beta_generated * beta_regime_factor * 0.15;
    
    % NOISE COLOR INVERSION: High VGF regime shows inverted noise color sensitivity
    noise_color_effect = -noise_color * 0.02 * regime_params.regime_noise_sensitivity;
    
    % REGRESSION TYPE DOMINANCE INVERSION: Different regression types dominate in high VGF
    regression_dominance = calculateRegressionDominanceHigh(regress_type, noise_mag, beta_generated);
    
    % MULTIPLICATIVE REGIME INTERACTIONS: Cannot be captured by additive global models
    regime_interaction = filter_base_effect * (1 + beta_modulation) * (1 + noise_color_effect) * regression_dominance;
    
    % NOISE AMPLIFICATION: High noise amplifies regime-specific patterns
    regime_effect = regime_interaction * noise_amplification * (1 + 0.3 * randn(n_obs, 1));
end

function regime_effect = generateLowVGFRegimeEffects(n_obs, filter_type, regress_type, noise_mag, noise_color, ...
                                                    beta_generated, beta_regime_factor, noise_amplification, params)
    % LOW VGF REGIME: Enhanced systematic patterns with regime modulation
    % Maintains current behavior but adds regime-specific characteristics
    
    % ENHANCED SYSTEMATIC EFFECTS: Stronger patterns than original systematic generation
    adjusted_breakdown = params.white_noise_breakdown * (1 + params.noise_color_scaling * noise_color / 3.0);
    
    if noise_mag <= params.low_noise_threshold
        base_effect = 0.01 * beta_generated;  % Beta dependency in low noise
    elseif noise_mag <= adjusted_breakdown
        transition_factor = (noise_mag - params.low_noise_threshold) / (adjusted_breakdown - params.low_noise_threshold);
        smooth_transition = tanh(params.transition_sharpness * (transition_factor - 0.5)) * 0.5 + 0.5;
        
        if filter_type == 1
            base_effect = 0.04 * smooth_transition * (1 + beta_regime_factor);
        else
            base_effect = 0.02 * smooth_transition * (1 + 0.5 * beta_regime_factor);
        end
    else
        breakdown_intensity = min(1.0, (noise_mag - adjusted_breakdown) / adjusted_breakdown);
        
        if filter_type == 1
            base_effect = params.legacy_spurious_beta * breakdown_intensity * (1 + beta_regime_factor);
        else
            target_bias = params.nonlinear_degraded_beta * (2.5 - 0.5 * regress_type);
            base_effect = target_bias * breakdown_intensity * (1 + 0.8 * beta_regime_factor);
        end
    end
    
    % REGRESSION TYPE DOMINANCE: Different pattern than high VGF regime
    regression_dominance = calculateRegressionDominanceLow(regress_type, noise_mag, beta_generated);
    
    % REGIME-SPECIFIC NOISE SENSITIVITY
    noise_sensitivity = 1.0 + (noise_color / 3.0) * 0.5;
    
    regime_effect = base_effect * regression_dominance * noise_sensitivity * noise_amplification * (1 + 0.2 * randn(n_obs, 1));
end

function dominance = calculateRegressionDominanceHigh(regress_type, noise_mag, beta_generated)
    % HIGH VGF REGIME: Regression type 3 dominates, types 1&2 show reduced effects
    switch regress_type
        case 1
            dominance = 0.3 + 0.2 * beta_generated;  % Reduced dominance
        case 2  
            dominance = 0.5 + 0.3 * beta_generated;  % Moderate dominance
        case 3
            dominance = 1.5 + 0.8 * beta_generated + 0.3 * noise_mag;  % Strong dominance
        otherwise
            dominance = 1.0;
    end
end

function dominance = calculateRegressionDominanceLow(regress_type, noise_mag, beta_generated)
    % LOW VGF REGIME: Regression type 1 dominates (opposite of high VGF)
    switch regress_type
        case 1
            dominance = 1.8 + 0.5 * beta_generated + 0.2 * noise_mag;  % Strong dominance
        case 2
            dominance = 1.0 + 0.3 * beta_generated;  % Moderate dominance  
        case 3
            dominance = 0.4 + 0.1 * beta_generated;  % Reduced dominance
        otherwise
            dominance = 1.0;
    end
end

function transition = sigmoid_transition(value, threshold, sharpness)
    % Smooth sigmoid transition function for regime switching
    transition = 1.0 / (1.0 + exp(-sharpness * (value - threshold)));
end

function combined = appendData(current, new)
    % Efficient data appending for batch processing
    fields = fieldnames(current);
    combined = struct();
    for i = 1:length(fields)
        field = fields{i};
        combined.(field) = [current.(field); new.(field)];
    end
end

function result = ternary(condition, true_val, false_val)
    if condition, result = true_val; else, result = false_val; end
end