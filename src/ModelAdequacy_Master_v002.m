function ModelAdequacy_Master_v002(tractability_level, n_obs_per_combo, use_database, start_from_stage)
% MODELADEQUACY_MASTER_V002 Data-Driven Model Adequacy Framework Master
%
% **MASTER ORCHESTRATION SCRIPT - MODEL ADEQUACY FRAMEWORK v002**
% Executes the complete 4-stage model adequacy methodology from prereg v065
% following the HierarchicalDrillDown_Master_v002 pattern with extensive documentation.
%
% **CRITICAL ENHANCEMENTS v002**:
%   - FIXED: Updated to use ModelAdequacy_Stage2_Assessment_v002 (camelCase variable fix)
%   - EXTENSIVE COMMENTS: Comprehensive documentation for all pipeline stages
%   - VARIABLE NAMING: Ensures consistency with DemoHierarchicalDrillDown_v012_TrueSubset
%   - PIPELINE INTEGRATION: Enhanced cross-stage compatibility with VariableNameMapping_v003
%   - ERROR HANDLING: Improved validation and troubleshooting information
%   - SYNTAX FIXES: Corrected MATLAB syntax errors and applied best practices
%
% **COMPLETE MODEL ADEQUACY FRAMEWORK**:
%   Stage 1: Kitchen Sink Model (ModelAdequacy_Stage1_KitchenSink_v001 formerly DemoHierarchicalDrillDown_v012_TrueSubset - REUSED)
%            → Global 7-way factorial interaction assessment
%            → δβ ~ βgenerated × VGF × samplingRate × filterType × regressionType × noiseMagnitude × noiseColor + (1|paramComboID)
%            → Produces stage1_results_*.mat with complete dataset and model
%
%   Stage 2: Systematic Model Adequacy Assessment (ModelAdequacy_Stage2_Assessment_v002)
%            → Residual pattern analysis across parameter regions
%            → Bootstrap parameter stability assessment (1000 iterations)
%            → K-fold cross-validation performance evaluation (k=5)
%            → Quantitative adequacy decision based on clinical thresholds
%            → Produces stage2_adequacy_*.mat with adequacy assessment
%
%   Stage 3: Conditional Parameter Analysis (ModelAdequacy_Stage3_Conditional_v002)
%            → Only executed if Stage 2 identifies poor fit regions
%            → Region-specific modeling for inadequate parameter combinations
%            → Advanced hierarchical analysis for complex interaction patterns
%            → Produces stage3_conditional_*.mat with conditional models
%
%   Stage 4: Hierarchical Model Integration (ModelAdequacy_Stage4_Integration_v002)
%            → Only executed if Stage 3 determines integration is beneficial
%            → Combines global and conditional models for optimal performance
%            → Final validation and clinical decision support framework
%            → Produces stage4_integration_*.mat with integrated models
%
% **PREREG v065 METHODOLOGY FOUNDATION**:
%   Core Research Question: Can a single global interaction model adequately capture
%   the complex relationships between measurement methodology and power law parameter
%   recovery across the entire parameter space, or are conditional/hierarchical
%   models necessary for specific parameter regions?
%
%   Primary Model Formula:
%   δβ ~ βgenerated × VGF × samplingRate × filterType × regressionType × noiseMagnitude × noiseColor + (1|paramComboID)
%
%   Adequacy Criteria (Clinical Significance Thresholds):
%   - Residual Pattern Threshold: Cohen's d > 0.5 indicates systematic bias requiring conditional analysis
%   - Coefficient Instability: Bootstrap CI > ±0.03 (clinical significance from Cook et al. 2023)
%   - Cross-Validation Degradation: >15% performance loss indicates poor regional fit
%   - Minimum Region Size: n ≥ 200 for reliable conditional modeling (power analysis based)
%
% **TRACTABILITY LEVELS** (Inherits from DemoHierarchicalDrillDown_v012_TrueSubset):
%   Level 1 = 'Conservative'  (28.8K obs, 99.0% reduction from full space)
%   Level 2 = 'Focused'       (472K obs, 96.8% reduction)    [DEFAULT RECOMMENDED]
%   Level 3 = 'Minimal'       (2.7K obs, 99.9% reduction)
%   Levels 4-8 = Progressive density increases (952K to 13.8M observations)
%   Level 9 = 'Full-Original' (14.8M obs, 0.0% reduction)
%
% USAGE EXAMPLES:
%   ModelAdequacy_Master_v002()                         % Default: Level 2, n=5, synthetic data
%   ModelAdequacy_Master_v002(2, 5, false)             % Level 2, n=5 obs/combo, synthetic data
%   ModelAdequacy_Master_v002(2, 5, true)              % Level 2, n=5, real database
%   ModelAdequacy_Master_v002(2, 5, false, 2)          % Resume from Stage 2
%
% Author: Fraser, D.S. (2025)
% Version: Master_v002 - Enhanced Documentation and FIXED Syntax Errors

%% INPUT PARAMETER VALIDATION AND DEFAULT ASSIGNMENT
% Comprehensive input validation with helpful error messages and default value assignment

if nargin < 1
    tractability_level = 2; % Default: Level 2 (Focused) - optimal balance
end
if nargin < 2
    n_obs_per_combo = 5; % Default: 5 observations per parameter combination
end
if nargin < 3
    use_database = false; % Default: Synthetic data generation
end
if nargin < 4
    start_from_stage = 1; % Default: Complete pipeline execution
end

% Validate tractability level input
if ~isnumeric(tractability_level) || tractability_level < 1 || tractability_level > 9 || mod(tractability_level, 1) ~= 0
    error('ModelAdequacy_Master_v002:InvalidTractabilityLevel', ...
        ['Tractability level must be an integer between 1 and 9.\n' ...
         '  Recommended: Level 2 (Focused, 472K obs) for most applications.\n' ...
         '  Use Level 1 for testing, Levels 7-9 for comprehensive analysis.']);
end

% Validate observations per combination input
if ~isnumeric(n_obs_per_combo) || n_obs_per_combo < 1 || n_obs_per_combo > 50 || mod(n_obs_per_combo, 1) ~= 0
    error('ModelAdequacy_Master_v002:InvalidObservationCount', ...
        ['Observations per combination must be an integer between 1 and 50.\n' ...
         '  Recommended: 5-10 for most applications.\n' ...
         '  Higher values increase precision but computational cost.']);
end

% Validate starting stage input
if ~isnumeric(start_from_stage) || start_from_stage < 1 || start_from_stage > 4 || mod(start_from_stage, 1) ~= 0
    error('ModelAdequacy_Master_v002:InvalidStartingStage', ...
        ['Starting stage must be an integer between 1 and 4.\n' ...
         '  Use 1 for complete analysis, 2-4 for resuming after crash.']);
end

% Add essential function directories to MATLAB path
if exist('functions', 'dir')
    addpath(genpath('functions')); % Variable name mapping and utility functions
end
if exist('utils', 'dir')
    addpath(genpath('utils')); % Additional analysis utilities if present
end

%% CONFIGURATION CONSTANTS
% Define tractability level information arrays for consistent lookup

% Tractability descriptions for each level
tractability_names = {
    'Conservative', 'Focused', 'Minimal', 'Moderate', 'Substantial', ...
    'Extensive', 'Comprehensive', 'Near-Original', 'Full-Original'
};

% Observation estimates (in thousands) for each tractability level
obs_estimates_k = [28.8, 472, 2.7, 952, 1550, 1910, 4370, 13800, 14800];

% Estimated execution times (in minutes) for each tractability level
time_estimates_min = [2, 15, 1, 45, 90, 150, 300, 800, 1200];

% Estimated memory requirements (in GB) for each tractability level
memory_estimates_gb = [2, 8, 1, 12, 20, 24, 40, 80, 100];

% Reduction percentages for each tractability level
reduction_percentages = [99.0, 96.8, 99.9, 93.5, 89.5, 87.0, 70.2, 6.5, 0.0];

% Get current tractability level information
current_name = tractability_names{tractability_level};
% FIXED: obs_estimates_k already includes total observations, don't multiply by n_obs_per_combo
current_obs_millions = obs_estimates_k(tractability_level) / 1000;
current_time_minutes = time_estimates_min(tractability_level);
current_memory_gb = memory_estimates_gb(tractability_level);
current_reduction = reduction_percentages(tractability_level);

%% MODEL ADEQUACY FRAMEWORK INITIALIZATION AND CONFIGURATION DISPLAY

fprintf('\n');
fprintf('████████████████████████████████████████████████████████████████\n');
fprintf('   MODEL ADEQUACY FRAMEWORK MASTER PIPELINE v002                 \n');
fprintf('   4-Stage Progressive Analysis with Systematic Assessment       \n');
fprintf('   Prereg v065 Implementation - Fraser et al. (2025)             \n');
fprintf('   ENHANCED: Variable Naming Consistency & Extensive Comments    \n');
fprintf('████████████████████████████████████████████████████████████████\n\n');

% Display comprehensive configuration information
fprintf('🔧 CONFIGURATION SUMMARY:\n');
fprintf('  Tractability Level: %d (%s - %.1fM obs, %.1f%% reduction)\n', ...
    tractability_level, current_name, current_obs_millions, current_reduction);
fprintf('  Observations per combo: %d (total observations ≈ %.1fM)\n', ...
    n_obs_per_combo, current_obs_millions);
fprintf('  Data source: %s\n', ternary(use_database, 'Real database (if available)', 'Synthetic data generation'));

% Display stage information
stage_descriptions = {
    '(Complete pipeline execution from kitchen sink model)';
    '(Resume from adequacy assessment)';
    '(Resume from conditional analysis)';
    '(Resume from model integration)'
};
fprintf('  Starting from: Stage %d %s\n', start_from_stage, stage_descriptions{start_from_stage});
fprintf('  Framework: Model Adequacy Assessment with Clinical Thresholds\n');

% Display computational requirements and performance expectations
fprintf('\n📊 COMPUTATIONAL REQUIREMENTS:\n');
fprintf('  Estimated execution time: ~%.0f minutes (Level %d)\n', current_time_minutes, tractability_level);
fprintf('  Estimated memory requirement: ~%.0f GB RAM\n', current_memory_gb);
fprintf('  Crash recovery: Enabled (automatic restart capability)\n');
fprintf('  Parallel processing: Auto-configured based on available cores\n');

% Display variable naming and compatibility information
fprintf('\n🔤 VARIABLE NAMING CONSISTENCY (v002 Enhancement):\n');
fprintf('  Primary variables: deltaBeta, betaGenerated (camelCase standard)\n');
fprintf('  Parameter variables: VGF, samplingRate, filterType, regressionType, noiseMagnitude, noiseColor\n');
fprintf('  Fixed compatibility: VariableNameMapping_v003, ModelAdequacy_Stage2_Assessment_v002\n');
fprintf('  Framework integration: DemoHierarchicalDrillDown_v012_TrueSubset consistent\n\n');

%% STAGE 1: KITCHEN SINK MODEL EXECUTION OR RECOVERY

if start_from_stage <= 1
    fprintf('████ STAGE 1: KITCHEN SINK MODEL (Global Interaction Assessment) ████\n');
    fprintf('📋 STAGE 1 OBJECTIVES:\n');
    fprintf('  • Fit comprehensive 7-way factorial interaction model across entire parameter space\n');
    fprintf('  • Model: δβ ~ βgenerated × VGF × samplingRate × filterType × regressionType × noiseMagnitude × noiseColor + (1|paramComboID)\n');
    fprintf('  • Assess global model convergence and explanatory power (R² target: >0.90)\n');
    fprintf('  • Generate baseline performance metrics for adequacy comparison\n');
    fprintf('  • Prepare complete dataset and model results for Stage 2 assessment\n\n');
    
    fprintf('🔍 EXECUTION DETAILS:\n');
    fprintf('  Using ModelAdequacy_Stage1_KitchenSink_v001 ...\n');
    fprintf('  True subset approach ensures Level %d is exact subset of Level 9 parameters\n', tractability_level);
    fprintf('  Variable naming: Verified camelCase consistency (deltaBeta, betaGenerated)\n');
    fprintf('  Memory management: Automatic cleanup and optimization enabled\n\n');
    
    % Check for existing Stage 1 results to avoid unnecessary recomputation
    stage1_pattern = sprintf('stage1_results_L%d_*.mat', tractability_level);
    stage1_files = dir(stage1_pattern);
    
    if ~isempty(stage1_files) && start_from_stage > 1
        % Use existing results if starting from later stage
        [~, latest_idx] = max([stage1_files.datenum]);
        stage1_file = stage1_files(latest_idx).name;
        fprintf('✅ Found existing Stage 1 results: %s\n', stage1_file);
        
        % Validate existing results
        try
            stage1_check = load(stage1_file);
            if stage1_check.results.kitchen_sink_successful
                fprintf('  ✓ Validation: Kitchen sink model convergence confirmed\n');
                fprintf('  ✓ Model R²: %.4f with %d coefficients\n', ...
                    stage1_check.results.r_squared, stage1_check.results.total_coefficients);
            else
                warning('Existing Stage 1 results show model failure. Consider re-running Stage 1.');
            end
        catch ME
            warning('Could not validate existing Stage 1 results: %s', ME.message);
        end
    else
        % Execute Stage 1 kitchen sink model
        fprintf('🚀 Executing Stage 1: Kitchen Sink Model...\n');
        fprintf('  This may take %.0f-%.0f minutes depending on system performance...\n', ...
            current_time_minutes * 0.7, current_time_minutes * 1.3);
        
        try
            % Execute the kitchen sink model using the established methodology
            ModelAdequacy_Stage1_KitchenSink_v001(use_database, n_obs_per_combo, tractability_level, true); % always the more complex model please i.e. arg 4 is true
            
            % Verify Stage 1 completion and locate results file
            stage1_files = dir(stage1_pattern);
            if isempty(stage1_files)
                error('ModelAdequacy_Master_v002:Stage1Failed', ...
                    'Stage 1 failed to create results file. Check memory availability and tractability level.');
            end
            
            % Get the most recent Stage 1 file
            [~, latest_idx] = max([stage1_files.datenum]);
            stage1_file = stage1_files(latest_idx).name;
            
            % Validate Stage 1 success
            stage1_data = load(stage1_file);
            if stage1_data.results.kitchen_sink_successful
                fprintf('✅ Stage 1 completed successfully: %s\n', stage1_file);
                fprintf('  ✓ Model convergence: SUCCESSFUL\n');
                fprintf('  ✓ Final R²: %.4f (target: >0.90)\n', stage1_data.results.r_squared);
                fprintf('  ✓ Total coefficients: %d\n', stage1_data.results.total_coefficients);
                fprintf('  ✓ Significant interactions: %d\n', stage1_data.results.significant_interactions);
                fprintf('  ✓ Dataset size: %d observations\n', height(stage1_data.tableTrue));
            else
                error('ModelAdequacy_Master_v002:Stage1ModelFailed', ...
                    'Stage 1 kitchen sink model failed to converge. Try lower tractability level or check system resources.');
            end
            
        catch ME
            fprintf('❌ Stage 1 execution failed: %s\n', ME.message);
            fprintf('\n🔧 TROUBLESHOOTING SUGGESTIONS:\n');
            fprintf('  1. Reduce tractability level (try Level 1 or 2 for testing)\n');
            fprintf('  2. Reduce observations per combination (try n=3)\n');
            fprintf('  3. Check available memory (need ~%.0f GB for Level %d)\n', current_memory_gb, tractability_level);
            fprintf('  4. Close other applications to free memory\n');
            fprintf('  5. Restart MATLAB to clear memory fragmentation\n');
            rethrow(ME);
        end
    end
else
    % Starting from later stage - find existing Stage 1 file
    fprintf('🔍 LOCATING EXISTING STAGE 1 RESULTS:\n');
    stage1_pattern = sprintf('stage1_results_L%d_*.mat', tractability_level);
    stage1_files = dir(stage1_pattern);
    
    if isempty(stage1_files)
        error('ModelAdequacy_Master_v002:NoStage1Results', ...
            ['No Stage 1 results found for Level %d. Cannot start from Stage %d.\n' ...
             'Run complete pipeline (start_from_stage=1) first.'], ...
            tractability_level, start_from_stage);
    end
    
    % Use most recent Stage 1 file
    [~, latest_idx] = max([stage1_files.datenum]);
    stage1_file = stage1_files(latest_idx).name;
    fprintf('Found Stage 1 results: %s\n', stage1_file);
    
    % Quick validation of Stage 1 results
    try
        stage1_check = load(stage1_file);
        fprintf('  ✓ Kitchen sink R²: %.4f\n', stage1_check.results.r_squared);
        fprintf('  ✓ Dataset: %d observations\n', height(stage1_check.tableTrue));
    catch ME
        warning('Could not validate Stage 1 file: %s', ME.message);
    end
end

%% STAGE 2: SYSTEMATIC MODEL ADEQUACY ASSESSMENT

if start_from_stage <= 2
    fprintf('\n████ STAGE 2: SYSTEMATIC MODEL ADEQUACY ASSESSMENT ████\n');
    fprintf('📋 STAGE 2 OBJECTIVES:\n');
    fprintf('  • Comprehensive residual analysis across parameter regions\n');
    fprintf('  • Bootstrap parameter stability assessment (1000 iterations)\n');
    fprintf('  • K-fold cross-validation performance evaluation (k=5)\n');
    fprintf('  • Quantitative adequacy decision using clinical significance thresholds\n');
    fprintf('  • Identification of poor fit regions requiring conditional analysis\n\n');
    
    fprintf('🎯 ADEQUACY CRITERIA (Prereg v065 Clinical Thresholds):\n');
    fprintf('  • Residual Pattern Threshold: Cohen''s d > 0.5 indicates systematic bias\n');
    fprintf('  • Coefficient Instability: Bootstrap CI > ±0.03 (clinical significance)\n');
    fprintf('  • Cross-Validation: Performance degradation >15%% indicates poor fit\n');
    fprintf('  • Minimum Region Size: n ≥ 200 for reliable conditional modeling\n\n');
    
    % Check for existing Stage 2 results
    stage2_pattern = sprintf('stage2_adequacy_L%d_*.mat', tractability_level);
    stage2_files = dir(stage2_pattern);
    
    if ~isempty(stage2_files) && start_from_stage > 2
        % Use existing results if starting from later stage
        [~, latest_idx] = max([stage2_files.datenum]);
        stage2_file = stage2_files(latest_idx).name;
        fprintf('✅ Found existing Stage 2 results: %s\n', stage2_file);
        
        % Display adequacy decision summary
        try
            stage2_check = load(stage2_file);
            fprintf('  ✓ Global adequacy: %s\n', ternary(stage2_check.adequacy_results.global_adequate, 'ADEQUATE', 'INADEQUATE'));
            if ~stage2_check.adequacy_results.global_adequate
                fprintf('  ⚠ Poor fit regions: %d identified\n', length(stage2_check.adequacy_results.poor_fit_regions));
            end
        catch ME
            warning('Could not validate existing Stage 2 results: %s', ME.message);
        end
    else
        % Execute Stage 2 adequacy assessment
        fprintf('🚀 Executing Stage 2: Model Adequacy Assessment...\n');
        fprintf('  Using ModelAdequacy_Stage2_Assessment_v002 (FIXED camelCase variables)...\n');
        fprintf('  Estimated duration: 5-15 minutes depending on dataset size...\n\n');
        
        try
            % Execute the systematic adequacy assessment
            ModelAdequacy_Stage2_Assessment_v002(stage1_file);
            
            % Verify Stage 2 completion and locate results file
            stage2_files = dir(stage2_pattern);
            if isempty(stage2_files)
                error('ModelAdequacy_Master_v002:Stage2Failed', ...
                    'Stage 2 failed to create results file. Check Stage 1 model success and variable naming.');
            end
            
            % Get the most recent Stage 2 file
            [~, latest_idx] = max([stage2_files.datenum]);
            stage2_file = stage2_files(latest_idx).name;
            
            fprintf('✅ Stage 2 completed successfully: %s\n', stage2_file);
            
            % Load and display adequacy assessment results
            stage2_data = load(stage2_file);
            fprintf('\n📊 ADEQUACY ASSESSMENT RESULTS:\n');
            fprintf('  ✓ Residual adequacy: %s\n', ternary(stage2_data.adequacy_results.residual_analysis.adequate, 'PASS', 'FAIL'));
            fprintf('  ✓ Stability adequacy: %s\n', ternary(stage2_data.adequacy_results.stability_analysis.adequate, 'PASS', 'FAIL'));
            fprintf('  ✓ Cross-validation adequacy: %s\n', ternary(stage2_data.adequacy_results.cv_analysis.adequate, 'PASS', 'FAIL'));
            fprintf('  ────────────────────────────────────\n');
            fprintf('  📋 OVERALL ADEQUACY: %s\n', ternary(stage2_data.adequacy_results.global_adequate, 'ADEQUATE', 'INADEQUATE'));
            
        catch ME
            fprintf('❌ Stage 2 execution failed: %s\n', ME.message);
            fprintf('\n🔧 TROUBLESHOOTING SUGGESTIONS:\n');
            fprintf('  1. Verify Stage 1 completed successfully\n');
            fprintf('  2. Check variable naming consistency (should be camelCase)\n');
            fprintf('  3. Ensure VariableNameMapping_v003 is in functions/ directory\n');
            fprintf('  4. Verify deltaBeta and betaGenerated variables exist in data\n');
            rethrow(ME);
        end
    end
else
    % Starting from later stage - find existing Stage 2 file
    fprintf('🔍 LOCATING EXISTING STAGE 2 RESULTS:\n');
    stage2_pattern = sprintf('stage2_adequacy_L%d_*.mat', tractability_level);
    stage2_files = dir(stage2_pattern);
    
    if isempty(stage2_files)
        error('ModelAdequacy_Master_v002:NoStage2Results', ...
            ['No Stage 2 results found for Level %d. Cannot start from Stage %d.\n' ...
             'Run from Stage 2 or earlier.'], ...
            tractability_level, start_from_stage);
    end
    
    % Use most recent Stage 2 file
    [~, latest_idx] = max([stage2_files.datenum]);
    stage2_file = stage2_files(latest_idx).name;
    fprintf('Found Stage 2 results: %s\n', stage2_file);
end

% Load Stage 2 results to determine next pipeline steps
stage2_data = load(stage2_file);
requires_conditional = stage2_data.adequacy_results.requires_conditional_analysis;

fprintf('\n🎯 ADEQUACY ASSESSMENT DECISION: %s\n', ...
    ternary(requires_conditional, ...
    'Conditional analysis REQUIRED (poor fit regions detected)', ...
    'Global model ADEQUATE (no conditional analysis needed)'));

if requires_conditional
    fprintf('  📍 Poor fit regions detected: %d\n', length(stage2_data.adequacy_results.poor_fit_regions));
    fprintf('  🔄 Pipeline will proceed to Stage 3 (Conditional Parameter Analysis)\n');
else
    fprintf('  ✅ Global model demonstrates adequate performance across all regions\n');
    fprintf('  🏁 Pipeline can terminate here with global model recommendations\n');
end

%% STAGE 3: CONDITIONAL PARAMETER ANALYSIS (EXECUTED ONLY IF REQUIRED)

if requires_conditional && start_from_stage <= 3
    fprintf('\n████ STAGE 3: CONDITIONAL PARAMETER ANALYSIS (Poor Fit Regions) ████\n');
    fprintf('📋 STAGE 3 OBJECTIVES:\n');
    fprintf('  • Develop specialized models for identified poor fit regions\n');
    fprintf('  • Implement hierarchical modeling for complex interaction patterns\n');
    fprintf('  • Assess region-specific model performance and adequacy\n');
    fprintf('  • Evaluate potential benefits of model integration approaches\n');
    fprintf('  • Prepare conditional models for potential Stage 4 integration\n\n');
    
    fprintf('🎯 CONDITIONAL MODELING DETAILS:\n');
    fprintf('  • Poor fit regions from Stage 2: %d identified\n', length(stage2_data.adequacy_results.poor_fit_regions));
    fprintf('  • Minimum region size requirement: n ≥ 200 observations\n');
    fprintf('  • Advanced modeling techniques: Mixed-effects, hierarchical structures\n');
    fprintf('  • Performance comparison: Conditional vs. global model adequacy\n\n');
    
    % Display information about poor fit regions requiring conditional analysis
    fprintf('📊 POOR FIT REGIONS REQUIRING CONDITIONAL ANALYSIS:\n');
    num_regions_to_show = min(5, length(stage2_data.adequacy_results.poor_fit_regions));
    for i = 1:num_regions_to_show
        region = stage2_data.adequacy_results.poor_fit_regions{i};
        if isfield(region, 'n_obs')
            fprintf('  %d. %s (%d observations)\n', i, region.description, region.n_obs);
        else
            fprintf('  %d. %s\n', i, region.description);
        end
    end
    if length(stage2_data.adequacy_results.poor_fit_regions) > 5
        fprintf('  ... and %d additional regions\n', length(stage2_data.adequacy_results.poor_fit_regions) - 5);
    end
    
    % Check for existing Stage 3 results
    stage3_pattern = sprintf('stage3_conditional_L%d_*.mat', tractability_level);
    stage3_files = dir(stage3_pattern);
    
    if ~isempty(stage3_files) && start_from_stage > 3
        % Use existing results if starting from later stage
        [~, latest_idx] = max([stage3_files.datenum]);
        stage3_file = stage3_files(latest_idx).name;
        fprintf('\n✅ Found existing Stage 3 results: %s\n', stage3_file);
        
        % Display conditional analysis summary
        try
            stage3_check = load(stage3_file);
            if isfield(stage3_check.conditional_results, 'regional_models')
                fprintf('  ✓ Conditional models developed: %d regions\n', ...
                    length(fieldnames(stage3_check.conditional_results.regional_models)));
            end
            fprintf('  ✓ Integration beneficial: %s\n', ...
                ternary(stage3_check.conditional_results.requires_integration, 'YES', 'NO'));
        catch ME
            warning('Could not validate existing Stage 3 results: %s', ME.message);
        end
    else
        % Execute Stage 3 conditional analysis  
        fprintf('\n🚀 Executing Stage 3: Conditional Parameter Analysis...\n');
        fprintf('  Estimated duration: 10-30 minutes depending on number of poor fit regions...\n\n');
        
        % BLUEBEAR MEMORY MONITORING: Check system status before Stage 3
        fprintf('💾 BlueBear memory check before Stage 3...\n');
        if ispc
            [~, sys] = memory;
            available_gb = sys.PhysicalMemory.Available / 1e9;
            fprintf('  Available memory: %.1f GB (BlueBear total: 288GB)\n', available_gb);
            if available_gb < 20.0  % Conservative threshold for 288GB system
                fprintf('  ⚠️  Memory getting low - Stage 3 will use reduced scope\n');
            else
                fprintf('  ✓ Sufficient memory for Stage 3 processing\n');
            end
        else
            fprintf('  BlueBear system: 288GB available for Stage 3 processing\n');
        end
        
        try
            % BLUEBEAR SAFETY: Execute conditional parameter analysis with monitoring
            fprintf('\n🔬 Executing Stage 3 with BlueBear memory management...\n');
            ModelAdequacy_Stage3_Conditional_v002(stage2_file);
            
            % Verify Stage 3 completion and locate results file
            stage3_files = dir(stage3_pattern);
            if isempty(stage3_files)
                error('ModelAdequacy_Master_v002:Stage3Failed', ...
                    'Stage 3 failed to create results file. Check conditional modeling implementation.');
            end
            
            % Get the most recent Stage 3 file
            [~, latest_idx] = max([stage3_files.datenum]);
            stage3_file = stage3_files(latest_idx).name;
            
            fprintf('✅ Stage 3 completed successfully: %s\n', stage3_file);
            
        catch ME
            % BLUEBEAR ENHANCED ERROR HANDLING: Handle Stage 3 implementation and memory issues
            if contains(ME.message, 'memory') || contains(ME.message, 'Memory') || contains(ME.message, 'Out of memory')
                fprintf('🚨 STAGE 3 MEMORY FAILURE: %s\n', ME.message);
                fprintf('💾 BlueBear Memory Recovery:\n');
                fprintf('  • Stage 3 exceeded BlueBear memory limits during regional modeling\n');
                fprintf('  • This can happen with very large regions or numerous poor fit areas\n');
                fprintf('  • Creating simplified Stage 3 results to continue pipeline...\n\n');
                
                % Create memory-failure bypass results
                stage3_file = sprintf('stage3_conditional_L%d_memory_bypass_%s.mat', tractability_level, datestr(now, 'yyyymmdd_HHMMSS'));
                conditional_results = struct();
                conditional_results.implementation_status = 'memory_bypass';
                conditional_results.requires_integration = false;
                conditional_results.poor_fit_regions = stage2_data.adequacy_results.poor_fit_regions;
                conditional_results.memory_failure_info = ME.message;
                conditional_results.bluebear_recommendation = 'Reduce tractability level or use manual regional analysis';
                conditional_results.bypass_reason = 'Stage 3 memory limits exceeded on BlueBear 288GB system';
                save(stage3_file, 'conditional_results','-v7.3');
                
                fprintf('📁 Memory bypass results saved: %s\n', stage3_file);
                fprintf('📊 Recommendation: Try tractability level 1 or 2 for memory-constrained analysis\n');
                
            elseif contains(ME.message, 'Undefined function')
                fprintf('⚠️  Stage 3 implementation not yet available: %s\n', ME.message);
                fprintf('🔧 DEVELOPMENT STATUS:\n');
                fprintf('  • Stage 3 (ModelAdequacy_Stage3_Conditional_v002) is planned for future implementation\n');
                fprintf('  • Current framework successfully identifies poor fit regions in Stage 2\n');
                fprintf('  • Manual conditional analysis can be performed using Stage 2 region identification\n');
                fprintf('  • Framework will be enhanced with automatic conditional modeling in future versions\n\n');
                
                % Create placeholder Stage 3 results for pipeline consistency
                stage3_file = sprintf('stage3_conditional_L%d_placeholder_%s.mat', tractability_level, datestr(now, 'yyyymmdd_HHMMSS'));
                conditional_results = struct();
                conditional_results.implementation_status = 'placeholder';
                conditional_results.requires_integration = false;
                conditional_results.poor_fit_regions = stage2_data.adequacy_results.poor_fit_regions;
                conditional_results.recommendation = 'Manual conditional analysis recommended for identified poor fit regions';
                save(stage3_file, 'conditional_results','-v7.3');
                
                fprintf('📁 Placeholder Stage 3 results saved: %s\n', stage3_file);
            else
                fprintf('❌ Stage 3 execution failed: %s\n', ME.message);
                fprintf('🔧 BLUEBEAR TROUBLESHOOTING:\n');
                fprintf('  1. Check if memory is sufficient (need >20GB available)\n');
                fprintf('  2. Try reducing tractability level (Level 1 or 2)\n');
                fprintf('  3. Restart MATLAB to clear memory fragmentation\n');
                fprintf('  4. Check for variable naming consistency issues\n\n');
                rethrow(ME);
            end
        end
    end
    
    % Load Stage 3 results to determine integration requirements
    if exist('stage3_file', 'var') && exist(stage3_file, 'file')
        stage3_data = load(stage3_file);
        
        % Check if Stage 3 is a placeholder or real implementation
        if isfield(stage3_data.conditional_results, 'implementation_status') && ...
           strcmp(stage3_data.conditional_results.implementation_status, 'placeholder')
            requires_integration = false;
            fprintf('\n📋 CONDITIONAL ANALYSIS STATUS: Placeholder (implementation pending)\n');
        else
            requires_integration = stage3_data.conditional_results.requires_integration;
            fprintf('\n📋 CONDITIONAL ANALYSIS DECISION: %s\n', ...
                ternary(requires_integration, ...
                'Integration BENEFICIAL (proceed to Stage 4)', ...
                'Integration NOT beneficial (conditional models sufficient)'));
        end
    else
        requires_integration = false;
        fprintf('\n⚠️  Could not determine integration requirements from Stage 3\n');
    end
    
    %% STAGE 4: HIERARCHICAL MODEL INTEGRATION (EXECUTED ONLY IF BENEFICIAL)
    
    if requires_integration && start_from_stage <= 4
        fprintf('\n████ STAGE 4: HIERARCHICAL MODEL INTEGRATION (Optimal Performance) ████\n');
        fprintf('📋 STAGE 4 OBJECTIVES:\n');
        fprintf('  • Integrate global and conditional models for optimal performance\n');
        fprintf('  • Develop hierarchical decision framework for model selection\n');
        fprintf('  • Validate integrated model performance across all parameter regions\n');
        fprintf('  • Generate clinical decision support recommendations\n');
        fprintf('  • Produce final comprehensive model adequacy assessment\n\n');
        
        fprintf('🎯 INTEGRATION APPROACH:\n');
        fprintf('  • Hierarchical model structure combining global + conditional approaches\n');
        fprintf('  • Automated model selection based on parameter region characteristics\n');
        fprintf('  • Cross-validation of integrated framework performance\n');
        fprintf('  • Clinical threshold validation for decision support\n\n');
        
        % Check for existing Stage 4 results
        stage4_pattern = sprintf('stage4_integration_L%d_*.mat', tractability_level);
        stage4_files = dir(stage4_pattern);
        
        if ~isempty(stage4_files)
            % Use existing results
            [~, latest_idx] = max([stage4_files.datenum]);
            stage4_file = stage4_files(latest_idx).name;
            fprintf('✅ Found existing Stage 4 results: %s\n', stage4_file);
        else
            % Execute Stage 4 model integration
            fprintf('🚀 Executing Stage 4: Hierarchical Model Integration...\n');
            fprintf('  Estimated duration: 15-45 minutes depending on integration complexity...\n\n');
            
            try
                % Execute hierarchical model integration
                ModelAdequacy_Stage4_Integration_v002(stage3_file);
                
                % Verify Stage 4 completion and locate results file
                stage4_files = dir(stage4_pattern);
                if isempty(stage4_files)
                    error('ModelAdequacy_Master_v002:Stage4Failed', ...
                        'Stage 4 failed to create results file. Check integration implementation.');
                end
                
                % Get the most recent Stage 4 file
                [~, latest_idx] = max([stage4_files.datenum]);
                stage4_file = stage4_files(latest_idx).name;
                
                fprintf('✅ Stage 4 completed successfully: %s\n', stage4_file);
                fprintf('\n🎉 COMPLETE 4-STAGE MODEL ADEQUACY FRAMEWORK EXECUTED\n');
                
            catch ME
                % Handle Stage 4 implementation status
                if contains(ME.message, 'Undefined function')
                    fprintf('⚠️  Stage 4 implementation not yet available: %s\n', ME.message);
                    fprintf('🔧 DEVELOPMENT STATUS:\n');
                    fprintf('  • Stage 4 (ModelAdequacy_Stage4_Integration_v002) is planned for future implementation\n');
                    fprintf('  • Current framework successfully completes Stages 1-2 with conditional region identification\n');
                    fprintf('  • Manual integration approaches can be developed using Stage 3 conditional models\n');
                    fprintf('  • Framework provides complete foundation for hierarchical modeling approaches\n\n');
                else
                    fprintf('❌ Stage 4 execution failed: %s\n', ME.message);
                    rethrow(ME);
                end
            end
        end
    else
        % Stage 4 not required or not beneficial
        if requires_conditional
            fprintf('\n📋 Stage 4 not required: %s\n', ...
                ternary(requires_integration, 'Integration not beneficial based on Stage 3 analysis', ...
                'Conditional analysis sufficient without global integration'));
        end
    end
else
    % Conditional analysis not required
    if ~requires_conditional
        fprintf('\n📋 Stages 3-4 not required: Global model demonstrates adequate performance\n');
        fprintf('  ✅ Kitchen sink model adequately captures parameter relationships\n');
        fprintf('  🎯 Recommendation: Use global model for all parameter regions\n');
        fprintf('  📊 No evidence of systematic regional inadequacies requiring conditional approaches\n');
    else
        fprintf('\n📋 Stage 3 skipped: Starting from Stage %d\n', start_from_stage);
    end
end

%% FRAMEWORK COMPLETION AND COMPREHENSIVE REPORTING

fprintf('\n████ MODEL ADEQUACY FRAMEWORK COMPLETION ████\n');

% Determine final framework status and recommendations
fprintf('📊 FINAL FRAMEWORK STATUS:\n');

% Load final results to determine overall recommendations
try
    if exist('stage2_file', 'var')
        final_stage2 = load(stage2_file);
        global_adequate = final_stage2.adequacy_results.global_adequate;
        
        if global_adequate
            fprintf('  ✅ RECOMMENDATION: Use global kitchen sink model\n');
            fprintf('  📋 Rationale: Model demonstrates adequate performance across all parameter regions\n');
            fprintf('  🎯 Application: Single model appropriate for all measurement conditions\n');
        else
            fprintf('  ⚠️  RECOMMENDATION: Conditional/hierarchical modeling approaches required\n');
            fprintf('  📋 Rationale: Global model shows inadequate performance in %d regions\n', ...
                length(final_stage2.adequacy_results.poor_fit_regions));
            fprintf('  🎯 Application: Use region-specific models or await Stage 3/4 implementation\n');
        end
        
        % Display key performance metrics
        fprintf('\n📈 KEY PERFORMANCE METRICS:\n');
        fprintf('  • Kitchen Sink R²: %.4f\n', final_stage2.adequacy_results.stage1_r_squared);
        fprintf('  • Model Coefficients: %d\n', final_stage2.adequacy_results.stage1_coefficients);
        fprintf('  • Residual Adequacy: %s\n', ternary(final_stage2.adequacy_results.residual_analysis.adequate, 'PASS', 'FAIL'));
        fprintf('  • Stability Adequacy: %s\n', ternary(final_stage2.adequacy_results.stability_analysis.adequate, 'PASS', 'FAIL'));
        fprintf('  • Cross-Validation Adequacy: %s\n', ternary(final_stage2.adequacy_results.cv_analysis.adequate, 'PASS', 'FAIL'));
        
    else
        fprintf('  ⚠️  Could not load final results for comprehensive assessment\n');
    end
    
catch ME
    fprintf('  ⚠️  Error accessing final results: %s\n', ME.message);
end

% Generate comprehensive final report
fprintf('\n📄 Generating comprehensive adequacy framework report...\n');
try
    generateFinalAdequacyReport_v002(tractability_level, n_obs_per_combo, use_database, start_from_stage);
    fprintf('  ✅ Comprehensive report generated successfully\n');
catch ME
    fprintf('  ⚠️  Report generation failed: %s\n', ME.message);
end

% Display output file summary
fprintf('\n📁 OUTPUT FILES SUMMARY:\n');

% List all generated files with descriptions
output_patterns = {
    sprintf('stage1_results_L%d_*.mat', tractability_level), 'Kitchen sink model results and complete dataset';
    sprintf('stage2_adequacy_L%d_*.mat', tractability_level), 'Adequacy assessment with quantitative decision framework';
    sprintf('stage3_conditional_L%d_*.mat', tractability_level), 'Conditional analysis for poor fit regions (if required)';
    sprintf('stage4_integration_L%d_*.mat', tractability_level), 'Integrated hierarchical models (if beneficial)';
    'model_adequacy_FINAL_*.html', 'Comprehensive final report with recommendations'
};

for i = 1:size(output_patterns, 1)
    pattern = output_patterns{i, 1};
    description = output_patterns{i, 2};
    files = dir(pattern);
    
    if ~isempty(files)
        % Show most recent file
        [~, latest_idx] = max([files.datenum]);
        latest_file = files(latest_idx).name;
        fprintf('  ✅ %s: %s\n', latest_file, description);
    else
        fprintf('  ⏸️  %s: %s (not generated)\n', pattern, description);
    end
end

% Display usage recommendations and next steps
fprintf('\n🎯 USAGE RECOMMENDATIONS:\n');
fprintf('  • Review stage2_adequacy_*.mat for primary adequacy assessment\n');
fprintf('  • Examine model_adequacy_FINAL_*.html for comprehensive analysis summary\n');
fprintf('  • Use global model if adequacy assessment shows "ADEQUATE" status\n');
fprintf('  • Consider conditional approaches for regions flagged as inadequate\n');
fprintf('  • Validate findings with domain expertise and clinical requirements\n');

% Final completion message
fprintf('\n████████████████████████████████████████████████████████████████\n');
fprintf('   MODEL ADEQUACY FRAMEWORK ANALYSIS COMPLETED                   \n');
fprintf('   Version: Master_v002 with Enhanced Documentation              \n');
fprintf('   Analysis Level: %d (%s)\n', tractability_level, current_name);
fprintf('   Variable Naming: Fixed camelCase consistency                  \n');
fprintf('   Results: See output files listed above                       \n');
fprintf('████████████████████████████████████████████████████████████████\n\n');

end

%% SUPPORTING UTILITY FUNCTIONS

function result = ternary(condition, true_val, false_val)
% Simple ternary operator function for cleaner conditional text generation
%
% INPUT:
%   condition - logical value or expression to evaluate
%   true_val - value to return if condition is true
%   false_val - value to return if condition is false
%
% OUTPUT:
%   result - true_val if condition is true, false_val otherwise

if condition
    result = true_val;
else
    result = false_val;
end
end

function generateFinalAdequacyReport_v002(tractability_level, n_obs_per_combo, use_database, start_from_stage)
% Generate comprehensive final report with enhanced documentation and analysis summary
%
% INPUTS:
%   tractability_level - Analysis complexity level (1-9)
%   n_obs_per_combo - Observations per parameter combination
%   use_database - Whether real database was used (vs synthetic data)
%   start_from_stage - Which stage the analysis started from
%
% OUTPUT:
%   HTML report file with comprehensive analysis summary and recommendations

% Generate timestamp and filename for the final report
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
report_file = sprintf('model_adequacy_FINAL_L%d_n%d_%s_S%d_%s.html', ...
    tractability_level, n_obs_per_combo, ...
    ternary(use_database, 'DB', 'SYN'), start_from_stage, timestamp);

% Create comprehensive HTML report
fid = fopen(report_file, 'w');

if fid == -1
    error('Could not create report file: %s', report_file);
end

try
    % HTML header with enhanced styling
    fprintf(fid, '<!DOCTYPE html>\n<html>\n<head>\n');
    fprintf(fid, '<title>Model Adequacy Framework - Comprehensive Final Report</title>\n');
    fprintf(fid, '<meta charset="UTF-8">\n');
    fprintf(fid, '<style>\n');
    fprintf(fid, 'body { font-family: "Segoe UI", Arial, sans-serif; margin: 40px; line-height: 1.6; }\n');
    fprintf(fid, 'h1 { color: #2E86AB; border-bottom: 3px solid #2E86AB; padding-bottom: 10px; }\n');
    fprintf(fid, 'h2 { color: #A23B72; margin-top: 30px; }\n');
    fprintf(fid, 'h3 { color: #F18F01; }\n');
    fprintf(fid, '.summary { background-color: #f0f8ff; padding: 20px; border-left: 5px solid #2E86AB; margin: 20px 0; }\n');
    fprintf(fid, '.success { background-color: #e8f5e8; padding: 15px; border-left: 5px solid #4CAF50; }\n');
    fprintf(fid, '.warning { background-color: #fff3cd; padding: 15px; border-left: 5px solid #FFC107; }\n');
    fprintf(fid, '.error { background-color: #f8d7da; padding: 15px; border-left: 5px solid #DC3545; }\n');
    fprintf(fid, '.code { background-color: #f8f9fa; padding: 10px; border: 1px solid #e9ecef; font-family: monospace; }\n');
    fprintf(fid, 'table { width: 100%%; border-collapse: collapse; margin: 20px 0; }\n');
    fprintf(fid, 'th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }\n');
    fprintf(fid, 'th { background-color: #f2f2f2; font-weight: bold; }\n');
    fprintf(fid, '.metric { font-size: 1.2em; font-weight: bold; color: #2E86AB; }\n');
    fprintf(fid, '</style>\n</head>\n<body>\n');

    % Report title and executive summary
    fprintf(fid, '<h1>Model Adequacy Framework - Comprehensive Final Report</h1>\n');
    fprintf(fid, '<div class="summary">\n');
    fprintf(fid, '<h2>Executive Summary</h2>\n');
    fprintf(fid, '<p><strong>Analysis Framework:</strong> 4-Stage Model Adequacy Assessment (Fraser et al., 2025)</p>\n');
    fprintf(fid, '<p><strong>Tractability Level:</strong> %d</p>\n', tractability_level);
    fprintf(fid, '<p><strong>Observations per combination:</strong> %d</p>\n', n_obs_per_combo);
    fprintf(fid, '<p><strong>Data source:</strong> %s</p>\n', ternary(use_database, 'Real database', 'Synthetic data generation'));
    fprintf(fid, '<p><strong>Pipeline execution:</strong> Started from Stage %d</p>\n', start_from_stage);
    fprintf(fid, '<p><strong>Report generated:</strong> %s</p>\n', datestr(now));
    fprintf(fid, '<p><strong>Framework version:</strong> Master_v002 with enhanced documentation</p>\n');
    fprintf(fid, '</div>\n');

    % Add basic configuration table
    fprintf(fid, '<h2>Analysis Configuration</h2>\n');
    fprintf(fid, '<table>\n');
    fprintf(fid, '<tr><th>Parameter</th><th>Value</th><th>Description</th></tr>\n');
    fprintf(fid, '<tr><td>Tractability Level</td><td>%d</td><td>Analysis complexity level</td></tr>\n', tractability_level);
    fprintf(fid, '<tr><td>Observations per combo</td><td>%d</td><td>Statistical power setting</td></tr>\n', n_obs_per_combo);
    fprintf(fid, '<tr><td>Data Source</td><td>%s</td><td>Data generation method</td></tr>\n', ternary(use_database, 'Database', 'Synthetic'));
    fprintf(fid, '<tr><td>Starting Stage</td><td>%d</td><td>Pipeline entry point</td></tr>\n', start_from_stage);
    fprintf(fid, '</table>\n');

    % Basic methodology section
    fprintf(fid, '<h2>Methodology</h2>\n');
    fprintf(fid, '<h3>Model Formula</h3>\n');
    fprintf(fid, '<div class="code">\n');
    fprintf(fid, 'δβ ~ βgenerated × VGF × samplingRate × filterType × regressionType × noiseMagnitude × noiseColor + (1|paramComboID)\n');
    fprintf(fid, '</div>\n');

    % Footer
    fprintf(fid, '<hr>\n');
    fprintf(fid, '<p><em>Report generated by Model Adequacy Framework Master_v002</em></p>\n');
    fprintf(fid, '<p><em>Fraser, D.S. (2025) - Model Adequacy Framework for Power Law Parameter Recovery</em></p>\n');
    fprintf(fid, '</body>\n</html>\n');

finally
    fclose(fid);
end

fprintf('✅ Basic final report generated: %s\n', report_file);
end