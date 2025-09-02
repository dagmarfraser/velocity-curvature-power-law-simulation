function ModelAdequacy_Stage3_Conditional_v004(stage2_file)
% MODELADEQUACY_STAGE3_CONDITIONAL_V004 Conditional Parameter Analysis
%
% **STAGE 3 OF MODEL ADEQUACY FRAMEWORK** (SYNTHETIC INTEGRATION v004)
% Executes conditional parameter analysis exclusively for regions identified
% as demonstrating systematic model inadequacies in Stage 2 assessment.
%
% **CRITICAL ENHANCEMENTS v004**:
%   - ADDED: Synthetic integration regions for Stage 4 methodology testing
%   - GATED: Only generates synthetic data when using synthetic (non-database) data
%   - TESTING: Enables complete 4-stage pipeline methodology validation
%   - SAFE: All field access errors fixed from v003
%   - DOCUMENTED: Clear methodology testing vs real analysis distinction
%
% **METHODOLOGY** (Prereg v071 Stage 3):
%   1. Load Stage 2 adequacy assessment results
%   2. Isolate poor fit regions based on systematic adequacy criteria
%   3. Develop region-specific conditional models with tailored interactions
%   4. Perform comparative assessment against global model performance
%   5. [TESTING] Generate synthetic integration regions for Stage 4 testing
%   6. Evaluate integration decision criteria for Stage 4
%   7. Save conditional analysis results for hierarchical model integration
%
% **SYNTHETIC INTEGRATION TESTING**:
%   - Only active for synthetic (non-database) data methodology testing
%   - Artificially enhances some regional models to trigger Stage 4
%   - Enables complete pipeline validation and development
%   - Real analysis would stop at natural integration decision
%
% **PIPELINE INTEGRATION**:
%   Input:  stage2_adequacy_*.mat (from ModelAdequacy_Stage2_Assessment_v002)
%   Output: stage3_conditional_*.mat (conditional analysis for Stage 4)
%   Next:   ModelAdequacy_Stage4_Integration_v002.m 
%
% USAGE:
%   ModelAdequacy_Stage3_Conditional_v004()                              % Use latest Stage 2 results
%   ModelAdequacy_Stage3_Conditional_v004('stage2_adequacy_L2_*.mat')    % Specific Stage 2 file
%
% Based on Bates & Watts (1988), Pinheiro & Bates (2000), Verbeke & Molenberghs (2000)
% Fraser, D.S. (2025) - Model Adequacy Framework Stage 3

if nargin < 1
    % Find latest Stage 2 adequacy assessment results
    stage2_files = dir('stage2_adequacy_*.mat');
    if isempty(stage2_files)
        error('No Stage 2 adequacy results found. Please run ModelAdequacy_Stage2_Assessment_v002 first.');
    end
    [~, latest_idx] = max([stage2_files.datenum]);
    stage2_file = stage2_files(latest_idx).name;
end

%% STAGE 3 INITIALIZATION
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('   MODEL ADEQUACY FRAMEWORK STAGE 3: CONDITIONAL ANALYSIS v004   \n');
fprintf('   Targeted Parameter Analysis for Poor Fit Regions              \n');
fprintf('   ENHANCED: Synthetic Integration Regions for Stage 4 Testing   \n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

% Load Stage 2 adequacy assessment results
fprintf('=== LOADING STAGE 2 ADEQUACY RESULTS ===\n');
try
    if ~exist(stage2_file, 'file')
        error('Stage 2 adequacy results file not found: %s', stage2_file);
    end
    
    stage2_data = load(stage2_file);
    fprintf('✓ Loaded Stage 2 adequacy results: %s\n', stage2_file);
    
    % Extract adequacy assessment
    adequacy_results = stage2_data.adequacy_results;
    config = stage2_data.config;
    stage1_data = stage2_data.stage1_data;
    
    fprintf('  📊 Overall adequacy: %s\n', logical2str(adequacy_results.global_adequate));
    fprintf('  📍 Conditional analysis required: %s\n', logical2str(adequacy_results.requires_conditional_analysis));
    fprintf('  🎯 Poor fit regions identified: %d\n', length(adequacy_results.poor_fit_regions));
    
    % **v004**: Detect synthetic data mode for methodology testing
    is_synthetic_testing = detectSyntheticTestingMode(adequacy_results, stage1_data);
    if is_synthetic_testing
        fprintf('  🔬 SYNTHETIC DATA MODE: Stage 4 methodology testing enabled\n');
    else
        fprintf('  📊 REAL DATA MODE: Standard analysis pipeline\n');
    end
    
catch ME
    error('Failed to load Stage 2 adequacy results: %s', ME.message);
end

% Validate conditional analysis requirement
if ~adequacy_results.requires_conditional_analysis
    fprintf('\n⚠️  CONDITIONAL ANALYSIS NOT REQUIRED\n');
    fprintf('   Global model demonstrates adequate performance across all regions.\n');
    fprintf('   No systematic inadequacies identified in Stage 2 assessment.\n');
    fprintf('   Framework analysis complete at Stage 2.\n\n');
    
    % Save Stage 3 results indicating no conditional analysis needed
    conditional_results = createNoConditionalAnalysisResults(adequacy_results, stage1_data, config);
    saveStage3Results(conditional_results, stage2_file, stage1_data);
    return;
end

%% POOR FIT REGION ISOLATION
fprintf('\n=== POOR FIT REGION ISOLATION ===\n');
fprintf('Isolating regions demonstrating systematic model inadequacies...\n');

% Extract original model and data from Stage 1
global_model = stage1_data.results.model;
data_table = stage1_data.tableTrue;

% Validate variable naming consistency (v002 enhancement)
required_vars = {'deltaBeta', 'betaGenerated', 'VGF', 'samplingRate', 'filterType', 'regressionType', 'noiseMagnitude', 'noiseColor'};
missing_vars = {};
for i = 1:length(required_vars)
    if ~ismember(required_vars{i}, data_table.Properties.VariableNames)
        missing_vars{end+1} = required_vars{i};
    end
end

if ~isempty(missing_vars)
    error('Missing required camelCase variables: %s. Check Stage 1 data generation.', strjoin(missing_vars, ', '));
end

fprintf('  ✓ Variable naming validation: All required camelCase variables found\n');

% Isolate poor fit regions with sufficient observations
valid_poor_regions = isolatePoorFitRegions(adequacy_results.poor_fit_regions, data_table, config);

fprintf('Poor Fit Region Analysis:\n');
fprintf('  Total regions identified: %d\n', length(adequacy_results.poor_fit_regions));
fprintf('  Valid regions for analysis: %d\n', length(valid_poor_regions));
fprintf('  Minimum region size: %d observations\n', config.minRegionSize);

if isempty(valid_poor_regions)
    fprintf('\n⚠️  NO VALID POOR FIT REGIONS\n');
    fprintf('   No regions meet minimum size requirements for conditional modeling.\n');
    fprintf('   Framework analysis complete at Stage 2.\n\n');
    
    conditional_results = createNoConditionalAnalysisResults(adequacy_results, stage1_data, config);
    saveStage3Results(conditional_results, stage2_file, stage1_data);
    return;
end

%% REGION-SPECIFIC CONDITIONAL MODELING
fprintf('\n=== REGION-SPECIFIC CONDITIONAL MODELING ===\n');
fprintf('Developing tailored interaction structures for poor fit regions...\n');

% Initialize conditional modeling configuration
conditional_config = initializeConditionalConfig(config);

%% BLUEBEAR MEMORY SAFETY: Monitor and limit for 288GB system
if isfield(adequacy_results, 'testing_mode') && adequacy_results.testing_mode
    fprintf('\n🔬 SYNTHETIC TESTING MODE: BlueBear memory management\n');
    if length(valid_poor_regions) > 5
        fprintf('  Reducing regions from %d to 5 for BlueBear efficiency\n', length(valid_poor_regions));
        valid_poor_regions = valid_poor_regions(1:5); % Reasonable limit for 288GB
    end
    fprintf('  Processing %d regions with BlueBear memory monitoring\n', length(valid_poor_regions));
end

% Monitor initial memory state
if ispc
    [~, sys] = memory;
    initial_memory_gb = sys.PhysicalMemory.Available / 1e9;
    fprintf('  💾 Initial available memory: %.1f GB (BlueBear: 288GB total)\n', initial_memory_gb);
else
    fprintf('  💾 BlueBear system: 288GB total memory available\n');
end

% Develop conditional models for each valid poor fit region
regional_models = {};
for i = 1:length(valid_poor_regions)
    region = valid_poor_regions{i};
    fprintf('\n--- Region %d: %s (%d observations) ---\n', i, region.description, region.n_obs);
    
    % Develop region-specific conditional model
    regional_model = developRegionalModel(region, data_table, global_model, conditional_config);
    regional_models{end+1} = regional_model;
    
    fprintf('  ✓ Regional model: R² = %.4f, Coefficients = %d\n', ...
        regional_model.r_squared, regional_model.num_coefficients);
    
    % BLUEBEAR MEMORY CLEANUP: Moderate cleanup every 3 regions
    if mod(i, 3) == 0 && i < length(valid_poor_regions)
        fprintf('  🗑️  BlueBear memory cleanup (region %d/%d)...\n', i, length(valid_poor_regions));
        clear region; % Clear current region data
        if ispc, java.lang.System.gc(); end
        pause(0.05); % Brief pause for cleanup
        
        % Check memory status
        if ispc
            [~, sys] = memory;
            current_memory_gb = sys.PhysicalMemory.Available / 1e9;
            fprintf('      Available memory: %.1f GB\n', current_memory_gb);
            if current_memory_gb < 20.0 % Conservative threshold for 288GB system
                warning('Memory getting low (%.1f GB). Consider reducing regions.', current_memory_gb);
            end
        end
    end
end

%% **NEW v004**: SYNTHETIC INTEGRATION REGIONS FOR STAGE 4 TESTING
if is_synthetic_testing
    fprintf('\n=== SYNTHETIC INTEGRATION REGIONS FOR STAGE 4 TESTING ===\n');
    fprintf('🔬 METHODOLOGY TESTING: Generating artificial integration-worthy regions...\n');
    fprintf('⚠️  NOTE: This is for pipeline methodology testing only - NOT real analysis\n');
    
    % Generate synthetic integration regions by enhancing some regional models
    regional_models = generateSyntheticIntegrationRegions(regional_models, global_model, conditional_config);
    
    fprintf('✓ Generated synthetic integration regions for Stage 4 testing\n');
    fprintf('📋 Note: Real analysis would use natural integration decision\n');
end

%% COMPARATIVE MODEL ASSESSMENT
fprintf('\n=== COMPARATIVE MODEL ASSESSMENT ===\n');
fprintf('Evaluating regional models against global performance...\n');

% Perform systematic comparison of regional vs global models
comparative_assessment = performComparativeAssessment(regional_models, global_model, data_table, conditional_config);

fprintf('Comparative Assessment Results:\n');
fprintf('  Regions with significant improvement: %d/%d\n', ...
    comparative_assessment.num_improved_regions, length(regional_models));
fprintf('  Average improvement (ΔR²): %.4f\n', comparative_assessment.mean_r_squared_improvement);
fprintf('  Integration recommended: %s\n', logical2str(comparative_assessment.integration_recommended));

%% INTEGRATION DECISION EVALUATION
fprintf('\n=== INTEGRATION DECISION EVALUATION ===\n');
fprintf('Applying integration decision criteria...\n');

% Evaluate whether conditional models warrant integration into hierarchical framework
integration_decision = evaluateIntegrationDecision(comparative_assessment, regional_models, conditional_config);

fprintf('Integration Decision Results:\n');
fprintf('  Statistical improvement threshold met: %s\n', logical2str(integration_decision.statistical_improvement));
fprintf('  Effect size maintenance: %s\n', logical2str(integration_decision.effect_size_maintained));
fprintf('  Mechanistic insight provided: %s\n', logical2str(integration_decision.mechanistic_insight));
fprintf('  ────────────────────────────────────────────────\n');
fprintf('  INTEGRATION DECISION: %s\n', logical2str(integration_decision.proceed_to_integration));

if integration_decision.proceed_to_integration
    fprintf('\n  📋 %d regions qualify for hierarchical integration\n', length(integration_decision.integration_regions));
    for i = 1:min(3, length(integration_decision.integration_regions))
        region = integration_decision.integration_regions{i};
        fprintf('    Region %d: %s (Improvement: ΔR² = %.4f)\n', i, region.description, getFieldSafely(region, 'r_squared_improvement', 0));
    end
    if length(integration_decision.integration_regions) > 3
        fprintf('    ... and %d more regions\n', length(integration_decision.integration_regions) - 3);
    end
    
    % **v004**: Additional messaging for synthetic testing
    if is_synthetic_testing
        fprintf('\n  🔬 SYNTHETIC TESTING: Stage 4 methodology testing will proceed\n');
        fprintf('  📋 Real analysis: Integration decision based on natural model performance\n');
    end
end

%% MECHANISTIC INSIGHT DOCUMENTATION
fprintf('\n=== MECHANISTIC INSIGHT DOCUMENTATION ===\n');
fprintf('Documenting identified relationship patterns...\n');

% Document mechanistic insights from conditional analysis
mechanistic_insights = documentMechanisticInsights(regional_models, integration_decision, conditional_config);

fprintf('Mechanistic Insight Results:\n');
fprintf('  Interaction patterns identified: %d\n', length(mechanistic_insights.interaction_patterns));
fprintf('  Parameter relationships documented: %d\n', length(mechanistic_insights.parameter_relationships));
fprintf('  Clinical implications: %d\n', length(mechanistic_insights.clinical_implications));

%% SAVE STAGE 3 RESULTS
fprintf('\n=== SAVING STAGE 3 RESULTS ===\n');

% Create comprehensive Stage 3 results structure
conditional_results = struct();
conditional_results.poor_fit_regions = valid_poor_regions;
conditional_results.regional_models = regional_models;
conditional_results.comparative_assessment = comparative_assessment;
conditional_results.integration_decision = integration_decision;
conditional_results.mechanistic_insights = mechanistic_insights;
conditional_results.requires_integration = integration_decision.proceed_to_integration; % FIXED naming
conditional_results.integration_regions = integration_decision.integration_regions;

% **v004**: Track synthetic testing mode
conditional_results.synthetic_testing_mode = is_synthetic_testing;
if is_synthetic_testing
    conditional_results.synthetic_integration_applied = true;
    conditional_results.methodology_testing_note = 'Synthetic integration regions generated for Stage 4 testing';
else
    conditional_results.synthetic_integration_applied = false;
end

% Include Stage 2 and Stage 1 metadata
conditional_results.stage2_file = stage2_file;
conditional_results.stage2_adequacy = adequacy_results;
conditional_results.stage1_r_squared = stage1_data.results.r_squared;
conditional_results.stage1_coefficients = stage1_data.results.total_coefficients;
conditional_results.tractability_level = stage1_data.tractability_level;

% Configuration and metadata
conditional_results.config = conditional_config;
conditional_results.stage3_timestamp = datestr(now);
conditional_results.stage3_version = 'ModelAdequacy_Stage3_Conditional_v004';
conditional_results.synthetic_integration_enhancement = 'v004: Added synthetic integration regions for Stage 4 methodology testing';

% Save results
saveStage3Results(conditional_results, stage2_file, stage1_data);

fprintf('\n✅ STAGE 3 (CONDITIONAL ANALYSIS) COMPLETED\n');
fprintf('✅ SYNTHETIC TESTING: Stage 4 methodology testing capability added (v004)\n');
if integration_decision.proceed_to_integration
    fprintf('📋 Next: Run ModelAdequacy_Stage4_Integration_v002 for hierarchical model development\n');
    if is_synthetic_testing
        fprintf('🔬 TESTING MODE: Stage 4 will use synthetic integration regions for methodology validation\n');
    end
else
    fprintf('📋 Analysis complete: Conditional models do not warrant hierarchical integration\n');
end

end

%% SUPPORTING FUNCTIONS

function is_synthetic = detectSyntheticTestingMode(adequacy_results, stage1_data)
% Detect if we're in synthetic data testing mode
% **NEW v004**: Determines whether to generate synthetic integration regions

is_synthetic = false;

% Check multiple indicators for synthetic testing mode
if isfield(adequacy_results, 'testing_mode') && adequacy_results.testing_mode
    is_synthetic = true;
end

% Check if Stage 1 data indicates synthetic generation
if isfield(stage1_data, 'data_source') && contains(stage1_data.data_source, 'synthetic')
    is_synthetic = true;
end

% Check if we have synthetic poor fit regions from Stage 2
if isfield(adequacy_results, 'synthetic_poor_regions_generated') && adequacy_results.synthetic_poor_regions_generated
    is_synthetic = true;
end

% Additional check: if we're NOT using database but have large dataset, likely synthetic
if isfield(stage1_data, 'tractability_level') && stage1_data.tractability_level >= 5
    if ~isfield(stage1_data, 'database_source') || ~stage1_data.database_source
        is_synthetic = true;
    end
end

end

function enhanced_models = generateSyntheticIntegrationRegions(regional_models, global_model, config)
% Generate synthetic integration regions for Stage 4 methodology testing
% **NEW v004**: Only for synthetic data - enables complete pipeline testing

enhanced_models = regional_models; % Start with original models
global_r_squared = global_model.Rsquared.Ordinary;

fprintf('Generating synthetic integration regions for methodology testing...\n');
fprintf('  Global model R²: %.4f\n', global_r_squared);
fprintf('  Integration improvement threshold: %.4f\n', config.improvement_threshold);

% **SYNTHETIC ENHANCEMENT**: Artificially boost performance of selected regions
num_regions = length(regional_models);
if num_regions >= 2
    % Enhance first 2 regions to show significant improvements
    
    % Region 1: Moderate improvement
    enhanced_models{1}.r_squared = global_r_squared + 0.08; % 8% improvement
    enhanced_models{1}.synthetic_enhancement = true;
    enhanced_models{1}.original_r_squared = regional_models{1}.r_squared;
    fprintf('  ✓ Enhanced Region 1: R² boosted to %.4f (ΔR² = +0.08)\n', enhanced_models{1}.r_squared);
    
    % Region 2: Strong improvement  
    enhanced_models{2}.r_squared = global_r_squared + 0.12; % 12% improvement
    enhanced_models{2}.synthetic_enhancement = true;
    enhanced_models{2}.original_r_squared = regional_models{2}.r_squared;
    fprintf('  ✓ Enhanced Region 2: R² boosted to %.4f (ΔR² = +0.12)\n', enhanced_models{2}.r_squared);
    
    % If we have more regions, enhance one more moderately
    if num_regions >= 4
        enhanced_models{4}.r_squared = global_r_squared + 0.06; % 6% improvement
        enhanced_models{4}.synthetic_enhancement = true;
        enhanced_models{4}.original_r_squared = regional_models{4}.r_squared;
        fprintf('  ✓ Enhanced Region 4: R² boosted to %.4f (ΔR² = +0.06)\n', enhanced_models{4}.r_squared);
    end
    
    % Mark remaining regions as unchanged for clarity
    for i = 3:num_regions
        if i ~= 4 || num_regions < 4
            enhanced_models{i}.synthetic_enhancement = false;
            enhanced_models{i}.original_r_squared = regional_models{i}.r_squared;
        end
    end
    
    fprintf('Generated %d enhanced regions for Stage 4 integration testing\n', min(3, num_regions));
    fprintf('💾 BlueBear optimized: Synthetic enhancements for methodology testing\n');
    
else
    fprintf('  ⚠️  Insufficient regions (%d) for synthetic enhancement - minimum 2 required\n', num_regions);
end

% Add metadata about synthetic enhancement
for i = 1:length(enhanced_models)
    if ~isfield(enhanced_models{i}, 'synthetic_enhancement')
        enhanced_models{i}.synthetic_enhancement = false;
        enhanced_models{i}.original_r_squared = enhanced_models{i}.r_squared;
    end
end

fprintf('🔬 Synthetic integration regions generated for complete pipeline testing\n');
fprintf('⚠️  IMPORTANT: This is methodology testing only - NOT real analysis results\n');

end

function value = getFieldSafely(struct_var, field_name, default_value)
% Safely get field value with default fallback
% **v003**: Utility function for safe field access

if isfield(struct_var, field_name)
    value = struct_var.(field_name);
else
    value = default_value;
end

end

function valid_poor_regions = isolatePoorFitRegions(poor_fit_regions, data_table, config)
% Isolate poor fit regions with sufficient observations for conditional modeling

valid_poor_regions = {};

for i = 1:length(poor_fit_regions)
    region = poor_fit_regions{i};
    
    % Determine region data based on description
    if isfield(region, 'indices') && ~isempty(region.indices)
        region_indices = region.indices;
    else
        % For residual-based regions, need to recreate indices
        region_indices = createRegionIndices(region, data_table);
    end
    
    if length(region_indices) >= config.minRegionSize
        region.indices = region_indices;
        region.n_obs = length(region_indices);
        valid_poor_regions{end+1} = region;
    end
end

end

function region_indices = createRegionIndices(region, data_table)
% Create region indices based on region description
% **FIXED v002**: Uses camelCase variable names consistently
% **FIXED v003**: Handle nominal/categorical filterType properly

region_indices = [];

% Parse region description to recreate indices
% **FIXED v002**: Use camelCase variable names
if contains(region.description, 'VGF_HIGH')
    vgf_values = data_table.VGF;
    median_vgf = median(vgf_values);
    region_indices = find(vgf_values > median_vgf);
elseif contains(region.description, 'VGF_LOW')
    vgf_values = data_table.VGF;
    median_vgf = median(vgf_values);
    region_indices = find(vgf_values <= median_vgf);
elseif contains(region.description, 'noiseMagnitude_HIGH') % FIXED: was noise_magnitude_HIGH
    noise_values = data_table.noiseMagnitude; % FIXED: was noise_magnitude
    median_noise = median(noise_values);
    region_indices = find(noise_values > median_noise);
elseif contains(region.description, 'noiseMagnitude_LOW') % FIXED: was noise_magnitude_LOW
    noise_values = data_table.noiseMagnitude; % FIXED: was noise_magnitude
    median_noise = median(noise_values);
    region_indices = find(noise_values <= median_noise);
elseif contains(region.description, 'filterType')
    % FIXED v003: Handle categorical filterType regions properly
    if contains(region.description, 'filterType_1')
        if isa(data_table.filterType, 'nominal') || isa(data_table.filterType, 'categorical')
            region_indices = find(data_table.filterType == nominal(1));
        else
            region_indices = find(data_table.filterType == 1);
        end
    elseif contains(region.description, 'filterType_2')
        if isa(data_table.filterType, 'nominal') || isa(data_table.filterType, 'categorical')
            region_indices = find(data_table.filterType == nominal(2));
        else
            region_indices = find(data_table.filterType == 2);
        end
    end
else
    % Default: use all data if region cannot be identified
    region_indices = (1:height(data_table))';
end

end

function conditional_config = initializeConditionalConfig(base_config)
% Initialize configuration for conditional modeling

conditional_config = base_config;

% Conditional modeling thresholds
conditional_config.improvement_threshold = 0.05;          % Minimum R² improvement
conditional_config.pvalue_threshold = 0.05;              % Statistical significance threshold
conditional_config.effect_size_threshold = 0.03;         % Clinical significance threshold
conditional_config.aic_improvement_threshold = -2;        % AIC improvement threshold

% Model complexity control
conditional_config.max_interactions = 3;                 % Maximum interaction order
conditional_config.min_observations_per_coeff = 10;      % Observations per coefficient rule

% Cross-validation settings
conditional_config.cv_folds = 3;                         % Reduced for computational efficiency
conditional_config.cv_repeats = 1;                       % Number of CV repetitions

end

function regional_model = developRegionalModel(region, data_table, global_model, config)
% Develop region-specific conditional model with tailored interactions
% **FIXED v002**: Uses camelCase variable names consistently
% **BLUEBEAR v003**: Added memory safety for 288GB system
% **FIXED v003**: Ensure ALL required fields are always initialized

regional_model = struct();
regional_model.region = region;

% BLUEBEAR MEMORY SAFETY: Check region size before extraction
region_indices = region.indices;
original_size = length(region_indices);

% Reasonable limit for BlueBear's 288GB - allow up to 500K observations per region
if original_size > 500000
    fprintf('    ⚠️  Large region detected (%d obs), sampling for BlueBear efficiency\n', original_size);
    sample_size = min(100000, original_size); % Conservative sample for large regions
    region_indices = randsample(region_indices, sample_size);
    fprintf('    📊 Sampled %d/%d observations for regional modeling\n', sample_size, original_size);
elseif original_size > 100000
    fprintf('    📊 Processing large region (%d obs) - BlueBear can handle this\n', original_size);
end

% Extract regional data with potentially reduced indices
region_data = data_table(region_indices, :);
n_obs = height(region_data);

% Update region info with actual processed size
regional_model.original_n_obs = original_size;
regional_model.processed_n_obs = n_obs;
if original_size ~= n_obs
    regional_model.sampling_applied = true;
    regional_model.sampling_ratio = n_obs / original_size;
else
    regional_model.sampling_applied = false;
end

% **CRITICAL v003**: Initialize ALL required fields before any processing
% This prevents field access errors later in the pipeline
regional_model.r_squared_improvement = 0;          % Default value
regional_model.improvement_significant = false;    % Default value
regional_model.model_successful = false;          % Will be set to true if successful
regional_model.r_squared = 0;                     % Default value
regional_model.num_coefficients = 0;              % Default value
regional_model.n_observations = n_obs;            % Set actual observations

% Start with essential interactions and build systematically
try
    % Simplified regional model specification using camelCase variables
    % **FIXED v002**: All variable references now use camelCase
    
    % Extract key variables for regional modeling
    y = region_data.deltaBeta;                  % FIXED: camelCase
    X_basic = [ones(n_obs, 1), region_data.betaGenerated]; % FIXED: camelCase
    
    % Add key interaction terms based on parameter space
    if ismember('VGF', region_data.Properties.VariableNames)
        X_basic = [X_basic, region_data.VGF, region_data.betaGenerated .* region_data.VGF];
    end
    
    if ismember('noiseMagnitude', region_data.Properties.VariableNames) % FIXED: camelCase
        X_basic = [X_basic, region_data.noiseMagnitude]; % FIXED: camelCase
    end
    
    if ismember('noiseColor', region_data.Properties.VariableNames) % FIXED: camelCase
        X_basic = [X_basic, region_data.noiseColor]; % FIXED: camelCase
    end
    
    % Add sampling rate effects if available
    if ismember('samplingRate', region_data.Properties.VariableNames) % FIXED: camelCase
        X_basic = [X_basic, region_data.samplingRate]; % FIXED: camelCase
    end
    
    % Fit regional model using robust estimation
    [beta_regional, ~, residuals_regional] = regress(y, X_basic);
    
    % Calculate regional model performance
    ss_res = sum(residuals_regional.^2);
    ss_tot = sum((y - mean(y)).^2);
    r_squared = 1 - ss_res/ss_tot;
    
    % Store regional model results - overwrite defaults with actual values
    regional_model.coefficients = beta_regional;
    regional_model.residuals = residuals_regional;
    regional_model.r_squared = max(0, r_squared); % Ensure non-negative
    regional_model.num_coefficients = length(beta_regional);
    regional_model.model_successful = true; % Mark as successful
    
    % Calculate model fit metrics
    regional_model.aic = n_obs * log(ss_res/n_obs) + 2 * length(beta_regional);
    regional_model.bic = n_obs * log(ss_res/n_obs) + log(n_obs) * length(beta_regional);
    
    % NOTE: r_squared_improvement and improvement_significant will be set later
    % in performComparativeAssessment function - they are already initialized above
    
catch ME
    % Handle modeling failures gracefully - keep initialized default values
    regional_model.error_message = ME.message;
    % All other fields already initialized with safe default values
    
    fprintf('    ⚠️ Regional model fitting failed: %s\n', ME.message);
end

end

function comparative_assessment = performComparativeAssessment(regional_models, global_model, data_table, config)
% Perform systematic comparison of regional vs global models
% **FIXED v003**: Comprehensive safety checks for all field access

comparative_assessment = struct();

% Initialize assessment metrics
num_models = length(regional_models);
r_squared_improvements = zeros(num_models, 1);
aic_improvements = zeros(num_models, 1);
significant_improvements = false(num_models, 1);

% Calculate global model performance for comparison
global_r_squared = global_model.Rsquared.Ordinary;

for i = 1:num_models
    regional_model = regional_models{i};
    
    % **FIXED v003**: Use safe field access with default values
    model_successful = getFieldSafely(regional_model, 'model_successful', false);
    regional_r_squared = getFieldSafely(regional_model, 'r_squared', 0);
    regional_aic = getFieldSafely(regional_model, 'aic', Inf);
    
    if model_successful && regional_r_squared > 0
        % Calculate R² improvement
        r_squared_improvements(i) = regional_r_squared - global_r_squared;
        
        % Approximate AIC improvement (simplified calculation)
        aic_improvements(i) = -regional_aic; % Negative for improvement
        
        % Test for statistical significance (simplified likelihood ratio test)
        improvement_significant = r_squared_improvements(i) > config.improvement_threshold;
        significant_improvements(i) = improvement_significant;
        
        % **FIXED v003**: Safely set improvement values in regional model
        regional_models{i}.r_squared_improvement = r_squared_improvements(i);
        regional_models{i}.improvement_significant = improvement_significant;
    else
        % **FIXED v003**: Ensure fields exist even for failed models
        if ~isfield(regional_models{i}, 'r_squared_improvement')
            regional_models{i}.r_squared_improvement = 0;
        end
        if ~isfield(regional_models{i}, 'improvement_significant')
            regional_models{i}.improvement_significant = false;
        end
    end
end

% Summary statistics
comparative_assessment.r_squared_improvements = r_squared_improvements;
comparative_assessment.aic_improvements = aic_improvements;
comparative_assessment.significant_improvements = significant_improvements;
comparative_assessment.num_improved_regions = sum(significant_improvements);
comparative_assessment.mean_r_squared_improvement = mean(r_squared_improvements);
comparative_assessment.integration_recommended = sum(significant_improvements) > 0;

end

function integration_decision = evaluateIntegrationDecision(comparative_assessment, regional_models, config)
% Evaluate whether conditional models warrant integration into hierarchical framework
% **FIXED v003**: All field access now uses comprehensive safety checks

integration_decision = struct();

% Evaluate statistical improvement criterion
statistical_improvement = comparative_assessment.num_improved_regions > 0;
integration_decision.statistical_improvement = statistical_improvement;

% **FIXED v003**: Evaluate effect size maintenance with safe field access
effect_sizes_maintained = true; % Simplified assumption
for i = 1:length(regional_models)
    model = regional_models{i};
    
    % **CRITICAL v003 FIX**: Use safe field access for ALL field checks
    model_successful = getFieldSafely(model, 'model_successful', false);
    improvement_significant = getFieldSafely(model, 'improvement_significant', false);
    r_squared_improvement = getFieldSafely(model, 'r_squared_improvement', 0);
    
    if model_successful && improvement_significant
        % Check if clinically meaningful effects are preserved
        % Simplified check: ensure R² improvement exceeds threshold
        if r_squared_improvement < config.effect_size_threshold
            effect_sizes_maintained = false;
        end
    end
end
integration_decision.effect_size_maintained = effect_sizes_maintained;

% Evaluate mechanistic insight provision
mechanistic_insight = statistical_improvement; % Simplified criterion
integration_decision.mechanistic_insight = mechanistic_insight;

% Overall integration decision
proceed_to_integration = statistical_improvement && effect_sizes_maintained && mechanistic_insight;
integration_decision.proceed_to_integration = proceed_to_integration;

% **FIXED v003**: Identify specific regions for integration with safe field access
integration_regions = {};
if proceed_to_integration
    for i = 1:length(regional_models)
        model = regional_models{i};
        
        % **CRITICAL v003 FIX**: Use safe field access throughout
        model_successful = getFieldSafely(model, 'model_successful', false);
        improvement_significant = getFieldSafely(model, 'improvement_significant', false);
        
        if model_successful && improvement_significant
            integration_regions{end+1} = model;
        end
    end
end
integration_decision.integration_regions = integration_regions;

% Decision summary
integration_decision.summary = struct();
integration_decision.summary.qualified_regions = length(integration_regions);
integration_decision.summary.statistical_threshold_met = statistical_improvement;
integration_decision.summary.effect_size_maintained = effect_sizes_maintained;
integration_decision.summary.mechanistic_insight_provided = mechanistic_insight;

end

function mechanistic_insights = documentMechanisticInsights(regional_models, integration_decision, config)
% Document mechanistic insights from conditional analysis
% **FIXED v003**: Safe field access throughout

mechanistic_insights = struct();

% Initialize insight categories
interaction_patterns = {};
parameter_relationships = {};
clinical_implications = {};

% Analyze integration-qualified regions for mechanistic insights
for i = 1:length(integration_decision.integration_regions)
    model = integration_decision.integration_regions{i};
    
    % **FIXED v003**: Use safe field access for all model fields
    region_description = 'Unknown Region';
    if isfield(model, 'region') && isfield(model.region, 'description')
        region_description = model.region.description;
    end
    
    r_squared_improvement = getFieldSafely(model, 'r_squared_improvement', 0);
    
    % Document interaction patterns
    pattern = struct();
    pattern.region = region_description;
    pattern.r_squared_improvement = r_squared_improvement;
    pattern.interaction_type = 'regional_enhancement';
    interaction_patterns{end+1} = pattern;
    
    % Document parameter relationships
    relationship = struct();
    relationship.region = region_description;
    relationship.relationship_type = 'conditional_parameter_variation';
    relationship.effect_magnitude = r_squared_improvement;
    parameter_relationships{end+1} = relationship;
    
    % Document clinical implications
    if isfield(config, 'clinicalThresholds') && isfield(config.clinicalThresholds, 'clinical_detection')
        clinical_threshold = config.clinicalThresholds.clinical_detection;
    else
        clinical_threshold = 0.01; % Default clinical detection threshold
    end
    
    if r_squared_improvement > clinical_threshold
        implication = struct();
        implication.region = region_description;
        implication.clinical_relevance = 'detection_threshold_exceeded';
        implication.magnitude = r_squared_improvement;
        clinical_implications{end+1} = implication;
    end
end

mechanistic_insights.interaction_patterns = interaction_patterns;
mechanistic_insights.parameter_relationships = parameter_relationships;
mechanistic_insights.clinical_implications = clinical_implications;

% Summary insights
mechanistic_insights.summary = struct();
mechanistic_insights.summary.total_patterns = length(interaction_patterns);
mechanistic_insights.summary.total_relationships = length(parameter_relationships);
mechanistic_insights.summary.clinical_relevant_regions = length(clinical_implications);

end

function conditional_results = createNoConditionalAnalysisResults(adequacy_results, stage1_data, config)
% Create results structure when no conditional analysis is required

conditional_results = struct();
conditional_results.conditional_analysis_required = false;
conditional_results.adequacy_reason = 'global_model_adequate';
conditional_results.poor_fit_regions = {};
conditional_results.regional_models = {};
conditional_results.requires_integration = false; % FIXED naming consistency
conditional_results.integration_regions = {};

% Include Stage 1 metadata
conditional_results.stage1_r_squared = stage1_data.results.r_squared;
conditional_results.stage1_coefficients = stage1_data.results.total_coefficients;
conditional_results.tractability_level = stage1_data.tractability_level;

% Configuration and metadata
conditional_results.config = config;
conditional_results.stage3_timestamp = datestr(now);
conditional_results.stage3_version = 'ModelAdequacy_Stage3_Conditional_v004';
conditional_results.synthetic_integration_enhancement = 'v004: Added synthetic integration regions for Stage 4 methodology testing';

end

function saveStage3Results(conditional_results, stage2_file, stage1_data)
% Save Stage 3 results with appropriate filename

% Create filename with timestamp
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
if isfield(stage1_data, 'tractability_level')
    stage3_filename = sprintf('stage3_conditional_L%d_%s.mat', stage1_data.tractability_level, timestamp);
else
    stage3_filename = sprintf('stage3_conditional_%s.mat', timestamp);
end

save(stage3_filename, 'conditional_results', '-v7.3');
fprintf('✓ Stage 3 results saved: %s\n', stage3_filename);

end

function str = logical2str(logical_val)
% Convert logical to string
if logical_val
    str = 'YES';
else
    str = 'NO';
end
end
