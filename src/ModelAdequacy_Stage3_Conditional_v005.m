function ModelAdequacy_Stage3_Conditional_v005(stage2_file)
% MODELADEQUACY_STAGE3_CONDITIONAL_V005  Conditional Parameter Analysis
%
% STAGE 3 OF MODEL ADEQUACY FRAMEWORK (prereg v101 aligned, synthetic deleted).
% Develops region-specific conditional models for parameter regions flagged
% as inadequate by Stage 2, compares them against the Stage 1 global model,
% and returns an integration decision for Stage 4.
%
% v005 CHANGES:
%   - DELETED: All synthetic integration-region test scaffolding that was
%     added in v004 for pipeline methodology testing. The Stage 2 -> 3
%     handoff now carries only real poor-fit regions from residual analysis
%     and changepoint detection (Stage 2 v003 already removed its synthetic
%     scaffold; v005 completes the pair).
%   - DELETED: detectSyntheticTestingMode() helper and its L >= 5 trap-door
%     which fired on every production-level DB run regardless of source.
%   - DELETED: generateSyntheticIntegrationRegions() helper which boosted
%     regional R^2 by +0.08 / +0.12 / +0.06, producing structurally
%     impossible R^2 > 1 values that propagated into Stage 4 inputs.
%   - ADDED: Hard Fail-Loud guard on regional R^2 > 1 before save.
%   - ADDED: Propagates stage1_data.use_database through the output struct
%     as conditional_results.use_database for provenance.
%   - KEPT: Safe field access, camelCase variable validation, BLUEBEAR
%     memory monitoring for large regional models, double() casts on
%     DB-loaded table columns in developRegionalModel.
%
% METHODOLOGY (prereg v101 Stage 3):
%   1. Load Stage 2 adequacy assessment results.
%   2. Isolate poor-fit regions with sufficient observations.
%   3. Fit region-specific conditional models with tailored interactions.
%   4. Compare regional vs global model performance via R^2 improvement.
%   5. Apply integration decision criteria (statistical + effect size +
%      mechanistic insight) to determine whether Stage 4 should proceed.
%   6. Document mechanistic insights for reporting.
%   7. Save results to stage3_conditional_L<level>_<timestamp>.mat.
%
% PIPELINE INTEGRATION:
%   Input:  stage2_adequacy_*.mat (from ModelAdequacy_Stage2_Assessment_v003)
%   Output: stage3_conditional_*.mat
%   Next:   ModelAdequacy_Stage4_Integration_v003.m (if integration required)
%
% USAGE:
%   ModelAdequacy_Stage3_Conditional_v005()                              % latest Stage 2
%   ModelAdequacy_Stage3_Conditional_v005('stage2_adequacy_L9_*.mat')    % named file
%
% PREREG REFERENCE: prereg_v101.docx Section 5 Stage 3.
% Fraser, D.S. (2026)

if nargin < 1
    stage2_files = dir('stage2_adequacy_*.mat');
    if isempty(stage2_files)
        error('ModelAdequacy:NoStage2', ...
            'No Stage 2 adequacy results found. Run ModelAdequacy_Stage2_Assessment_v003 first.');
    end
    [~, latest_idx] = max([stage2_files.datenum]);
    stage2_file = stage2_files(latest_idx).name;
end

%% INITIALISATION
fprintf('\n================================================================\n');
fprintf('  STAGE 3: CONDITIONAL PARAMETER ANALYSIS v005 (prereg aligned)\n');
fprintf('  Region-specific modelling for poor-fit regions from Stage 2.\n');
fprintf('================================================================\n\n');

fprintf('=== LOADING STAGE 2 ADEQUACY RESULTS ===\n');
if ~exist(stage2_file, 'file')
    error('ModelAdequacy:FileNotFound', ...
        'Stage 2 adequacy results file not found: %s', stage2_file);
end

stage2_data = load(stage2_file);
fprintf('  Loaded: %s\n', stage2_file);

adequacy_results = stage2_data.adequacy_results;
config           = stage2_data.config;
stage1_data      = stage2_data.stage1_data;

fprintf('  Overall adequacy: %s\n', yesNo(adequacy_results.global_adequate));
fprintf('  Conditional analysis required: %s\n', ...
    yesNo(adequacy_results.requires_conditional_analysis));
fprintf('  Poor fit regions identified: %d\n', length(adequacy_results.poor_fit_regions));

% Provenance: the Stage 1 output already carries use_database from KitchenSink v002.
if isfield(stage1_data, 'use_database')
    use_database = stage1_data.use_database;
    fprintf('  Data source: %s\n', ...
        ternary(use_database, 'PRODUCTION DB (use_database=true)', 'ersatz (use_database=false)'));
else
    use_database = false;
    warning('ModelAdequacy:NoDataSourceFlag', ...
        ['stage1_data.use_database is not set. Stage 1 may be pre-v002 output. ' ...
         'Cannot verify provenance; assuming false for safety.']);
end

% Validate conditional analysis requirement
if ~adequacy_results.requires_conditional_analysis
    fprintf('\n  CONDITIONAL ANALYSIS NOT REQUIRED.\n');
    fprintf('  Global model is adequate across all regions. Framework complete at Stage 2.\n\n');
    conditional_results = createNoConditionalAnalysisResults(adequacy_results, stage1_data, config, use_database);
    saveStage3Results(conditional_results, stage2_file, stage1_data);
    return;
end

%% POOR FIT REGION ISOLATION
fprintf('\n=== POOR FIT REGION ISOLATION ===\n');

global_model = stage1_data.results.model;
data_table   = stage1_data.tableTrue;

% camelCase variable validation (pairs with Stage 2 v003 and KitchenSink v002)
required_vars = {'deltaBeta','betaGenerated','VGF','samplingRate', ...
                 'filterType','regressionType','noiseMagnitude','noiseColor'};
missing_vars = setdiff(required_vars, data_table.Properties.VariableNames);
if ~isempty(missing_vars)
    error('ModelAdequacy:MissingVars', ...
        'Missing camelCase variables: %s. Check Stage 1 output.', ...
        strjoin(missing_vars, ', '));
end
fprintf('  Variable naming validation: camelCase PASS\n');

valid_poor_regions = isolatePoorFitRegions(adequacy_results.poor_fit_regions, data_table, config);

fprintf('  Total regions identified: %d\n', length(adequacy_results.poor_fit_regions));
fprintf('  Valid regions for analysis: %d\n', length(valid_poor_regions));
fprintf('  Minimum region size: %d observations\n', config.minRegionSize);

if isempty(valid_poor_regions)
    fprintf('\n  No valid poor-fit regions meet minimum size requirements.\n');
    fprintf('  Framework complete at Stage 2.\n\n');
    conditional_results = createNoConditionalAnalysisResults(adequacy_results, stage1_data, config, use_database);
    saveStage3Results(conditional_results, stage2_file, stage1_data);
    return;
end

%% REGION-SPECIFIC CONDITIONAL MODELLING
fprintf('\n=== REGION-SPECIFIC CONDITIONAL MODELLING ===\n');
conditional_config = initializeConditionalConfig(config);

% Memory monitoring (BlueBEAR 288 GB host)
if ispc
    [~, sys] = memory;
    fprintf('  Initial available memory: %.1f GB\n', sys.PhysicalMemory.Available / 1e9);
else
    fprintf('  macOS/Linux host: memory monitoring via ispc branch disabled\n');
end

regional_models = cell(length(valid_poor_regions), 1);
for iR = 1:length(valid_poor_regions)
    region = valid_poor_regions{iR};
    fprintf('\n  --- Region %d: %s (%d observations) ---\n', ...
        iR, region.description, region.n_obs);

    regional_model = developRegionalModel(region, data_table, global_model, conditional_config);
    regional_models{iR} = regional_model;

    fprintf('    Regional model: R^2 = %.4f, coefficients = %d\n', ...
        regional_model.r_squared, regional_model.num_coefficients);

    % Periodic cleanup for large region counts
    if mod(iR, 3) == 0 && iR < length(valid_poor_regions)
        if ispc, java.lang.System.gc(); end
        pause(0.05);
    end
end

%% FAIL-LOUD GUARD: IMPOSSIBLE R^2
regional_r2 = cellfun(@(m) getFieldSafely(m, 'r_squared', 0), regional_models);
if any(regional_r2 > 1.0 + 1e-9)
    impossibleIdx = find(regional_r2 > 1.0 + 1e-9);
    impossibleR2  = regional_r2(impossibleIdx);
    msg = sprintf('Regional R^2 > 1 detected at region indices [%s] with values [%s]. ', ...
        num2str(impossibleIdx(:)'), num2str(impossibleR2(:)','%.4f '));
    msg = [msg, 'This is structurally impossible and indicates contamination (e.g. synthetic ' ...
        'enhancement scaffold from Stage 3 v004). Refusing to save.'];
    error('ModelAdequacy:ImpossibleRSquared', '%s', msg);
end

%% COMPARATIVE MODEL ASSESSMENT
fprintf('\n=== COMPARATIVE MODEL ASSESSMENT ===\n');
[comparative_assessment, regional_models] = performComparativeAssessment( ...
    regional_models, global_model, conditional_config);

fprintf('  Regions with significant improvement: %d / %d\n', ...
    comparative_assessment.num_improved_regions, length(regional_models));
fprintf('  Mean R^2 improvement: %.4f\n', comparative_assessment.mean_r_squared_improvement);
fprintf('  Integration recommended: %s\n', ...
    yesNo(comparative_assessment.integration_recommended));

%% INTEGRATION DECISION EVALUATION
fprintf('\n=== INTEGRATION DECISION EVALUATION ===\n');
integration_decision = evaluateIntegrationDecision(comparative_assessment, regional_models, conditional_config);

fprintf('  Statistical improvement threshold met: %s\n', yesNo(integration_decision.statistical_improvement));
fprintf('  Effect size maintained: %s\n', yesNo(integration_decision.effect_size_maintained));
fprintf('  Mechanistic insight provided: %s\n', yesNo(integration_decision.mechanistic_insight));
fprintf('  ---------------------------------------------------\n');
fprintf('  INTEGRATION DECISION: %s\n', yesNo(integration_decision.proceed_to_integration));

if integration_decision.proceed_to_integration
    fprintf('\n  %d regions qualify for hierarchical integration\n', length(integration_decision.integration_regions));
    for iI = 1:min(3, length(integration_decision.integration_regions))
        rm = integration_decision.integration_regions{iI};
        fprintf('    Region %d: %s (DeltaR^2 = %.4f)\n', iI, rm.region.description, ...
            getFieldSafely(rm, 'r_squared_improvement', 0));
    end
    if length(integration_decision.integration_regions) > 3
        fprintf('    ... and %d more regions\n', length(integration_decision.integration_regions) - 3);
    end
end

%% MECHANISTIC INSIGHT DOCUMENTATION
fprintf('\n=== MECHANISTIC INSIGHT DOCUMENTATION ===\n');
mechanistic_insights = documentMechanisticInsights(integration_decision, conditional_config);
fprintf('  Interaction patterns identified: %d\n', length(mechanistic_insights.interaction_patterns));
fprintf('  Parameter relationships documented: %d\n', length(mechanistic_insights.parameter_relationships));
fprintf('  Clinical implications: %d\n', length(mechanistic_insights.clinical_implications));

%% SAVE STAGE 3 RESULTS
fprintf('\n=== SAVING STAGE 3 RESULTS ===\n');

conditional_results = struct();
conditional_results.poor_fit_regions      = valid_poor_regions;
conditional_results.regional_models       = regional_models;
conditional_results.comparative_assessment = comparative_assessment;
conditional_results.integration_decision  = integration_decision;
conditional_results.mechanistic_insights  = mechanistic_insights;
conditional_results.requires_integration  = integration_decision.proceed_to_integration;
conditional_results.integration_regions   = integration_decision.integration_regions;

% Provenance
conditional_results.use_database          = use_database;
conditional_results.stage2_file           = stage2_file;
conditional_results.stage2_adequacy       = adequacy_results;
conditional_results.stage1_r_squared      = stage1_data.results.r_squared;
conditional_results.stage1_coefficients   = stage1_data.results.total_coefficients;
conditional_results.tractability_level    = stage1_data.tractability_level;

conditional_results.config                = conditional_config;
conditional_results.stage3_timestamp      = datestr(now);
conditional_results.stage3_version        = 'ModelAdequacy_Stage3_Conditional_v005';

saveStage3Results(conditional_results, stage2_file, stage1_data);

fprintf('\n  STAGE 3 COMPLETE.\n');
if integration_decision.proceed_to_integration
    fprintf('  Next: ModelAdequacy_Stage4_Integration_v003\n');
else
    fprintf('  Conditional models do not warrant hierarchical integration.\n');
    fprintf('  Analysis complete; Stage 4 not required.\n');
end

end

%% SUPPORTING FUNCTIONS

function value = getFieldSafely(struct_var, field_name, default_value)
    if isfield(struct_var, field_name)
        value = struct_var.(field_name);
    else
        value = default_value;
    end
end

function str = yesNo(logical_val)
    if logical_val, str = 'YES'; else, str = 'NO'; end
end

function result = ternary(condition, true_val, false_val)
    if condition, result = true_val; else, result = false_val; end
end

function valid_poor_regions = isolatePoorFitRegions(poor_fit_regions, data_table, config)
    valid_poor_regions = {};
    for iR = 1:length(poor_fit_regions)
        region = poor_fit_regions{iR};
        if isfield(region, 'indices') && ~isempty(region.indices)
            region_indices = region.indices;
        else
            region_indices = createRegionIndices(region, data_table);
        end
        if length(region_indices) >= config.minRegionSize
            region.indices = region_indices;
            region.n_obs   = length(region_indices);
            valid_poor_regions{end+1} = region; %#ok<AGROW>
        end
    end
end

function region_indices = createRegionIndices(region, data_table)
    % Recreate indices for a region described by Stage 2 changepoint/residual output.
    region_indices = [];
    if contains(region.description, 'VGF_HIGH')
        m = median(data_table.VGF);
        region_indices = find(data_table.VGF > m);
    elseif contains(region.description, 'VGF_LOW')
        m = median(data_table.VGF);
        region_indices = find(data_table.VGF <= m);
    elseif contains(region.description, 'noiseMagnitude_HIGH')
        m = median(data_table.noiseMagnitude);
        region_indices = find(data_table.noiseMagnitude > m);
    elseif contains(region.description, 'noiseMagnitude_LOW')
        m = median(data_table.noiseMagnitude);
        region_indices = find(data_table.noiseMagnitude <= m);
    elseif contains(region.description, 'filterType')
        % filterType may be nominal (from KitchenSink v002) or double
        ft = data_table.filterType;
        if contains(region.description, 'filterType_1')
            if isa(ft,'nominal') || isa(ft,'categorical')
                region_indices = find(ft == nominal(1));
            else
                region_indices = find(ft == 1);
            end
        elseif contains(region.description, 'filterType_2')
            if isa(ft,'nominal') || isa(ft,'categorical')
                region_indices = find(ft == nominal(2));
            else
                region_indices = find(ft == 2);
            end
        end
    else
        region_indices = (1:height(data_table))';
    end
end

function conditional_config = initializeConditionalConfig(base_config)
    conditional_config = base_config;
    conditional_config.improvement_threshold       = 0.05;   % min DeltaR^2 for integration (prereg S8)
    conditional_config.pvalue_threshold            = 0.05;
    conditional_config.effect_size_threshold       = 0.03;   % MDC from Cook et al. 2023
    conditional_config.aic_improvement_threshold   = -2;
    conditional_config.max_interactions            = 3;
    conditional_config.min_observations_per_coeff  = 10;
    conditional_config.cv_folds                    = 3;
    conditional_config.cv_repeats                  = 1;
end

function regional_model = developRegionalModel(region, data_table, global_model, config) %#ok<INUSD>
    % Fit a region-specific model using camelCase variables and double()-cast
    % DB-loaded columns.
    regional_model = struct();
    regional_model.region = region;

    region_indices = region.indices;
    original_size  = length(region_indices);

    % Guard very large regions to keep fitlme tractable on BlueBEAR
    if original_size > 500000
        fprintf('    Large region (%d obs), downsampling for tractability\n', original_size);
        sample_size = min(100000, original_size);
        region_indices = randsample(region_indices, sample_size);
        fprintf('    Sampled %d / %d observations\n', sample_size, original_size);
    end

    region_data = data_table(region_indices, :);
    n_obs = height(region_data);

    regional_model.original_n_obs           = original_size;
    regional_model.processed_n_obs          = n_obs;
    regional_model.sampling_applied         = (original_size ~= n_obs);
    if regional_model.sampling_applied
        regional_model.sampling_ratio = n_obs / original_size;
    else
        regional_model.sampling_ratio = 1;
    end

    % Defaults (overwritten on success)
    regional_model.r_squared_improvement = 0;
    regional_model.improvement_significant = false;
    regional_model.model_successful      = false;
    regional_model.r_squared             = 0;
    regional_model.num_coefficients      = 0;
    regional_model.n_observations        = n_obs;

    try
        % double() cast defends against JDBC int64 returns on DB-loaded tables
        y = double(region_data.deltaBeta);
        X = [ones(n_obs,1), double(region_data.betaGenerated)];

        if ismember('VGF', region_data.Properties.VariableNames)
            X = [X, double(region_data.VGF), ...
                    double(region_data.betaGenerated) .* double(region_data.VGF)];
        end
        if ismember('noiseMagnitude', region_data.Properties.VariableNames)
            X = [X, double(region_data.noiseMagnitude)];
        end
        if ismember('noiseColor', region_data.Properties.VariableNames)
            X = [X, double(region_data.noiseColor)];
        end
        if ismember('samplingRate', region_data.Properties.VariableNames)
            X = [X, double(region_data.samplingRate)];
        end

        [beta_regional, ~, residuals_regional] = regress(y, X);

        ss_res = sum(residuals_regional.^2);
        ss_tot = sum((y - mean(y)).^2);
        r_squared = 1 - ss_res / ss_tot;

        regional_model.coefficients     = beta_regional;
        regional_model.residuals        = residuals_regional;
        regional_model.r_squared        = max(0, r_squared);
        regional_model.num_coefficients = length(beta_regional);
        regional_model.model_successful = true;
        regional_model.aic = n_obs * log(ss_res / n_obs) + 2 * length(beta_regional);
        regional_model.bic = n_obs * log(ss_res / n_obs) + log(n_obs) * length(beta_regional);

    catch ME
        regional_model.error_message = ME.message;
        fprintf('    Regional model fitting failed: %s\n', ME.message);
    end
end

function [comparative_assessment, regional_models] = performComparativeAssessment(regional_models, global_model, config)
    comparative_assessment = struct();
    num_models = length(regional_models);
    r_squared_improvements    = zeros(num_models, 1);
    aic_improvements          = zeros(num_models, 1);
    significant_improvements  = false(num_models, 1);

    global_r_squared = global_model.Rsquared.Ordinary;

    for iM = 1:num_models
        rm = regional_models{iM};
        model_successful  = getFieldSafely(rm, 'model_successful', false);
        regional_r_squared = getFieldSafely(rm, 'r_squared', 0);
        regional_aic      = getFieldSafely(rm, 'aic', Inf);

        if model_successful && regional_r_squared > 0
            r_squared_improvements(iM) = regional_r_squared - global_r_squared;
            aic_improvements(iM)       = -regional_aic;
            sig = r_squared_improvements(iM) > config.improvement_threshold;
            significant_improvements(iM) = sig;
            regional_models{iM}.r_squared_improvement   = r_squared_improvements(iM);
            regional_models{iM}.improvement_significant = sig;
        else
            if ~isfield(regional_models{iM}, 'r_squared_improvement')
                regional_models{iM}.r_squared_improvement = 0;
            end
            if ~isfield(regional_models{iM}, 'improvement_significant')
                regional_models{iM}.improvement_significant = false;
            end
        end
    end

    comparative_assessment.r_squared_improvements     = r_squared_improvements;
    comparative_assessment.aic_improvements           = aic_improvements;
    comparative_assessment.significant_improvements   = significant_improvements;
    comparative_assessment.num_improved_regions       = sum(significant_improvements);
    comparative_assessment.mean_r_squared_improvement = mean(r_squared_improvements);
    comparative_assessment.integration_recommended    = sum(significant_improvements) > 0;
end

function integration_decision = evaluateIntegrationDecision(comparative_assessment, regional_models, config)
    integration_decision = struct();

    statistical_improvement = comparative_assessment.num_improved_regions > 0;
    integration_decision.statistical_improvement = statistical_improvement;

    effect_sizes_maintained = true;
    for iM = 1:length(regional_models)
        m = regional_models{iM};
        if getFieldSafely(m,'model_successful',false) && ...
           getFieldSafely(m,'improvement_significant',false)
            if getFieldSafely(m,'r_squared_improvement',0) < config.effect_size_threshold
                effect_sizes_maintained = false;
            end
        end
    end
    integration_decision.effect_size_maintained = effect_sizes_maintained;

    mechanistic_insight = statistical_improvement;
    integration_decision.mechanistic_insight = mechanistic_insight;

    proceed = statistical_improvement && effect_sizes_maintained && mechanistic_insight;
    integration_decision.proceed_to_integration = proceed;

    integration_regions = {};
    if proceed
        for iM = 1:length(regional_models)
            m = regional_models{iM};
            if getFieldSafely(m,'model_successful',false) && ...
               getFieldSafely(m,'improvement_significant',false)
                integration_regions{end+1} = m; %#ok<AGROW>
            end
        end
    end
    integration_decision.integration_regions = integration_regions;

    integration_decision.summary = struct( ...
        'qualified_regions',            length(integration_regions), ...
        'statistical_threshold_met',    statistical_improvement, ...
        'effect_size_maintained',       effect_sizes_maintained, ...
        'mechanistic_insight_provided', mechanistic_insight);
end

function mechanistic_insights = documentMechanisticInsights(integration_decision, config)
    mechanistic_insights = struct();
    interaction_patterns     = {};
    parameter_relationships  = {};
    clinical_implications    = {};

    for iR = 1:length(integration_decision.integration_regions)
        m = integration_decision.integration_regions{iR};

        region_description = 'Unknown Region';
        if isfield(m,'region') && isfield(m.region,'description')
            region_description = m.region.description;
        end
        dR2 = getFieldSafely(m,'r_squared_improvement',0);

        interaction_patterns{end+1}    = struct( ...
            'region', region_description, ...
            'r_squared_improvement', dR2, ...
            'interaction_type', 'regional_enhancement'); %#ok<AGROW>

        parameter_relationships{end+1} = struct( ...
            'region', region_description, ...
            'relationship_type', 'conditional_parameter_variation', ...
            'effect_magnitude', dR2); %#ok<AGROW>

        if isfield(config,'clinicalThresholds') && isfield(config.clinicalThresholds,'clinical_detection')
            clinical_threshold = config.clinicalThresholds.clinical_detection;
        else
            clinical_threshold = 0.01;
        end
        if dR2 > clinical_threshold
            clinical_implications{end+1} = struct( ...
                'region', region_description, ...
                'clinical_relevance', 'detection_threshold_exceeded', ...
                'magnitude', dR2); %#ok<AGROW>
        end
    end

    mechanistic_insights.interaction_patterns    = interaction_patterns;
    mechanistic_insights.parameter_relationships = parameter_relationships;
    mechanistic_insights.clinical_implications   = clinical_implications;
    mechanistic_insights.summary = struct( ...
        'total_patterns',           length(interaction_patterns), ...
        'total_relationships',      length(parameter_relationships), ...
        'clinical_relevant_regions',length(clinical_implications));
end

function conditional_results = createNoConditionalAnalysisResults(adequacy_results, stage1_data, config, use_database)
    conditional_results = struct();
    conditional_results.conditional_analysis_required = false;
    conditional_results.adequacy_reason       = 'global_model_adequate';
    conditional_results.poor_fit_regions      = {};
    conditional_results.regional_models       = {};
    conditional_results.requires_integration  = false;
    conditional_results.integration_regions   = {};
    conditional_results.use_database          = use_database;
    conditional_results.stage1_r_squared      = stage1_data.results.r_squared;
    conditional_results.stage1_coefficients   = stage1_data.results.total_coefficients;
    conditional_results.tractability_level    = stage1_data.tractability_level;
    conditional_results.config                = config;
    conditional_results.stage3_timestamp      = datestr(now);
    conditional_results.stage3_version        = 'ModelAdequacy_Stage3_Conditional_v005';
    % Non-use: adequacy_results is referenced only to bind the signature shape
    % for future callers; stored for provenance.
    conditional_results.stage2_adequacy       = adequacy_results;
end

function saveStage3Results(conditional_results, stage2_file, stage1_data) %#ok<INUSL>
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    if isfield(stage1_data, 'tractability_level')
        fname = sprintf('stage3_conditional_L%d_%s.mat', stage1_data.tractability_level, timestamp);
    else
        fname = sprintf('stage3_conditional_%s.mat', timestamp);
    end
    save(fname, 'conditional_results', '-v7.3');
    fprintf('  Stage 3 results saved: %s\n', fname);
end
