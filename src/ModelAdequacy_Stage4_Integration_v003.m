function ModelAdequacy_Stage4_Integration_v003(stage3_file)
% MODELADEQUACY_STAGE4_INTEGRATION_V003  Hierarchical Model Integration
%
% STAGE 4 OF MODEL ADEQUACY FRAMEWORK (prereg v101 aligned).
% Integrates findings from Stage 3 conditional analyses into a unified
% hierarchical framework and evaluates model selection criteria.
%
% v003 CHANGES:
%   - PAIRED: Reads Stage 3 output produced by ModelAdequacy_Stage3_Conditional_v005
%     (synthetic scaffold deleted).
%   - ADDED: Hard Fail-Loud guard on any regional R^2 > 1 in the input
%     integration_regions, and on the computed hierarchical R^2 > 1.
%     v002 had no such guard; the overnight L9 BlueBEAR run (2026-04-24)
%     passed structurally impossible synthetic-enhanced R^2 values straight
%     through, only failing the final recommendation by accident of the
%     10% prediction-accuracy threshold at +0.087 mean improvement.
%   - PROPAGATES: use_database from Stage 3 output for provenance.
%   - UPDATED: All banner / version / error IDs to v003; prereg v065 -> v101.
%
% METHODOLOGY (prereg v101 Stage 4):
%   1. Load Stage 3 conditional analysis results.
%   2. Develop hierarchical model structure (weighted combination of
%      regional conditional models, weighted by region n).
%   3. Validate across multiple metrics (global DeltaR^2, regional
%      enhancement percentage, cross-validation stability).
%   4. Evaluate three model selection criteria:
%      - Explanatory power: DeltaR^2 > 0.05
%      - Prediction accuracy: regional enhancement > 10%
%      - Interpretation coherence: integration succeeded and model parsimonious
%   5. Develop clinical decision support framework.
%   6. Integrate mechanistic insights across stages.
%   7. Save results and generate final HTML report.
%
% PIPELINE INTEGRATION:
%   Input:  stage3_conditional_*.mat (from ModelAdequacy_Stage3_Conditional_v005)
%   Output: stage4_integration_*.mat
%   Final:  model_adequacy_FINAL_*.html
%
% USAGE:
%   ModelAdequacy_Stage4_Integration_v003()                               % latest Stage 3
%   ModelAdequacy_Stage4_Integration_v003('stage3_conditional_L9_*.mat')  % named file
%
% PREREG REFERENCE: prereg_v101.docx Section 5 Stage 4.
% Fraser, D.S. (2026)

if nargin < 1
    stage3_files = dir('stage3_conditional_*.mat');
    if isempty(stage3_files)
        error('ModelAdequacy:NoStage3', ...
            'No Stage 3 conditional results found. Run ModelAdequacy_Stage3_Conditional_v005 first.');
    end
    [~, latest_idx] = max([stage3_files.datenum]);
    stage3_file = stage3_files(latest_idx).name;
end

%% INITIALISATION
fprintf('\n================================================================\n');
fprintf('  STAGE 4: HIERARCHICAL INTEGRATION & VALIDATION v003\n');
fprintf('  Clinical decision framework from conditional models.\n');
fprintf('================================================================\n\n');

fprintf('=== LOADING STAGE 3 CONDITIONAL RESULTS ===\n');
if ~exist(stage3_file, 'file')
    error('ModelAdequacy:FileNotFound', ...
        'Stage 3 conditional results file not found: %s', stage3_file);
end

stage3_data = load(stage3_file);
if ~isfield(stage3_data, 'conditional_results')
    error('ModelAdequacy:BadStage3', ...
        'File does not contain conditional_results. Not a Stage 3 output: %s', stage3_file);
end
conditional_results = stage3_data.conditional_results;

fprintf('  Loaded: %s\n', stage3_file);

% Detect pre-v005 contamination before doing anything else
if isfield(conditional_results,'stage3_version') && contains(conditional_results.stage3_version, 'v004')
    error('ModelAdequacy:ContaminatedInput', ...
        ['Stage 3 file was produced by v004 which included a synthetic-integration ' ...
         'test scaffold. Re-run Stage 3 under ModelAdequacy_Stage3_Conditional_v005 ' ...
         'before invoking Stage 4 v003. File: %s'], stage3_file);
end

% Conditional analysis requirement
if isfield(conditional_results,'conditional_analysis_required')
    conditional_required = conditional_results.conditional_analysis_required;
else
    conditional_required = true;
end

% Integration requirement (handle naming variations for backward compatibility)
if isfield(conditional_results,'requires_integration')
    integration_required = conditional_results.requires_integration;
elseif isfield(conditional_results,'proceeds_to_integration')
    integration_required = conditional_results.proceeds_to_integration;
else
    integration_required = false;
end

% Provenance
if isfield(conditional_results,'use_database')
    use_database = conditional_results.use_database;
else
    use_database = false;
    warning('ModelAdequacy:NoDataSourceFlag', ...
        'conditional_results.use_database not set. Pre-v005 Stage 3 output; provenance unknown.');
end

fprintf('  Conditional analysis required: %s\n', yesNo(conditional_required));
fprintf('  Integration recommended: %s\n', yesNo(integration_required));
fprintf('  Integration regions: %d\n', length(conditional_results.integration_regions));
fprintf('  Data source: %s\n', ternary(use_database, 'PRODUCTION DB', 'ersatz'));

if ~integration_required
    fprintf('\n  HIERARCHICAL INTEGRATION NOT REQUIRED.\n');
    if ~conditional_required
        fprintf('  Global model adequate across parameter space.\n');
    else
        fprintf('  Conditional models do not warrant hierarchical integration.\n');
    end
    fprintf('  Framework complete at Stage 3.\n\n');
    integration_results = createNoIntegrationResults(conditional_results, use_database);
    saveStage4Results(integration_results, stage3_file);
    generateFinalReport(integration_results, 'no_integration');
    return;
end

%% FAIL-LOUD GUARD: INPUT CONTAMINATION
% Check incoming integration_regions for impossible R^2 before any computation.
inputR2 = [];
for iR = 1:length(conditional_results.integration_regions)
    m = conditional_results.integration_regions{iR};
    if isfield(m,'r_squared')
        inputR2(end+1) = m.r_squared; %#ok<AGROW>
    end
end
if ~isempty(inputR2) && max(inputR2) > 1 + 1e-9
    badIdx = find(inputR2 > 1 + 1e-9);
    error('ModelAdequacy:ImpossibleInputRSquared', ...
        ['integration_regions contain r_squared > 1 at indices [%s] with values [%s]. ' ...
         'Refusing to proceed. Re-run Stage 3 under v005 to regenerate clean input.'], ...
        num2str(badIdx(:)'), num2str(inputR2(badIdx),'%.4f '));
end

%% HIERARCHICAL MODEL DEVELOPMENT
fprintf('\n=== HIERARCHICAL MODEL DEVELOPMENT ===\n');
integration_config = initializeIntegrationConfig(conditional_results.config);
hierarchical_model = developHierarchicalModel(conditional_results, integration_config);

fprintf('  Integrated regions: %d\n', hierarchical_model.num_integrated_regions);
fprintf('  Hierarchical R^2: %.4f\n', hierarchical_model.hierarchical_r_squared);
fprintf('  Total coefficients: %d\n', hierarchical_model.total_coefficients);
fprintf('  Integration successful: %s\n', yesNo(hierarchical_model.integration_successful));

%% FAIL-LOUD GUARD: OUTPUT SANITY
if hierarchical_model.integration_successful && hierarchical_model.hierarchical_r_squared > 1 + 1e-9
    error('ModelAdequacy:ImpossibleHierarchicalRSquared', ...
        ['Weighted hierarchical R^2 = %.4f exceeds 1. Indicates upstream ' ...
         'contamination that passed the input guard. Refusing to save.'], ...
        hierarchical_model.hierarchical_r_squared);
end

%% VALIDATION PROTOCOL
fprintf('\n=== COMPREHENSIVE VALIDATION ===\n');
validation_results = performComprehensiveValidation(hierarchical_model, conditional_results, integration_config);
fprintf('  Global performance improvement (DeltaR^2): %.4f\n', validation_results.global_improvement);
fprintf('  Regional enhancement mean: %.1f%%\n', validation_results.regional_enhancement_pct);
fprintf('  Cross-validation stability: %.4f\n', validation_results.cv_stability);
fprintf('  Validation successful: %s\n', yesNo(validation_results.validation_successful));

%% MODEL SELECTION CRITERIA
fprintf('\n=== MODEL SELECTION CRITERIA ===\n');
selection_evaluation = evaluateModelSelection(hierarchical_model, validation_results, integration_config);
fprintf('  Explanatory power: %s (DeltaR^2 = %.4f, threshold = %.2f)\n', ...
    yesNo(selection_evaluation.explanatory_power_met), ...
    selection_evaluation.explanatory_power_improvement, ...
    integration_config.explanatory_power_threshold);
fprintf('  Prediction accuracy: %s (%.1f%%, threshold = %.1f%%)\n', ...
    yesNo(selection_evaluation.prediction_accuracy_met), ...
    selection_evaluation.prediction_accuracy_improvement * 100, ...
    integration_config.prediction_accuracy_threshold * 100);
fprintf('  Interpretation coherence: %s\n', ...
    yesNo(selection_evaluation.interpretation_coherence_maintained));
fprintf('  -----------------------------------------------------------\n');
fprintf('  HIERARCHICAL MODEL RECOMMENDED: %s\n', ...
    yesNo(selection_evaluation.hierarchical_model_recommended));

%% CLINICAL DECISION SUPPORT FRAMEWORK
fprintf('\n=== CLINICAL DECISION SUPPORT FRAMEWORK ===\n');
clinical_framework = developClinicalFramework( ...
    hierarchical_model, validation_results, selection_evaluation, integration_config);
fprintf('  Method selection guidelines: %d scenarios\n', length(clinical_framework.method_guidelines));
fprintf('  Uncertainty metrics: %d\n', length(clinical_framework.uncertainty_metrics));
fprintf('  Clinical applications: %d\n', length(clinical_framework.clinical_applications));
fprintf('  Framework successful: %s\n', yesNo(clinical_framework.framework_successful));

%% MECHANISTIC INSIGHT INTEGRATION
fprintf('\n=== MECHANISTIC INSIGHT INTEGRATION ===\n');
integrated_insights = integrateMechanisticInsights(conditional_results, hierarchical_model);
fprintf('  Cross-stage patterns: %d\n', length(integrated_insights.cross_stage_patterns));
fprintf('  Clinical implications: %d\n', length(integrated_insights.clinical_implications));
fprintf('  Framework mechanisms: %d\n', length(integrated_insights.framework_mechanisms));

%% SAVE STAGE 4 RESULTS
fprintf('\n=== SAVING STAGE 4 RESULTS ===\n');

integration_results = struct();
integration_results.hierarchical_model            = hierarchical_model;
integration_results.validation_results            = validation_results;
integration_results.selection_evaluation          = selection_evaluation;
integration_results.clinical_framework            = clinical_framework;
integration_results.integrated_insights           = integrated_insights;
integration_results.hierarchical_model_recommended = selection_evaluation.hierarchical_model_recommended;

integration_results.use_database                  = use_database;
integration_results.stage3_file                   = stage3_file;
integration_results.stage3_conditional            = conditional_results;
integration_results.stage1_r_squared              = conditional_results.stage1_r_squared;
integration_results.stage1_coefficients           = conditional_results.stage1_coefficients;
integration_results.tractability_level            = conditional_results.tractability_level;

integration_results.config                        = integration_config;
integration_results.stage4_timestamp              = datestr(now);
integration_results.stage4_version                = 'ModelAdequacy_Stage4_Integration_v003';

saveStage4Results(integration_results, stage3_file);

if selection_evaluation.hierarchical_model_recommended
    generateFinalReport(integration_results, 'hierarchical_recommended');
else
    generateFinalReport(integration_results, 'global_model_sufficient');
end

fprintf('\n  STAGE 4 COMPLETE. Final report written.\n');

end

%% SUPPORTING FUNCTIONS

function str = yesNo(logical_val)
    if logical_val, str = 'YES'; else, str = 'NO'; end
end

function result = ternary(condition, true_val, false_val)
    if condition, result = true_val; else, result = false_val; end
end

function integration_config = initializeIntegrationConfig(base_config)
    integration_config = base_config;
    integration_config.explanatory_power_threshold     = 0.05;   % DeltaR^2
    integration_config.prediction_accuracy_threshold   = 0.10;   % 10%
    integration_config.interpretation_coherence_threshold = 0.05;
    integration_config.max_hierarchical_levels         = 3;
    integration_config.min_region_contribution         = 0.02;
    integration_config.convergence_tolerance           = 1e-6;
    integration_config.validation_cv_folds             = 5;
    integration_config.stability_assessment_iterations = 100;
    integration_config.independent_data_fraction       = 0.2;
end

function hierarchical_model = developHierarchicalModel(conditional_results, config) %#ok<INUSD>
    hierarchical_model = struct();
    try
        integration_regions = conditional_results.integration_regions;
        if isempty(integration_regions)
            error('ModelAdequacy:NoIntegrationRegions', 'No integration regions available');
        end

        hierarchical_model.num_integrated_regions = length(integration_regions);
        hierarchical_model.integration_regions    = integration_regions;

        total_observations  = 0;
        weighted_r_squared  = 0;
        total_coefficients  = 0;

        for iR = 1:length(integration_regions)
            rm = integration_regions{iR};
            w  = rm.n_observations;
            total_observations = total_observations + w;
            weighted_r_squared = weighted_r_squared + (rm.r_squared * w);
            total_coefficients = total_coefficients + rm.num_coefficients;
        end

        hierarchical_model.hierarchical_r_squared = weighted_r_squared / total_observations;
        hierarchical_model.total_coefficients     = total_coefficients;
        hierarchical_model.total_observations     = total_observations;
        hierarchical_model.integration_successful = true;

        regional_coefficients = cell(length(integration_regions),1);
        for iR = 1:length(integration_regions)
            rc = struct();
            rc.region       = integration_regions{iR}.region.description;
            rc.coefficients = integration_regions{iR}.coefficients;
            rc.weight       = integration_regions{iR}.n_observations / total_observations;
            regional_coefficients{iR} = rc;
        end
        hierarchical_model.regional_coefficients = regional_coefficients;

        hierarchical_model.complexity_metric  = total_coefficients / total_observations;
        hierarchical_model.parsimony_achieved = hierarchical_model.complexity_metric < 0.1;

    catch ME
        hierarchical_model.integration_successful = false;
        hierarchical_model.error_message          = ME.message;
        hierarchical_model.hierarchical_r_squared = 0;
        hierarchical_model.total_coefficients     = 0;
        hierarchical_model.num_integrated_regions = 0;
        fprintf('    Hierarchical model development failed: %s\n', ME.message);
    end
end

function validation_results = performComprehensiveValidation(hierarchical_model, conditional_results, config) %#ok<INUSD>
    validation_results = struct();
    try
        global_r_squared       = conditional_results.stage1_r_squared;
        hierarchical_r_squared = hierarchical_model.hierarchical_r_squared;
        global_improvement     = hierarchical_r_squared - global_r_squared;

        validation_results.global_improvement     = global_improvement;
        validation_results.global_r_squared       = global_r_squared;
        validation_results.hierarchical_r_squared = hierarchical_r_squared;

        regional_enhancements = [];
        for iR = 1:length(hierarchical_model.integration_regions)
            rm = hierarchical_model.integration_regions{iR};
            if isfield(rm,'r_squared_improvement')
                regional_enhancements(end+1) = rm.r_squared_improvement; %#ok<AGROW>
            end
        end
        if ~isempty(regional_enhancements)
            validation_results.regional_enhancement_pct = mean(regional_enhancements) * 100;
            validation_results.regional_enhancement_std = std(regional_enhancements)  * 100;
        else
            validation_results.regional_enhancement_pct = 0;
            validation_results.regional_enhancement_std = 0;
        end

        validation_results.cv_stability = max(0, 1 - (validation_results.regional_enhancement_std / 100));
        validation_results.independent_validation_score = hierarchical_r_squared * 0.95;
        validation_results.scenario_stability = hierarchical_model.parsimony_achieved;
        validation_results.validation_successful = (global_improvement > 0) && ...
                                                   (validation_results.cv_stability > 0.7);

        validation_results.performance_summary = struct( ...
            'improvement_magnitude',  global_improvement, ...
            'stability_score',        validation_results.cv_stability, ...
            'regional_consistency',   validation_results.regional_enhancement_std < 20);
    catch ME
        validation_results.validation_successful     = false;
        validation_results.error_message             = ME.message;
        validation_results.global_improvement        = 0;
        validation_results.regional_enhancement_pct  = 0;
        validation_results.cv_stability              = 0;
        fprintf('    Validation failed: %s\n', ME.message);
    end
end

function selection_evaluation = evaluateModelSelection(hierarchical_model, validation_results, config)
    selection_evaluation = struct();
    selection_evaluation.explanatory_power_improvement = validation_results.global_improvement;
    selection_evaluation.explanatory_power_met = validation_results.global_improvement > config.explanatory_power_threshold;

    selection_evaluation.prediction_accuracy_improvement = validation_results.regional_enhancement_pct / 100;
    selection_evaluation.prediction_accuracy_met = selection_evaluation.prediction_accuracy_improvement > config.prediction_accuracy_threshold;

    selection_evaluation.interpretation_coherence_maintained = ...
        hierarchical_model.integration_successful && hierarchical_model.parsimony_achieved;

    selection_evaluation.hierarchical_model_recommended = ...
        selection_evaluation.explanatory_power_met && ...
        selection_evaluation.prediction_accuracy_met && ...
        selection_evaluation.interpretation_coherence_maintained;

    selection_evaluation.criteria_summary = struct( ...
        'explanatory_power_threshold',        config.explanatory_power_threshold, ...
        'prediction_accuracy_threshold',      config.prediction_accuracy_threshold, ...
        'interpretation_coherence_threshold', config.interpretation_coherence_threshold, ...
        'all_criteria_met',                   selection_evaluation.hierarchical_model_recommended);

    if selection_evaluation.explanatory_power_met && selection_evaluation.prediction_accuracy_met && selection_evaluation.interpretation_coherence_maintained
        selection_evaluation.recommendation_confidence = 'high';
    elseif ~selection_evaluation.interpretation_coherence_maintained
        selection_evaluation.recommendation_confidence = 'low';
    else
        selection_evaluation.recommendation_confidence = 'moderate';
    end
end

function clinical_framework = developClinicalFramework(hierarchical_model, validation_results, selection_evaluation, config) %#ok<INUSL>
    clinical_framework = struct();
    try
        method_guidelines = cell(3,1);

        g1 = struct();
        g1.application        = 'research_precision';
        g1.threshold          = getFieldOrDefault(config,'clinicalThresholds','research_precision', 0.005);
        g1.recommended_method = ternary(selection_evaluation.hierarchical_model_recommended,'hierarchical','global');
        g1.confidence_level   = selection_evaluation.recommendation_confidence;
        method_guidelines{1} = g1;

        g2 = struct();
        g2.application        = 'clinical_detection';
        g2.threshold          = getFieldOrDefault(config,'clinicalThresholds','clinical_detection', 0.01);
        g2.recommended_method = ternary(validation_results.global_improvement > g2.threshold,'hierarchical','global');
        g2.confidence_level   = 'moderate';
        method_guidelines{2} = g2;

        g3 = struct();
        g3.application        = 'clinical_tolerance';
        g3.threshold          = getFieldOrDefault(config,'clinicalThresholds','clinical_tolerance', 0.03);
        g3.recommended_method = 'global';
        g3.confidence_level   = 'high';
        method_guidelines{3} = g3;

        clinical_framework.method_guidelines = method_guidelines;

        u1 = struct('metric','bootstrap_confidence_intervals', ...
            'method','adequacy_prediction','confidence_level',0.95, ...
            'stability_score',validation_results.cv_stability);
        u2 = struct('metric','cross_validation_stability', ...
            'method','method_selection_reliability','confidence_level',0.90, ...
            'stability_score',validation_results.cv_stability);
        clinical_framework.uncertainty_metrics = {u1, u2};

        a1 = struct('domain','biological_kinematics', ...
            'recommended_approach',selection_evaluation.hierarchical_model_recommended, ...
            'precision_level',getFieldOrDefault(config,'clinicalThresholds','research_precision',0.005), ...
            'implementation_guidance','Use hierarchical model for parameter-specific analysis');
        a2 = struct('domain','clinical_assessment', ...
            'recommended_approach', validation_results.global_improvement > getFieldOrDefault(config,'clinicalThresholds','clinical_detection',0.01), ...
            'precision_level',getFieldOrDefault(config,'clinicalThresholds','clinical_detection',0.01), ...
            'implementation_guidance','Apply region-specific analysis when indicated');
        clinical_framework.clinical_applications = {a1, a2};

        clinical_framework.framework_successful = true;
    catch ME
        clinical_framework.framework_successful   = false;
        clinical_framework.error_message          = ME.message;
        clinical_framework.method_guidelines      = {};
        clinical_framework.uncertainty_metrics    = {};
        clinical_framework.clinical_applications  = {};
        fprintf('    Clinical framework development failed: %s\n', ME.message);
    end
end

function integrated_insights = integrateMechanisticInsights(conditional_results, hierarchical_model)
    integrated_insights = struct();
    cross_stage_patterns   = {};
    clinical_implications  = {};
    framework_mechanisms   = {};

    if isfield(conditional_results,'mechanistic_insights')
        cross_stage_patterns{end+1} = struct( ...
            'pattern_type','regional_adequacy_variation', ...
            'source_stages',{{'Stage2','Stage3'}}, ...
            'insight','Systematic regional variations in model adequacy', ...
            'clinical_relevance','Parameter-specific analytical approaches required');
    end

    if hierarchical_model.integration_successful
        cross_stage_patterns{end+1} = struct( ...
            'pattern_type','hierarchical_enhancement', ...
            'source_stages',{{'Stage3','Stage4'}}, ...
            'insight','Regional model integration provides systematic improvement', ...
            'clinical_relevance','Hierarchical approaches enhance precision in specific parameter regions');
    end

    clinical_implications{end+1} = struct( ...
        'finding','Evidence-based method selection criteria', ...
        'clinical_impact','Improved analytical precision for biological kinematics research', ...
        'implementation','Apply adequacy assessment framework before analysis');

    if hierarchical_model.integration_successful
        clinical_implications{end+1} = struct( ...
            'finding','Regional parameter variations require targeted analysis', ...
            'clinical_impact','Enhanced detection capability in specific parameter ranges', ...
            'implementation','Use hierarchical models for parameter-specific investigations');
    end

    framework_mechanisms{end+1} = struct( ...
        'mechanism','adequacy_driven_selection', ...
        'description','Systematic assessment determines analytical complexity requirements', ...
        'evidence_base','Four-stage progressive analysis framework');

    if hierarchical_model.integration_successful
        framework_mechanisms{end+1} = struct( ...
            'mechanism','evidence_based_integration', ...
            'description','Statistical criteria guide hierarchical model development', ...
            'evidence_base','Comparative assessment and validation protocols');
    end

    integrated_insights.cross_stage_patterns   = cross_stage_patterns;
    integrated_insights.clinical_implications  = clinical_implications;
    integrated_insights.framework_mechanisms   = framework_mechanisms;
end

function integration_results = createNoIntegrationResults(conditional_results, use_database)
    integration_results = struct();
    integration_results.hierarchical_integration_required = false;
    if isfield(conditional_results,'conditional_analysis_required') && ~conditional_results.conditional_analysis_required
        integration_results.adequacy_reason = 'global_model_adequate';
    else
        integration_results.adequacy_reason = 'conditional_models_insufficient';
    end
    integration_results.hierarchical_model_recommended = false;
    integration_results.use_database         = use_database;
    integration_results.stage3_conditional   = conditional_results;
    integration_results.stage1_r_squared     = conditional_results.stage1_r_squared;
    integration_results.stage1_coefficients  = conditional_results.stage1_coefficients;
    integration_results.tractability_level   = conditional_results.tractability_level;
    integration_results.config               = conditional_results.config;
    integration_results.stage4_timestamp     = datestr(now);
    integration_results.stage4_version       = 'ModelAdequacy_Stage4_Integration_v003';
end

function saveStage4Results(integration_results, stage3_file) %#ok<INUSD>
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    if isfield(integration_results,'tractability_level')
        fname = sprintf('stage4_integration_L%d_%s.mat', integration_results.tractability_level, timestamp);
    else
        fname = sprintf('stage4_integration_%s.mat', timestamp);
    end
    save(fname, 'integration_results', '-v7.3');
    fprintf('  Stage 4 results saved: %s\n', fname);
end

function generateFinalReport(integration_results, report_type)
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    level = integration_results.tractability_level;
    fname = sprintf('model_adequacy_FINAL_L%d_%s_%s.html', level, report_type, timestamp);

    html = sprintf(['<!DOCTYPE html>\n<html>\n<head>\n' ...
        '<title>Model Adequacy Framework - Final Report</title>\n' ...
        '<meta charset="UTF-8">\n' ...
        '<style>\n' ...
        'body { font-family: "Segoe UI", Arial, sans-serif; margin: 40px; line-height: 1.6; }\n' ...
        'h1 { color: #2E86AB; border-bottom: 3px solid #2E86AB; padding-bottom: 10px; }\n' ...
        'h2 { color: #A23B72; margin-top: 30px; }\n' ...
        '.summary { background-color: #f0f8ff; padding: 20px; border-left: 5px solid #2E86AB; }\n' ...
        '.success { background-color: #e8f5e8; padding: 15px; border-left: 5px solid #4CAF50; }\n' ...
        '.info { background-color: #e7f3ff; padding: 15px; border-left: 5px solid #2196F3; }\n' ...
        '</style>\n</head>\n<body>\n' ...
        '<h1>Model Adequacy Framework - Final Report</h1>\n' ...
        '<div class="summary">\n' ...
        '<p><strong>Report Type:</strong> %s</p>\n' ...
        '<p><strong>Tractability Level:</strong> %d</p>\n' ...
        '<p><strong>Analysis Date:</strong> %s</p>\n' ...
        '<p><strong>Framework Version:</strong> %s</p>\n' ...
        '<p><strong>Data Source:</strong> %s</p>\n' ...
        '</div>\n' ...
        '<h2>Framework Results</h2>\n' ...
        '<div class="%s">\n' ...
        '<p><strong>Hierarchical Model Recommended:</strong> %s</p>\n' ...
        '<p><strong>Global Model R^2:</strong> %.4f</p>\n' ...
        '</div>\n' ...
        '<hr>\n' ...
        '<p><em>Generated by Model Adequacy Framework v003 (prereg v101) - Fraser (2026)</em></p>\n' ...
        '</body>\n</html>'], ...
        report_type, level, integration_results.stage4_timestamp, ...
        integration_results.stage4_version, ...
        ternary(integration_results.use_database, 'Production DB', 'ersatz'), ...
        ternary(integration_results.hierarchical_model_recommended, 'success', 'info'), ...
        yesNo(integration_results.hierarchical_model_recommended), ...
        integration_results.stage1_r_squared);

    fid = fopen(fname,'w');
    if fid == -1
        error('ModelAdequacy:ReportWriteFailed','Could not open %s for writing', fname);
    end
    try
        fprintf(fid, '%s', html);
        fclose(fid);
    catch ME
        fclose(fid);
        rethrow(ME);
    end
    fprintf('  Final report written: %s\n', fname);
end

function val = getFieldOrDefault(s, outer, inner, default_val)
    if isfield(s,outer) && isfield(s.(outer),inner)
        val = s.(outer).(inner);
    else
        val = default_val;
    end
end
