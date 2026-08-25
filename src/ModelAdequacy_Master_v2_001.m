function ModelAdequacy_Master_v2_001(tractability_level, n_obs_per_combo, use_database, start_from_stage)
% MODELADEQUACY_MASTER_V004  Master pipeline for Model Adequacy Framework.
%
% 4-Stage progressive analysis pipeline (prereg v101 Section 5) orchestrating
% Kitchen Sink global LMM through hierarchical integration and clinical
% decision support framework generation.
%
% v004 CHANGES:
%   - CALLS: ModelAdequacy_Stage3_Conditional_v005 (synthetic scaffold deleted)
%   - CALLS: ModelAdequacy_Stage4_Integration_v003 (R^2 > 1 guard, v005 input)
%   - PROPAGATES: use_database through Stage 3 / Stage 4 outputs for provenance
%   - DROPPED: "Undefined function" placeholder branches for Stages 3 and 4
%     (they exist and are called unconditionally when required)
%   - DROPPED: memory_bypass fallback that created placeholder Stage 3 .mat
%     files when fitlme ran out of memory - this was for ersatz-L9 on Windows
%     and is obsolete now that L9 runs cleanly on BlueBEAR
%   - FIXED: use_regime_based_generation at line ~242 was hardcoded true;
%     now set to ~use_database so DB runs pass false (regime-based generation
%     is an ersatz-only concept)
%   - UPDATED: All banner / error ID / version strings from v002/v065 to
%     v004 / v101
%
% PIPELINE:
%   Stage 1: Kitchen Sink Model (ModelAdequacy_Stage1_KitchenSink_v2_001)
%            -> stage1_results_*.mat (global 7-way factorial LMM)
%   Stage 2: Model Adequacy Assessment (ModelAdequacy_Stage2_Assessment_v003)
%            -> stage2_adequacy_*.mat (3 prereg criteria)
%   Stage 3: Conditional Parameter Analysis (ModelAdequacy_Stage3_Conditional_v005)
%            -> stage3_conditional_*.mat (only if Stage 2 flags poor regions)
%   Stage 4: Hierarchical Integration (ModelAdequacy_Stage4_Integration_v003)
%            -> stage4_integration_*.mat (only if Stage 3 flags integration)
%   Report:  model_adequacy_FINAL_*.html
%
% TRACTABILITY LEVELS (KitchenSink v2_001 loads from DB when use_database=true):
%   1 Conservative  (28.8K)   2 Focused      (472K)    3 Minimal       (2.7K)
%   4 Moderate      (952K)    5 Substantial  (1.55M)   6 Extensive     (1.91M)
%   7 Comprehensive (4.37M)   8 Near-Orig    (13.8M)   9 Full-Original (14.8M)
%
% USAGE:
%   ModelAdequacy_Master_v2_001()                       % defaults: L2, n=5, ersatz
%   ModelAdequacy_Master_v2_001(2, 5, false)            % L2 on ersatz data
%   ModelAdequacy_Master_v2_001(2, 5, true)             % L2 on production DB
%   ModelAdequacy_Master_v2_001(9, 5, true)             % L9 on production DB (BlueBEAR)
%   ModelAdequacy_Master_v2_001(9, 5, true, 2)          % resume from Stage 2
%
% PREREG REFERENCE: prereg_v101.docx Section 5.
% Fraser, D.S. (2026)

%% INPUT VALIDATION
if nargin < 1, tractability_level = 2; end
if nargin < 2, n_obs_per_combo    = 5; end
if nargin < 3, use_database       = false; end
if nargin < 4, start_from_stage   = 1; end

if ~isnumeric(tractability_level) || tractability_level < 1 || tractability_level > 9 || mod(tractability_level,1) ~= 0
    error('ModelAdequacy_Master_v2_001:InvalidTractabilityLevel', ...
        'Tractability level must be an integer 1-9 (got %g).', tractability_level);
end
if ~isnumeric(n_obs_per_combo) || n_obs_per_combo < 1 || n_obs_per_combo > 50 || mod(n_obs_per_combo,1) ~= 0
    error('ModelAdequacy_Master_v2_001:InvalidObservationCount', ...
        'Observations per combo must be an integer 1-50 (got %g).', n_obs_per_combo);
end
if ~isnumeric(start_from_stage) || start_from_stage < 1 || start_from_stage > 4 || mod(start_from_stage,1) ~= 0
    error('ModelAdequacy_Master_v2_001:InvalidStartingStage', ...
        'Starting stage must be an integer 1-4 (got %g).', start_from_stage);
end

% Path anchoring
masterDir_MA = fileparts(mfilename('fullpath'));
if exist(fullfile(masterDir_MA, 'functions'), 'dir')
    addpath(genpath(fullfile(masterDir_MA, 'functions')));
end
if exist(fullfile(masterDir_MA, 'utils'), 'dir')
    addpath(genpath(fullfile(masterDir_MA, 'utils')));
end

%% TRACTABILITY METADATA
tract_names = {'Conservative','Focused','Minimal','Moderate','Substantial', ...
               'Extensive','Comprehensive','Near-Original','Full-Original'};
obs_k       = [28.8, 472, 2.7, 952, 1550, 1910, 4370, 13800, 14800];
time_min    = [2, 15, 1, 45, 90, 150, 300, 800, 1200];
mem_gb      = [2, 8, 1, 12, 20, 24, 40, 80, 100];
red_pct     = [99.0, 96.8, 99.9, 93.5, 89.5, 87.0, 70.2, 6.5, 0.0];

current_name       = tract_names{tractability_level};
current_obs_m      = obs_k(tractability_level) / 1000;
current_time       = time_min(tractability_level);
current_mem        = mem_gb(tractability_level);
current_reduction  = red_pct(tractability_level);

%% BANNER
fprintf('\n');
fprintf('================================================================\n');
fprintf('  MODEL ADEQUACY FRAMEWORK MASTER PIPELINE v004 (prereg v101)\n');
fprintf('  4-Stage Progressive Analysis with Systematic Assessment\n');
fprintf('  Fraser (2026) - frozen prereg v101, DB v057\n');
fprintf('================================================================\n\n');

fprintf('Configuration:\n');
fprintf('  Tractability: Level %d (%s - %.2fM obs, %.1f%% reduction)\n', ...
    tractability_level, current_name, current_obs_m, current_reduction);
fprintf('  Observations per combo: %d\n', n_obs_per_combo);
fprintf('  Data source: %s\n', ternary(use_database, ...
    'PRODUCTION DB (powerlaw_multiverse_v057.db)', 'ersatz (synthetic heuristic)'));

stage_descriptions = { ...
    '(full pipeline from Stage 1 Kitchen Sink)'; ...
    '(resume from Stage 2 Adequacy Assessment)'; ...
    '(resume from Stage 3 Conditional Analysis)'; ...
    '(resume from Stage 4 Integration)'};
fprintf('  Starting from: Stage %d %s\n', start_from_stage, stage_descriptions{start_from_stage});

fprintf('\nComputational estimate:\n');
fprintf('  Execution time: ~%.0f minutes\n', current_time);
fprintf('  Memory: ~%.0f GB RAM\n', current_mem);
fprintf('  Crash recovery: automatic via Stage 1 checkpoints\n\n');

%% STAGE 1: KITCHEN SINK MODEL
if start_from_stage <= 1
    fprintf('==== STAGE 1: KITCHEN SINK MODEL ====\n');
    fprintf('  Model: deltaBeta ~ betaGenerated * VGF * samplingRate * filterType\n');
    fprintf('                    * regressionType * noiseMagnitude * noiseColor\n');
    fprintf('                    + (1|paramComboID)\n\n');

    stage1_pattern = sprintf('stage1_results_L%d_*.mat', tractability_level);
    stage1_files = dir(stage1_pattern);

    if ~isempty(stage1_files) && start_from_stage > 1
        [~, i] = max([stage1_files.datenum]);
        stage1_file = stage1_files(i).name;
        fprintf('  Using existing Stage 1 file: %s\n', stage1_file);
    else
        fprintf('  Executing Stage 1...\n');
        % Regime-based generation is an ersatz-only concept (makes the synthetic
        % data harder for Stage 2 to pass trivially). DB runs pass false.
        use_regime = ~use_database;
        try
            ModelAdequacy_Stage1_KitchenSink_v2_001(use_database, n_obs_per_combo, tractability_level, use_regime);
        catch ME
            fprintf('  Stage 1 FAILED: %s\n', ME.message);
            rethrow(ME);
        end

        stage1_files = dir(stage1_pattern);
        if isempty(stage1_files)
            error('ModelAdequacy_Master_v2_001:Stage1NoOutput', ...
                'Stage 1 did not produce a results file matching %s', stage1_pattern);
        end
        [~, i] = max([stage1_files.datenum]);
        stage1_file = stage1_files(i).name;
    end
else
    stage1_pattern = sprintf('stage1_results_L%d_*.mat', tractability_level);
    stage1_files = dir(stage1_pattern);
    if isempty(stage1_files)
        error('ModelAdequacy_Master_v2_001:NoStage1', ...
            'Cannot start from Stage %d: no Stage 1 results for Level %d.', ...
            start_from_stage, tractability_level);
    end
    [~, i] = max([stage1_files.datenum]);
    stage1_file = stage1_files(i).name;
    fprintf('  Locating existing Stage 1: %s\n', stage1_file);
end

% Validate Stage 1 carries the use_database flag we need for provenance
stage1_check = load(stage1_file);
if ~isfield(stage1_check,'use_database')
    warning('ModelAdequacy_Master_v2_001:LegacyStage1', ...
        ['Stage 1 file %s does not contain use_database flag. ' ...
         'Likely a pre-v002 KitchenSink output. Provenance will be degraded.'], ...
        stage1_file);
else
    if stage1_check.use_database ~= use_database
        warning('ModelAdequacy_Master_v2_001:DataSourceMismatch', ...
            ['Stage 1 file was produced with use_database=%s but Master_v004 was ' ...
             'invoked with use_database=%s. Proceeding with Stage 1''s value.'], ...
            bool2str(stage1_check.use_database), bool2str(use_database));
    end
end

%% STAGE 2: ADEQUACY ASSESSMENT
if start_from_stage <= 2
    fprintf('\n==== STAGE 2: MODEL ADEQUACY ASSESSMENT ====\n');
    fprintf('  Criterion 1: Wald CI > +/-0.03 (coefficient stability)\n');
    fprintf('  Criterion 2: Interaction-level Cohen''s d > 0.5\n');
    fprintf('  Criterion 3: Changepoint analysis of regional R^2\n\n');

    stage2_pattern = sprintf('stage2_adequacy_L%d_*.mat', tractability_level);
    stage2_files = dir(stage2_pattern);

    if ~isempty(stage2_files) && start_from_stage > 2
        [~, i] = max([stage2_files.datenum]);
        stage2_file = stage2_files(i).name;
        fprintf('  Using existing Stage 2 file: %s\n', stage2_file);
    else
        fprintf('  Executing Stage 2...\n');
        try
            ModelAdequacy_Stage2_Assessment_v003(stage1_file);
        catch ME
            fprintf('  Stage 2 FAILED: %s\n', ME.message);
            rethrow(ME);
        end

        stage2_files = dir(stage2_pattern);
        if isempty(stage2_files)
            error('ModelAdequacy_Master_v2_001:Stage2NoOutput', ...
                'Stage 2 did not produce a results file matching %s', stage2_pattern);
        end
        [~, i] = max([stage2_files.datenum]);
        stage2_file = stage2_files(i).name;
    end
else
    stage2_pattern = sprintf('stage2_adequacy_L%d_*.mat', tractability_level);
    stage2_files = dir(stage2_pattern);
    if isempty(stage2_files)
        error('ModelAdequacy_Master_v2_001:NoStage2', ...
            'Cannot start from Stage %d: no Stage 2 results for Level %d.', ...
            start_from_stage, tractability_level);
    end
    [~, i] = max([stage2_files.datenum]);
    stage2_file = stage2_files(i).name;
    fprintf('  Locating existing Stage 2: %s\n', stage2_file);
end

% Decision point: does Stage 2 want conditional analysis?
stage2_data = load(stage2_file);
requires_conditional = stage2_data.adequacy_results.requires_conditional_analysis;

fprintf('\n  Stage 2 verdict: %s\n', ternary(requires_conditional, ...
    'CONDITIONAL ANALYSIS REQUIRED', 'global model ADEQUATE'));
if requires_conditional
    fprintf('  Poor-fit regions: %d\n', length(stage2_data.adequacy_results.poor_fit_regions));
end

%% STAGE 3: CONDITIONAL PARAMETER ANALYSIS
stage3_file = '';
requires_integration = false;

if requires_conditional && start_from_stage <= 3
    fprintf('\n==== STAGE 3: CONDITIONAL PARAMETER ANALYSIS ====\n');
    fprintf('  Poor-fit regions to model: %d\n\n', length(stage2_data.adequacy_results.poor_fit_regions));

    stage3_pattern = sprintf('stage3_conditional_L%d_*.mat', tractability_level);
    stage3_files = dir(stage3_pattern);

    if ~isempty(stage3_files) && start_from_stage > 3
        [~, i] = max([stage3_files.datenum]);
        stage3_file = stage3_files(i).name;
        fprintf('  Using existing Stage 3 file: %s\n', stage3_file);
    else
        fprintf('  Executing Stage 3...\n');
        try
            ModelAdequacy_Stage3_Conditional_v005(stage2_file);
        catch ME
            fprintf('  Stage 3 FAILED: %s\n', ME.message);
            rethrow(ME);
        end

        stage3_files = dir(stage3_pattern);
        if isempty(stage3_files)
            error('ModelAdequacy_Master_v2_001:Stage3NoOutput', ...
                'Stage 3 did not produce a results file matching %s', stage3_pattern);
        end
        [~, i] = max([stage3_files.datenum]);
        stage3_file = stage3_files(i).name;
    end

    stage3_data = load(stage3_file);
    requires_integration = stage3_data.conditional_results.requires_integration;
    fprintf('\n  Stage 3 verdict: %s\n', ternary(requires_integration, ...
        'INTEGRATION BENEFICIAL', 'integration NOT beneficial'));

elseif requires_conditional && start_from_stage > 3
    % Resuming from Stage 4; locate existing Stage 3
    stage3_pattern = sprintf('stage3_conditional_L%d_*.mat', tractability_level);
    stage3_files = dir(stage3_pattern);
    if isempty(stage3_files)
        error('ModelAdequacy_Master_v2_001:NoStage3', ...
            'Cannot start from Stage %d: no Stage 3 results for Level %d.', ...
            start_from_stage, tractability_level);
    end
    [~, i] = max([stage3_files.datenum]);
    stage3_file = stage3_files(i).name;
    stage3_data = load(stage3_file);
    requires_integration = stage3_data.conditional_results.requires_integration;
end

%% STAGE 4: HIERARCHICAL INTEGRATION
if requires_conditional && requires_integration && start_from_stage <= 4
    fprintf('\n==== STAGE 4: HIERARCHICAL INTEGRATION ====\n\n');
    try
        ModelAdequacy_Stage4_Integration_v003(stage3_file);
    catch ME
        fprintf('  Stage 4 FAILED: %s\n', ME.message);
        rethrow(ME);
    end
    fprintf('\n  Stage 4 complete.\n');
elseif requires_conditional && ~requires_integration
    fprintf('\n  Stages 4 skipped: Stage 3 concluded integration is not warranted.\n');
elseif ~requires_conditional
    fprintf('\n  Stages 3-4 skipped: Stage 2 concluded global model is adequate.\n');
end

%% FINAL SUMMARY
fprintf('\n==== FRAMEWORK COMPLETION ====\n');

try
    final_s2 = load(stage2_file);
    fprintf('  Kitchen Sink R^2: %.4f\n', final_s2.adequacy_results.stage1_r_squared);
    fprintf('  Model coefficients: %d\n', final_s2.adequacy_results.stage1_coefficients);
    fprintf('  Criterion 1 (stability): %s\n', ynFromAdequate(final_s2.adequacy_results, 'stability_analysis'));
    fprintf('  Criterion 2 (residuals): %s\n', ynFromAdequate(final_s2.adequacy_results, 'residual_analysis'));
    fprintf('  Criterion 3 (CV):        %s\n', ynFromAdequate(final_s2.adequacy_results, 'cv_analysis'));
    fprintf('  Overall: %s\n', ternary(final_s2.adequacy_results.global_adequate, 'ADEQUATE', 'INADEQUATE'));
catch
    fprintf('  (Could not reload Stage 2 for summary.)\n');
end

fprintf('\n  Output files (latest per stage):\n');
outputPatterns = { ...
    sprintf('stage1_results_L%d_*.mat',    tractability_level), 'Stage 1'; ...
    sprintf('stage2_adequacy_L%d_*.mat',   tractability_level), 'Stage 2'; ...
    sprintf('stage3_conditional_L%d_*.mat',tractability_level), 'Stage 3'; ...
    sprintf('stage4_integration_L%d_*.mat',tractability_level), 'Stage 4'; ...
    'model_adequacy_FINAL_*.html',                              'Final HTML'};
for iP = 1:size(outputPatterns,1)
    p = outputPatterns{iP,1};
    lbl = outputPatterns{iP,2};
    files = dir(p);
    % Exclude renamed contaminated outputs from the "latest" selection so the
    % summary doesn't misreport an obsolete file as authoritative.
    if ~isempty(files)
        files = files(~contains({files.name}, 'CONTAMINATED'));
    end
    if ~isempty(files)
        [~, i] = max([files.datenum]);
        fprintf('    %s: %s\n', lbl, files(i).name);
    else
        fprintf('    %s: (not generated)\n', lbl);
    end
end

fprintf('\n================================================================\n');
fprintf('  MODEL ADEQUACY PIPELINE v004 COMPLETE\n');
fprintf('  Level %d (%s) | Data: %s\n', tractability_level, current_name, ...
    ternary(use_database,'DB','ersatz'));
fprintf('================================================================\n\n');

end

%% SUPPORTING UTILITY FUNCTIONS

function result = ternary(condition, true_val, false_val)
    if condition, result = true_val; else, result = false_val; end
end

function s = bool2str(b)
    if b, s = 'true'; else, s = 'false'; end
end

function s = ynFromAdequate(adequacy_results, fieldname)
    if isfield(adequacy_results, fieldname) && isfield(adequacy_results.(fieldname),'adequate')
        if adequacy_results.(fieldname).adequate, s = 'PASS'; else, s = 'FAIL'; end
    else
        s = 'UNKNOWN';
    end
end
