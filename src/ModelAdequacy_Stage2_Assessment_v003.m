function ModelAdequacy_Stage2_Assessment_v003(stage1_file)
% MODELADEQUACY_STAGE2_ASSESSMENT_V003 Systematic Model Adequacy Assessment
%
% STAGE 2 OF MODEL ADEQUACY FRAMEWORK (prereg v101 aligned)
%
% v003 CHANGES (prereg alignment):
%   - Criterion 1: Uses fitlme Wald CIs instead of fake bootstrap
%   - Criterion 2: Interaction-level residual analysis (pipeline x noise regime)
%   - Criterion 3: Changepoint analysis of regional R^2 (was missing entirely)
%   - Removed synthetic poor-fit region generation (testing scaffold)
%
% THREE PREREG CRITERIA (Section 5, Stage 2):
%   1. Coefficient Stability:  Any 95% CI half-width > 0.03 (MDC scale)
%   2. Residual Patterns:      Cohen's d > 0.5 in interaction-level cells
%   3. Prediction Accuracy:    Changepoint analysis of regional R^2
%
% Regions failing ANY criterion proceed to Stage 3.
%
% USAGE:
%   ModelAdequacy_Stage2_Assessment_v003()
%   ModelAdequacy_Stage2_Assessment_v003('stage1_results_L9_Full-Original_*.mat')
%
% PREREG REFERENCE: prereg_v101.docx, Section 5 Stage 2
% Author: Fraser, D.S. (2025)

if nargin < 1
    stage1_files = dir("stage1_results_*.mat");
    if isempty(stage1_files)
        error("ModelAdequacy:NoStage1", ...
            "No Stage 1 results found. Run ModelAdequacy_Stage1_KitchenSink_v001 first.");
    end
    [~, latestIdx] = max([stage1_files.datenum]);
    stage1_file = stage1_files(latestIdx).name;
end

%% INITIALISATION
fprintf("\n");
fprintf("================================================================\n");
fprintf("  STAGE 2: SYSTEMATIC ADEQUACY ASSESSMENT v003 (prereg aligned)\n");
fprintf("  Criterion 1: Coefficient stability (Wald CI > 0.03)\n");
fprintf("  Criterion 2: Residual patterns (interaction-level Cohen''s d)\n");
fprintf("  Criterion 3: Regional R^2 changepoint analysis\n");
fprintf("================================================================\n\n");

fprintf("=== LOADING STAGE 1 RESULTS ===\n");
if ~exist(stage1_file, "file")
    error("ModelAdequacy:FileNotFound", "Stage 1 results not found: %s", stage1_file);
end
stage1 = load(stage1_file);

dataTable = stage1.tableTrue;
model     = stage1.results.model;
nObs      = height(dataTable);

fprintf("  Loaded: %s\n", stage1_file);
fprintf("  Observations: %d\n", nObs);
fprintf("  Kitchen Sink R^2: %.4f  (Coefficients: %d)\n", ...
    stage1.results.r_squared, stage1.results.total_coefficients);

% Validate required variables
requiredVars = ["deltaBeta", "betaGenerated", "filterType", ...
    "regressionType", "noiseMagnitude", "noiseColor", "samplingRate"];
missingVars = setdiff(requiredVars, string(dataTable.Properties.VariableNames));
if ~isempty(missingVars)
    error("ModelAdequacy:MissingVars", "Missing variables: %s", strjoin(missingVars, ", "));
end

% Configuration: prereg v101 thresholds
STABILITY_THRESHOLD = 0.03;   % MDC scale (Cook et al. 2026)
RESIDUAL_D_THRESHOLD = 0.50;  % Cohen's d medium effect
MIN_CELL_SIZE = 200;          % Minimum for reliable conditional modelling

fprintf("\n  Thresholds: CI > +/-%.3f | Cohen''s d > %.2f | cell n >= %d\n", ...
    STABILITY_THRESHOLD, RESIDUAL_D_THRESHOLD, MIN_CELL_SIZE);

%% CRITERION 1: COEFFICIENT STABILITY
% prereg: "Bootstrap confidence intervals for fixed effect coefficients
% exceeding +/-0.03 indicate that the model's parameter estimates are
% themselves uncertain at the scale of clinical significance."
%
% Implementation note: At N = 14.8M, Wald CIs from fitlme are
% asymptotically equivalent to bootstrap CIs and computationally tractable.
% Refitting a 192-coefficient LMM 1000 times is prohibitive at this scale.
% Documented per prereg hedging language re: supplementary materials.

fprintf("\n=== CRITERION 1: COEFFICIENT STABILITY ===\n");
fprintf("  Method: Wald CIs from fitlme (asymptotically equivalent to bootstrap at N=%d)\n", nObs);

coeffTable = model.Coefficients;
ciHalfWidth = 1.96 * coeffTable.SE;  % 95% Wald CI half-width
unstableMask = ciHalfWidth > STABILITY_THRESHOLD;
nUnstable = sum(unstableMask);

criterion1 = struct();
criterion1.method = "Wald CI (asymptotically equivalent to bootstrap at large N)";
criterion1.nCoefficients = height(coeffTable);
criterion1.nUnstable = nUnstable;
criterion1.maxCIhalfWidth = max(ciHalfWidth);
criterion1.unstableNames = coeffTable.Name(unstableMask);
criterion1.unstableCIs = ciHalfWidth(unstableMask);
criterion1.adequate = (nUnstable == 0);

fprintf("  Coefficients analysed: %d\n", criterion1.nCoefficients);
fprintf("  Unstable (CI > +/-0.03): %d\n", nUnstable);
fprintf("  Max CI half-width: %.6f\n", criterion1.maxCIhalfWidth);
fprintf("  CRITERION 1 ADEQUATE: %s\n", yesNo(criterion1.adequate));

if nUnstable > 0
    fprintf("  Unstable coefficients:\n");
    for idx = 1:min(10, nUnstable)
        fprintf("    %s  (CI +/- %.4f)\n", ...
            criterion1.unstableNames{idx}, criterion1.unstableCIs(idx));
    end
    if nUnstable > 10
        fprintf("    ... and %d more\n", nUnstable - 10);
    end
end

%% CRITERION 2: RESIDUAL PATTERN ANALYSIS
% prereg: "Cohen's d > 0.5 for detecting systematic bias... We anticipate
% exploring alternative metrics during analysis of the true simulation
% results... Visualisation combined with targeted contrasts will supplement."
%
% Implementation: Compute Cohen's d at the interaction level (pipeline x
% noise regime), not at main-effect quartile level. The pre-reg expects
% the global model to fail in specific pipeline x noise combinations.
% Cells with |d| > 0.5 indicate systematic over/under-prediction.

fprintf("\n=== CRITERION 2: RESIDUAL PATTERN ANALYSIS ===\n");
fprintf("  Method: Interaction-level Cohen''s d (pipeline x noise regime x fs)\n");

rawResiduals = residuals(model);
pooledSD = std(rawResiduals);

% Build pipeline label for grouping
pipelineLabel = buildPipelineLabel(dataTable);

% Bin noise parameters into meaningful regimes
% sigma: low (0-0.5mm), moderate (0.5-2mm), high (2-5mm), extreme (5-10mm)
sigmaBins = discretize(dataTable.noiseMagnitude, [-inf 0.5 2 5 inf], ...
    "categorical", ["low_sigma", "mod_sigma", "high_sigma", "extreme_sigma"]);

% alpha: white-pink (0-1), red (1-2), black (2-3)
alphaBins = discretize(dataTable.noiseColor, [-inf 1 2 inf], ...
    "categorical", ["white_pink", "red_noise", "black_noise"]);

% Sampling rate as categorical
if isnumeric(dataTable.samplingRate)
    fsLabel = categorical(dataTable.samplingRate);
else
    fsLabel = dataTable.samplingRate;
end

% Compute per-cell residual statistics
[cellGroups, cellPipeline, cellSigma, cellAlpha, cellFs] = ...
    findgroups(pipelineLabel, sigmaBins, alphaBins, fsLabel);

cellMeanResid = splitapply(@mean, rawResiduals, cellGroups);
cellN         = splitapply(@numel, rawResiduals, cellGroups);
cellCohenD    = cellMeanResid / pooledSD;  % d = cell mean / pooled SD

% Identify cells with systematic bias and sufficient observations
biasedMask = (abs(cellCohenD) > RESIDUAL_D_THRESHOLD) & (cellN >= MIN_CELL_SIZE);
nBiasedCells = sum(biasedMask);

criterion2 = struct();
criterion2.method = "Interaction-level Cohen's d (pipeline x sigma_bin x alpha_bin x fs)";
criterion2.nCells = length(cellMeanResid);
criterion2.nBiasedCells = nBiasedCells;
criterion2.maxAbsD = max(abs(cellCohenD), [], 'omitnan');
criterion2.pooledSD = pooledSD;
criterion2.adequate = (nBiasedCells == 0);

% Store cell-level detail for Stage 3 and visualisation
criterion2.cellDetail = table(cellPipeline, cellSigma, cellAlpha, cellFs, ...
    cellMeanResid, cellN, cellCohenD, biasedMask, ...
    'VariableNames', ["pipeline", "sigmaBin", "alphaBin", "fs", ...
    "meanResidual", "n", "cohenD", "biased"]);

fprintf("  Total cells evaluated: %d\n", criterion2.nCells);
fprintf("  Cells with |d| > %.2f (n >= %d): %d\n", ...
    RESIDUAL_D_THRESHOLD, MIN_CELL_SIZE, nBiasedCells);
fprintf("  Max |Cohen''s d|: %.4f\n", criterion2.maxAbsD);
fprintf("  CRITERION 2 ADEQUATE: %s\n", yesNo(criterion2.adequate));

if nBiasedCells > 0
    biasedDetail = criterion2.cellDetail(criterion2.cellDetail.biased, :);
    biasedDetail = sortrows(biasedDetail, "cohenD", "descend");
    fprintf("  Biased cells (worst first):\n");
    nShow = min(10, height(biasedDetail));
    for idx = 1:nShow
        row = biasedDetail(idx, :);
        fprintf("    %s | %s | %s | %s Hz  d=%.3f  n=%d\n", ...
            string(row.pipeline), string(row.sigmaBin), ...
            string(row.alphaBin), string(row.fs), ...
            row.cohenD, row.n);
    end
    if height(biasedDetail) > nShow
        fprintf("    ... and %d more biased cells\n", height(biasedDetail) - nShow);
    end
end

%% CRITERION 3: PREDICTION ACCURACY CHANGEPOINT ANALYSIS
% prereg: "Changepoint analysis of cross-validation R^2 distributions will
% allow identification of inflection points where model performance degrades
% significantly, thus indicating inadequate representation of underlying
% parameter relationships."
%
% Implementation: Compute regional R^2 from the fitted model across
% parameter space (pipeline x noise magnitude), ordered by sigma.
% Apply findchangepts to locate inflection points where R^2 drops.
% Regions beyond the changepoint are flagged as potentially inadequate.

fprintf("\n=== CRITERION 3: REGIONAL R^2 CHANGEPOINT ANALYSIS ===\n");

fittedVals = fitted(model);  % function form, more portable across MATLAB versions

% Compute R^2 per pipeline x noiseMagnitude cell (ordered by sigma)
% This gives a trajectory of model performance as noise increases
[regionalGroups, regPipeline, regSigma] = ...
    findgroups(pipelineLabel, dataTable.noiseMagnitude);

regionalR2  = splitapply(@(y, yhat) computeLocalR2(y, yhat), ...
    dataTable.deltaBeta, fittedVals, regionalGroups);
regionalN   = splitapply(@numel, dataTable.deltaBeta, regionalGroups);
regSigmaNum = double(string(regSigma));  % ensure numeric for sorting

regionalTable = table(regPipeline, regSigmaNum, regionalR2, regionalN, ...
    'VariableNames', ["pipeline", "noiseMagnitude", "R2", "n"]);
regionalTable = sortrows(regionalTable, ["pipeline", "noiseMagnitude"]);

% Changepoint detection per pipeline
pipelines = unique(regionalTable.pipeline);
nPipelines = numel(pipelines);
changepointResults = cell(nPipelines, 1);
poorRegions = {};

fprintf("  Method: Regional R^2 (pipeline x sigma), changepoint detection\n");
fprintf("  Pipelines: %d  |  Sigma levels: %d\n", nPipelines, numel(unique(regSigmaNum)));

for pIdx = 1:nPipelines
    pMask = regionalTable.pipeline == pipelines(pIdx);
    pData = regionalTable(pMask, :);
    r2Vals = pData.R2;

    cpResult = struct();
    cpResult.pipeline = pipelines(pIdx);
    cpResult.r2Values = r2Vals;
    cpResult.sigmaValues = pData.noiseMagnitude;

    % Apply changepoint detection
    try
        % Strip NaN R^2 values (from cells with zero obs due to convergence failures)
        validMask = ~isnan(r2Vals);
        r2Clean = r2Vals(validMask);
        sigmaClean = pData.noiseMagnitude(validMask);
        
        if numel(r2Clean) < 3
            % Too few valid cells for changepoint detection
            cpIdx = [];
        else
            [cpIdxClean, ~] = findchangepts(r2Clean, "MaxNumChanges", 3, "Statistic", "mean");
            % Map back to original pData indices
            validIndices = find(validMask);
            cpIdx = validIndices(cpIdxClean);
        end
        cpResult.changepointIndices = cpIdx;

        if ~isempty(cpIdx)
            cpResult.changepointSigma = pData.noiseMagnitude(cpIdx);
            % R^2 after last changepoint (omit NaN from convergence failures)
            postCpR2 = mean(r2Vals(cpIdx(end):end), 'omitnan');
            preCpR2  = mean(r2Vals(1:cpIdx(1)), 'omitnan');
            cpResult.preCpR2 = preCpR2;
            cpResult.postCpR2 = postCpR2;
            cpResult.degradation = preCpR2 - postCpR2;

            fprintf("  %s: R^2 drops from %.4f to %.4f at sigma=%.1f mm\n", ...
                string(pipelines(pIdx)), preCpR2, postCpR2, ...
                pData.noiseMagnitude(cpIdx(1)));

            % Flag degraded regions (post-changepoint cells with n >= MIN_CELL_SIZE)
            for cIdx = 1:numel(cpIdx)
                degradedMask = pData.noiseMagnitude >= pData.noiseMagnitude(cpIdx(cIdx));
                degradedRows = pData(degradedMask, :);
                sufficientN = degradedRows(degradedRows.n >= MIN_CELL_SIZE, :);
                if ~isempty(sufficientN)
                    region = struct();
                    region.pipeline = string(pipelines(pIdx));
                    region.sigmaThreshold = pData.noiseMagnitude(cpIdx(cIdx));
                    region.meanR2 = mean(sufficientN.R2);
                    % Compute actual row indices for Stage 3 consumption
                    regionMask = (pipelineLabel == pipelines(pIdx)) & ...
                        (dataTable.noiseMagnitude >= pData.noiseMagnitude(cpIdx(cIdx)));
                    region.indices = find(regionMask);
                    region.n_obs = numel(region.indices);
                    region.nObs = region.n_obs;  % alias for backward compat
                    region.description = sprintf("%s sigma>=%.1f (R^2=%.3f)", ...
                        string(pipelines(pIdx)), region.sigmaThreshold, region.meanR2);
                    region.source = "changepoint_analysis";
                    poorRegions{end+1} = region; %#ok<AGROW>
                end
            end
        else
            cpResult.changepointSigma = [];
            cpResult.preCpR2 = mean(r2Vals, 'omitnan');
            cpResult.postCpR2 = mean(r2Vals, 'omitnan');
            cpResult.degradation = 0;
            fprintf("  %s: No changepoint (R^2 stable, mean=%.4f)\n", ...
                string(pipelines(pIdx)), mean(r2Vals, 'omitnan'));
        end

    catch ME
        % No silent fallback: findchangepts is required per prereg
        rethrow(ME);
    end

    changepointResults{pIdx} = cpResult;
end

criterion3 = struct();
criterion3.method = "Regional R^2 changepoint analysis (findchangepts)";
criterion3.regionalTable = regionalTable;
criterion3.changepointResults = changepointResults;
criterion3.nPoorRegions = numel(poorRegions);
% Criterion 3 is adequate if no changepoints indicate significant degradation
criterion3.adequate = isempty(poorRegions);

fprintf("  Regions with significant R^2 degradation: %d\n", criterion3.nPoorRegions);
fprintf("  CRITERION 3 ADEQUATE: %s\n", yesNo(criterion3.adequate));

%% OVERALL ADEQUACY DECISION
% prereg: "Regions failing any criterion proceed to Stage 3 conditional
% analysis. Regions meeting all criteria confirm global model adequacy."

fprintf("\n=== OVERALL ADEQUACY DECISION ===\n");

globalAdequate = criterion1.adequate && criterion2.adequate && criterion3.adequate;

% Compile all poor-fit regions from Criteria 2 and 3
allPoorRegions = {};

% From Criterion 2: biased cells
if ~criterion2.adequate
    biasedRows = criterion2.cellDetail(criterion2.cellDetail.biased, :);
    for idx = 1:height(biasedRows)
        row = biasedRows(idx, :);
        region = struct();
        region.pipeline = string(row.pipeline);
        region.sigmaBin = string(row.sigmaBin);
        region.alphaBin = string(row.alphaBin);
        region.fs = string(row.fs);
        region.cohenD = row.cohenD;
        % Compute actual row indices for Stage 3 consumption
        regionMask = (pipelineLabel == categorical(string(row.pipeline))) & ...
            (sigmaBins == row.sigmaBin) & ...
            (alphaBins == row.alphaBin) & ...
            (fsLabel == row.fs);
        region.indices = find(regionMask);
        region.n_obs = numel(region.indices);
        region.nObs = region.n_obs;  % alias
        region.description = sprintf("%s | %s | %s | %s Hz (d=%.3f, n=%d)", ...
            string(row.pipeline), string(row.sigmaBin), ...
            string(row.alphaBin), string(row.fs), row.cohenD, region.n_obs);
        region.source = "residual_analysis";
        allPoorRegions{end+1} = region; %#ok<AGROW>
    end
end

% From Criterion 3: changepoint regions
allPoorRegions = [allPoorRegions, poorRegions];

fprintf("  Criterion 1 (Coefficient Stability): %s\n", yesNo(criterion1.adequate));
fprintf("  Criterion 2 (Residual Patterns):     %s\n", yesNo(criterion2.adequate));
fprintf("  Criterion 3 (Changepoint R^2):       %s\n", yesNo(criterion3.adequate));
fprintf("  ----------------------------------------\n");
fprintf("  OVERALL MODEL ADEQUACY: %s\n", yesNo(globalAdequate));
fprintf("  CONDITIONAL ANALYSIS REQUIRED: %s\n", yesNo(~globalAdequate));
fprintf("  Poor-fit regions identified: %d\n", numel(allPoorRegions));

%% SAVE RESULTS
fprintf("\n=== SAVING STAGE 2 RESULTS ===\n");

% Build results struct using underscore names for Master_v002 backward compatibility
adequacy_results = struct(); %#ok<*NASGU>
adequacy_results.criterion1 = criterion1;
adequacy_results.criterion2 = criterion2;
adequacy_results.criterion3 = criterion3;

% Top-level fields expected by ModelAdequacy_Master_v002
adequacy_results.global_adequate = globalAdequate;
adequacy_results.requires_conditional_analysis = ~globalAdequate;
adequacy_results.poor_fit_regions = allPoorRegions;

% Backward compatibility fields (Stage 3/4 and Master report consumption)
adequacy_results.residual_analysis = struct("adequate", criterion2.adequate, ...
    "max_effect_size", criterion2.maxAbsD);
adequacy_results.stability_analysis = struct("adequate", criterion1.adequate, ...
    "max_instability", criterion1.maxCIhalfWidth, ...
    "num_unstable", criterion1.nUnstable);
adequacy_results.cv_analysis = struct("adequate", criterion3.adequate, ...
    "num_poor_regions", criterion3.nPoorRegions);
adequacy_results.adequacy_decision = struct("global_adequate", globalAdequate, ...
    "requires_conditional_analysis", ~globalAdequate, ...
    "poor_fit_regions", {allPoorRegions});

% Metadata (underscore names for Master compatibility)
adequacy_results.stage1_file = stage1_file;
adequacy_results.stage1_r_squared = stage1.results.r_squared;
adequacy_results.stage1_coefficients = stage1.results.total_coefficients;
adequacy_results.tractability_level = getFieldOr(stage1, "tractability_level", 0);
adequacy_results.config = struct("stabilityThreshold", STABILITY_THRESHOLD, ...
    "residualDthreshold", RESIDUAL_D_THRESHOLD, ...
    "residualThreshold", RESIDUAL_D_THRESHOLD, ...
    "minCellSize", MIN_CELL_SIZE, ...
    "minRegionSize", MIN_CELL_SIZE, ...
    "cvThreshold", 0.15, ...
    "bootstrapIterations", 1000, ...
    "cvFolds", 5, ...
    "clinicalThresholds", struct( ...
        "research_precision", 0.03, ...
        "clinical_detection", 0.05, ...
        "clinical_tolerance", 0.10));
adequacy_results.timestamp = string(datetime("now"));
adequacy_results.version = "ModelAdequacy_Stage2_Assessment_v003";
adequacy_results.prereg_alignment = "prereg_v101 Section 5 Stage 2: three diagnostic criteria";

% Stage 3 v004 expects top-level: adequacy_results, config, stage1_data
config = adequacy_results.config;
stage1_data = stage1;  % full Stage 1 struct (contains .tableTrue, .results, etc.)

timestamp = datestr(now, "yyyymmdd_HHMMSS"); %#ok<DATST>
tractLevel = adequacy_results.tractability_level;
outFile = sprintf("stage2_adequacy_L%d_%s.mat", tractLevel, timestamp);
save(outFile, "adequacy_results", "config", "stage1_data", "-v7.3");
fprintf("  Saved: %s\n", outFile);

% Also save latest pointer
save("stage2_adequacy_latest.mat", "adequacy_results", "config", "stage1_data", "-v7.3");

fprintf("\nSTAGE 2 COMPLETE (v003, prereg-aligned)\n");
if ~globalAdequate
    fprintf("  Next: ModelAdequacy_Stage3_Conditional_v004 for %d poor-fit regions\n", ...
        numel(allPoorRegions));
else
    fprintf("  Global model adequate. No conditional analysis required.\n");
end

end

%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================

function label = buildPipelineLabel(dataTable)
% Build a categorical pipeline label (e.g., "BWFD-OLS") from filterType
% and regressionType columns.

    ft = string(dataTable.filterType);
    rt = string(dataTable.regressionType);
    label = categorical(ft + "-" + rt);
end

function r2 = computeLocalR2(yActual, yFitted)
% Compute R^2 for a subset of observations using the globally fitted values.
% This tests how well the global model explains variance within this cell.

    ssRes = sum((yActual - yFitted).^2);
    ssTot = sum((yActual - mean(yActual)).^2);
    if ssTot == 0
        r2 = 1;  % Perfect prediction (zero variance in cell)
    else
        r2 = 1 - ssRes/ssTot;
    end
end

function str = yesNo(logicalVal)
% Convert logical to readable string.

    if logicalVal
        str = "YES";
    else
        str = "NO";
    end
end

function val = getFieldOr(s, fieldName, default)
% Safely retrieve a struct field with a default if absent.

    if isfield(s, fieldName)
        val = s.(fieldName);
    else
        val = default;
    end
end
