% checkInvertibilityForEmpirical_v2_003.m  v2_002 with Fraser on the N=48 corpus and the corrected iPad scale (2026-11-04).
% ONLY CHANGES: Fraser's registry entry reads constellationFraser_v002.mat, noiseCharacterisation_fraserN48.mat and
% importDB_fraser_v002, with sigmaToMM = 1/ipadScaleM4_v001() (v2_002: 1/10.41793, 2x too small; it enters the VGF
% conversion s^(1-b)); a missing constellation or noise file now stops the run (v2_002 warned and skipped the dataset);
% outputs invertibilityCoverage_v2_003.*. Other datasets unchanged. v2_002's header follows.
% checkInvertibilityForEmpirical_v2_002.m  Adds Fraser to the v2_001 registry
%
% Identical mechanism to v2_001 (per-trial rising-branch coverage against
% monotonicSegments_v2_002.mat -- see v2_001 header for full method detail).
% Fraser (Analysis_2026) is added as a seventh registry entry, replacing
% Pilot's role as the confirmatory test-retest dataset. Note Pilot itself
% was never in this registry -- v2_001's six datasets (Zarandi, Cook CTRL/ASD,
% Dhieb, Hickman PLAC/HALO) are what produced the "35/36 FAIL" R6 statistic;
% Pilot's SEM-adequacy (R5) and local-monotonicity (R6) numbers instead come
% from a direct grid-centroid lookup and loopClosureSlopeMargin_v001.m
% respectively -- neither of which goes through this script. Fraser is being
% added here at Dagmar's explicit direction (2026-08-06): Fraser is a
% canonical/core dataset in its own right, on the same footing as Cook and
% Hickman, not merely inheriting whatever path Pilot happened to take.
%
% v2_001 is left untouched -- its existing invertibilityCoverage_v2_001.mat
% and _summary.txt are cited directly in R6 ("35/36 FAIL"); this version
% produces separately-named v2_002 outputs so that citation stays valid
% until the paper's R6 prose is deliberately updated to the 42-cell (7
% dataset) picture.
%
% Fraser, D.S. (2026)
% See also: checkInvertibilityForEmpirical_v2_001, constellationFraser_v001
%
% Per-dataset, per-pipeline empirical-trial coverage of the invertibility
% maps built in monotonicSegments_v2_001.mat. Implements the go/no-go gate
% for Step 3 (constellationCCC_v002).
%
% For each empirical trial:
%   1. Snap (alpha, sigma_mm, fs, VGF_mm-equiv) to simulation grid via
%      checkTrialInvertibility_v001.
%   2. Test for each pipeline whether the trial's beta_rec sits inside the
%      rising-branch image, descending-branch image, both, or neither at
%      the matched grid coord.
%
% Per-(dataset x pipeline) coverage rates feed the gate:
%   Pass         : rising-branch coverage >= 95%  -> proceed Step 3, low-beta default
%   Conditional  : 90-95%                         -> proceed Step 3 with prominent disclosure
%   Fail         : < 90%                          -> reconsider rescue strategies (TODO Sec 5)
%
% Inputs (registry-driven, mirrors runConstellationCCC_v001):
%   - src/monotonicSegments_v2_001.mat                          (Step 1 output)
%   - src/constellation<Dataset>_v001.mat                     (betaCanon, vgfAll, canonLabels)
%   - src/noiseCharacterisation_<dataset>.mat                 (bioResults)
%   - functions/importDB_<dataset>_v001/v002 importer handles (per-trial fs, subjectID)
%
% Outputs:
%   src/invertibilityCoverage_v2_002.mat
%     .perTrial         table, one row per (dataset, trial, pipeline) -- long form
%     .perDataPipe      table, one row per (dataset, pipeline)        -- summary
%     .perDataset       table, one row per dataset                    -- pipeline-mean
%     .registry         dataset registry used
%     .opts             options bundle (BetaStrategy etc.)
%   src/invertibilityCoverage_v2_002_summary.txt   human-readable gate verdict
%
% Refs: TODO_CCC_Invertibility.md Sec 4 Step 2 + Step 1 status block;
%       claude.md Key findings #15 (inverse-lookup pivot), #17 (bend-back
%       is dominant regime, 70-79% beta-ambiguous coverage in simulation grid).
%
% Fraser, D.S. (2026)

clearvars

srcDir = fileparts(mfilename("fullpath"));
addpath(genpath(fullfile(srcDir, "functions")));
cd(srcDir);

%% Options (matched to constellationCCC_v001 defaults for like-for-like)
opts = struct();
opts.BetaStrategy   = "conditional_median";   % how to estimate per-pipeline beta_rec
                                              % NOTE: for Step 2 we use the FULL [nTri x 6]
                                              % betaCanon directly; opts.BetaStrategy is
                                              % retained only for documentation parity
                                              % with constellationCCC_v001.
opts.VGFStrategy    = "conditional_median";
opts.VGFConversion  = @(b, s) s^(1 - b);      % dimensional analysis (see constellationCCC_v001)
opts.AlphaSource    = "ira_alphaMean";
opts.GatePass       = 0.95;
opts.GateConditional= 0.90;

%% Load Step 1 output
monoFile = fullfile(srcDir, "monotonicSegments_v2_002.mat");
if ~isfile(monoFile)
    error("checkInvertibility:MissingStep1", ...
        "Step 1 output not found: %s\nRun buildMonotonicSegments_v001 first.", monoFile);
end
fprintf("=== checkInvertibilityForEmpirical_v2_003 ===\n");
fprintf("Loading %s ...\n", monoFile);
monoBundle = load(monoFile, "params", "pipeOrder", ...
    "alphaGrid", "sigmaGrid", "fsGrid", "betaGenGrid", "VGFGrid", ...
    "betaSeg", "VGFSeg");
nP = numel(monoBundle.pipeOrder);
fprintf("  %d pipelines, %d alpha, %d sigma, %d fs, %d VGF grid points\n", ...
    nP, numel(monoBundle.alphaGrid), numel(monoBundle.sigmaGrid), ...
    numel(monoBundle.fsGrid), numel(monoBundle.VGFGrid));

%% Dataset registry (mirrors runConstellationCCC_v001)
registry = { ...
    struct( ...
        "name",             "Zarandi", ...
        "constellationMat", fullfile(srcDir, "constellationZarandi_v001.mat"), ...
        "noiseMat",         fullfile(srcDir, "noiseCharacterisation_zarandi.mat"), ...
        "importerHandle",   @importDB_zarandi_v001, ...
        "importerArgs",     {{ "Verbose", false }}, ...
        "sigmaToMM",        10.0); ...
    struct( ...
        "name",             "Cook_CTRL", ...
        "constellationMat", fullfile(srcDir, "constellationCook_v001.mat"), ...
        "noiseMat",         fullfile(srcDir, "noiseCharacterisation_cook.mat"), ...
        "importerHandle",   @importDB_cook_v002, ...
        "importerArgs",     {{ "Group", "CTRL", "Tasks", 7, "Verbose", false }}, ...
        "sigmaToMM",        0.248); ...
    struct( ...
        "name",             "Cook_ASD", ...
        "constellationMat", fullfile(srcDir, "constellationCookASD_v001.mat"), ...
        "noiseMat",         fullfile(srcDir, "noiseCharacterisation_cookASD.mat"), ...
        "importerHandle",   @importDB_cook_v002, ...
        "importerArgs",     {{ "Group", "ASD", "Tasks", 7, "Verbose", false }}, ...
        "sigmaToMM",        0.248); ...
    struct( ...
        "name",             "Dhieb", ...
        "constellationMat", fullfile(srcDir, "constellationDhieb_v001.mat"), ...
        "noiseMat",         fullfile(srcDir, "noiseCharacterisation_dhieb.mat"), ...
        "importerHandle",   @importDB_dhieb_v001, ...
        "importerArgs",     {{}}, ...
        "sigmaToMM",        0.1478); ...
    struct( ...
        "name",             "Hickman_PLAC", ...
        "constellationMat", fullfile(srcDir, "constellationHickmanPLAC_v001.mat"), ...
        "noiseMat",         fullfile(srcDir, "noiseCharacterisation_hickmanPLAC.mat"), ...
        "importerHandle",   @importDB_hickman_v002, ...
        "importerArgs",     {{ "Study", 2, "Group", "PLAC", "Shapes", 3, "Verbose", false }}, ...
        "sigmaToMM",        0.248); ...
    struct( ...
        "name",             "Hickman_HALO", ...
        "constellationMat", fullfile(srcDir, "constellationHickmanHALO_v001.mat"), ...
        "noiseMat",         fullfile(srcDir, "noiseCharacterisation_hickmanHALO.mat"), ...
        "importerHandle",   @importDB_hickman_v002, ...
        "importerArgs",     {{ "Study", 2, "Group", "HALO", "Shapes", 3, "Verbose", false }}, ...
        "sigmaToMM",        0.248); ...
    struct( ...
        "name",             "Fraser", ...
        "constellationMat", fullfile(srcDir, "constellationFraser_v002.mat"), ...          % v2_003
        "noiseMat",         fullfile(srcDir, "noiseCharacterisation_fraserN48.mat"), ...  % v2_003
        "importerHandle",   @importDB_fraser_v002, ...                                    % v2_003
        "importerArgs",     {{ "Verbose", false }}, ...
        "sigmaToMM",        1/ipadScaleM4_v001()); ...                                    % v2_003: UIKit points
    };
nSets = numel(registry);

%% Per-dataset loop
allLong       = cell(nSets, 1);          % long-form per (dataset, trial, pipeline)
perDataPipeC  = cell(nSets, 1);          % per (dataset, pipeline) summary
datasetNames  = strings(nSets, 1);

for k = 1:nSets
    spec = registry{k};
    fprintf("\n--- Dataset %d/%d: %s ---\n", k, nSets, spec.name);
    datasetNames(k) = spec.name;

    % --- Load constellation + noise ---
    if ~isfile(spec.constellationMat)
        error("checkInvertibility:MissingFile", "%s", "Constellation .mat not found: " + spec.constellationMat);   % v2_003: fail loud
    end
    if ~isfile(spec.noiseMat)
        error("checkInvertibility:MissingFile", "%s", "Noise .mat not found: " + spec.noiseMat);   % v2_003: fail loud
    end
    emp = load(spec.constellationMat, "betaCanon", "vgfAll", "canonLabels");
    nz  = load(spec.noiseMat, "bioResults");
    bio = nz.bioResults;
    nTri = size(emp.betaCanon, 1);

    % --- Pipeline order check ---
    if ~isequal(string(emp.canonLabels(:))', monoBundle.pipeOrder(:)')
        error("checkInvertibility:PipelineOrder", ...
            "%s canonLabels [%s] != monoBundle.pipeOrder [%s]", ...
            spec.name, strjoin(string(emp.canonLabels), ","), ...
            strjoin(monoBundle.pipeOrder, ","));
    end

    % --- Re-import for per-trial fs and subjectID ---
    trials = spec.importerHandle(spec.importerArgs{:});
    if numel(trials) ~= nTri
        error("checkInvertibility:LengthMismatch", ...
            "%s: importer returned %d trials but constellation has %d.", ...
            spec.name, numel(trials), nTri);
    end
    subIDs   = arrayfun(@(t) string(t.subjectID), trials).';
    fsPerTri = arrayfun(@(t) double(t.fs),       trials).';

    % --- Empirical alpha & sigma per trial ---
    if ~ismember(opts.AlphaSource, bio.Properties.VariableNames)
        available = strjoin(string(bio.Properties.VariableNames), ", ");
        error("checkInvertibility:NoAlphaCol", ...
            "AlphaSource '%s' not in bioResults for %s. Available: %s", ...
            opts.AlphaSource, spec.name, available);
    end
    alphaEmp   = double(bio.(opts.AlphaSource));
    sigmaEmpMM = double(bio.sigmaMean) * spec.sigmaToMM;

    % --- VGF: median across pipelines, then beta-aware unit conversion ---
    vgfCanon = canonicaliseVGF(emp.vgfAll);
    switch opts.VGFStrategy
        case "conditional_median"
            vgfTrueEst = median(vgfCanon, 2, "omitnan");
        case "conditional_sgIRLS"
            vgfTrueEst = vgfCanon(:, 6);
        otherwise
            error("checkInvertibility:BadVGFStrategy", ...
                "Unknown VGFStrategy '%s'.", opts.VGFStrategy);
    end
    switch opts.BetaStrategy
        case "conditional_median"
            betaTrueEst = median(emp.betaCanon, 2, "omitnan");
        case "conditional_sgIRLS"
            betaTrueEst = emp.betaCanon(:, 6);
        otherwise
            error("checkInvertibility:BadBetaStrategy", ...
                "Unknown BetaStrategy '%s'.", opts.BetaStrategy);
    end
    vgfTrueScaled = NaN(nTri, 1);
    for t = 1:nTri
        if ~isnan(vgfTrueEst(t)) && ~isnan(betaTrueEst(t))
            vgfTrueScaled(t) = vgfTrueEst(t) * opts.VGFConversion(betaTrueEst(t), spec.sigmaToMM);
        end
    end

    % --- Per-trial loop ---
    inRiseAll  = false(nTri, nP);
    inDescAll  = false(nTri, nP);
    inAmbAll   = false(nTri, nP);
    inNeitherAll = false(nTri, nP);
    snapWarn   = strings(nTri, 1);
    coordValid = false(nTri, 1);

    for t = 1:nTri
        result = checkTrialInvertibility_v001( ...
            alphaEmp(t), sigmaEmpMM(t), fsPerTri(t), vgfTrueScaled(t), ...
            emp.betaCanon(t, :), monoBundle);

        coordValid(t) = ~any(isnan([alphaEmp(t), sigmaEmpMM(t), ...
                                    fsPerTri(t), vgfTrueScaled(t)]));
        inRiseAll(t, :)    = result.inRise;
        inDescAll(t, :)    = result.inDesc;
        inAmbAll(t, :)     = result.inAmbiguous;
        inNeitherAll(t, :) = result.inNeither;
        snapWarn(t)        = result.snapWarnings;
    end

    % --- Long-form table for downstream re-aggregation ---
    longRows = cell(nTri * nP, 1);
    rowIdx = 0;
    for t = 1:nTri
        for p = 1:nP
            rowIdx = rowIdx + 1;
            longRows{rowIdx} = struct( ...
                "dataset",     string(spec.name), ...
                "trialIdx",    t, ...
                "subjectID",   subIDs(t), ...
                "pipeline",    monoBundle.pipeOrder(p), ...
                "alphaEmp",    alphaEmp(t), ...
                "sigmaEmpMM",  sigmaEmpMM(t), ...
                "fsEmp",       fsPerTri(t), ...
                "vgfScaled",   vgfTrueScaled(t), ...
                "betaRec",     emp.betaCanon(t, p), ...
                "coordValid",  coordValid(t), ...
                "inRise",      inRiseAll(t, p), ...
                "inDesc",      inDescAll(t, p), ...
                "inAmbiguous", inAmbAll(t, p), ...
                "inNeither",   inNeitherAll(t, p), ...
                "snapWarn",    snapWarn(t));
        end
    end
    allLong{k} = struct2table([longRows{:}]);

    % --- Per-pipeline summary for this dataset ---
    pipeRows = cell(nP, 1);
    for p = 1:nP
        valid = coordValid & ~isnan(emp.betaCanon(:, p));
        nValid = sum(valid);
        if nValid == 0
            pipeRows{p} = struct( ...
                "dataset",          string(spec.name), ...
                "pipeline",         monoBundle.pipeOrder(p), ...
                "nTrials",          nTri, ...
                "nValid",           0, ...
                "nInRise",          0, ...
                "nInDesc",          0, ...
                "nInAmbiguous",     0, ...
                "nInNeither",       0, ...
                "coverageRise",     NaN, ...
                "coverageDesc",     NaN, ...
                "coverageAmbiguous",NaN, ...
                "coverageNeither",  NaN, ...
                "verdict",          "no-data");
            continue
        end
        nInRise = sum(inRiseAll(valid, p));
        nInDesc = sum(inDescAll(valid, p));
        nInAmb  = sum(inAmbAll(valid, p));
        nInNeit = sum(inNeitherAll(valid, p));
        covRise = nInRise / nValid;
        if      covRise >= opts.GatePass,        verdict = "PASS";
        elseif  covRise >= opts.GateConditional, verdict = "CONDITIONAL";
        else                                     verdict = "FAIL";
        end
        pipeRows{p} = struct( ...
            "dataset",          string(spec.name), ...
            "pipeline",         monoBundle.pipeOrder(p), ...
            "nTrials",          nTri, ...
            "nValid",           nValid, ...
            "nInRise",          nInRise, ...
            "nInDesc",          nInDesc, ...
            "nInAmbiguous",     nInAmb, ...
            "nInNeither",       nInNeit, ...
            "coverageRise",     covRise, ...
            "coverageDesc",     nInDesc / nValid, ...
            "coverageAmbiguous",nInAmb  / nValid, ...
            "coverageNeither",  nInNeit / nValid, ...
            "verdict",          verdict);
    end
    perDataPipeC{k} = struct2table([pipeRows{:}]);

    % --- Per-dataset console snapshot ---
    fprintf("  N trials: %d   coord-valid: %d\n", nTri, sum(coordValid));
    fprintf("  Per-pipeline rising-branch coverage:\n");
    for p = 1:nP
        r = perDataPipeC{k};
        if r.nValid(p) == 0
            fprintf("    %-12s  no valid trials\n", monoBundle.pipeOrder(p));
        else
            fprintf("    %-12s  %5.1f%%  [%s]\n", ...
                monoBundle.pipeOrder(p), 100*r.coverageRise(p), r.verdict(p));
        end
    end
end

%% Cross-dataset aggregation
perTrial    = vertcat(allLong{~cellfun(@isempty, allLong)});
perDataPipe = vertcat(perDataPipeC{~cellfun(@isempty, perDataPipeC)});

% Per-dataset (pipeline-mean coverage)
perDataset = table();
[G, gName] = findgroups(perDataPipe.dataset);
perDataset.dataset           = gName;
perDataset.meanCoverageRise  = splitapply(@(x) mean(x, "omitnan"), perDataPipe.coverageRise, G);
perDataset.minCoverageRise   = splitapply(@(x) min(x, [], "omitnan"), perDataPipe.coverageRise, G);
perDataset.meanCoverageDesc  = splitapply(@(x) mean(x, "omitnan"), perDataPipe.coverageDesc, G);
perDataset.meanCoverageAmbig = splitapply(@(x) mean(x, "omitnan"), perDataPipe.coverageAmbiguous, G);

%% Save
outMat = fullfile(srcDir, "invertibilityCoverage_v2_003.mat");
save(outMat, "perTrial", "perDataPipe", "perDataset", "registry", "opts", "monoFile", "-v7");
d = dir(outMat);
fprintf("\nSaved %s (%.1f MB)\n", outMat, d.bytes/2^20);

%% Human-readable summary
sumPath = fullfile(srcDir, "invertibilityCoverage_v2_003_summary.txt");
writeSummary(sumPath, perDataPipe, perDataset, opts);
fprintf("Saved %s\n", sumPath);
fprintf("=== checkInvertibilityForEmpirical_v2_003 COMPLETE ===\n");


%% ========================  local functions  ============================

function vgfCanon = canonicaliseVGF(vgfAll)
% Reproduce constellationCCC_v001's canonicalisation. Cook-style nDeriv=2
% maps canonIdx=[1,2]; Zarandi-style nDeriv=3 maps canonIdx=[1,3].
    nTri = size(vgfAll, 1);
    nDeriv = size(vgfAll, 2);
    nReg = size(vgfAll, 3);
    if nReg ~= 3
        error("checkInvertibility:VGFAllShape", ...
            "Expected vgfAll [nTri x nDeriv x 3]; got [%d x %d x %d].", nTri, nDeriv, nReg);
    end
    switch nDeriv
        case 2,  canonIdx = [1, 2];
        case 3,  canonIdx = [1, 3];
        otherwise
            error("checkInvertibility:VGFAllShape", ...
                "Unrecognised vgfAll nDeriv = %d. Expected 2 or 3.", nDeriv);
    end
    vgfCanon = NaN(nTri, 6);
    col = 0;
    for d = canonIdx
        for r = 1:3
            col = col + 1;
            vgfCanon(:, col) = vgfAll(:, d, r);
        end
    end
end

function writeSummary(path, perDataPipe, perDataset, opts)
    fid = fopen(path, "w");
    if fid < 0
        warning("checkInvertibility:NoWrite", "Could not open %s for write", path);
        return
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, "checkInvertibilityForEmpirical_v2_003 summary\n");
    fprintf(fid, "Generated: %s\n\n", string(datetime("now")));
    fprintf(fid, "Gate thresholds: PASS >= %.0f%%, CONDITIONAL >= %.0f%%, else FAIL\n\n", ...
        100*opts.GatePass, 100*opts.GateConditional);

    fprintf(fid, "PER (DATASET, PIPELINE) COVERAGE\n");
    fprintf(fid, "%-14s %-12s %6s %6s %8s %8s %8s %8s  %s\n", ...
        "Dataset", "Pipeline", "nTri", "nVal", "covRise", "covDesc", "covAmb", "covNeit", "Verdict");
    fprintf(fid, "%s\n", repmat('-', 1, 95));
    for r = 1:height(perDataPipe)
        fprintf(fid, "%-14s %-12s %6d %6d %7.1f%% %7.1f%% %7.1f%% %7.1f%%  %s\n", ...
            perDataPipe.dataset(r), perDataPipe.pipeline(r), ...
            perDataPipe.nTrials(r), perDataPipe.nValid(r), ...
            100*perDataPipe.coverageRise(r), 100*perDataPipe.coverageDesc(r), ...
            100*perDataPipe.coverageAmbiguous(r), 100*perDataPipe.coverageNeither(r), ...
            perDataPipe.verdict(r));
    end

    fprintf(fid, "\nPER DATASET (mean across pipelines)\n");
    fprintf(fid, "%-14s %12s %12s %12s %12s\n", ...
        "Dataset", "meanCovRise", "minCovRise", "meanCovDesc", "meanCovAmb");
    fprintf(fid, "%s\n", repmat('-', 1, 70));
    for r = 1:height(perDataset)
        fprintf(fid, "%-14s %11.1f%% %11.1f%% %11.1f%% %11.1f%%\n", ...
            perDataset.dataset(r), ...
            100*perDataset.meanCoverageRise(r), ...
            100*perDataset.minCoverageRise(r), ...
            100*perDataset.meanCoverageDesc(r), ...
            100*perDataset.meanCoverageAmbig(r));
    end

    % Verdict summary
    fprintf(fid, "\nGATE OUTCOME (count of (dataset, pipeline) cells per verdict)\n");
    [G, vNames] = findgroups(perDataPipe.verdict);
    counts = splitapply(@numel, perDataPipe.verdict, G);
    for k = 1:numel(vNames)
        fprintf(fid, "  %-12s  %d\n", vNames(k), counts(k));
    end
end
