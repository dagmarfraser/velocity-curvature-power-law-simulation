% extractSimpleEffects_v003.m
% Six planned simple effects from the Stage 1 kitchen-sink LMM (prereg v101 §8.1).
%
% v003 vs v002: fixes S4 noise-conditional snap — centroid (alpha, sigma)
% values are z-scored using results.predictor_scaling before snapping against
% the (now z-scored) noiseColor/noiseMagnitude grid. Previously all three
% drawing datasets snapped to the same cell because empirical original-scale
% values were compared against z-scored grid values.
%
% v002 vs v001: applies predictor z-scoring to tableTrue before predict()
% when the stage1 mat was built with Stage1_KitchenSink_v2_001.
% Computes estimated marginal means and pairwise simple effects for the six
% prereg-specified kinematic-derivation × regression contrasts, both grand-mean
% marginalised and conditional on four empirical noise centroids (Finding #6).
%
% Runs against any stage1_results_*.mat. Validated on L2 (667K obs); designed
% for L9 (14.5M obs). Run this script the moment stage1_results_latest.mat
% is updated after BlueBEAR L9 sync.
%
% OUTPUT
%   src/simpleEffects_L<N>_<timestamp>.mat   — full results struct
%   results/simpleEffects_L<N>_summary.txt  — human-readable table for paper
%
% PIPELINE CODING (treatment/dummy coding, BWFD-OLS is reference)
%   filterType    : 2 = BWFD (reference), 6 = SG
%   regressionType: 3 = OLS  (reference), 4 = LMLS, 5 = IRLS
%
% PREREG REFERENCE: prereg_v101.docx §8.1 (Anticipated Results)
% Fraser, D.S. (2026)

clear; clc;
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'functions')));

rng(42, 'twister');   % reproducible Monte Carlo sampling

% Open diary log immediately so all console output is captured
srcDir_  = fileparts(mfilename('fullpath'));
diaryFile = fullfile(srcDir_, sprintf('simpleEffects_log_%s.txt', ...
    datestr(now, 'yyyymmdd_HHMMSS'))); %#ok<DATST>
diary(diaryFile);
diary on
fprintf('Diary: %s\n\n', diaryFile);

%% S1: LOAD AND VALIDATE
fprintf('\n================================================================\n');
fprintf('  EXTRACT SIMPLE EFFECTS v003  (prereg v101 §8.1)\n');
fprintf('================================================================\n\n');

% Prefer L9 mat; fall back to most-recent
srcDir = fileparts(mfilename('fullpath'));
mats   = dir(fullfile(srcDir, 'stage1_results_*.mat'));
if isempty(mats)
    error('extractSimpleEffects:NoStage1', '%s', ...
        'No stage1_results_*.mat found in src/. Sync from BlueBEAR first.');
end
l9Mask = contains({mats.name}, '_L9_');
if any(l9Mask)
    [~, idx] = max([mats(l9Mask).datenum]);
    tmp     = mats(l9Mask);
    matFile = fullfile(srcDir, tmp(idx).name);
else
    [~, idx] = max([mats.datenum]);
    matFile  = fullfile(srcDir, mats(idx).name);
end
fprintf('Loading: %s\n', matFile);
s = load(matFile);

% --- field validation (Fail Loud) ---
if ~isfield(s, 'results') || ~isfield(s.results, 'model')
    error('extractSimpleEffects:BadStruct', '%s', ...
        'stage1 mat missing results.model — re-run Stage 1.');
end
if ~isfield(s, 'tableTrue')
    error('extractSimpleEffects:BadStruct', '%s', ...
        'stage1 mat missing tableTrue — re-run Stage 1.');
end
model    = s.results.model;
tbl      = s.tableTrue;
nObs     = height(tbl);
nCoef    = s.results.total_coefficients;
r2       = s.results.r_squared;
tractLvl = s.tractability_level;

requiredCols = {'deltaBeta','betaGenerated','VGF','samplingRate', ...
                'filterType','regressionType','noiseMagnitude','noiseColor'};
missing = setdiff(requiredCols, tbl.Properties.VariableNames);
if ~isempty(missing)
    error('extractSimpleEffects:MissingCols', '%s', ...
        ['tableTrue missing columns: ' strjoin(missing, ', ')]);
end

if nCoef ~= 192
    warning('extractSimpleEffects:WrongLevel', '%s', ...
        sprintf(['Expected 192 coefficients (L9); got %d (L%d). ' ...
                 'Results valid but not full-space.'], nCoef, tractLvl));
end

fprintf('  Tractability level : L%d\n',   tractLvl);
fprintf('  Observations       : %d\n',    nObs);
fprintf('  R-squared          : %.4f\n',  r2);
fprintf('  Coefficients       : %d\n\n',  nCoef);

%% S1b: APPLY PREDICTOR SCALING (v2_001 stage1 mats only)
% Stage1_KitchenSink_v2_001 z-scores continuous predictors inside
% performKitchenSink (pass-by-value), so tableTrue in the mat is original
% scale but the model was fitted on z-scored data.  Apply the stored
% scaling before any predict() call.  Fail loud if scaling is present but
% a required field is missing.
if isfield(s.results, 'predictor_scaling')
    sc = s.results.predictor_scaling;
    contVars = fieldnames(sc);
    fprintf('Applying predictor z-scoring from results.predictor_scaling:\n');
    for k = 1:numel(contVars)
        v = contVars{k};
        if ~ismember(v, tbl.Properties.VariableNames)
            error('extractSimpleEffects:ScalingColMissing', '%s', ...
                sprintf('Scaling field "%s" not found in tableTrue columns.', v));
        end
        mu = sc.(v).mean;
        sg = sc.(v).std;
        if sg < eps
            error('extractSimpleEffects:ZeroStd', '%s', ...
                sprintf('predictor_scaling.%s.std = 0 — cannot z-score.', v));
        end
        tbl.(v) = (tbl.(v) - mu) / sg;
        fprintf('  %-18s  mu=%.4g  std=%.4g\n', v, mu, sg);
    end
    fprintf('  tableTrue z-scored. predict() will use scaled values.\n\n');
else
    fprintf('  predictor_scaling not found — tableTrue used as-is (v001 mat).\n\n');
end

%% S2: OVERALL filterType × regressionType INTERACTION
% Omnibus F-test for the two-way interaction (prereg §7 primary tests).
fprintf('--- S2: Overall filterType x regressionType interaction ---\n');

try
    anovaTable  = anova(model);
    termNames   = string(anovaTable.Term);
    % The pure two-way term contains both predictor names, exactly one colon,
    % and no other predictor names.
    twoWayMask  = contains(termNames, 'filterType') & ...
                  contains(termNames, 'regressionType') & ...
                  (count(termNames, ':') == 1) & ...
                  ~contains(termNames, 'betaGenerated') & ...
                  ~contains(termNames, 'VGF') & ...
                  ~contains(termNames, 'samplingRate') & ...
                  ~contains(termNames, 'noiseMagnitude') & ...
                  ~contains(termNames, 'noiseColor');
    if any(twoWayMask)
        fRow = anovaTable(twoWayMask, :);
        intF   = fRow.FStat;
        intP   = fRow.pValue;
        intDF1 = fRow.DF1;
        intDF2 = fRow.DF2;
    else
        warning('extractSimpleEffects:TermNotFound', '%s', ...
            'Two-way filterType:regressionType term not found in anova table.');
        [intF, intP, intDF1, intDF2] = deal(NaN, NaN, NaN, NaN);
    end
catch ME
    warning('extractSimpleEffects:ANOVAFailed', '%s', ME.message);
    anovaTable = table();
    [intF, intP, intDF1, intDF2] = deal(NaN, NaN, NaN, NaN);
end

interactionTest.F   = intF;
interactionTest.p   = intP;
interactionTest.df1 = intDF1;
interactionTest.df2 = intDF2;
interactionTest.sig = ~isnan(intP) && intP < 0.05;
fprintf('  F(%d,%d) = %.2f, p = %.2e, significant = %s\n\n', ...
    intDF1, intDF2, intF, intP, tf2str(interactionTest.sig));

%% S3: SIX §8.1 CRITICAL SIMPLE EFFECTS — GRAND-MEAN MARGINALISED
% Monte Carlo marginalisation over the empirical covariate distribution.
% For each pair of pipelines, override filterType/regressionType on the
% same sampled rows and predict. Elementwise (predA - predB) gives a
% within-observation simple effect; its mean and SE are the estimates.
fprintf('--- S3: Six §8.1 simple effects (grand-mean marginalised) ---\n');

N_MC   = min(10000, nObs);
idxMC  = sampleStratified(tbl, N_MC);
tblMC  = tbl(idxMC, :);

% Pipeline codes {filterType, regressionType}
PIPE = struct( ...
    'BWFD_OLS',  [2, 3], ...
    'BWFD_LMLS', [2, 4], ...
    'BWFD_IRLS', [2, 5], ...
    'SG_OLS',    [6, 3], ...
    'SG_LMLS',   [6, 4], ...
    'SG_IRLS',   [6, 5]);

fprintf('  Predicting %d observations x 6 pipelines...\n', N_MC);
preds  = computeAllPredictions(model, tblMC, PIPE);
obsSD  = computeObservedSD(tbl, PIPE);

% Six contrasts: {label, pipeA, pipeB, expectedSign, rationale}
%   expectedSign: +1 = A > B expected; -1 = A < B expected; 0 = no prediction
%   Note: deltaBeta = betaGenerated - betaRecovered, so a "better" pipeline
%   has deltaBeta closer to zero. Positive delta means pipeA further from truth.
contrasts = { ...
    'BWFD-OLS vs SG-OLS',    'BWFD_OLS',  'SG_OLS',   +1, 'SG >> BWFD with OLS (large expected advantage)'; ...
    'BWFD-LMLS vs SG-LMLS',  'BWFD_LMLS', 'SG_LMLS',  +1, 'SG > BWFD with LMLS (moderate expected advantage)'; ...
    'BWFD-IRLS vs SG-IRLS',  'BWFD_IRLS', 'SG_IRLS',  +1, 'SG >= BWFD with IRLS (small expected advantage)'; ...
    'BWFD: OLS vs IRLS',     'BWFD_OLS',  'BWFD_IRLS', +1, 'IRLS rescues BWFD artefacts (expected large)'; ...
    'SG: OLS vs IRLS',       'SG_OLS',    'SG_IRLS',   +1, 'IRLS advantage smaller with SG (expected small)'; ...
    'SG-OLS vs BWFD-IRLS',   'SG_OLS',    'BWFD_IRLS', 0,  'Compensation test: kinematic vs regression dominance' ...
    };

nContrasts = size(contrasts, 1);
grandSE(nContrasts) = struct();

fprintf('\n  %-28s  %8s  %12s  %6s  %8s\n', 'Comparison', 'delta', '95% CI', 'd', 'Flag');
fprintf('  %s\n', repmat('-', 1, 74));

for k = 1:nContrasts
    label       = contrasts{k, 1};
    nameA       = contrasts{k, 2};
    nameB       = contrasts{k, 3};
    expectedSgn = contrasts{k, 4};

    diff  = preds.(nameA) - preds.(nameB);
    delta = mean(diff);
    se    = std(diff) / sqrt(N_MC);
    ciLo  = delta - 1.96 * se;
    ciHi  = delta + 1.96 * se;
    dCoh  = delta / sqrt(0.5 * (obsSD.(nameA)^2 + obsSD.(nameB)^2));
    flag  = interpretFlag(delta, ciLo, ciHi, expectedSgn);

    grandSE(k).label       = label;
    grandSE(k).pipeA       = nameA;
    grandSE(k).pipeB       = nameB;
    grandSE(k).delta       = delta;
    grandSE(k).ciLo        = ciLo;
    grandSE(k).ciHi        = ciHi;
    grandSE(k).se          = se;
    grandSE(k).d           = dCoh;
    grandSE(k).expectedSgn = expectedSgn;
    grandSE(k).flag        = flag;
    grandSE(k).rationale   = contrasts{k, 5};

    fprintf('  %-28s  %+8.4f  [%+6.4f %+6.4f]  %+5.2f  %s\n', ...
        label, delta, ciLo, ciHi, dCoh, flag);
end
fprintf('\n');

%% S4: NOISE-CONDITIONAL SIMPLE EFFECTS AT EMPIRICAL DATASET CENTROIDS
% Six contrasts repeated at the four empirical noise centroids (Finding #6
% median values). Nearest-grid-point snap on (noiseColor, noiseMagnitude).
fprintf('--- S4: Noise-conditional simple effects ---\n\n');

% [alpha_median, sigma_median_mm, dataset_label] — Finding #6 median table
CENTROIDS = { ...
    2.531, 4.440, 'Zarandi'; ...
    3.649, 6.470, 'Cook CTRL'; ...
    4.600, 5.677, 'Hickman PLAC'; ...
    2.867, 1.490, 'Dagenais' ...
    };
nCentroids = size(CENTROIDS, 1);

alphaGrid = unique(tbl.noiseColor);      % z-scored values (if v2_001 mat)
sigmaGrid = unique(tbl.noiseMagnitude);  % z-scored values (if v2_001 mat)

% Pre-extract scaling for centroid z-scoring (empty if v001 mat)
if isfield(s.results, 'predictor_scaling')
    sc_nc = s.results.predictor_scaling.noiseColor;
    sc_nm = s.results.predictor_scaling.noiseMagnitude;
    centroidsAreScaled = true;
    fprintf('  Centroid snap: empirical (alpha, sigma) will be z-scored before grid match.\n\n');
else
    centroidsAreScaled = false;
    fprintf('  Centroid snap: original scale (v001 mat).\n\n');
end

noiseSE(nCentroids * nContrasts) = struct();
entry = 0;

fprintf('  %-14s  %-28s  %8s  %12s  %6s  %5s\n', ...
    'Dataset', 'Comparison', 'delta', '95% CI', 'd', 'Flag');
fprintf('  %s\n', repmat('-', 1, 84));

for c = 1:nCentroids
    alphaEmp = CENTROIDS{c, 1};
    sigmaEmp = CENTROIDS{c, 2};
    dsLabel  = CENTROIDS{c, 3};

    % Convert empirical centroid to z-scored space if needed
    if centroidsAreScaled
        alphaSnap_query = (alphaEmp - sc_nc.mean) / sc_nc.std;
        sigmaSnap_query = (sigmaEmp - sc_nm.mean) / sc_nm.std;
    else
        alphaSnap_query = alphaEmp;
        sigmaSnap_query = sigmaEmp;
    end

    % Nearest-grid-point snap (in z-scored space)
    [~, iA] = min(abs(alphaGrid - alphaSnap_query));
    [~, iS] = min(abs(sigmaGrid - sigmaSnap_query));
    alphaSnap = alphaGrid(iA);
    sigmaSnap = sigmaGrid(iS);

    % Snap distance reported in original units for readability
    if centroidsAreScaled
        alphaSnap_orig = alphaSnap * sc_nc.std + sc_nc.mean;
        sigmaSnap_orig = sigmaSnap * sc_nm.std + sc_nm.mean;
    else
        alphaSnap_orig = alphaSnap;
        sigmaSnap_orig = sigmaSnap;
    end
    snapDist = sqrt((alphaEmp - alphaSnap_orig)^2 + (sigmaEmp - sigmaSnap_orig)^2);

    mask   = (tbl.noiseColor == alphaSnap) & (tbl.noiseMagnitude == sigmaSnap);
    nCell  = sum(mask);
    if nCell < 30
        warning('extractSimpleEffects:ThinCell', '%s', ...
            sprintf('%s: only %d rows at snapped centroid — estimates unreliable.', ...
            dsLabel, nCell));
    end
    tblCell   = tbl(mask, :);
    predsCell = computeAllPredictions(model, tblCell, PIPE);
    sdCell    = computeObservedSD(tblCell, PIPE);

    fprintf('  %s  (alpha %.2f->%.2f, sigma %.2f->%.2f mm, snap=%.3f, n=%d)\n', ...
        dsLabel, alphaEmp, alphaSnap_orig, sigmaEmp, sigmaSnap_orig, snapDist, nCell);

    for k = 1:nContrasts
        nameA = contrasts{k, 2};
        nameB = contrasts{k, 3};
        diff  = predsCell.(nameA) - predsCell.(nameB);
        delta = mean(diff);
        se    = std(diff) / max(sqrt(nCell), 1);
        ciLo  = delta - 1.96 * se;
        ciHi  = delta + 1.96 * se;
        sdA   = sdCell.(nameA);
        sdB   = sdCell.(nameB);
        dCoh  = delta / sqrt(0.5 * (sdA^2 + sdB^2));
        flag  = interpretFlag(delta, ciLo, ciHi, contrasts{k, 4});

        entry = entry + 1;
        noiseSE(entry).dataset   = dsLabel;
        noiseSE(entry).alphaEmp  = alphaEmp;
        noiseSE(entry).sigmaEmp  = sigmaEmp;
        noiseSE(entry).alphaSnap = alphaSnap_orig;
        noiseSE(entry).sigmaSnap = sigmaSnap_orig;
        noiseSE(entry).snapDist  = snapDist;
        noiseSE(entry).nCell     = nCell;
        noiseSE(entry).label     = contrasts{k, 1};
        noiseSE(entry).pipeA     = nameA;
        noiseSE(entry).pipeB     = nameB;
        noiseSE(entry).delta     = delta;
        noiseSE(entry).ciLo      = ciLo;
        noiseSE(entry).ciHi      = ciHi;
        noiseSE(entry).d         = dCoh;
        noiseSE(entry).flag      = flag;

        fprintf('  %-14s  %-28s  %+8.4f  [%+6.4f %+6.4f]  %+5.2f  %s\n', ...
            dsLabel, contrasts{k, 1}, delta, ciLo, ciHi, dCoh, flag);
    end
    fprintf('\n');
end
noiseSE = noiseSE(1:entry);   % trim pre-allocation

%% S5: WITHIN-FAMILY PIPELINE RANKINGS
% Marginal means per pipeline and ordinal ranking within BWFD and SG.
fprintf('--- S5: Within-family pipeline rankings ---\n');

piNames   = fieldnames(PIPE);
margMeans = struct();
for p = 1:numel(piNames)
    margMeans.(piNames{p}) = mean(preds.(piNames{p}));
    fprintf('  %-12s  mean(deltaBeta) = %+.4f\n', piNames{p}, margMeans.(piNames{p}));
end

fprintf('\n  BWFD ranking (smaller |deltaBeta| = less bias):\n');
bwfdNames = {'BWFD_OLS','BWFD_LMLS','BWFD_IRLS'};
bwfdM = cellfun(@(n) abs(margMeans.(n)), bwfdNames);
[~, bOrd] = sort(bwfdM);
for r = 1:3
    fprintf('    %d. %s  (|delta|=%.4f)\n', r, bwfdNames{bOrd(r)}, bwfdM(bOrd(r)));
end

fprintf('\n  SG ranking:\n');
sgNames = {'SG_OLS','SG_LMLS','SG_IRLS'};
sgM = cellfun(@(n) abs(margMeans.(n)), sgNames);
[~, sOrd] = sort(sgM);
for r = 1:3
    fprintf('    %d. %s  (|delta|=%.4f)\n', r, sgNames{sOrd(r)}, sgM(sOrd(r)));
end
fprintf('\n');

%% S6: SAVE AND WRITE PLAIN-TEXT SUMMARY
fprintf('--- S6: Saving results ---\n');

results.metadata.matFile       = matFile;
results.metadata.tractLvl      = tractLvl;
results.metadata.nObs          = nObs;
results.metadata.r2            = r2;
results.metadata.nCoef         = nCoef;
results.metadata.N_MC          = N_MC;
results.metadata.timestamp     = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<DATST>
results.metadata.preregSection = 'v101 §8.1';
results.metadata.note          = ['deltaBeta = betaGenerated - betaRecovered. ' ...
                                  'Positive delta means pipeline over-estimates bias.'];
results.interactionTest        = interactionTest;
results.grandSE                = grandSE;
results.noiseSE                = noiseSE;
results.marginalMeans          = margMeans;
results.pipelineCoding         = PIPE;
results.empiricalCentroids     = CENTROIDS;

timestamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<DATST>
outMat    = fullfile(srcDir, sprintf('simpleEffects_L%d_%s.mat', tractLvl, timestamp));
save(outMat, 'results', '-v7');
fprintf('  Mat:  %s\n', outMat);

% Plain-text summary
resDir = fullfile(fileparts(srcDir), 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
outTxt = fullfile(resDir, sprintf('simpleEffects_L%d_summary.txt', tractLvl));
writeSummary(outTxt, results);
fprintf('  Text: %s\n', outTxt);

fprintf('\nextractSimpleEffects_v003 COMPLETE.\n\n');
diary off

%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================

function mask = pipeMask(col, val)
% Compare a nominal column against a numeric value by converting via string.
    mask = strcmp(cellstr(char(col)), num2str(val));
end

function preds = computeAllPredictions(model, tblIn, PIPE)
% Predict deltaBeta under each pipeline for all rows in tblIn.
% Returns struct with one field per pipeline (column vectors).
    names = fieldnames(PIPE);
    preds = struct();
    for p = 1:numel(names)
        code           = PIPE.(names{p});
        tblTmp         = tblIn;
        tblTmp.filterType(:)     = nominal(code(1));
        tblTmp.regressionType(:) = nominal(code(2));
        preds.(names{p}) = predict(model, tblTmp, 'Conditional', false);
    end
end

function sdMap = computeObservedSD(tblIn, PIPE)
% Observed SD of deltaBeta per pipeline (for Cohen's d denominator).
    names = fieldnames(PIPE);
    sdMap = struct();
    for p = 1:numel(names)
        code = PIPE.(names{p});
        mask = pipeMask(tblIn.filterType, code(1)) & pipeMask(tblIn.regressionType, code(2));
        vals = tblIn.deltaBeta(mask);
        sdMap.(names{p}) = ternary(numel(vals) >= 2, std(vals), NaN);
    end
end

function idx = sampleStratified(tbl, N_target)
% Sample up to N_target rows stratified over the 6 pipelines × 4 sigma bins.
    nPerStratum = floor(N_target / 24);
    sigmaEdges  = [0, 0.5, 2, 5, Inf];
    idx = [];
    for ft = [2, 6]
        for rt = [3, 4, 5]
            for b = 1:4
                mask = pipeMask(tbl.filterType, ft) & pipeMask(tbl.regressionType, rt) & ...
                       (tbl.noiseMagnitude >= sigmaEdges(b)) & ...
                       (tbl.noiseMagnitude <  sigmaEdges(b+1));
                cand = find(mask);
                if isempty(cand), continue; end
                take = cand(randperm(numel(cand), min(nPerStratum, numel(cand))));
                idx  = [idx; take(:)]; %#ok<AGROW>
            end
        end
    end
    % Top up to N_target
    if numel(idx) < N_target
        remaining = setdiff((1:height(tbl))', idx);
        extra = remaining(randperm(numel(remaining), ...
                    min(N_target - numel(idx), numel(remaining))));
        idx = [idx; extra(:)];
    end
end

function flag = interpretFlag(delta, ciLo, ciHi, expectedSgn)
% Classify the simple effect against the §8.1 directional prediction.
%   expectedSgn: +1 (A > B), -1 (A < B), 0 (no prediction)
    if ciLo <= 0 && ciHi >= 0
        flag = 'NULL';     return;
    end
    if expectedSgn == 0
        flag = 'NO_PRED';  return;
    end
    if sign(delta) == expectedSgn
        flag = 'EXPECTED';
    else
        flag = 'UNEXPECTED';
    end
end

function v = ternary(cond, a, b)
% Inline ternary: v = cond ? a : b
    if cond, v = a; else, v = b; end
end

function str = tf2str(val)
    if val, str = 'YES'; else, str = 'NO'; end
end

function writeSummary(outFile, res)
% Write plain-text summary suitable for pasting into supplementary materials.
    fid = fopen(outFile, 'w');
    if fid < 0
        warning('extractSimpleEffects:WriteError', '%s', ...
            ['Cannot open for writing: ' outFile]);
        return;
    end
    fprintf(fid, 'Simple Effects Summary — prereg v101 §8.1\n');
    fprintf(fid, 'Generated : %s\n', res.metadata.timestamp);
    fprintf(fid, 'Source    : L%d  N=%d  R2=%.4f  nCoef=%d\n\n', ...
        res.metadata.tractLvl, res.metadata.nObs, ...
        res.metadata.r2, res.metadata.nCoef);
    fprintf(fid, 'Note: %s\n\n', res.metadata.note);

    fprintf(fid, 'INTERACTION TEST\n');
    fprintf(fid, '  filterType x regressionType: F(%d,%d) = %.2f, p = %.2e, sig = %s\n\n', ...
        res.interactionTest.df1, res.interactionTest.df2, ...
        res.interactionTest.F,   res.interactionTest.p, ...
        tf2str(res.interactionTest.sig));

    fprintf(fid, 'GRAND-MEAN SIMPLE EFFECTS  (N_MC = %d, Monte Carlo marginalisation)\n', ...
        res.metadata.N_MC);
    fprintf(fid, '%-28s  %8s  %12s  %6s  %8s  Rationale\n', ...
        'Comparison', 'delta', '95% CI', 'd', 'Flag');
    fprintf(fid, '%s\n', repmat('-', 1, 100));
    for k = 1:numel(res.grandSE)
        e = res.grandSE(k);
        fprintf(fid, '%-28s  %+8.4f  [%+6.4f %+6.4f]  %+5.2f  %-10s  %s\n', ...
            e.label, e.delta, e.ciLo, e.ciHi, e.d, e.flag, e.rationale);
    end

    fprintf(fid, '\nNOISE-CONDITIONAL SIMPLE EFFECTS\n');
    fprintf(fid, '  Snapped to nearest (alpha, sigma) grid point.\n\n');
    curDs = '';
    for k = 1:numel(res.noiseSE)
        e = res.noiseSE(k);
        if ~strcmp(e.dataset, curDs)
            curDs = e.dataset;
            fprintf(fid, '  %s  alpha %.2f->%.2f  sigma %.2f->%.2f mm  snap=%.3f  n=%d\n', ...
                e.dataset, e.alphaEmp, e.alphaSnap, e.sigmaEmp, e.sigmaSnap, ...
                e.snapDist, e.nCell);
        end
        fprintf(fid, '    %-28s  %+8.4f  [%+6.4f %+6.4f]  %+5.2f  %s\n', ...
            e.label, e.delta, e.ciLo, e.ciHi, e.d, e.flag);
    end

    fprintf(fid, '\nPIPELINE CODING\n');
    fprintf(fid, '  filterType    : 2 = BWFD (reference), 6 = SG\n');
    fprintf(fid, '  regressionType: 3 = OLS (reference),  4 = LMLS, 5 = IRLS\n');
    fprintf(fid, '  Treatment (dummy) coding; BWFD-OLS is the model baseline.\n');
    fprintf(fid, '\nFLAG LEGEND\n');
    fprintf(fid, '  EXPECTED   : sign matches §8.1 prediction, CI excludes zero\n');
    fprintf(fid, '  UNEXPECTED : sign opposes §8.1 prediction, CI excludes zero\n');
    fprintf(fid, '  NULL       : 95%% CI spans zero\n');
    fprintf(fid, '  NO_PRED    : no directional prediction registered in §8.1\n');
    fclose(fid);
end
