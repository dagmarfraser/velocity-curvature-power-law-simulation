function results = runLoopClosureFftnoise_v028(datasetName, varargin)
% runLoopClosureFftnoise_v028  v027 + SegmentRule option (auto | tol | budget | detect; endToEndSegment_v006). GENERATED: see makeRunLoopClosureV028_v001.m. Edit the generator, not this file.
% runLoopClosureFftnoise_v028  v026 + NoiseScale option (end-matched noise variance multiplier, for the budget sensitivity). GENERATED: see makeRunLoopClosureV027_v001.m. Edit the generator, not this file.
% runLoopClosureFftnoise_v028  v025 + v024's window bug fixed + SegmentBudget option (endToEndSegment_v005). GENERATED: see makeRunLoopClosureV026_v001.m. Edit the generator, not this file.
% runLoopClosureFftnoise_v028  v024 + NHarmonics option (template harmonic count; integer or 'fc'). GENERATED: see makeRunLoopClosureV025_v001.m. Edit the generator, not this file.
% runLoopClosureFftnoise_v028  v023 + the production template's uncovered (zero) tail excluded from the end-matched segment search. GENERATED: see makeRunLoopClosureV024_v001.m. Edit the generator, not this file.
% runLoopClosureFftnoise_v028  v022 + NoiseModel 'shapedXuSeg' (the shapeXu sensitivity arm on the same segments, length, windows and seeds as 'endMatched'). GENERATED: see makeRunLoopClosureV023_v001.m. Edit the generator, not this file.
% runLoopClosureFftnoise_v028  v021 + NoiseModel 'endMatched' (fairShare common-length segments, shortened synthetic, beta_obs on a same-length central window with the whole trial as sensitivity). GENERATED: see makeRunLoopClosureV022_v001.m. Edit the generator, not this file.
% runLoopClosureFftnoise_v028  v020 + CarrierAlpha option (irasa | pwlw | irasaRecordBand) for the forward-map noise carrier. GENERATED: see makeRunLoopClosureV021_v001.m. Edit the generator, not this file.
% runLoopClosureFftnoise_v028  v019 + Template option (endAnchored | production), dataset FraserN48 (importDB_fraser_v002, 5.209 pt/mm), counted EBR failures. GENERATED: see makeRunLoopClosureV020_v001.m. Edit the generator, not this file.
% runLoopClosureFftnoise_v019  v017 + END-ANCHORED template (one final window ending at the last sample). GENERATED: see makeRunLoopClosureV019_v001.m. Edit the generator, not this file.
% runLoopClosureFftnoise_v017  v016 + checkpointed chunks, a wall-clock deadline and an atomic save.
%
% GENERATED, not hand-written: see makeRunLoopClosureV017_v001.m. Edit the generator, not this file.
%
% Options 'ChunkSize' (0), 'Deadline' (NaT), 'OutTag' (''). With ChunkSize > 0 the trials run in chunks, each saved
% to <output>_chunks/ as it finishes; a restart loads them and continues. Each trial seeds its own stream
% (rng(RngSeed + index)) so a resumed run equals an uninterrupted one. Deadline pauses cleanly between chunks.
% The final save is atomic. At ChunkSize 0 the numerics are v016's. v016's own header follows unchanged.
% runLoopClosureFftnoise_v016  v015 + options for sampling rate and filters (D-D16 stage S2).
%
% GENERATED, not hand-written: see makeRunLoopClosureV016_v001.m. Edit the generator, not this file.
%
% THE ONLY CHANGES FROM v015: options 'Subsample' (1), 'BWCutoff' (10 Hz) and 'SGRef' (17), beside v015's
% 'EdgeClip' (20). Clean subsampling (x(1:Subsample:end), no filter) is applied to the selected trials and
% everything then runs at FS_RUN = FS_DS / Subsample. Non-default options are added to the output and log names and all
% are saved in config. r.alphaIRA and r.sigma still come from the 240 Hz noise file (labels only); the
% forward map re-estimates its noise at FS_RUN, as v015 does at FS_DS.
% ACCEPTANCE TEST: at the defaults v016 is v015 (runD16HPC_v001 checks it on Zarandi first).
% See docs/SPEC_D16_SamplingRateCutoff_v001.md. v015's own header follows unchanged.
% runLoopClosureFftnoise_v015  v013 + persist the per-trial forward map.
%
% GENERATED, not hand-written: see makeRunLoopClosureV015_v001.m. Edit the
% generator, not this file, or the next regeneration will overwrite you.
%
% THE ONLY CHANGE FROM v013 is in packResult_local: two extra fields,
%   betaRecCurveMed  [nPipelines x N_BETA]  median over reps
%   betaRecCurveSD   [nPipelines x N_BETA]  std over reps
% These are the exact two reductions invertBeta_local already computes
% (v013 lines 582-583) and then throws away. v013 persisted only
% betaRecSlice, a [N_REPS x 6] replicate cloud at ONE grid point, which
% is not a forward map and cannot predict a held-out pipeline's betaObs.
%
% WHY: the strong, out-of-sample form of leave-one-pipeline-out -- infer
% beta_gen* from five pipelines, predict the sixth's betaObs through its
% own forward map -- is the residue of GPT SOL's circularity critique
% (section 2) and the manuscript's own most substantive open item. The
% scoped form (Finding #199) compares inverted values only and cannot
% answer it.
%
% ACCEPTANCE TEST: v015 changes no numerics, so betaGenStar, invertStatus,
% ciLo/ciHi and betaObs must match v013's EXACTLY on any dataset run with
% the same seed and config. runLoopClosureAll_v015.m checks this and
% refuses to proceed if it fails.
%
% v013's own header follows unchanged.
% runLoopClosureFftnoise_v013  Loop-closure inversion -- v012 + one CI bug fix.
%
% WHY THIS VERSION EXISTS (2026-08-16): v012's own invertBeta_local had a
% real bug in its CI computation, found via a downstream convergence-
% check script's own output, not by inspection. The CI is built from two
% candidate bounds (ciA from the 5th percentile curve, ciB from the 95th),
% each of which correctly returns NaN when the observed beta falls
% outside what that particular percentile curve spans -- a genuine "this
% side of the interval cannot be confidently bounded" case, not an error.
% But the next line, `ciLo = min(ciA, ciB); ciHi = max(ciA, ciB)`, relies
% on MATLAB's two-argument max/min, which SILENTLY IGNORE NaN when
% comparing two scalars (max(NaN, 5) returns 5, not NaN). So whenever
% either bound was genuinely undefined, the code silently collapsed to
% the single defined value for BOTH ciLo and ciHi, reporting a false,
% deceptively precise ZERO-WIDTH confidence interval instead of correctly
% reporting the interval as undefined.
%
% CAUGHT how, concretely: `checkNRepsConvergence_HPC_v001.m` (a
% convergence-check script built the same session, testing whether
% N_REPS=20 is empirically sufficient) duplicated this exact CI logic
% from v012 (matching this project's own duplication-over-refactor
% convention, since v012's own invertBeta_local is a private local
% function). Running it across all 7 canonical datasets, Cook ASD came
% back with ciWidthMean=ciWidthSD=0.0000 across all 50 independent
% subsamples, at THREE separate subset sizes (40, 60, 100) simultaneously
% -- not a plausible coincidence of genuine convergence, a clear signature
% of every single subsample independently hitting the same silent-NaN-
% collapse path. Confirmed as this project's own long-standing pattern,
% not confined to the new script: v012's invertBeta_local (below,
% unchanged except for this ONE fix) and the standalone tool's own
% `processTrialLoopClosure_v001.m` (fixed the same session, see
% docs/TODO_ToolBuild_v001.md) both carry the identical bug, inherited
% from a common lineage.
%
% THE FIX, the ONE change from v012, isolated to invertBeta_local: if
% either ciA or ciB is NaN, both ciLo and ciHi are now explicitly set to
% NaN, rather than letting max/min silently paper over the missing bound.
% Everything else in this file, including betaGenStar itself, invertStatus
% classification, and the gate mechanism (findBothBranches_v008), is
% IDENTICAL to v012 -- this bug affects CI width reporting only, not
% point estimates or PASS/CONDITIONAL/FAIL-style classification (which
% never used CI width in the first place, see
% checkInvertibilityForEmpirical_v2_002.m's own rise-coverage-only gate).
%
% RE-RUN IMPLICATIONS, stated plainly rather than assumed: any existing
% v012 output whose downstream use depends on a specific CI width or
% range value (not just the point estimate or invertStatus) may be
% understated wherever this exact edge case was hit -- how often that
% happened across the existing 7-dataset v012 corpus is NOT yet known,
% this file's own existence does not retroactively audit that. v012's
% own existing .mat outputs are left untouched (this project's own
% version-lineage-immutability convention); v013 writes its own,
% separately-named outputs (loopClosureResults_..._v013.mat), so nothing
% about v012's existing citable results is silently altered by this file
% existing.
%
% v012's own original header, preserved below as historical record
% (unchanged, still accurate as a description of v012's own work):
%
% runLoopClosureFftnoise_v012  Loop-closure inversion -- R11: noise-
% tolerant local-monotonicity gate (findBothBranches_v008), fixing the
% grid-density sensitivity found and fully resolved this session
% (Findings #149-#158), plus v011's local-monotonicity structure capture.
%
% WHY THIS VERSION EXISTS (2026-08-11): v010/v011's gate (findBothBranches_v002)
% classified per-trial invertibility correctly given a fixed simulation
% grid, but was highly sensitive to that grid's own density (N_BETA) --
% 72.2% of a real 36-cell benchmark flipped classification purely from
% resampling the same curve at different resolutions (Finding #149).
% Root cause, traced to a specific step: the strict-monotonicity TRIM
% (walking forward on the raw curve, hard-stopping the instant ANY
% decrease occurred, however tiny) had zero noise tolerance (Finding
% #153) -- a genuine but noise-scale reversal at fine resolution could
% truncate an otherwise-robust segment by a razor-thin margin, while a
% coarser grid's sample points simply skipped over the same location.
%
% Two other candidate mechanisms were tested and explicitly ruled out
% before finding the real one: width-based smoothing/segment-floor fixes
% (findBothBranches_v003/v004) helped substantially (72.2%->44.4%,
% Finding #151) but did not fully resolve it; noise-calibrated
% slope-DETECTION thresholds (findBothBranches_v005/v006/v007) were
% actively counterproductive -- any finite-difference slope's estimation-
% uncertainty formula scales as 1/dx by construction, making detection
% MORE permissive at coarse resolution, the opposite of what is needed
% (Finding #152/#153, v007 reached 86.1%, worse than the original
% problem). A literature review (Finding #144/#154) and cross-validation
% against R's MonotonicityTest package (Hall & Heckman 2000, Finding
% #157) grounded the eventual fix, which turned out to need only a
% MINIMAL change once correctly targeted (Finding #158): apply a noise-
% calibrated TOLERANCE to the trim step specifically (a decrease only
% ends a segment if it exceeds k*sigma_y_local, not if merely negative),
% leaving detection (width-based smoothing/segment-floor, flat slopeTol)
% completely unchanged. Result: 0/36 flip (0.0%) on the same benchmark,
% confirmed non-degenerate (real coverage, not a collapse) and
% independently validated against R.
%
% v010/v011 are left unmodified -- their existing outputs (the confirmed-
% unchanged 40/2/0 headline table, Finding #147) remain valid at their
% own resolution. v012 (and now v013) supersede both for any per-trial
% classification work going forward -- v002/v001's gates are known-
% sensitive to grid density in a way v008 is not.
%
% NEW PARAMETER vs v011: MonoTrimToleranceK (default 1.0) -- the noise-
% tolerance multiplier for the trim step. MonoSmoothWindow/MonoMinSegLength
% are RENAMED MonoSmoothWidth/MonoMinSegWidth (beta_gen-unit widths, not
% grid-point counts, matching v003/v004's own already-established fix).
%
% Fraser, D.S. (2026)  v013
% See also: runLoopClosureFftnoise_v012 (identical except the CI fix
%           above), runLoopClosureFftnoise_v011 (local-monotonicity
%           capture, unchanged here), findBothBranches_v008 (the fix),
%           docs/FINDINGS_REFERENCE.md Findings #149-#158 (the full
%           investigation), docs/TODO_NoiseCalibratedInvertibilityGate_v001.md

%% ---------- Parse inputs --------------------------------------------------
p = inputParser;
addRequired(p,  'datasetName',    @ischar);
addParameter(p, 'TrialSelection', 'per_subject', ...
    @(x) ismember(x, {'per_subject','all'}));
addParameter(p, 'TrialID',        '',    @ischar);
addParameter(p, 'UseParfor',      false, @islogical);
addParameter(p, 'N_BETA',         25,    @isscalar);
addParameter(p, 'BetaRange',      [0 0.75], @(x) numel(x)==2);
addParameter(p, 'N_REPS',         20,    @isscalar);
addParameter(p, 'N_R7',           10,    @isscalar);
addParameter(p, 'RngSeed',        1729,  @isscalar);
addParameter(p, 'EdgeClip',       20,    @isscalar);
addParameter(p, 'Subsample',      1,     @(x) isscalar(x) && x >= 1 && x == round(x));   % v016: clean stride, no filter
addParameter(p, 'BWCutoff',       10,    @(x) isscalar(x) && x > 0);   % v016: Butterworth cutoff, Hz
addParameter(p, 'SGRef',          17,    @(x) isscalar(x) && x > 0);   % v016: SG reference frame at 100 Hz
addParameter(p, 'ChunkSize',      0,     @(x) isscalar(x) && x >= 0 && x == round(x));   % v017: trials per checkpoint chunk, 0 = off
addParameter(p, 'Deadline',       NaT,   @(x) isa(x, 'datetime') && isscalar(x));   % v017: pause between chunks if the next would pass it
addParameter(p, 'Template',       'endAnchored', @(x) ismember(x, {'endAnchored','production'}));   % v020
addParameter(p, 'SegmentRule',    'auto', @(x) ismember(string(x), ["auto","tol","budget","detect"]));   % v028
addParameter(p, 'NoiseScale',     1, @(x) isnumeric(x) && isscalar(x) && x > 0);   % v027
addParameter(p, 'SegmentBudget',  [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0 && x < 1));   % v026
addParameter(p, 'NHarmonics',     4, @(x) (isnumeric(x) && isscalar(x) && x >= 1 && x == round(x)) || strcmp(string(x), "fc"));   % v025
addParameter(p, 'CarrierAlpha',   'irasa', @(x) ismember(x, {'irasa','pwlw','irasaRecordBand'}));   % v021
addParameter(p, 'OutTag',         '',    @(x) ischar(x) || isstring(x));   % v017: smoke tests only
addParameter(p, 'MaxChunks',      Inf,   @(x) isscalar(x) && x >= 0);   % v017: pause after this many NEW chunks (smoke tests)
addParameter(p, 'ForceRerun',     false, @islogical);
addParameter(p, 'NoiseModel',     'shaped_xu', ...
    @(x) ismember(x, {'auto','fftnoise','xu','shaped_xu','bootstrap','endMatched','shapedXuSeg'}));   % v022: endMatched; v023: shapedXuSeg
addParameter(p, 'MonoSlopeTol',       0.05,   @isscalar);
addParameter(p, 'MonoSmoothWidth',    0.0625, @isscalar);
addParameter(p, 'MonoMinSegWidth',    0.0625, @isscalar);
addParameter(p, 'MonoTrimToleranceK', 1.0,    @isscalar);
parse(p, datasetName, varargin{:});
opt = p.Results;

DATASET         = opt.datasetName;
TRIAL_SELECTION = opt.TrialSelection;
TRIAL_ID_OVR    = opt.TrialID;
USE_PARFOR      = opt.UseParfor;
N_BETA          = opt.N_BETA;
BETA_RANGE      = opt.BetaRange;
N_REPS          = opt.N_REPS;
N_R7            = opt.N_R7;
RNG_SEED        = opt.RngSeed;
EDGE_CLIP       = opt.EdgeClip;
SUBSAMPLE       = opt.Subsample;  % v016
BW_CUTOFF       = opt.BWCutoff;   % v016
SG_REF          = opt.SGRef;      % v016
CHUNK           = opt.ChunkSize;  % v017
DEADLINE        = opt.Deadline;   % v017
MAXCHUNKS       = opt.MaxChunks;  % v017
NOISE_MODEL_OPT = opt.NoiseModel;

monoParams = struct('slopeTol', opt.MonoSlopeTol, ...
    'smoothWidth', opt.MonoSmoothWidth, 'minSegWidth', opt.MonoMinSegWidth, ...
    'trimToleranceK', opt.MonoTrimToleranceK);

%% ---------- Parfor guard --------------------------------------------------
if USE_PARFOR && isempty(gcp('nocreate'))
    isHPC = ~isempty(getenv('SLURM_JOB_ID')) || feature('numcores') > 16;
    if isHPC
        nW = feature('numcores'); pType = 'Processes';
    else
        nW = min(max(1, feature('numcores')-1), 16); pType = 'Processes';
    end
    fprintf('No pool -- creating %s pool with %d workers...\n', pType, nW);
    parpool(pType, nW);
end

%% ---------- Paths ---------------------------------------------------------
srcDir = fileparts(mfilename('fullpath'));
addpath(srcDir);
addpath(genpath(fullfile(srcDir, 'functions')));
addpath(genpath(fullfile(srcDir, 'req')));

%% ---------- Pipeline config -----------------------------------------------
derivConfigs   = struct('label',  {"BWFD","SG"}, ...
                        'filterType',   {2, 6}, ...
                        'filterParams', {[2 BW_CUTOFF 1],[4 SG_REF]});   % v016: v015 literals were [2 10 1],[4 17]
regressConfigs = struct('label', {"OLS","LMLS","IRLS"}, 'type', {3,4,5});
pipelineLabels = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];
nDeriv         = numel(derivConfigs);
nRegress       = numel(regressConfigs);
nPipelines     = nDeriv * nRegress;
betaGenVec     = linspace(BETA_RANGE(1), BETA_RANGE(2), N_BETA);

%% ---------- Dataset config ------------------------------------------------
switch DATASET
    case 'Pilot'
        matFile  = fullfile(srcDir, 'noiseCharacterisation_pilot.mat');
        sigToMM  = 1 / 9.73;
        FS_DS    = 240;
        importFn = @() importDB_pilot_v001(Verbose=false);
    case 'Cook CTRL'
        matFile  = fullfile(srcDir, 'noiseCharacterisation_cook.mat');
        sigToMM  = 0.248;  FS_DS = 133;
        importFn = @() importDB_cook_v002(Group="CTRL", Tasks=7, Verbose=false);
    case 'Cook ASD'
        matFile  = fullfile(srcDir, 'noiseCharacterisation_cookASD.mat');
        sigToMM  = 0.248;  FS_DS = 133;
        importFn = @() importDB_cook_v002(Group="ASD",  Tasks=7, Verbose=false);
    case 'Hickman PLAC'
        matFile  = fullfile(srcDir, 'noiseCharacterisation_hickmanPLAC.mat');
        sigToMM  = 0.248;  FS_DS = 133;
        importFn = @() importDB_hickman_v003(Study=2, Group="PLAC", Verbose=false);
    case 'Hickman HALO'
        matFile  = fullfile(srcDir, 'noiseCharacterisation_hickmanHALO.mat');
        sigToMM  = 0.248;  FS_DS = 133;
        importFn = @() importDB_hickman_v003(Study=2, Group="HALO", Verbose=false);
    case 'Zarandi'
        matFile  = fullfile(srcDir, 'noiseCharacterisation_zarandi.mat');
        sigToMM  = 10.0;   FS_DS = 100;
        importFn = @() importDB_zarandi_v001(Verbose=false);
    case 'Dhieb'
        matFile  = fullfile(srcDir, 'noiseCharacterisation_dhieb.mat');
        sigToMM  = 0.1478; FS_DS = 100;
        importFn = @() importDB_dhieb_v001(Verbose=false);
    case 'Fraser'
        matFile  = fullfile(srcDir, 'noiseCharacterisation_fraser.mat');
        sigToMM  = 1 / 5.208965;
        FS_DS    = 240;
        importFn = @() importDB_fraser_v001(Verbose=false);
    case 'FraserN48'   % v020
        matFile  = fullfile(srcDir, 'noiseCharacterisation_fraserN48.mat');
        sigToMM  = 1 / 5.209;   % UIKit points per mm, M4 13" iPad (Step 1 audit, 2026-11-03)
        FS_DS    = 240;
        importFn = @() importDB_fraser_v002(Verbose=false);
    otherwise
        error('runLoopClosureFftnoise_v028:UnknownDataset', '%s', ...
            ['Unknown DATASET "' DATASET '"']);
end

%% ---------- v016: run rate and option tag ---------------------------------
FS_RUN  = FS_DS / SUBSAMPLE;
OPT_TAG = "";
if SUBSAMPLE ~= 1,  OPT_TAG = OPT_TAG + "_sub" + SUBSAMPLE; end
if BW_CUTOFF ~= 10, OPT_TAG = OPT_TAG + "_bw" + strrep(string(BW_CUTOFF), ".", "p"); end
if SG_REF ~= 17,    OPT_TAG = OPT_TAG + "_sg" + strrep(string(SG_REF), ".", "p"); end
END_ANCHORED = strcmp(opt.Template, 'endAnchored');   % v020
if ~END_ANCHORED, OPT_TAG = OPT_TAG + "_tplProd"; end   % v020: v015's production windows
SEG_BUDGET = opt.SegmentBudget;  if ~isempty(SEG_BUDGET), OPT_TAG = OPT_TAG + "_bud" + string(round(100 * SEG_BUDGET, 4)); end   % v026
SEG_RULE = string(opt.SegmentRule);  if ~ismember(SEG_RULE, ["auto", "budget"]), OPT_TAG = OPT_TAG + "_seg" + SEG_RULE; end   % v028
NOISE_SCALE = opt.NoiseScale;  if NOISE_SCALE ~= 1, OPT_TAG = OPT_TAG + "_ns" + string(NOISE_SCALE); end   % v027
NH_OPT = opt.NHarmonics;  if ~(isnumeric(NH_OPT) && NH_OPT == 4), OPT_TAG = OPT_TAG + "_nH" + string(NH_OPT); end   % v025
CARRIER_ALPHA = char(opt.CarrierAlpha);   % v021
if strcmp(CARRIER_ALPHA, 'pwlw'), OPT_TAG = OPT_TAG + "_carPWLW"; elseif strcmp(CARRIER_ALPHA, 'irasaRecordBand'), OPT_TAG = OPT_TAG + "_carRB"; end   % v021
if EDGE_CLIP ~= 20, OPT_TAG = OPT_TAG + "_clip" + EDGE_CLIP; end
OPT_TAG = OPT_TAG + string(opt.OutTag);   % v017: smoke tests only

%% ---------- Load noise characterisation -----------------------------------
bio = load(matFile, 'bioResults').bioResults;
if ~ismember('ira_alphaMean', bio.Properties.VariableNames)
    error('runLoopClosureFftnoise_v028:NoIRASA', '%s', ...
        'bioResults lacks ira_alphaMean -- run recharacteriseBiologicalNoise_v001 first.');
end
alphaVec    = bio.ira_alphaMean;
sigmaVec_mm = bio.sigmaMean * sigToMM;
NOISE_MODEL = NOISE_MODEL_OPT;

%% ---------- Output filename -----------------------------------------------
tag = strrep(DATASET, ' ', '_');
if ~isempty(TRIAL_ID_OVR)
    outFile = fullfile(srcDir, sprintf('loopClosureResults_%s_%s_%s_v028.mat', ...
        tag, TRIAL_ID_OVR, NOISE_MODEL));
else
    outFile = fullfile(srcDir, sprintf('loopClosureResults_%s_%s_%s_v028.mat', ...
        tag, TRIAL_SELECTION, NOISE_MODEL));
end

outFile = char(replace(string(outFile), "_v028.mat", OPT_TAG + "_v028.mat"));   % v016: no-op at the defaults

if isfile(outFile) && ~opt.ForceRerun
    fprintf('[v028] Output exists -- skipping %s.\n', outFile);
    fprintf('  Pass ForceRerun=true to overwrite.\n');
    if nargout > 0
        S = load(outFile, 'results'); results = S.results;
    end
    return
end

rng(RNG_SEED, 'twister');

%% ---------- Import trials -------------------------------------------------
fprintf('=== runLoopClosureFftnoise_v028: %s | %s | model=%s | parfor=%d ===\n', ...
    DATASET, TRIAL_SELECTION, NOISE_MODEL, USE_PARFOR);
fprintf('  v016 options: subsample stride %d (run at %g Hz), BW cutoff %g Hz, SG reference frame %g, edge clip %d, tag "%s"\n', ...
    SUBSAMPLE, FS_RUN, BW_CUTOFF, SG_REF, EDGE_CLIP, OPT_TAG);
fprintf('Importing trials...\n');
trials   = importFn();
trialIDs = string({trials.trialID}');
fprintf('  %d trials imported.\n', numel(trials));

%% ---------- Trial selection -----------------------------------------------
if ~isempty(TRIAL_ID_OVR)
    tidMatch = find(string(bio.trialID) == string(TRIAL_ID_OVR), 1);
    if isempty(tidMatch)
        error('runLoopClosureFftnoise_v028:TrialNotFound', '%s', ...
            ['TrialID "' TRIAL_ID_OVR '" not found in bioResults.']);
    end
    selIdx = tidMatch;
    fprintf('  TrialID override: %s\n', TRIAL_ID_OVR);
else
    switch TRIAL_SELECTION
        case 'per_subject'
            subjIDs = unique(bio.subjectID);
            nSel_   = numel(subjIDs);
            selIdx  = zeros(nSel_, 1);
            for s = 1:nSel_
                mask  = bio.subjectID == subjIDs(s);
                aV    = alphaVec(mask);
                sV    = sigmaVec_mm(mask);
                valid = ~isnan(aV) & ~isnan(sV);
                if ~any(valid)
                    selIdx(s) = find(mask, 1);
                    continue
                end
                muA = median(aV(valid));
                muS = median(sV(valid));
                d   = (aV - muA).^2 / max(var(aV(valid)), eps) + ...
                      (sV - muS).^2 / max(var(sV(valid)), eps);
                d(~valid) = Inf;
                bioRows = find(mask);
                [~, best] = min(d);
                selIdx(s) = bioRows(best);
            end
            fprintf('  per_subject: %d subjects.\n', nSel_);
        case 'all'
            selIdx = (1:height(bio))';
            fprintf('  all: %d trials.\n', numel(selIdx));
    end
end

%% ---------- Match bioResults -> imported structs --------------------------
nSel        = numel(selIdx);
selTrialIDs = bio.trialID(selIdx);
matchIdx    = zeros(nSel, 1);
for s = 1:nSel
    m = find(trialIDs == selTrialIDs(s), 1);
    if isempty(m)
        warning('runLoopClosureFftnoise_v028:NoMatch', '%s', ...
            ['Trial ' char(selTrialIDs(s)) ' not found -- skipping.']);
    else
        matchIdx(s) = m;
    end
end
validMask = matchIdx > 0;
selIdx    = selIdx(validMask);
matchIdx  = matchIdx(validMask);
nSel      = numel(selIdx);
fprintf('  %d trials matched.\n', nSel);

alphaSelected = alphaVec(selIdx);
sigmaSelected = sigmaVec_mm(selIdx);

%% ---------- Slice trials for parfor (broadcast-variable safety, R9 pattern) -
trialsSel = trials(matchIdx);
for k = 1:numel(trialsSel), trialsSel(k).cellKey = selIdx(k); end   % v022: fixed identity for the endMatched substreams
clear trials
if SUBSAMPLE > 1   % v016: clean subsampling, every SUBSAMPLE-th sample from the first, no filter
    for k = 1:numel(trialsSel)
        trialsSel(k).x = trialsSel(k).x(1:SUBSAMPLE:end);
        trialsSel(k).y = trialsSel(k).y(1:SUBSAMPLE:end);
        if isfield(trialsSel, 't'),  trialsSel(k).t  = (0:numel(trialsSel(k).x) - 1)' / FS_RUN; end
        if isfield(trialsSel, 'fs'), trialsSel(k).fs = FS_RUN; end
    end
    fprintf('  v016: %d trials subsampled (stride %d, no filter) to %g Hz.\n', numel(trialsSel), SUBSAMPLE, FS_RUN);
end

%% ---------- Log setup -----------------------------------------------------
logsDir = fullfile(srcDir, 'logs');
if ~isfolder(logsDir), mkdir(logsDir); end
logFile = fullfile(logsDir, sprintf('loopClosureRun_%s_%s_%s.log', ...
    [tag char(OPT_TAG)], TRIAL_SELECTION, datestr(now, 'yyyymmdd_HHMMSS')));
fprintf('Log: %s\n\n', logFile);

%% ---------- Timing probe --------------------------------------------------
nTime  = min(2, nSel) * (CHUNK == 0);   % v017: chunk mode has no serial timing probe
tTimer = tic;
for k = 1:nTime
    processOneTrial_local(trialsSel(k), FS_RUN, sigToMM, betaGenVec, ...
        N_REPS, N_R7, EDGE_CLIP, derivConfigs, regressConfigs, ...
        nDeriv, nRegress, nPipelines, NOISE_MODEL, END_ANCHORED, CARRIER_ALPHA, RNG_SEED, NH_OPT, SEG_BUDGET, NOISE_SCALE, SEG_RULE);
end
secsPerTrial = toc(tTimer) / max(nTime, 1);
estMins      = secsPerTrial * nSel / 60;
if USE_PARFOR
    nWorkers   = gcp('nocreate').NumWorkers;
    fprintf('  %.1fs/trial -> %.0f min serial / ~%.0f min (%d workers)\n\n', ...
        secsPerTrial, estMins, estMins/nWorkers, nWorkers);
else
    nWorkers = 1;
    fprintf('  %.1fs/trial -> est. %.0f min total\n\n', secsPerTrial, estMins);
end

%% ---------- Main loop -----------------------------------------------------
resultsCells = cell(nSel, 1);
logStrings   = cell(nSel, 1);
tRun         = tic;

if CHUNK > 0   % v017: checkpointed chunks
    if ~USE_PARFOR, error('runLoopClosureFftnoise_v028:ChunkSerial', '%s', 'ChunkSize needs UseParfor=true.'); end
    ckDir = [outFile(1:end-4) '_chunks'];
    if opt.ForceRerun && isfolder(ckDir), rmdir(ckDir, 's'); end
    if ~isfolder(ckDir), mkdir(ckDir); end
    sig = struct('dataset', DATASET, 'model', NOISE_MODEL, 'nBeta', N_BETA, 'betaRange', BETA_RANGE, ...
        'nReps', N_REPS, 'nR7', N_R7, 'seed', RNG_SEED, 'edgeClip', EDGE_CLIP, 'sub', SUBSAMPLE, ...
        'bw', BW_CUTOFF, 'sgRef', SG_REF, 'chunk', CHUNK, 'nSel', nSel, 'trialIDs', {string({trialsSel.trialID})});
    nChunks   = ceil(nSel / CHUNK);
    chunkSecs = [];  newDone = 0;
    for c = 1:nChunks
        idx    = (c - 1) * CHUNK + 1 : min(c * CHUNK, nSel);
        ckFile = fullfile(ckDir, sprintf('chunk_%03d.mat', c));
        if isfile(ckFile)
            K = load(ckFile);   % a corrupt file errors here: fail loud
            if ~isequal(K.sig, sig) || ~isequal(K.idx, idx)
                error('runLoopClosureFftnoise_v028:ChunkMismatch', '%s', ...
                    ['Checkpoint ' ckFile ' was written with different options or trials. Move the _chunks folder away, or pass ForceRerun=true.']);
            end
            resultsCells(idx) = K.rc;  logStrings(idx) = K.lg;  chunkSecs(end + 1) = K.secs; %#ok<AGROW>
            fprintf('  chunk %d/%d (trials %d-%d): loaded from checkpoint (took %.1f min)\n', c, nChunks, idx(1), idx(end), K.secs / 60);
            continue
        end
        if newDone >= MAXCHUNKS || (~isnat(DEADLINE) && datetime('now') + seconds(1.15 * max([chunkSecs 0])) > DEADLINE)
            fprintf('\nPAUSED before chunk %d/%d (%s; deadline %s; MaxChunks %g). %d of %d chunks are saved in %s. Run again to continue.\n', ...
                c, nChunks, char(datetime('now')), char(DEADLINE), MAXCHUNKS, c - 1, nChunks, ckDir);
            results = [];
            return
        end
        pl = gcp('nocreate');
        if isempty(pl), error('runLoopClosureFftnoise_v028:NoPool', '%s', 'The parallel pool is gone; restart the launcher (finished chunks are kept).'); end
        nC = numel(idx);  rcC = cell(nC, 1);  lgC = cell(nC, 1);
        trialsC = trialsSel(idx);  aC = alphaSelected(idx);  sC = sigmaSelected(idx);
        t1 = tic;
        parfor j = 1:nC
            rng(RNG_SEED + idx(j), 'twister');   % v017: per-trial stream
            tr  = trialsC(j);
            res = processOneTrial_local(tr, FS_RUN, sigToMM, betaGenVec, ...
                N_REPS, N_R7, EDGE_CLIP, derivConfigs, regressConfigs, ...
                nDeriv, nRegress, nPipelines, NOISE_MODEL, END_ANCHORED, CARRIER_ALPHA, RNG_SEED, NH_OPT, SEG_BUDGET, NOISE_SCALE, SEG_RULE);
            [betaGenStar, ciLo, ciHi, betaGenStarMed, nearestIdx, invertInfo] = ...
                invertBeta_local(res, betaGenVec, nPipelines, monoParams); %#ok<PFBNS>
            betaRecSlice = extractSlice_local(res, nearestIdx);
            r = packResult_local(tr, res, betaRecSlice, betaGenStar, ...
                betaGenStarMed, ciLo, ciHi, aC(j), sC(j), NOISE_MODEL, invertInfo); %#ok<PFBNS>
            r = addWholeTrial_local(r, res, betaGenVec, nPipelines, monoParams, NOISE_MODEL);   % v022
            rcC{j} = r;
            lgC{j} = sprintf('%4d  %-28s  alpha=%6.3f  gen*_med=%+.4f  gen*_sg=%+.4f  status_sg=%s', ...
                idx(j), r.trialID, r.alphaIRA, betaGenStarMed, betaGenStar(6), r.invertStatus(6));
        end
        secs = toc(t1);  rc = rcC;  lg = lgC;
        tmpFile = fullfile(ckDir, sprintf('chunk_%03d_tmp.mat', c));
        save(tmpFile, 'rc', 'lg', 'idx', 'sig', 'secs', '-v7.3');
        movefile(tmpFile, ckFile);
        resultsCells(idx) = rcC;  logStrings(idx) = lgC;  chunkSecs(end + 1) = secs;  newDone = newDone + 1; %#ok<AGROW>
        fprintf('  chunk %d/%d (trials %d-%d) done in %.1f min | %s | ETA for the rest at this pace %.0f min\n', ...
            c, nChunks, idx(1), idx(end), secs / 60, char(datetime('now')), mean(chunkSecs) * (nChunks - c) / 60);
    end

elseif USE_PARFOR
    dq             = parallel.pool.DataQueue;
    reportInterval = max(floor(nSel/20), 1);
    afterEach(dq, @(s) reportProgress_local(s, nSel, tRun, reportInterval));

    parfor s = 1:nSel
        tr  = trialsSel(s);
        res = processOneTrial_local(tr, FS_RUN, sigToMM, betaGenVec, ...
            N_REPS, N_R7, EDGE_CLIP, derivConfigs, regressConfigs, ...
            nDeriv, nRegress, nPipelines, NOISE_MODEL, END_ANCHORED, CARRIER_ALPHA, RNG_SEED, NH_OPT, SEG_BUDGET, NOISE_SCALE, SEG_RULE);

        [betaGenStar, ciLo, ciHi, betaGenStarMed, nearestIdx, invertInfo] = ...
            invertBeta_local(res, betaGenVec, nPipelines, monoParams); %#ok<PFBNS>

        betaRecSlice = extractSlice_local(res, nearestIdx);
        r = packResult_local(tr, res, betaRecSlice, betaGenStar, ...
            betaGenStarMed, ciLo, ciHi, ...
            alphaSelected(s), sigmaSelected(s), NOISE_MODEL, invertInfo); %#ok<PFBNS>
            r = addWholeTrial_local(r, res, betaGenVec, nPipelines, monoParams, NOISE_MODEL);   % v022
        resultsCells{s} = r;
        logStrings{s}   = sprintf('%4d  %-28s  alpha=%6.3f  gen*_med=%+.4f  gen*_sg=%+.4f  status_sg=%s', ...
            s, r.trialID, r.alphaIRA, betaGenStarMed, betaGenStar(6), r.invertStatus(6));
        send(dq, s); %#ok<PFBNS>
    end

else
    for s = 1:nSel
        tr  = trialsSel(s);
        res = processOneTrial_local(tr, FS_RUN, sigToMM, betaGenVec, ...
            N_REPS, N_R7, EDGE_CLIP, derivConfigs, regressConfigs, ...
            nDeriv, nRegress, nPipelines, NOISE_MODEL, END_ANCHORED, CARRIER_ALPHA, RNG_SEED, NH_OPT, SEG_BUDGET, NOISE_SCALE, SEG_RULE);

        [betaGenStar, ciLo, ciHi, betaGenStarMed, nearestIdx, invertInfo] = ...
            invertBeta_local(res, betaGenVec, nPipelines, monoParams);

        betaRecSlice = extractSlice_local(res, nearestIdx);
        r = packResult_local(tr, res, betaRecSlice, betaGenStar, ...
            betaGenStarMed, ciLo, ciHi, ...
            alphaSelected(s), sigmaSelected(s), NOISE_MODEL, invertInfo);
        r = addWholeTrial_local(r, res, betaGenVec, nPipelines, monoParams, NOISE_MODEL);   % v022
        resultsCells{s} = r;
        elapsed = toc(tRun);
        eta     = elapsed / s * (nSel - s);
        fprintf('[%4d/%4d] %-28s  gen*_med=%+.4f  gen*(SG-IRLS)=%+.4f  status(SG-IRLS)=%-9s  ETA %.0f min\n', ...
            s, nSel, r.trialID, betaGenStarMed, betaGenStar(6), r.invertStatus(6), eta/60);
        logStrings{s} = sprintf('%4d  %-28s  alpha=%6.3f  gen*_med=%+.4f  gen*_sg=%+.4f  status_sg=%s', ...
            s, r.trialID, r.alphaIRA, betaGenStarMed, betaGenStar(6), r.invertStatus(6));
    end
end

results  = [resultsCells{:}];
totalMin = toc(tRun) / 60;

%% ---------- R2: loop-closure CCC (post-hoc) ------------------------------
loopCCC = computeLoopCCC_local(results, nPipelines);
fprintf('\nLoop-closure CCC (constellation-centre, N_trials=%d): %.4f\n', ...
    nSel, loopCCC);

%% ---------- R8: per-pipeline invertibility coverage (post-hoc) -----------
fprintf('\nR8 per-pipeline invertibility status (N_trials=%d):\n', nSel);
fprintf(' %-12s  %6s  %6s  %10s  %8s  %8s\n', ...
    'Pipeline','rise','desc','ambiguous','neither','no_obs');
statusCounts = zeros(nPipelines, 5);
for pp = 1:nPipelines
    st = arrayfun(@(r) r.invertStatus(pp), results);
    statusCounts(pp,1) = sum(st == "rise");
    statusCounts(pp,2) = sum(st == "desc");
    statusCounts(pp,3) = sum(st == "ambiguous");
    statusCounts(pp,4) = sum(st == "neither");
    statusCounts(pp,5) = sum(st == "no_beta_obs");
    fprintf(' %-12s  %6d  %6d  %10d  %8d  %8d\n', pipelineLabels(pp), statusCounts(pp,:));
end

%% ---------- Write log -----------------------------------------------------
fLog = fopen(logFile, 'w');
fprintf(fLog, 'runLoopClosureFftnoise_v028  |  %s  |  %s  |  model=%s  |  parfor=%d\n', ...
    DATASET, TRIAL_SELECTION, NOISE_MODEL, USE_PARFOR);
fprintf(fLog, 'N=%d  N_BETA=%d  BETA=[%.2f %.2f]  N_REPS=%d  N_R7=%d  SEED=%d\n', ...
    nSel, N_BETA, BETA_RANGE(1), BETA_RANGE(2), N_REPS, N_R7, RNG_SEED);
fprintf(fLog, 'MonoSlopeTol=%.3f  MonoSmoothWidth=%.4f  MonoMinSegWidth=%.4f  MonoTrimToleranceK=%.2f\n', ...
    monoParams.slopeTol, monoParams.smoothWidth, monoParams.minSegWidth, monoParams.trimToleranceK);
fprintf(fLog, 'Loop-closure CCC: %.4f\n', loopCCC);
fprintf(fLog, 'Start: %s\n\n', datestr(now));
for s = 1:nSel
    if ~isempty(logStrings{s}), fprintf(fLog, '%s\n', logStrings{s}); end
end
fprintf(fLog, '\nEnd: %s  Total: %.1f min\n', datestr(now), totalMin);
fclose(fLog);

%% ---------- Save ----------------------------------------------------------
config = struct('DATASET', DATASET, 'TRIAL_SELECTION', TRIAL_SELECTION, ...
    'TRIAL_ID_OVR', TRIAL_ID_OVR, 'NOISE_MODEL', NOISE_MODEL, ...
    'NOISE_MODEL_OPT', NOISE_MODEL_OPT, ...
    'N_BETA', N_BETA, 'BETA_RANGE', BETA_RANGE, ...
    'N_REPS', N_REPS, 'N_R7', N_R7, ...
    'RNG_SEED', RNG_SEED, 'USE_PARFOR', USE_PARFOR, 'N_WORKERS', nWorkers, ...
    'monoParams', monoParams, 'SUBSAMPLE', SUBSAMPLE, 'FS_DS', FS_DS, 'FS_RUN', FS_RUN, ...
    'BW_CUTOFF', BW_CUTOFF, 'SG_REF', SG_REF, 'EDGE_CLIP', EDGE_CLIP, 'OPT_TAG', OPT_TAG, 'CHUNK_SIZE', CHUNK, 'END_ANCHORED_TEMPLATE', END_ANCHORED, 'TEMPLATE', opt.Template, 'SIG_TO_MM', sigToMM, 'CARRIER_ALPHA', CARRIER_ALPHA, 'V022_ENDMATCHED', "endToEndSegment_v004 commonLength Tol 1; endMatchedNoise_v001 Tol 0.3 MaxDraw 20; betaObs = central window of length L, betaObsWhole = whole trial; substream cellKey*1e5 + bi*1e3 + rep", 'NHARMONICS', NH_OPT, 'SEGMENT_BUDGET', SEG_BUDGET, 'NOISE_SCALE', NOISE_SCALE, 'SEGMENT_RULE', SEG_RULE);   % v028   % v027   % v026   % v022, v025   % v016, v017, v019, v020
runDate = datestr(now);
tmpOut = [outFile(1:end-4) '_tmp.mat'];   % v017: write under a temporary name, then rename (atomic)
save(tmpOut, 'results', 'pipelineLabels', 'betaGenVec', ...
    'config', 'runDate', 'loopCCC', '-v7.3');
movefile(tmpOut, outFile);
fprintf('Saved: %s\n', outFile);
if CHUNK > 0   % v017: the checkpoints are redundant once the output is written and readable
    w = whos('-file', outFile);
    if ~ismember('results', {w.name}), error('runLoopClosureFftnoise_v028:SaveCheck', '%s', ['Saved file has no results variable: ' outFile]); end
    rmdir(ckDir, 's');
end

%% ---------- Console summary -----------------------------------------------
fprintf('\n=== SUMMARY: %s | %s ===\n', DATASET, NOISE_MODEL);
fprintf(' %-12s  %8s  %9s  %9s  %8s  %8s  %6s\n', ...
    'Pipeline','med(obs)','med(gen*)','mn(gen*)','SD(gen*)','compress','N_conv');
for pp = 1:nPipelines
    bObs  = arrayfun(@(r) r.betaObs(pp),     results);
    bStar = arrayfun(@(r) r.betaGenStar(pp), results);
    fprintf(' %-12s  %8.4f  %9.4f  %9.4f  %8.4f  %8.4f  %6d\n', pipelineLabels(pp), ...
        median(bObs,'omitnan'), median(bStar,'omitnan'), ...
        mean(bStar,'omitnan'),  std(bStar, 0, 'omitnan'), ...
        median(bObs,'omitnan') - median(bStar,'omitnan'), sum(isfinite(bStar)));
end
bgsm = arrayfun(@(r) r.betaGenStarMed, results);
fprintf('\n betaGenStarMed  med=%.4f  mean=%.4f  SD=%.4f  N_fin=%d\n', ...
    median(bgsm,'omitnan'), mean(bgsm,'omitnan'), std(bgsm, 0, 'omitnan'), sum(isfinite(bgsm)));
fprintf(' Loop-closure CCC:                      %.4f\n', loopCCC);

if strcmp(NOISE_MODEL, 'shaped_xu')
    gapMaj = arrayfun(@(r) r.surrogateAlphaGapMaj, results);
    gapMin = arrayfun(@(r) r.surrogateAlphaGapMin, results);
    fprintf(' R7 surrogate alpha gap (maj/min):      %.3f / %.3f (mean |gap|)\n', ...
        mean(gapMaj, 'omitnan'), mean(gapMin, 'omitnan'));
end
fprintf('\nDone at %s  (%.1f min)\nLog: %s\n', ...
    datestr(now,'HH:MM:SS'), totalMin, logFile);

end % main function

%% ==========================================================================
%% LOCAL FUNCTIONS
%% ==========================================================================

function reportProgress_local(s, nSel, tRun, reportInterval)
    if mod(s, reportInterval) == 0 || s == nSel
        fprintf('  parfor progress: %d/%d  (%.0f min)\n', s, nSel, toc(tRun)/60);
    end
end

% --------------------------------------------------------------------------
function slice = extractSlice_local(res, nearestIdx)
    if all(isnan(res.betaRec(:)))
        slice = NaN(size(res.betaRec, 1), size(res.betaRec, 2));
    else
        slice = squeeze(res.betaRec(:, :, nearestIdx));
    end
end

% --------------------------------------------------------------------------
function r = packResult_local(tr, res, betaRecSlice, betaGenStar, ...
        betaGenStarMed, ciLo, ciHi, alphaIRA, sigmaMM, noiseModel, invertInfo)
    r.trialID           = char(tr.trialID);
    r.subjectID         = char(tr.subjectID);
    r.betaObs           = res.betaObs_vec;
    r.vgfObs            = res.vgfObs_vec;
    r.betaGenStar       = betaGenStar;
    r.betaGenStarMed    = betaGenStarMed;
    r.betaRecSlice      = betaRecSlice;
    % v015: persist the forward map the inversion is performed against.
    % These are the SAME two reductions invertBeta_local consumes (median
    % and std over reps), so any out-of-sample prediction made from them
    % uses the identical curve the inversion used, not a reconstruction.
    % [nPipelines x N_BETA]. v013 computed and discarded these.
    r.betaRecCurveMed   = squeeze(median(res.betaRec, 1, 'omitnan'));
    r.betaRecCurveSD    = squeeze(std(res.betaRec, 0, 1, 'omitnan'));
    r.ciLo              = ciLo;
    r.ciHi              = ciHi;
    r.smokeA            = res.smokeA;
    r.nEbrErr           = res.nEbrErr;     % v020
    r.ebrErrMsg         = res.ebrErrMsg;   % v020
    r.f0                = res.f0;
    r.a_mm              = res.a_mm;
    r.b_mm              = res.b_mm;
    r.theta             = res.theta;
    r.M                 = res.M;
    if isfield(res, 'nHarm'), r.nHarm = res.nHarm; else, r.nHarm = NaN; end   % v025
    r.alphaIRA          = alphaIRA;
    r.sigmaMM           = sigmaMM;
    r.alphaMaj          = res.alphaMaj;
    r.alphaMin          = res.alphaMin;
    if isfield(res, 'carrierAlphaMaj'), r.carrierAlphaMaj = res.carrierAlphaMaj; r.carrierAlphaMin = res.carrierAlphaMin; else, r.carrierAlphaMaj = NaN; r.carrierAlphaMin = NaN; end   % v021
    r.surrogateAlphaGapMaj = res.surrogateAlphaGapMaj;
    r.surrogateAlphaGapMin = res.surrogateAlphaGapMin;
    r.noiseModel        = noiseModel;
    r.invertStatus       = [invertInfo.status];
    r.riseInvertible      = [invertInfo.riseInvertible];
    r.descInvertible       = [invertInfo.descInvertible];
    r.nRiseRuns          = [invertInfo.nRiseRuns];
    r.nDescRuns          = [invertInfo.nDescRuns];
    r.riseRunLengths     = {invertInfo.riseRunLengths};
    r.descRunLengths     = {invertInfo.descRunLengths};
    r.matchedRunLength   = [invertInfo.matchedRunLength];
    r.monotonicCoverageFrac = [invertInfo.monotonicCoverageFrac];
end

% --------------------------------------------------------------------------
function frac = computeCoverageFrac_local(both, xGrid)
% Fraction of the swept beta_gen RANGE (not grid point count) covered by
% the UNION of every qualifying rise+desc run, independent of any one
% trial's betaObs. Answers "how locally monotonic is this curve overall",
% not "was this one value recoverable". Overlapping/adjacent run extents
% are merged, not double-counted.
    allRuns = [both.rise, both.desc];
    if isempty(allRuns) || isempty(xGrid)
        frac = 0;
        return
    end
    totalRange = max(xGrid) - min(xGrid);
    if totalRange <= 0
        frac = 0;
        return
    end

    intervals = zeros(numel(allRuns), 2);
    for k = 1:numel(allRuns)
        intervals(k, :) = [min(allRuns(k).x), max(allRuns(k).x)];
    end
    intervals = sortrows(intervals, 1);

    merged  = intervals(1, :);
    covered = 0;
    for k = 2:size(intervals, 1)
        if intervals(k, 1) <= merged(end, 2)
            merged(end, 2) = max(merged(end, 2), intervals(k, 2));
        else
            covered = covered + (merged(end, 2) - merged(end, 1));
            merged  = [merged; intervals(k, :)]; %#ok<AGROW>
        end
    end
    covered = covered + (merged(end, 2) - merged(end, 1));

    frac = covered / totalRange;
end

% --------------------------------------------------------------------------
function [betaGenStar, ciLo, ciHi, betaGenStarMed, nearestIdx, invertInfo] = ...
        invertBeta_local(res, betaGenVec, nPipelines, monoParams)
% R11 (v012, unchanged in v013): noise-tolerant local-monotonicity gate.
%
% Calls findBothBranches_v008 (grid-invariant width-based
% smoothing/segment-floor, PLUS a noise-calibrated tolerance on the
% strict-monotonicity trim step -- see that function's own docstring,
% and docs/FINDINGS_REFERENCE.md Finding #156, for the fix itself).
% Detection still returns ALL qualifying runs per direction (v002's own
% contribution, R9, unchanged in spirit) -- v008 changes ONLY how a
% candidate run gets trimmed once found, not which points are grouped
% into one. Membership is still checked across the FULL combined set of
% rise+desc runs -- a pipeline whose observed beta falls inside more
% than one qualifying run TOTAL (same direction or opposite) is
% genuinely ambiguous at this trial -- left NaN (fail loud) rather than
% guessing which run is "correct".
%
% invertInfo(pp).status one of:
%   "rise"        - inverted via a rising branch (exactly one run covered it)
%   "desc"        - inverted via a descending branch (exactly one run covered it)
%   "ambiguous"   - covered by more than one qualifying run total; left NaN
%   "neither"     - covered by no qualifying run; left NaN
%   "no_beta_obs" - betaObs itself was NaN (pipeline regression failed)

    betaGenStar = NaN(nPipelines, 1);
    ciLo        = NaN(nPipelines, 1);
    ciHi        = NaN(nPipelines, 1);
    invertInfo(nPipelines, 1) = struct('status', "no_beta_obs", ...
        'riseInvertible', false, 'descInvertible', false, ...
        'nRiseRuns', 0, 'nDescRuns', 0, ...
        'riseRunLengths', [], 'descRunLengths', [], ...
        'matchedRunLength', NaN, 'monotonicCoverageFrac', 0);
    for pp = 1:nPipelines
        invertInfo(pp) = struct('status', "no_beta_obs", ...
            'riseInvertible', false, 'descInvertible', false, ...
            'nRiseRuns', 0, 'nDescRuns', 0, ...
            'riseRunLengths', [], 'descRunLengths', [], ...
            'matchedRunLength', NaN, 'monotonicCoverageFrac', 0);
    end

    for pp = 1:nPipelines
        bObs = res.betaObs_vec(pp);
        if isnan(bObs), continue; end

        medC = squeeze(median(res.betaRec(:, pp, :), 1, 'omitnan'))';
        sigmaY = squeeze(std(res.betaRec(:, pp, :), 0, 1, 'omitnan'))';
        finMask = isfinite(medC) & isfinite(sigmaY) & isfinite(betaGenVec);
        if sum(finMask) < 2
            invertInfo(pp).status = "neither";
            continue
        end

        both = findBothBranches_v008(betaGenVec(finMask), medC(finMask), sigmaY(finMask), monoParams);
        invertInfo(pp).riseInvertible = ~isempty(both.rise);
        invertInfo(pp).descInvertible = ~isempty(both.desc);

        % Local-monotonicity structure, recorded here rather than discarded
        % once classification is done -- see this file's own header for why.
        % Computed regardless of this trial's own betaObs/classification
        % outcome, so
        % it answers "what does this curve's monotonicity structure look
        % like" independent of whether THIS trial's value was recoverable.
        invertInfo(pp).nRiseRuns    = numel(both.rise);
        invertInfo(pp).nDescRuns    = numel(both.desc);
        invertInfo(pp).riseRunLengths = arrayfun(@(seg) seg.n, both.rise);
        invertInfo(pp).descRunLengths = arrayfun(@(seg) seg.n, both.desc);
        invertInfo(pp).monotonicCoverageFrac = ...
            computeCoverageFrac_local(both, betaGenVec(finMask));

        % Combine rise+desc runs into one list, tag each with its own
        % direction so the CI step below knows which .x/.y to trust.
        allRuns = [both.rise, both.desc];
        nRiseRuns = numel(both.rise);
        covering = false(numel(allRuns), 1);
        for k = 1:numel(allRuns)
            covering(k) = bObs >= min(allRuns(k).y) - eps && bObs <= max(allRuns(k).y) + eps;
        end
        nCover = sum(covering);

        if nCover == 0
            invertInfo(pp).status = "neither";
            continue
        elseif nCover > 1
            invertInfo(pp).status = "ambiguous";
            continue
        end

        matchIdx = find(covering);
        branch = allRuns(matchIdx);
        invertInfo(pp).matchedRunLength = branch.n;
        if matchIdx <= nRiseRuns
            invertInfo(pp).status = "rise";
        else
            invertInfo(pp).status = "desc";
        end

        betaGenStar(pp) = interp1(branch.y, branch.x, bObs, 'linear', 'extrap');

        % CI: same approach as v009 -- restrict the 5th/95th percentile
        % curves to the SAME grid points as the matched branch, sort+interp
        % within that sub-range. Direction-agnostic min/max fix (v009's own
        % 2026-08-09 correction) carried forward unchanged.
        branchMask = ismember(betaGenVec, branch.x);
        p5c  = squeeze(prctile(res.betaRec(:, pp, :),  5, 1))';
        p95c = squeeze(prctile(res.betaRec(:, pp, :), 95, 1))';
        gBranch = betaGenVec(branchMask);

        [p5s,  o5]  = sort(p5c(branchMask),  'ascend');
        [p95s, o95] = sort(p95c(branchMask), 'ascend');
        ciA = safeInterp1_local(p5s,  gBranch(o5),  bObs);
        ciB = safeInterp1_local(p95s, gBranch(o95), bObs);
        if isnan(ciA) || isnan(ciB)
            % v013's ONE change from v012 (Finding pending, 2026-08-16):
            % genuinely undefined, not zero-width. MATLAB's max/min
            % silently ignore NaN when comparing two scalars
            % (max(NaN,5) returns 5, not NaN), which was collapsing an
            % undefined bound into a false, deceptively precise
            % zero-width CI. Found via a downstream convergence-check
            % script's own output (Cook ASD dataset, three separate
            % N_REPS levels each showing ciWidthMean=ciWidthSD=0.0000
            % across all 50 independent subsamples simultaneously -- not
            % a coincidence of genuine convergence). See this file's own
            % header for the full account and re-run implications.
            ciLo(pp) = NaN;
            ciHi(pp) = NaN;
        else
            ciLo(pp) = min(ciA, ciB);
            ciHi(pp) = max(ciA, ciB);
        end
    end

    finite_      = betaGenStar(isfinite(betaGenStar));
    if isempty(finite_)
        betaGenStarMed = NaN;
        nearestIdx     = 1;
    else
        betaGenStarMed = median(finite_);
        [~, nearestIdx] = min(abs(betaGenVec - betaGenStarMed));
    end
end

% --------------------------------------------------------------------------
function ccc = computeLoopCCC_local(results, nPipelines)
    nTrials = numel(results);
    maxPairs = nTrials * nPipelines;
    allObs  = NaN(maxPairs, 1);
    allPred = NaN(maxPairs, 1);
    cnt = 0;
    for s = 1:nTrials
        r = results(s);
        if all(isnan(r.betaRecSlice(:))), continue; end
        for pp = 1:nPipelines
            bObs  = r.betaObs(pp);
            bPred = median(r.betaRecSlice(:, pp), 'omitnan');
            if isfinite(bObs) && isfinite(bPred)
                cnt = cnt + 1;
                allObs(cnt)  = bObs;
                allPred(cnt) = bPred;
            end
        end
    end
    ccc = linCCC_local(allObs(1:cnt), allPred(1:cnt));
end

% --------------------------------------------------------------------------
function ccc = linCCC_local(x, y)
    valid = isfinite(x) & isfinite(y);
    x = x(valid); y = y(valid);
    if numel(x) < 3, ccc = NaN; return; end
    mx  = mean(x); my  = mean(y);
    sx2 = var(x, 0); sy2 = var(y, 0);
    sxy = mean((x - mx) .* (y - my));
    denom = sx2 + sy2 + (mx - my)^2;
    if denom < eps, ccc = NaN; return; end
    ccc = 2*sxy / denom;
end

% --------------------------------------------------------------------------
function res = processOneTrial_local(tr, FS, sigToMM, betaGenVec, N_REPS, N_R7, ...
        EDGE_CLIP, derivConfigs, regressConfigs, nDeriv, nRegress, nPipelines, noiseModel, endAnchored, carrierAlpha, rngSeed, nHopt, segBudget, noiseScale, segRule)
% UNCHANGED from runner versions v008/v009 (the OLD grid-based/first
% per-trial runner FILES -- not to be confused with findBothBranches_v008,
% the segment-finder function invertBeta_local calls below, an unrelated
% same-numbered file) -- generates the per-trial forward map. Only the
% inversion of this map's output differs in v012/v013 (see invertBeta_local).

    N_BETA  = numel(betaGenVec);
    needIRA = ismember(noiseModel, {'xu', 'shaped_xu'});
    doR7    = strcmp(noiseModel, 'shaped_xu');
    res.nEbrErr = 0;  res.ebrErrMsg = "";   % v020: regressDataEBR errors, counted (were silent)

    betaObs = NaN(nDeriv, nRegress);
    vgfObs  = NaN(nDeriv, nRegress);
    for dIdx = 1:nDeriv
        [dx, dy] = differentiateKinematicsEBR(tr.x, tr.y, ...
            derivConfigs(dIdx).filterType, derivConfigs(dIdx).filterParams, FS);
        vx = dx(EDGE_CLIP:end-EDGE_CLIP, 2);
        vy = dy(EDGE_CLIP:end-EDGE_CLIP, 2);
        ax = dx(EDGE_CLIP:end-EDGE_CLIP, 3);
        ay = dy(EDGE_CLIP:end-EDGE_CLIP, 3);
        sp = hypot(vx, vy);
        kp = curvatureKinematicEBR(vx, vy, ax, ay);
        lm = [1, -1/3];
        for rIdx = 1:nRegress
            try
                [b, v] = regressDataEBR(sp, kp, regressConfigs(rIdx).type, lm, 0, 0);
                betaObs(dIdx, rIdx) = b;
                vgfObs(dIdx, rIdx)  = v;
                if rIdx == 1, lm = [v, b]; end
            catch ME   % v020: counted, not silent
                res.nEbrErr = res.nEbrErr + 1;  res.ebrErrMsg = string(ME.identifier) + ": " + string(ME.message);
            end
        end
    end
    res.betaObs_vec = betaObs(:)';
    res.vgfObs_vec  = vgfObs(:)';

    res.f0 = estimateF0_local(double(tr.x), double(tr.y), FS);
    if isnumeric(nHopt), nHarm = nHopt; else, nHarm = max(1, floor(derivConfigs(1).filterParams(2) / res.f0)); end   % v025
    res.nHarm = nHarm;
    [xResid, yResid, xFit, yFit, tsClip, tsTail] = templateSubtract_local( ...   % v024: + tail   % v022: + clip
        double(tr.x), double(tr.y), FS, res.f0, nHarm, 4, endAnchored);   % v025: nHarm (was 4)
    M = numel(xResid); res.M = M;

    nanTemplate = NaN(N_REPS, nPipelines, N_BETA);
    nanSlimVec  = NaN(1, nPipelines);
    if M < 4*(EDGE_CLIP + 50)
        res.betaRec            = nanTemplate;
        res.smokeA             = nanSlimVec;
        res.a_mm               = NaN; res.b_mm  = NaN; res.theta = NaN;
        res.alphaMaj           = NaN; res.alphaMin = NaN;
        res.surrogateAlphaGapMaj = NaN;
        res.surrogateAlphaGapMin = NaN;
        return
    end

    xMM = double(tr.x)*sigToMM; xMM = xMM - mean(xMM);
    yMM = double(tr.y)*sigToMM; yMM = yMM - mean(yMM);
    [V, D]  = eig(cov(xMM, yMM));
    lams    = diag(D); [lams, ord] = sort(lams, 'descend'); V = V(:, ord);
    a_mm    = sqrt(2*max(lams(1), 0));
    b_mm    = sqrt(2*max(lams(2), 0));
    theta   = atan2(V(2,1), V(1,1));
    a_nat   = a_mm / sigToMM;
    b_nat   = b_mm / sigToMM;
    res.a_mm = a_mm; res.b_mm = b_mm; res.theta = theta;

    resid_maj = xResid*cos(theta) + yResid*sin(theta);
    resid_min = -xResid*sin(theta) + yResid*cos(theta);
    fftMaj    = fft(resid_maj(:));
    fftMin    = fft(resid_min(:));

    alphaMaj = NaN; alphaMin = NaN;
    if needIRA
        IRA_FLOW  = 1.0;
        IRA_FHIGH = min(20.0, FS/2 - 1);
        IRA_HSET  = 1.1:0.05:1.9;
        try
            alphaMaj = iraAlphaSigma_v001(resid_maj(:), FS, IRA_FLOW, IRA_FHIGH, IRA_HSET);
            alphaMin = iraAlphaSigma_v001(resid_min(:), FS, IRA_FLOW, IRA_FHIGH, IRA_HSET);
        catch ME
            warning('processOneTrial_local:iraFail', '%s', ...
                ['IRASA per-axis alpha failed: ' ME.message]);
        end
        if ~isfinite(alphaMaj) || ~isfinite(alphaMin)
            res.betaRec  = nanTemplate;
            res.smokeA   = nanSlimVec;
            res.alphaMaj = alphaMaj; res.alphaMin = alphaMin;
            res.surrogateAlphaGapMaj = NaN;
            res.surrogateAlphaGapMin = NaN;
            return
        end
    end
    res.alphaMaj = alphaMaj; res.alphaMin = alphaMin;
    cMaj = alphaMaj;  cMin = alphaMin;   % v021: carrier alpha for the forward-map noise (default = v020's IRASA)
    if needIRA && ~strcmp(carrierAlpha, 'irasa')
        cMaj = NaN;  cMin = NaN;
        try
            switch carrierAlpha
                case 'pwlw'
                    cMaj = alphaPWLW_v001(resid_maj(:), FS, [1.0, min(20.0, FS/2 - 1)]);
                    cMin = alphaPWLW_v001(resid_min(:), FS, [1.0, min(20.0, FS/2 - 1)]);
                case 'irasaRecordBand'
                    fLoRB = max(1.0, 20 / (numel(resid_maj) / FS));
                    cMaj = iraAlphaSigma_v001(resid_maj(:), FS, fLoRB, min(20.0, FS/2 - 1), 1.1:0.05:1.9);
                    cMin = iraAlphaSigma_v001(resid_min(:), FS, fLoRB, min(20.0, FS/2 - 1), 1.1:0.05:1.9);
            end
        catch ME
            warning('processOneTrial_local:carrierFail', '%s', ['Carrier alpha (' carrierAlpha ') failed: ' ME.message]);
        end
        if ~isfinite(cMaj) || ~isfinite(cMin)
            res.betaRec  = nanTemplate;  res.smokeA = nanSlimVec;
            res.carrierAlphaMaj = cMaj;  res.carrierAlphaMin = cMin;
            res.surrogateAlphaGapMaj = NaN;  res.surrogateAlphaGapMin = NaN;
            return
        end
    end
    res.carrierAlphaMaj = cMaj;  res.carrierAlphaMin = cMin;   % v021

    surrogateAlphaGapMaj = NaN;
    surrogateAlphaGapMin = NaN;
    if doR7 && isfinite(alphaMaj) && isfinite(alphaMin)
        IRA_FLOW  = 1.0;
        IRA_FHIGH = min(20.0, FS/2 - 1);
        IRA_HSET  = 1.1:0.05:1.9;
        gapsMaj = NaN(N_R7, 1);
        gapsMin = NaN(N_R7, 1);
        for k = 1:N_R7
            try
                [nMaj_r7, nMin_r7] = generateLoopClosureNoise_v003('shaped_xu', ...
                    resid_maj(:), resid_min(:), FS, alphaMaj, alphaMin);
                aMaj_r7 = iraAlphaSigma_v001(nMaj_r7, FS, IRA_FLOW, IRA_FHIGH, IRA_HSET);
                aMin_r7 = iraAlphaSigma_v001(nMin_r7, FS, IRA_FLOW, IRA_FHIGH, IRA_HSET);
                gapsMaj(k) = abs(aMaj_r7 - alphaMaj);
                gapsMin(k) = abs(aMin_r7 - alphaMin);
            catch ME
                warning('processOneTrial_local:r7Fail', '%s', ...
                    sprintf('R7 IRASA check (rep %d): %s', k, ME.message));
            end
        end
        surrogateAlphaGapMaj = mean(gapsMaj, 'omitnan');
        surrogateAlphaGapMin = mean(gapsMin, 'omitnan');
    end
    res.surrogateAlphaGapMaj = surrogateAlphaGapMaj;
    res.surrogateAlphaGapMin = surrogateAlphaGapMin;

    MTPL = M;   % v022: template length (the whole residual unless endMatched)
    if ismember(noiseModel, {'endMatched', 'shapedXuSeg'})   % v022; v023: shapedXuSeg too
        res.emL = NaN;  res.emLfrac = NaN;  res.emSeg = NaN(2, 2);  res.emAlphaSeg = [NaN NaN];  res.emStatus = "ok";  res.nNoiseRejected = 0;
        res.betaObsWhole_vec = res.betaObs_vec;  res.vgfObsWhole_vec = res.vgfObs_vec;
        try
            if N_REPS >= 1000 || N_BETA >= 100, error('runLoopClosureFftnoise_v028:emSubstream', 'substream layout needs N_REPS < 1000 and N_BETA < 100'); end
            Mc = M - tsTail;  res.emTail = tsTail;   % v024: search only the covered residual
            [sI0, sI1, emInfo] = endToEndSegment_v006([resid_maj(1:Mc), resid_min(1:Mc)], Rule=segRule, Tol=1, Budget=segBudget);   % v028: v006 (auto = v005)   % v026: v005 (= v004 when segBudget is empty)
            res.emRule = emInfo.rule;
            Lsyn = emInfo.L;
            if Lsyn < 4*(EDGE_CLIP + 50), error('runLoopClosureFftnoise_v028:emShort', 'end-matched length %d is below the runner minimum %d', Lsyn, 4*(EDGE_CLIP + 50)); end
            segMaj = resid_maj(sI0(1):sI1(1));  segMin = resid_min(sI0(2):sI1(2));  segMaj = segMaj(:);  segMin = segMin(:);
            aSeg = [alphaPWLW_v001(segMaj, FS, [1.0, min(20.0, FS/2 - 1)]), alphaPWLW_v001(segMin, FS, [1.0, min(20.0, FS/2 - 1)])];
            if ~all(isfinite(aSeg)), error('runLoopClosureFftnoise_v028:emAlpha', 'segment alpha not finite (%g, %g)', aSeg(1), aSeg(2)); end
            c0 = floor((Mc - Lsyn) / 2) + 1;  wA = c0 + tsClip - 1;  wB = wA + Lsyn - 1;   % v026: restored (v024 had commented out wA, wB); central window of the covered span, length Lsyn
            [res.betaObs_vec, res.vgfObs_vec, nErrW, msgW] = betaObsOn_local(double(tr.x(wA:wB)), double(tr.y(wA:wB)), FS, EDGE_CLIP, derivConfigs, regressConfigs, nDeriv, nRegress);
            res.nEbrErr = res.nEbrErr + nErrW;  if nErrW > 0, res.ebrErrMsg = msgW; end
            res.emL = Lsyn;  res.emLfrac = emInfo.Lfrac;  res.emSeg = [sI0(:), sI1(:)];  res.emAlphaSeg = aSeg;  res.emWindow = [wA wB];
            MTPL = Lsyn;
        catch ME
            res.emStatus = "FAIL: " + string(ME.message);
            res.betaObs_vec(:) = NaN;  res.vgfObs_vec(:) = NaN;
            res.betaRec = nanTemplate;  res.smokeA = nanSlimVec;
            return
        end
    end

    nPhi  = 10000;
    phi   = linspace(0, 2*pi, nPhi)';
    dsdp  = sqrt((a_nat.*sin(phi)).^2 + (b_nat.*cos(phi)).^2);
    den   = (a_nat^2.*sin(phi).^2 + b_nat^2.*cos(phi).^2).^1.5;
    kp_phi = (a_nat*b_nat) ./ max(den, eps);
    tVec  = (0:MTPL-1)' / FS;   % v022: MTPL = M unless endMatched
    templates = zeros(MTPL, 2, N_BETA);   % v022
    for bi = 1:N_BETA
        wt   = dsdp .* kp_phi.^betaGenVec(bi);
        cumT = cumsum(wt); cumT = cumT / cumT(end);
        tN   = min(mod(tVec*res.f0, 1), 1 - 1e-9);
        phiA = interp1(cumT, phi, tN, 'linear', 'extrap');
        xE   = a_nat*cos(phiA);
        yE   = b_nat*sin(phiA);
        templates(:, 1, bi) = xE*cos(theta) - yE*sin(theta);
        templates(:, 2, bi) = xE*sin(theta) + yE*cos(theta);
    end

    betaRec = NaN(N_REPS, nPipelines, N_BETA);
    for bi = 1:N_BETA
        xTpl = templates(:, 1, bi);
        yTpl = templates(:, 2, bi);
        for rep = 1:N_REPS
            switch noiseModel
                case 'fftnoise'
                    [nMaj, nMin] = generateLoopClosureNoise_v003('fftnoise', ...
                        fftMaj, fftMin, FS);
                case 'xu'
                    [nMaj, nMin] = generateLoopClosureNoise_v003('xu', ...
                        resid_maj(:), resid_min(:), FS, cMaj, cMin);   % v021: carrier alpha
                case 'shaped_xu'
                    [nMaj, nMin] = generateLoopClosureNoise_v003('shaped_xu', ...
                        resid_maj(:), resid_min(:), FS, cMaj, cMin);   % v021: carrier alpha
                case 'endMatched'   % v022: per-cell substream (README_SeedPractice), each axis's own end-matched segment
                    rsC = RandStream("threefry4x64_20", "Seed", rngSeed);  rsC.Substream = tr.cellKey * 1e5 + bi * 1e3 + rep;
                    rsOld = RandStream.setGlobalStream(rsC);
                    [nMaj, nMin] = endMatchedNoise_v001(segMaj, segMin, FS, aSeg(1), aSeg(2));
                    nMaj = nMaj * sqrt(noiseScale);  nMin = nMin * sqrt(noiseScale);   % v027: variance x noiseScale (1 = unchanged)
                    RandStream.setGlobalStream(rsOld);
                    if any(isnan(nMaj)) || any(isnan(nMin)), res.nNoiseRejected = res.nNoiseRejected + 1; end
                case 'shapedXuSeg'   % v023: production shaped_xu on each axis's segment, carrier at the segment's PWLW alpha, same substreams
                    rsC = RandStream("threefry4x64_20", "Seed", rngSeed);  rsC.Substream = tr.cellKey * 1e5 + bi * 1e3 + rep;
                    rsOld = RandStream.setGlobalStream(rsC);
                    [nMaj, nMin] = generateLoopClosureNoise_v003('shaped_xu', segMaj, segMin, FS, aSeg(1), aSeg(2));
                    RandStream.setGlobalStream(rsOld);
                    if any(isnan(nMaj)) || any(isnan(nMin)), res.nNoiseRejected = res.nNoiseRejected + 1; end
                case 'bootstrap'
                    [nMaj, nMin] = generateLoopClosureNoise_v003('bootstrap', ...
                        resid_maj(:), resid_min(:), FS);
            end
            if any(isnan(nMaj)) || any(isnan(nMin))
                continue
            end
            xSyn = xTpl + nMaj*cos(theta) - nMin*sin(theta);
            ySyn = yTpl + nMaj*sin(theta) + nMin*cos(theta);
            for dIdx = 1:nDeriv
                [dx, dy] = differentiateKinematicsEBR(xSyn, ySyn, ...
                    derivConfigs(dIdx).filterType, ...
                    derivConfigs(dIdx).filterParams, FS);
                vx = dx(EDGE_CLIP:end-EDGE_CLIP, 2);
                vy = dy(EDGE_CLIP:end-EDGE_CLIP, 2);
                ax = dx(EDGE_CLIP:end-EDGE_CLIP, 3);
                ay = dy(EDGE_CLIP:end-EDGE_CLIP, 3);
                sp = hypot(vx, vy);
                kp = curvatureKinematicEBR(vx, vy, ax, ay);
                lm = [1, -1/3];
                for rIdx = 1:nRegress
                    ppCM = (rIdx-1)*nDeriv + dIdx;
                    try
                        [b, v] = regressDataEBR(sp, kp, ...
                            regressConfigs(rIdx).type, lm, 0, 0);
                        betaRec(rep, ppCM, bi) = b;
                        if rIdx == 1, lm = [v, b]; end
                    catch ME   % v020: counted, not silent
                        res.nEbrErr = res.nEbrErr + 1;  res.ebrErrMsg = string(ME.identifier) + ": " + string(ME.message);
                    end
                end
            end
        end
    end
    res.betaRec = betaRec;

    nMaj_det = real(ifft(fftMaj));
    nMin_det = real(ifft(fftMin));
    xDet     = xFit + nMaj_det*cos(theta) - nMin_det*sin(theta);
    yDet     = yFit + nMaj_det*sin(theta) + nMin_det*cos(theta);
    smokeA   = NaN(1, nPipelines);
    for dIdx = 1:nDeriv
        [dx, dy] = differentiateKinematicsEBR(xDet, yDet, ...
            derivConfigs(dIdx).filterType, derivConfigs(dIdx).filterParams, FS);
        vx = dx(EDGE_CLIP:end-EDGE_CLIP, 2);
        vy = dy(EDGE_CLIP:end-EDGE_CLIP, 2);
        ax = dx(EDGE_CLIP:end-EDGE_CLIP, 3);
        ay = dy(EDGE_CLIP:end-EDGE_CLIP, 3);
        sp = hypot(vx, vy);
        kp = curvatureKinematicEBR(vx, vy, ax, ay);
        lm = [1, -1/3];
        for rIdx = 1:nRegress
            ppCM = (rIdx-1)*nDeriv + dIdx;
            try
                [b, ~] = regressDataEBR(sp, kp, ...
                    regressConfigs(rIdx).type, lm, 0, 0);
                smokeA(ppCM) = b - res.betaObs_vec(ppCM);
            catch ME   % v020: counted, not silent
                res.nEbrErr = res.nEbrErr + 1;  res.ebrErrMsg = string(ME.identifier) + ": " + string(ME.message);
            end
            if rIdx == 1
                lm = [res.vgfObs_vec(ppCM), res.betaObs_vec(ppCM)];
            end
        end
    end
    if ismember(noiseModel, {'endMatched', 'shapedXuSeg'}), smokeA(:) = NaN; end   % v023: shapedXuSeg too;   % v022: the determinism check is whole-trial fftnoise; not comparable
    res.smokeA = smokeA;
end

% --------------------------------------------------------------------------
function [bObs, vObs, nErr, msg] = betaObsOn_local(x, y, FS, EDGE_CLIP, derivConfigs, regressConfigs, nDeriv, nRegress)
% v022: the runner's beta_obs block (processOneTrial_local, top) on a given stretch of trial; ppCM = (rIdx-1)*nDeriv + dIdx.
    betaObs = NaN(nDeriv, nRegress);  vgfObs = NaN(nDeriv, nRegress);  nErr = 0;  msg = "";
    for dIdx = 1:nDeriv
        [dx, dy] = differentiateKinematicsEBR(x, y, derivConfigs(dIdx).filterType, derivConfigs(dIdx).filterParams, FS);
        vx = dx(EDGE_CLIP:end-EDGE_CLIP, 2);  vy = dy(EDGE_CLIP:end-EDGE_CLIP, 2);
        ax = dx(EDGE_CLIP:end-EDGE_CLIP, 3);  ay = dy(EDGE_CLIP:end-EDGE_CLIP, 3);
        sp = hypot(vx, vy);  kp = curvatureKinematicEBR(vx, vy, ax, ay);  lm = [1, -1/3];
        for rIdx = 1:nRegress
            try
                [b, v] = regressDataEBR(sp, kp, regressConfigs(rIdx).type, lm, 0, 0);
                betaObs(dIdx, rIdx) = b;  vgfObs(dIdx, rIdx) = v;
                if rIdx == 1, lm = [v, b]; end
            catch ME
                nErr = nErr + 1;  msg = string(ME.identifier) + ": " + string(ME.message);
            end
        end
    end
    bObs = betaObs(:)';  vObs = vgfObs(:)';
end

% --------------------------------------------------------------------------
function r = addWholeTrial_local(r, res, betaGenVec, nPipelines, monoParams, noiseModel)
% v022: for endMatched, the record of the end-matched noise and a second inversion from the whole-trial beta_obs against the
% same forward map. Every endMatched result gets the same fields (NaN where the trial stopped early), so the struct array joins.
    if ~ismember(noiseModel, {'endMatched', 'shapedXuSeg'}), return, end   % v023: shapedXuSeg too
    nanP = NaN(nPipelines, 1);
    r.betaObsBasis = "window";
    if ~isfield(res, 'emStatus'), r.betaObs(:) = NaN;  r.vgfObs(:) = NaN; end   % never reached the end-matched step: no window beta_obs
    fl = {'emL', NaN; 'emLfrac', NaN; 'emSeg', NaN(2, 2); 'emAlphaSeg', [NaN NaN]; 'emWindow', [NaN NaN]; 'emTail', NaN; 'emRule', ""; 'emStatus', "not reached (residual below the runner minimum or alpha failure)"; 'nNoiseRejected', NaN; 'betaObsWhole_vec', res.betaObs_vec; 'vgfObsWhole_vec', res.vgfObs_vec};
    for k = 1:size(fl, 1)
        if isfield(res, fl{k, 1}), v = res.(fl{k, 1}); else, v = fl{k, 2}; end
        r.(strrep(fl{k, 1}, '_vec', '')) = v;
    end
    r.betaGenStarWhole = nanP;  r.betaGenStarMedWhole = NaN;  r.ciLoWhole = nanP;  r.ciHiWhole = nanP;  r.invertStatusWhole = repmat("no_beta_obs", 1, nPipelines);
    if isfield(res, 'betaObsWhole_vec') && isfield(res, 'betaRec') && any(isfinite(res.betaRec(:)))
        resW = res;  resW.betaObs_vec = res.betaObsWhole_vec;
        [bW, loW, hiW, bWmed, ~, infoW] = invertBeta_local(resW, betaGenVec, nPipelines, monoParams);
        r.betaGenStarWhole = bW;  r.betaGenStarMedWhole = bWmed;  r.ciLoWhole = loW;  r.ciHiWhole = hiW;  r.invertStatusWhole = [infoW.status];
    end
end

% --------------------------------------------------------------------------
function f0 = estimateF0_local(x, y, fs)
    N    = numel(x); nfft = 2^nextpow2(4*N);
    Xf   = abs(fft(detrend(x(:), 'linear'), nfft));
    Yf   = abs(fft(detrend(y(:), 'linear'), nfft));
    fAx  = (0:nfft-1)' * fs / nfft;
    band = fAx > 0.1 & fAx < 5; idx = find(band);
    [~, pkX] = max(Xf(band)); [~, pkY] = max(Yf(band));
    f0x = fAx(idx(pkX)); f0y = fAx(idx(pkY));
    if abs(f0x - f0y) / max(f0x, 0.01) < 0.2
        f0 = mean([f0x, f0y]);
    elseif max(Xf(band)) > max(Yf(band))
        f0 = f0x;
    else
        f0 = f0y;
    end
end

% --------------------------------------------------------------------------
function [resX, resY, fitX, fitY, cs, tailN] = templateSubtract_local(x, y, fs, f0, nH, nCW, endAnchored)   % v024: + tailN   % v022: + cs
    N   = numel(x);
    win = round(nCW/f0*fs);
    hop = round(win/2);
    win = min(win, N);
    win = max(win, round(2/f0*fs));
    resX = zeros(N,1); resY = zeros(N,1); ww = zeros(N,1);
    starts = 1:hop:(N - win + 1);   % the production windows (identical to: s = 1; while s + win - 1 <= N, s = s + hop)
    if win > N, starts = []; end
    if endAnchored && ~isempty(starts) && starts(end) + win - 1 < N, starts(end+1) = N - win + 1; end   % v019: END-ANCHORED final window; v020: only when Template = 'endAnchored'
    for s = starts
        e   = s + win - 1;
        tw  = (0:win-1)'/fs;
        han = 0.5*(1 - cos(2*pi*(0:win-1)'/(win-1)));
        D   = [ones(win,1), tw];
        for h = 1:nH
            D(:,end+1) = cos(2*pi*h*f0*tw); %#ok<AGROW>
            D(:,end+1) = sin(2*pi*h*f0*tw); %#ok<AGROW>
        end
        Dw        = D .* han;
        bX        = Dw \ (x(s:e) .* han);
        bY        = Dw \ (y(s:e) .* han);
        resX(s:e) = resX(s:e) + (x(s:e) - D*bX) .* han;
        resY(s:e) = resY(s:e) + (y(s:e) - D*bY) .* han;
        ww(s:e)   = ww(s:e) + han;
    end
    ok = ww > 0;
    resX(ok) = resX(ok) ./ ww(ok);
    resY(ok) = resY(ok) ./ ww(ok);
    cl   = max(round(win/4), 1);
    cs   = cl;
    ce   = min(max(N - cl, cs + 100), N);
    lastEnd = 0;  if ~isempty(starts), lastEnd = starts(end) + win - 1; end   % v024: end of the last overlap-add window
    tailN = max(0, ce - lastEnd);   % v024: samples at the end of the kept residual that no window reached (exactly zero)
    fitX = x(cs:ce) - resX(cs:ce);
    fitY = y(cs:ce) - resY(cs:ce);
    resX = resX(cs:ce);
    resY = resY(cs:ce);
end

% --------------------------------------------------------------------------
function v = safeInterp1_local(xV, yV, xi)
    finMask = isfinite(xV) & isfinite(yV);
    if sum(finMask) < 2, v = NaN; return; end
    xVf = xV(finMask); yVf = yV(finMask);
    if xi < min(xVf) - eps || xi > max(xVf) + eps
        v = NaN;
    else
        v = interp1(xVf, yVf, xi, 'linear', 'extrap');
    end
end
