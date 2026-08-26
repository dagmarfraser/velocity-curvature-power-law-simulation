%RUNCONSTELLATIONRSA_HPC_V007 Per-dataset and pooled cross-dataset RSA
% (Mantel), full-N, BlueBEAR HPC version, Fraser included. NPerm raised
% to 10000 for the real, citable result.
%
% v007 changes from v006 (2026-08-07, Fraser D.S.):
%   v006 confirmed the full pipeline end-to-end at NPerm=100: no memory
%   issue (confirmed clean via sacct-informed 16-worker pooled step),
%   steady live DataQueue progress throughout, 52.2s pooled step total.
%   Only two changes here, following the exact discipline this project's
%   own history established (docs/TODO_ConstellationRSA_v001.md) --
%   confirm clean at low NPerm, THEN raise, never the reverse:
%     1. CONFIG.NPerm: 100 -> 10000.
%     2. CONFIG.ProgressEvery: 5 -> 500 (at NPerm=10000, printing every 5th
%        permutation would produce 2000 log lines; 500 gives ~20 update
%        lines across the run, enough to confirm it is alive without
%        flooding the diary).
%   CONFIG.OutMat also changed so v006's confirmed-working NPerm=100
%   result is not overwritten.
%
% Expected runtime (linear extrapolation from v006's clean NPerm=100
% timings): per-dataset loop (72 workers, largest cost from Pilot+Fraser)
% roughly 30-35 min; pooled step (16 workers, deliberately fewer -- see
% v006 header for why) roughly 87 min. Total approximately 2 hours.
% Request walltime accordingly (recommend >=3h for buffer).
%
% Everything else (72-worker per-dataset pool unchanged, separate
% 16-worker pool for the pooled step, single/uint16 broadcast shrink with
% doubles actually cleared, consistency check against mantelTest's own
% serial statistic, no full RDM matrices saved for the pooled result) is
% unchanged from v006 -- see that file's header for the full incident
% history (v003-v005) that led here.
%
% PORTABILITY FIX, 2026-08-26 (Fraser D.S., same session as the clean-clone
% Part G verification pass): CONFIG.SrcDir/MantelToolboxDir/ConstellationMetricsMat/
% LogDir/OutMat were hardcoded to this account's own BlueBEAR RDS absolute
% path, unusable by anyone else who clones this repo onto their own HPC.
% Replaced with fileparts(mfilename("fullpath")) (i.e. wherever src/ actually
% lives) plus relative fullfile() joins; MantelToolboxDir dropped entirely --
% callers are expected to addpath the published mantel toolbox
% (github.com/dagmarfraser/mantel) themselves, same assumption
% buildConstellationRDM_v002.m already makes. PooledWorkers now capped by
% feature("numcores") rather than assuming 16 are always available. No
% change to CONFIG.Datasets, CONFIG.NPerm, or any other computational
% parameter -- this is a path/portability fix only, and does not alter the
% already-computed, already-cited constellationRSA_HPC_v007.mat result.
% Edited in place (not v008) so every existing citation of "v007" stays
% valid.
%
% FLAGGED, NOT FIXED, same pass: CONFIG.Datasets below pools "Pilot" AND
% "Fraser" together, despite Pilot's 2026-08-07 retirement (Fraser is meant
% to supersede Pilot project-wide, not be pooled alongside it). This dataset
% list predates that decision reaching this file. NOT changed here -- this
% pass is a path fix only. Finding #138's pooled mantelR (+0.2165) was
% computed with Pilot included; whether to re-run pooled without Pilot is
% Dagmar's call, not made unilaterally here.
%
% Fraser, D.S. (2026) v007 (HPC variant)

%% -- CONFIG ------------------------------------------------------------
CONFIG = struct();
CONFIG.SrcDir                  = fileparts(mfilename("fullpath"));
CONFIG.ConstellationMetricsMat = fullfile(CONFIG.SrcDir, "constellationMetrics_v003.mat");
CONFIG.LogDir      = fullfile(CONFIG.SrcDir, "logs");
CONFIG.NoiseModel  = "shaped_xu";
CONFIG.Method      = "spearman";
CONFIG.NPerm       = 10000;    % RAISED from v006's confirmed-clean 100 -- see header
CONFIG.RngSeed     = 1729;
CONFIG.NWorkers       = feature("numcores");  % per-dataset loop
CONFIG.PooledWorkers  = min(16, feature("numcores"));  % pooled step ONLY -- capped 2026-08-26 so this doesn't request more workers than a smaller allocation/local machine actually has
CONFIG.Datasets    = ["Pilot", "Fraser", "Zarandi", "Cook CTRL", "Cook ASD", ...
                      "Hickman PLAC", "Hickman HALO"];  % Dhieb excluded (EPP contamination). FLAGGED 2026-08-26: Pilot+Fraser both pooled here, predates Pilot's retirement -- see header note, not changed this pass.
CONFIG.SaveMat     = true;
CONFIG.OutMat      = fullfile(CONFIG.SrcDir, "constellationRSA_HPC_v007.mat");
CONFIG.ProgressEvery = 500;    % RAISED from v006's 5 -- see header

%% -- Diary -----------------------------------------------------------------
if ~isfolder(CONFIG.LogDir), mkdir(CONFIG.LogDir); end
logFile = fullfile(CONFIG.LogDir, sprintf("runConstellationRSA_HPC_%s.log", ...
    string(datetime("now", "Format", "yyyyMMdd_HHmmss"))));
diary(logFile);
diaryCleanup = onCleanup(@() diary("off"));
fprintf("Diary: %s\n", logFile);
fprintf("Started: %s\n\n", string(datetime("now")));

%% -- Path setup --------------------------------------------------------------
addpath(CONFIG.SrcDir);
if exist("mantelTest", "file") ~= 2
    error("runConstellationRSA_HPC_v007:mantelTestNotFound", "%s", ...
        "mantelTest.m not found on path -- addpath the mantel toolbox " + ...
        "(github.com/dagmarfraser/mantel) before running this script.");
end
if exist("buildConstellationRDM_v002", "file") ~= 2
    error("runConstellationRSA_HPC_v007:builderNotFound", "%s", ...
        "buildConstellationRDM_v002.m not found on path (expected in src/).");
end
if ~isfile(CONFIG.ConstellationMetricsMat)
    error("runConstellationRSA_HPC_v007:missingMat", "%s", ...
        sprintf("%s not found.", CONFIG.ConstellationMetricsMat));
end

delete(gcp("nocreate"));
pool = parpool("local", CONFIG.NWorkers);
fprintf("Parpool started: %d workers (per-dataset loop)\n\n", pool.NumWorkers);

%% -- Load canonical per-trial constellation vectors -------------------------
S = load(CONFIG.ConstellationMetricsMat, "results");
allDatasets  = S.results.byNoiseModel.(CONFIG.NoiseModel).perDataset;
datasetNames = string({allDatasets.dataset});

%% -- Per-dataset RSA (UNCHANGED method from v003-v006 -- proven working) ----
fprintf("=== Per-dataset RSA (Mantel test on trial x trial RDMs) ===\n");
fprintf("Method=%s NPerm=%d UseParfor=true (full N, no cap)\n\n", CONFIG.Method, CONFIG.NPerm);

perDataset = struct("dataset", {}, "nTrials", {}, "mantelR", {}, ...
    "pValue", {}, "cvPred", {}, "cvObs", {});
poolData = struct("dataset", {}, "deltaBetaPred", {}, "deltaBetaObs", {}, ...
    "finiteMask", {});

for d = 1:numel(CONFIG.Datasets)
    dsName = CONFIG.Datasets(d);
    matchI = find(datasetNames == dsName, 1);
    if isempty(matchI)
        error("runConstellationRSA_HPC_v007:datasetNotFound", "%s", ...
            sprintf("Dataset '%s' not found in constellationMetrics_v003.mat.", dsName));
    end
    rec    = allDatasets(matchI);
    dbPred = rec.deltaBetaPred;
    dbObs  = rec.deltaBetaObs;
    finiteMask = isfinite(dbObs(:, 1)) & isfinite(dbPred(:, 1));

    poolData(end+1) = struct("dataset", dsName, "deltaBetaPred", dbPred, ...
        "deltaBetaObs", dbObs, "finiteMask", finiteMask); %#ok<AGROW>

    fprintf("Starting %s (n=%d)...\n", dsName, sum(finiteMask));
    tic;
    res = buildConstellationRDM_v002(dbPred, dbObs, finiteMask, ...
        Method=CONFIG.Method, NPerm=CONFIG.NPerm, Exact="off", ...
        UseParfor=true, Label=dsName);
    tElapsed = toc;

    perDataset(end+1) = struct("dataset", dsName, "nTrials", res.nTrials, ...
        "mantelR", res.mantelR, "pValue", res.pValue, ...
        "cvPred", res.cvPred, "cvObs", res.cvObs); %#ok<AGROW>

    fprintf("  %-15s n=%5d  mantelR=%+.4f  p=%.4f  (cvPred=%.3f cvObs=%.3f)  [%.1fs]\n\n", ...
        dsName, res.nTrials, res.mantelR, res.pValue, res.cvPred, res.cvObs, tElapsed);
end

%% -- Build pooled arrays -----------------------------------------------------
fprintf("=== Pooled cross-dataset RSA (Fraser included, Dhieb excluded, full N throughout) ===\n");
pooledPred = [];
pooledObs  = [];
pooledDatasetLabel = strings(0, 1);
for d = 1:numel(poolData)
    dsName = poolData(d).dataset;
    fm     = poolData(d).finiteMask;
    idx    = find(fm);
    fprintf("  %-15s full N=%d included (no cap)\n", dsName, numel(idx));
    pooledPred = [pooledPred; poolData(d).deltaBetaPred(idx, :)]; %#ok<AGROW>
    pooledObs  = [pooledObs;  poolData(d).deltaBetaObs(idx, :)];  %#ok<AGROW>
    pooledDatasetLabel = [pooledDatasetLabel; repmat(dsName, numel(idx), 1)]; %#ok<AGROW>
end
nPooled = size(pooledPred, 1);
fprintf("  Total pooled n = %d\n\n", nPooled);

if nPooled > 65535
    error("runConstellationRSA_HPC_v007:tooLargeForUint16", "%s", ...
        sprintf("nPooled=%d exceeds uint16 range (65535).", nPooled));
end

%% -- Smaller, fresh pool for the pooled step ONLY ---------------------------
fprintf("=== Starting SEPARATE %d-worker pool for the pooled step ===\n", CONFIG.PooledWorkers);
delete(gcp("nocreate"));
tic;
pool = parpool("local", CONFIG.PooledWorkers);
fprintf("  Pooled-step pool started: %d workers [%.1fs]\n\n", pool.NumWorkers, toc);

%% -- Build pooled RDMs (full double precision, needed for consistency check)
fprintf("Building pooled RDMs (double precision)...\n");
tic;
RDM_pred = squareEuclideanRDM_v007_local(pooledPred);
RDM_obs  = squareEuclideanRDM_v007_local(pooledObs);
fprintf("  Both RDMs built: %.2fs\n\n", toc);

triIdx = triu(true(nPooled), 1);
cvPred = std(RDM_pred(triIdx)) / mean(RDM_pred(triIdx));
cvObs  = std(RDM_obs(triIdx)) / mean(RDM_obs(triIdx));

[mantelR_toolbox, ~, ~] = mantelTest(RDM_pred, RDM_obs, ...
    Method=CONFIG.Method, NPerm=1, Exact="off");
fprintf("Reference mantelR (mantelTest, serial, NPerm=1, double): %.6f\n\n", mantelR_toolbox);

lowerMask = tril(true(nPooled), -1);
[iLowRaw, jLowRaw] = find(lowerMask);
v1 = RDM_pred(lowerMask);
v2 = RDM_obs(lowerMask);
v1rank_double = rankVector_v007_local(v1);
C0 = corrcoef(v1rank_double, rankVector_v007_local(v2));
mantelR = C0(1, 2);
absR = abs(mantelR);
fprintf("Observed mantelR (this script's own computation, double): %.6f\n", mantelR);

if abs(mantelR - mantelR_toolbox) > 1e-9
    error("runConstellationRSA_HPC_v007:consistencyCheckFailed", "%s", ...
        sprintf("Observed mantelR (%.10f) disagrees with mantelTest's own " + ...
        "serial mantelR (%.10f) by more than 1e-9.", mantelR, mantelR_toolbox));
end
fprintf("Consistency check passed (double precision): |%.10f - %.10f| < 1e-9\n\n", ...
    mantelR, mantelR_toolbox);

%% -- Build shrunk broadcast copies, THEN ACTUALLY CLEAR THE DOUBLES ---------
fprintf("=== Shrinking broadcast copies, then clearing doubles ===\n");
iLow = uint16(iLowRaw);
jLow = uint16(jLowRaw);
RDM_obs_bcast = single(RDM_obs);
v1rank = single(v1rank_double);
clear RDM_pred RDM_obs v1 v2 v1rank_double iLowRaw jLowRaw
bPerWorker = (numel(iLow)*2 + numel(jLow)*2 + numel(RDM_obs_bcast)*4 + numel(v1rank)*4) / 1e6;
fprintf("  Doubles cleared. Estimated broadcast per worker: %.1f MB (x %d workers = %.1f GB total)\n\n", ...
    bPerWorker, CONFIG.PooledWorkers, bPerWorker*CONFIG.PooledWorkers/1000);

%% -- Pooled permutation loop, 16 workers, live DataQueue progress -----------
fprintf("=== Pooled permutation loop, NPerm=%d, %d workers, live progress every %d ===\n", ...
    CONFIG.NPerm, CONFIG.PooledWorkers, CONFIG.ProgressEvery);

q = parallel.pool.DataQueue;
startTime = tic;
progEvery = CONFIG.ProgressEvery;
afterEach(q, @(k) reportProgress_v007_local(k, progEvery, startTime));

fprintf("Entering pooled parfor now -- WATCH FOR LIVE OUTPUT:\n");
NPerm = CONFIG.NPerm;
hits = zeros(NPerm, 1);
parfor k = 1:NPerm
    permIdx = uint16(randperm(nPooled)); %#ok<PFBNS>
    linIdx = sub2ind([nPooled, nPooled], double(permIdx(iLow)), double(permIdx(jLow))); %#ok<PFBNS>
    v2perm = RDM_obs_bcast(linIdx); %#ok<PFBNS>
    if std(double(v2perm)) == 0
        hits(k) = 0;
    else
        v2p = rankVector_v007_local(double(v2perm));
        Ck = corrcoef(double(v1rank), v2p); %#ok<PFBNS>
        hits(k) = abs(Ck(1,2)) >= absR;
    end
    send(q, k);
end
tPooled = toc(startTime);
fprintf("ALL %d POOLED PERMUTATIONS DONE. Total: %.2fs\n", NPerm, tPooled);

nGE = sum(hits);
pValue = (nGE + 1) / (NPerm + 1);

resPooled = struct();
resPooled.label      = "Pooled (7 datasets, Fraser included, Dhieb excluded, full N)";
resPooled.nTrials     = nPooled;
resPooled.mantelR     = mantelR;
resPooled.pValue      = pValue;
resPooled.nPermUsed   = NPerm;
resPooled.cvPred      = cvPred;
resPooled.cvObs       = cvObs;

fprintf("  Pooled: n=%d  mantelR=%+.4f  p=%.4f  (cvPred=%.3f cvObs=%.3f)  [%.1fs]\n", ...
    resPooled.nTrials, resPooled.mantelR, resPooled.pValue, ...
    resPooled.cvPred, resPooled.cvObs, tPooled);

%% -- Save --------------------------------------------------------------------
results = struct();
results.perDataset         = perDataset;
results.pooled              = resPooled;
results.pooledDatasetLabel  = pooledDatasetLabel;
results.config              = CONFIG;
results.runDate             = string(datetime("now"));
results.note = "HPC v007: NPerm raised to 10000 (from v006's confirmed-clean " + ...
    "100) -- this is the real, citable result. Same memory-safe structure as " + ...
    "v006 (16-worker pooled step, doubles cleared, no full RDM matrices saved " + ...
    "for the pooled result -- rebuild from constellationMetrics_v003.mat if " + ...
    "needed).";

if CONFIG.SaveMat
    save(CONFIG.OutMat, "results", "-v7.3");
    fprintf("\nSaved: %s\n", CONFIG.OutMat);
end

fprintf("Finished: %s\n", string(datetime("now")));
delete(gcp("nocreate"));
% diary closes automatically via diaryCleanup (onCleanup)

%% ============================ LOCAL HELPERS ================================

function D = squareEuclideanRDM_v007_local(X)
n = size(X, 1);
D = zeros(n, n);
for i = 1:n-1
    diffs = X(i+1:n, :) - X(i, :);
    d = sqrt(sum(diffs .^ 2, 2));
    D(i+1:n, i) = d;
    D(i, i+1:n) = d';
end
end

function r = rankVector_v007_local(x)
x = x(:);
n = numel(x);
[sortedX, sortIdx] = sort(x);
isNewGroup = [true; diff(sortedX) ~= 0];
groupID = cumsum(isNewGroup);
positions = (1:n)';
meanRankPerGroup = accumarray(groupID, positions) ./ accumarray(groupID, 1);
r = zeros(n, 1);
r(sortIdx) = meanRankPerGroup(groupID);
end

function reportProgress_v007_local(k, every, startTime)
if mod(k, every) == 0 || k == 1
    fprintf("  permutation %d done at t=%.2fs\n", k, toc(startTime));
end
end
