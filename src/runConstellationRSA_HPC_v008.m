%RUNCONSTELLATIONRSA_HPC_V008 Per-dataset and pooled cross-dataset RSA
% (Mantel), full-N, HPC version -- Pilot-retirement rerun.
%
% v008 changes from v007 (2026-08-26, Fraser D.S.):
%   v007's own header flagged (same session, portability-fix pass) that
%   CONFIG.Datasets pooled "Pilot" AND "Fraser" together, despite Pilot's
%   2026-08-07 retirement decision (Fraser supersedes Pilot project-wide,
%   is not meant to be pooled alongside it) -- that dataset list predated
%   the retirement decision reaching this file. This version is the actual
%   fix: "Pilot" removed from CONFIG.Datasets. Six datasets remain (Fraser,
%   Zarandi, Cook CTRL, Cook ASD, Hickman PLAC, Hickman HALO); Dhieb stays
%   excluded on its own separate, unrelated grounds (EPP contamination,
%   unchanged from v007/v006).
%
%   CONFIG.OutMat changed to its own new filename so v007's already-cited
%   Pilot-included result (Finding #138) is not overwritten -- that mat
%   remains on disk as the historical record, per this project's own
%   "dated correction notes, not silent rewrites" convention. Nothing else
%   changes: same NPerm=10000, same worker-pool structure, same
%   consistency check against mantelTest's own serial statistic. Local
%   helper functions renamed _v007_local -> _v008_local for clarity only,
%   identical logic.
%
%   Expected runtime: v007's own header logged ~30-35 min (per-dataset,
%   72 workers) + ~87 min (pooled, 16 workers) at nPooled=7893. Removing
%   Pilot (4237 of those 7893 trials) drops nPooled to ~3656 -- the pooled
%   step's cost scales with the lower-triangle size (~n^2/2), so that
%   stage alone should fall to very roughly ~19 min (3656/7893)^2 of
%   v007's 87 min. The per-dataset loop also loses its single largest job
%   (Pilot's own n=4237 RDM was bigger than Fraser's n=2729), so total
%   wall time is likely well under an hour -- this is arithmetic, not a
%   measurement, so still request a walltime buffer (recommend >=1.5h)
%   rather than trusting the estimate exactly.
%
% Everything else (relative-path CONFIG via fileparts(mfilename("fullpath")),
% no hardcoded MantelToolboxDir, PooledWorkers capped by feature("numcores"),
% 72-worker per-dataset pool / separate pooled-step pool, single/uint16
% broadcast shrink with doubles actually cleared, no full RDM matrices
% saved for the pooled result) is unchanged from v007 -- see that file's
% header for the full v003-v007 incident history that led here.
%
% Fraser, D.S. (2026) v008 (HPC variant, Pilot-retirement rerun)

%% -- CONFIG ------------------------------------------------------------
CONFIG = struct();
CONFIG.SrcDir                  = fileparts(mfilename("fullpath"));
CONFIG.ConstellationMetricsMat = fullfile(CONFIG.SrcDir, "constellationMetrics_v003.mat");
CONFIG.LogDir      = fullfile(CONFIG.SrcDir, "logs");
CONFIG.NoiseModel  = "shaped_xu";
CONFIG.Method      = "spearman";
CONFIG.NPerm       = 10000;
CONFIG.RngSeed     = 1729;
CONFIG.NWorkers       = feature("numcores");  % per-dataset loop
CONFIG.PooledWorkers  = min(16, feature("numcores"));  % pooled step ONLY
CONFIG.Datasets    = ["Fraser", "Zarandi", "Cook CTRL", "Cook ASD", ...
                      "Hickman PLAC", "Hickman HALO"];  % Pilot REMOVED 2026-08-26 (retired 2026-08-07, was erroneously still pooled in v007). Dhieb excluded (EPP contamination, unchanged from v007).
CONFIG.SaveMat     = true;
CONFIG.OutMat      = fullfile(CONFIG.SrcDir, "constellationRSA_HPC_v008.mat");
CONFIG.ProgressEvery = 500;

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
    error("runConstellationRSA_HPC_v008:mantelTestNotFound", "%s", ...
        "mantelTest.m not found on path -- addpath the mantel toolbox " + ...
        "(github.com/dagmarfraser/mantel) before running this script.");
end
if exist("buildConstellationRDM_v002", "file") ~= 2
    error("runConstellationRSA_HPC_v008:builderNotFound", "%s", ...
        "buildConstellationRDM_v002.m not found on path (expected in src/).");
end
if ~isfile(CONFIG.ConstellationMetricsMat)
    error("runConstellationRSA_HPC_v008:missingMat", "%s", ...
        sprintf("%s not found.", CONFIG.ConstellationMetricsMat));
end

delete(gcp("nocreate"));
pool = parpool("local", CONFIG.NWorkers);
fprintf("Parpool started: %d workers (per-dataset loop)\n\n", pool.NumWorkers);

%% -- Load canonical per-trial constellation vectors -------------------------
S = load(CONFIG.ConstellationMetricsMat, "results");
allDatasets  = S.results.byNoiseModel.(CONFIG.NoiseModel).perDataset;
datasetNames = string({allDatasets.dataset});

%% -- Per-dataset RSA (UNCHANGED method from v003-v007 -- proven working) ----
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
        error("runConstellationRSA_HPC_v008:datasetNotFound", "%s", ...
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
fprintf("=== Pooled cross-dataset RSA (Pilot removed, Dhieb excluded, full N throughout) ===\n");
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
    error("runConstellationRSA_HPC_v008:tooLargeForUint16", "%s", ...
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
RDM_pred = squareEuclideanRDM_v008_local(pooledPred);
RDM_obs  = squareEuclideanRDM_v008_local(pooledObs);
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
v1rank_double = rankVector_v008_local(v1);
C0 = corrcoef(v1rank_double, rankVector_v008_local(v2));
mantelR = C0(1, 2);
absR = abs(mantelR);
fprintf("Observed mantelR (this script's own computation, double): %.6f\n", mantelR);

if abs(mantelR - mantelR_toolbox) > 1e-9
    error("runConstellationRSA_HPC_v008:consistencyCheckFailed", "%s", ...
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

%% -- Pooled permutation loop, live DataQueue progress -----------------------
fprintf("=== Pooled permutation loop, NPerm=%d, %d workers, live progress every %d ===\n", ...
    CONFIG.NPerm, CONFIG.PooledWorkers, CONFIG.ProgressEvery);

q = parallel.pool.DataQueue;
startTime = tic;
progEvery = CONFIG.ProgressEvery;
afterEach(q, @(k) reportProgress_v008_local(k, progEvery, startTime));

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
        v2p = rankVector_v008_local(double(v2perm));
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
resPooled.label      = "Pooled (6 datasets, Pilot removed per its 2026-08-07 retirement, Dhieb excluded, full N)";
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
results.note = "HPC v008: Pilot removed from CONFIG.Datasets (retired 2026-08-07, " + ...
    "was erroneously still pooled alongside Fraser in v007 -- see v007's own " + ...
    "flagged header note). NPerm=10000 unchanged from v007. Same memory-safe " + ...
    "structure (separate pooled-step pool, doubles cleared, no full RDM " + ...
    "matrices saved for the pooled result -- rebuild from " + ...
    "constellationMetrics_v003.mat if needed). v007's own Pilot-included " + ...
    "result (Finding #138) is retained on disk, not overwritten -- this is " + ...
    "a new, separate finding, not a correction in place.";

if CONFIG.SaveMat
    save(CONFIG.OutMat, "results", "-v7.3");
    fprintf("\nSaved: %s\n", CONFIG.OutMat);
end

fprintf("Finished: %s\n", string(datetime("now")));
delete(gcp("nocreate"));
% diary closes automatically via diaryCleanup (onCleanup)

%% ============================ LOCAL HELPERS ================================

function D = squareEuclideanRDM_v008_local(X)
n = size(X, 1);
D = zeros(n, n);
for i = 1:n-1
    diffs = X(i+1:n, :) - X(i, :);
    d = sqrt(sum(diffs .^ 2, 2));
    D(i+1:n, i) = d;
    D(i, i+1:n) = d';
end
end

function r = rankVector_v008_local(x)
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

function reportProgress_v008_local(k, every, startTime)
if mod(k, every) == 0 || k == 1
    fprintf("  permutation %d done at t=%.2fs\n", k, toc(startTime));
end
end
