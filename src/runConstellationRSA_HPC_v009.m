%RUNCONSTELLATIONRSA_HPC_V009 Per-dataset and pooled cross-dataset RSA (Mantel), full-N, HPC, on the v015 basis.
%
% v009 changes from v009 (2026-10-28, Fraser D.S.; SESSION_LOG Addendum 57):
%   v009 read constellationMetrics_v003.mat (v007/v009 corpora, N_REPS = 20). v009 reads the saved v015 files
%   (loopClosureResults_<dataset>_all_shaped_xu_v015.mat, gated, N_REPS = 200) directly from the folder this file is in, and rebuilds the 15 predicted and
%   observed inter-pipeline contrasts per trial exactly as constellationMetrics_v005 / v004 do (predicted = median of betaRecSlice; a trial is finite if all six
%   observed and predicted betas and betaGenStarMed are finite). finiteMask is the v009 definition (finite first contrast in both), i.e. every inverted trial, not
%   only non-degenerate ones. Fail-loud checks: file presence, N_REPS = 200, N_BETA = 25, seed 1729. Nothing else changes: Method = spearman, NPerm = 10000, seed 1729,
%   the same six datasets (Dhieb excluded, EPP/measurement contamination), the same memory-safe pooled step. Output constellationRSA_HPC_v009.mat (v009 retained).
%   Run on BlueBEAR from the folder holding the v015 files and the mantel toolbox on the path; needs a walltime of at least 1.5 h (v009's arithmetic; Fraser's n is now about 2,738).
%
% Fraser, D.S. (2026) v009 (HPC variant, v015 basis)

%% -- CONFIG ------------------------------------------------------------
CONFIG = struct();
CONFIG.SrcDir                  = fileparts(mfilename("fullpath"));
CONFIG.DataDir     = CONFIG.SrcDir;   % folder holding loopClosureResults_<dataset>_all_shaped_xu_v015.mat
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
CONFIG.OutMat      = fullfile(CONFIG.SrcDir, "constellationRSA_HPC_v009.mat");
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
% mantelTest is only an independent cross-check of this script's own Mantel statistic (the toolbox was built from
% these bespoke functions). If it is absent the check is SKIPPED, loudly, and recorded in the result; it is not an error.
% The toolbox ships with the project: src/mantelToolbox (RDS copy) or src/req/dagmarfraser-mantel-*/toolbox (local).
mtDirs = [string(fullfile(CONFIG.SrcDir, "mantelToolbox")); string(fullfile({dir(fullfile(CONFIG.SrcDir, "req", "dagmarfraser-mantel-*", "toolbox")).folder}, {dir(fullfile(CONFIG.SrcDir, "req", "dagmarfraser-mantel-*", "toolbox")).name}))'];
for mtDir = mtDirs(:)'
    if exist("mantelTest", "file") ~= 2 && isfile(fullfile(mtDir, "mantelTest.m"))
        addpath(genpath(mtDir));  fprintf("mantel toolbox added from: %s\n", mtDir);
    end
end
haveMantelToolbox = exist("mantelTest", "file") == 2;
if ~haveMantelToolbox
    fprintf(2, "!!! WARNING: mantelTest.m not on path -- the consistency check against the mantel toolbox will be SKIPPED (recorded as mantelToolboxCheck = ""SKIPPED"").\n\n");
end
if exist("buildConstellationRDM_v002", "file") ~= 2
    error("runConstellationRSA_HPC_v009:builderNotFound", "%s", ...
        "buildConstellationRDM_v002.m not found on path (expected in src/).");
end

delete(gcp("nocreate"));
pool = parpool("local", CONFIG.NWorkers);
fprintf("Parpool started: %d workers (per-dataset loop)\n\n", pool.NumWorkers);

%% -- Load canonical per-trial constellation vectors -------------------------
allDatasets  = loadV015Constellation_v009_local(CONFIG.DataDir, CONFIG.Datasets);
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
        error("runConstellationRSA_HPC_v009:datasetNotFound", "%s", ...
            sprintf("Dataset '%s' not found among the v015 files loaded.", dsName));
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
    error("runConstellationRSA_HPC_v009:tooLargeForUint16", "%s", ...
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
RDM_pred = squareEuclideanRDM_v009_local(pooledPred);
RDM_obs  = squareEuclideanRDM_v009_local(pooledObs);
fprintf("  Both RDMs built: %.2fs\n\n", toc);

triIdx = triu(true(nPooled), 1);
cvPred = std(RDM_pred(triIdx)) / mean(RDM_pred(triIdx));
cvObs  = std(RDM_obs(triIdx)) / mean(RDM_obs(triIdx));

if haveMantelToolbox
    [mantelR_toolbox, ~, ~] = mantelTest(RDM_pred, RDM_obs, ...
        Method=CONFIG.Method, NPerm=1, Exact="off");
    fprintf("Reference mantelR (mantelTest, serial, NPerm=1, double): %.6f\n\n", mantelR_toolbox);
end

lowerMask = tril(true(nPooled), -1);
[iLowRaw, jLowRaw] = find(lowerMask);
v1 = RDM_pred(lowerMask);
v2 = RDM_obs(lowerMask);
v1rank_double = rankVector_v009_local(v1);
C0 = corrcoef(v1rank_double, rankVector_v009_local(v2));
mantelR = C0(1, 2);
absR = abs(mantelR);
fprintf("Observed mantelR (this script's own computation, double): %.6f\n", mantelR);

if haveMantelToolbox
    if abs(mantelR - mantelR_toolbox) > 1e-9
        error("runConstellationRSA_HPC_v009:consistencyCheckFailed", "%s", ...
            sprintf("Observed mantelR (%.10f) disagrees with mantelTest's own " + ...
            "serial mantelR (%.10f) by more than 1e-9.", mantelR, mantelR_toolbox));
    end
    fprintf("Consistency check passed (double precision): |%.10f - %.10f| < 1e-9\n\n", ...
        mantelR, mantelR_toolbox);
    mantelToolboxCheck = "PASSED";
else
    fprintf(2, "!!! Consistency check against mantelTest SKIPPED (toolbox not on path).\n\n");
    mantelToolboxCheck = "SKIPPED";
end

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
afterEach(q, @(k) reportProgress_v009_local(k, progEvery, startTime));

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
        v2p = rankVector_v009_local(double(v2perm));
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
resPooled.mantelToolboxCheck = mantelToolboxCheck;
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
results.note = "HPC v009: v015 basis (gated, N_REPS = 200), contrasts rebuilt from the saved v015 files; same method as v009 (Spearman, NPerm 10000, six datasets, Dhieb excluded). v009 result retained, not overwritten.";

if CONFIG.SaveMat
    save(CONFIG.OutMat, "results", "-v7.3");
    fprintf("\nSaved: %s\n", CONFIG.OutMat);
end

fprintf("Finished: %s\n", string(datetime("now")));
delete(gcp("nocreate"));
% diary closes automatically via diaryCleanup (onCleanup)

%% ============================ LOCAL HELPERS ================================

function D = squareEuclideanRDM_v009_local(X)
n = size(X, 1);
D = zeros(n, n);
for i = 1:n-1
    diffs = X(i+1:n, :) - X(i, :);
    d = sqrt(sum(diffs .^ 2, 2));
    D(i+1:n, i) = d;
    D(i, i+1:n) = d';
end
end

function r = rankVector_v009_local(x)
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

function reportProgress_v009_local(k, every, startTime)
if mod(k, every) == 0 || k == 1
    fprintf("  permutation %d done at t=%.2fs\n", k, toc(startTime));
end
end


function ds = loadV015Constellation_v009_local(dataDir, datasets)
%LOADV015CONSTELLATION_V009_LOCAL  Per-dataset 15-contrast predicted / observed matrices from the saved v015 shaped_xu files (as constellationMetrics_v004/v005).
    nPipe = 6;  [pP, pQ] = pairs_v009_local(nPipe);  ds = struct("dataset", {}, "deltaBetaPred", {}, "deltaBetaObs", {});
    for d = datasets
        f = fullfile(dataDir, "loopClosureResults_" + strrep(d, " ", "_") + "_all_shaped_xu_v015.mat");
        if ~isfile(f), error("runConstellationRSA_HPC_v009:missingMat", "No file: %s", f); end
        S0 = load(f, "results", "config");  c = S0.config;
        if c.N_REPS ~= 200 || c.N_BETA ~= 25 || c.RNG_SEED ~= 1729, error("runConstellationRSA_HPC_v009:config", "%s: N_REPS/N_BETA/seed differ from 200/25/1729", f); end
        r = S0.results;  N = numel(r);  dO = NaN(N, numel(pP));  dP = dO;
        for t = 1:N
            bo = r(t).betaObs(:)';  sl = r(t).betaRecSlice;
            if numel(bo) ~= nPipe || any(~isfinite(bo)) || isempty(sl) || size(sl, 2) ~= nPipe || ~isfinite(r(t).betaGenStarMed), continue; end
            bp = median(sl, 1, "omitnan");  if any(~isfinite(bp)), continue; end
            dO(t, :) = bo(pP) - bo(pQ);  dP(t, :) = bp(pP) - bp(pQ);
        end
        ds(end+1) = struct("dataset", d, "deltaBetaPred", dP, "deltaBetaObs", dO); %#ok<AGROW>
    end
end

function [pP, pQ] = pairs_v009_local(n)
    pP = [];  pQ = [];
    for p = 2:n, for q = 1:p-1, pP(end+1) = p; pQ(end+1) = q; end, end %#ok<AGROW>
end
