% buildMonotonicSegments_v2_002.m  Step 1 of TODO_CCC_Invertibility.md
%
% v2_002 vs v2_001: parfor parallelisation of both sweep loops.
%   Each iteration is independent (read-only grid access, pure findBothBranches).
%   Results collected in cell arrays then assembled serially into betaSeg/VGFSeg.
%   HPC detection mirrors Toolchain_caller_v058 / loopClosureMasterHPC_v001:
%     isHPC = SLURM_JOB_ID set OR numcores > 16 (catches interactive BlueBEAR).
%   HPC  : ProcessPool, all cores (~72 on BlueBEAR, ~3-5 min per sweep).
%   Local: ThreadPool,  numcores-1  (desktop stays responsive).
%
% Input:   src/perCoordinateSEM_v2_001.mat   (coordTable)
% Output:  src/monotonicSegments_v2_002.mat
%
% Fraser, D.S. (2026)

clearvars

srcDir  = fileparts(mfilename("fullpath"));
addpath(genpath(fullfile(srcDir, "functions")));

inFile  = fullfile(srcDir, "perCoordinateSEM_v2_001.mat");
outFile = fullfile(srcDir, "monotonicSegments_v2_002.mat");

%% Parameters
params = struct();
params.slopeTol      = 0.05;
params.smoothWindow  = 3;
params.minSegLength  = 3;
params.pipeOrder     = ["BWFD-OLS", "BWFD-LMLS", "BWFD-IRLS", ...
                        "SG-OLS",   "SG-LMLS",   "SG-IRLS"];
params.dateGenerated = datetime("now");
params.srcFile       = inFile;

%% Load
fprintf("=== buildMonotonicSegments_v2_002 ===\n");
if ~isfile(inFile)
    error("buildMonotonicSegments:MissingInput", "%s", ...
        "Input not found: " + inFile + newline + ...
        "Run computePerCoordinateSEM_v2_001 first.");
end
fprintf("Loading %s ...\n", inFile);
S = load(inFile, "coordTable");
T = S.coordTable;
for col = ["meanVGFrec", "meanBetaRec"]
    if ~ismember(col, T.Properties.VariableNames)
        error("buildMonotonicSegments:MissingColumn", "%s", ...
            "coordTable lacks " + col + ". Re-run computePerCoordinateSEM_v2_001.");
    end
end
fprintf("  %d coordinate rows.\n", height(T));

%% Grid axes
alphaGrid   = unique(T.alpha);
sigmaGrid   = unique(T.sigma);
fsGrid      = unique(T.fs);
betaGenGrid = unique(T.betaGen);
VGFGrid     = unique(T.VGF);
pipeOrder   = params.pipeOrder;

nA = numel(alphaGrid);  nS = numel(sigmaGrid);  nF = numel(fsGrid);
nB = numel(betaGenGrid); nV = numel(VGFGrid);   nP = numel(pipeOrder);

fprintf("Grid: alpha=%d  sigma=%d  fs=%d  betaGen=%d  VGF=%d  pipe=%d\n", ...
    nA, nS, nF, nB, nV, nP);

%% Fill 6D grids
pipeStr = string(T.pipeline);
[~, aIdx] = ismember(T.alpha,   alphaGrid);
[~, sIdx] = ismember(T.sigma,   sigmaGrid);
[~, fIdx] = ismember(T.fs,      fsGrid);
[~, bIdx] = ismember(T.betaGen, betaGenGrid);
[~, vIdx] = ismember(T.VGF,     VGFGrid);
[~, pIdx] = ismember(pipeStr,   pipeOrder);

valid = aIdx>0 & sIdx>0 & fIdx>0 & bIdx>0 & vIdx>0 & pIdx>0;
if any(~valid)
    warning("buildMonotonicSegments:UnknownPipeline", "%s", ...
        sprintf("%d rows skipped (unknown pipeline): %s", ...
            sum(~valid), strjoin(unique(pipeStr(~valid)), ", ")));
end

gridSize = [nA, nS, nF, nB, nV, nP];
betaGrid = NaN(gridSize);
VGFgrid  = NaN(gridSize);
linIdx   = sub2ind(gridSize, ...
    aIdx(valid), sIdx(valid), fIdx(valid), bIdx(valid), vIdx(valid), pIdx(valid));
betaGrid(linIdx) = T.meanBetaRec(valid);
VGFgrid(linIdx)  = T.meanVGFrec(valid);
fprintf("  beta_rec grid filled: %.2f%%\n", 100*nnz(~isnan(betaGrid))/numel(betaGrid));
fprintf("  VGF_rec  grid filled: %.2f%%\n", 100*nnz(~isnan(VGFgrid))/numel(VGFgrid));

%% Parallel pool setup
isHPC     = ~isempty(getenv("SLURM_JOB_ID")) || feature("numcores") > 16;
nCores    = feature("numcores");
if isHPC
    nWorkers = nCores;
    poolType = "Processes";
    fprintf("\nHPC detected (%d cores) — ProcessPool, %d workers.\n", nCores, nWorkers);
else
    nWorkers = max(1, nCores - 1);
    poolType = "Threads";
    fprintf("\nLocal machine (%d cores) — ThreadPool, %d workers.\n", nCores, nWorkers);
end

existingPool = gcp("nocreate");
if isempty(existingPool)
    parpool(poolType, nWorkers);
elseif existingPool.NumWorkers ~= nWorkers
    delete(existingPool);
    parpool(poolType, nWorkers);
else
    fprintf("Reusing existing pool (%d workers).\n", existingPool.NumWorkers);
end

%% Beta sweep  (P x A x S x F x V)  ----------------------------------------
fprintf("\n--- beta sweep: %d slices ---\n", nP*nA*nS*nF*nV);
betaSegDims  = [nP, nA, nS, nF, nV];
totalBeta    = prod(betaSegDims);

% Collect per-slice results in parallel
rawBeta = cell(totalBeta, 1);
dq      = parallel.pool.DataQueue;
count   = 0;                                                    %#ok<NASGU>
step    = max(1, round(totalBeta / 20));
afterEach(dq, @(~) reportProgress());

t0 = tic;
parfor flat = 1:totalBeta
    [p, a, s, f, v] = ind2sub(betaSegDims, flat);
    yRaw = squeeze(betaGrid(a, s, f, :, v, p));                 % [nB x 1]
    rawBeta{flat} = findBothBranches(betaGenGrid, yRaw, params);
    send(dq, flat);
end
fprintf("\n  beta sweep: %.1f s\n", toc(t0));

% Serial assembly into struct
betaSeg = initSegStore(betaSegDims);
for flat = 1:totalBeta
    [p, a, s, f, v] = ind2sub(betaSegDims, flat);
    betaSeg = storeSegment(betaSeg, [p, a, s, f, v], rawBeta{flat});
end
clear rawBeta

%% VGF sweep   (P x A x S x F x B)  -----------------------------------------
fprintf("\n--- VGF sweep: %d slices ---\n", nP*nA*nS*nF*nB);
VGFSegDims  = [nP, nA, nS, nF, nB];
totalVGF    = prod(VGFSegDims);

rawVGF = cell(totalVGF, 1);
dq2    = parallel.pool.DataQueue;
afterEach(dq2, @(~) reportProgress());
count  = 0;                                                     %#ok<NASGU>

t0 = tic;
parfor flat = 1:totalVGF
    [p, a, s, f, b] = ind2sub(VGFSegDims, flat);
    yRaw = squeeze(VGFgrid(a, s, f, b, :, p));                  % [nV x 1]
    rawVGF{flat} = findBothBranches(VGFGrid, yRaw, params);
    send(dq2, flat);
end
fprintf("\n  VGF  sweep: %.1f s\n", toc(t0));

VGFSeg = initSegStore(VGFSegDims);
for flat = 1:totalVGF
    [p, a, s, f, b] = ind2sub(VGFSegDims, flat);
    VGFSeg = storeSegment(VGFSeg, [p, a, s, f, b], rawVGF{flat});
end
clear rawVGF

%% Summary
printInvertibilitySummary("beta", "rising",     betaSeg.rise.invertible,  pipeOrder, betaSegDims);
printInvertibilitySummary("beta", "descending", betaSeg.desc.invertible,  pipeOrder, betaSegDims);
printInvertibilitySummary("beta", "ambiguous",  betaSeg.ambiguous,        pipeOrder, betaSegDims);
printInvertibilitySummary("VGF",  "rising",     VGFSeg.rise.invertible,   pipeOrder, VGFSegDims);
printInvertibilitySummary("VGF",  "descending", VGFSeg.desc.invertible,   pipeOrder, VGFSegDims);
printInvertibilitySummary("VGF",  "ambiguous",  VGFSeg.ambiguous,         pipeOrder, VGFSegDims);

%% Save
fprintf("\nSaving %s ...\n", outFile);
save(outFile, "params", "pipeOrder", ...
    "alphaGrid", "sigmaGrid", "fsGrid", "betaGenGrid", "VGFGrid", ...
    "betaSeg", "VGFSeg", "-v7");
d = dir(outFile);
fprintf("  %.1f MB\n", d.bytes/2^20);
fprintf("=== buildMonotonicSegments_v2_002 COMPLETE ===\n");


%% ===========================  local functions  ============================

function reportProgress()
    % Called from DataQueue afterEach — prints a dot every ~5% of slices.
    % Uses a persistent counter shared across calls in the same session.
    persistent n step_
    if isempty(n),     n = 0;     end
    if isempty(step_), step_ = 1; end
    n = n + 1;
    if mod(n, step_) == 0
        fprintf(".");
    end
end

function printInvertibilitySummary(param, branch, invertMat, pipeOrder, dims)
    fprintf("\n--- %s %s invertibility ---\n", param, branch);
    for p = 1:numel(pipeOrder)
        slice = squeeze(invertMat(p, :, :, :, :));
        k = nnz(slice);
        n = numel(slice);
        fprintf("  %-12s: %7d / %7d  (%5.1f%%)\n", pipeOrder(p), k, n, 100*k/n);
    end
end

function seg = initSegStore(dims)
    seg.rise      = initBranchStore(dims);
    seg.desc      = initBranchStore(dims);
    seg.ambiguous = false(dims);
end

function b = initBranchStore(dims)
    b.segXcell   = cell(dims);
    b.segYcell   = cell(dims);
    b.segRange   = NaN([dims, 2]);
    b.segImage   = NaN([dims, 2]);
    b.nPoints    = zeros(dims);
    b.meanSlope  = NaN(dims);
    b.invertible = false(dims);
end

function store = storeSegment(store, idxVec, both)
    store.rise = storeBranch(store.rise, idxVec, both.rise);
    store.desc = storeBranch(store.desc, idxVec, both.desc);
    subs = num2cell(idxVec);
    store.ambiguous(subs{:}) = both.rise.invertible && both.desc.invertible;
end

function branch = storeBranch(branch, idxVec, seg)
    subs = num2cell(idxVec);
    branch.segXcell{subs{:}}   = seg.x;
    branch.segYcell{subs{:}}   = seg.y;
    branch.nPoints(subs{:})    = seg.n;
    branch.meanSlope(subs{:})  = seg.meanSlope;
    branch.invertible(subs{:}) = seg.invertible;
    if seg.n >= 1
        branch.segRange(subs{:}, 1) = min(seg.x);
        branch.segRange(subs{:}, 2) = max(seg.x);
        branch.segImage(subs{:}, 1) = min(seg.y);
        branch.segImage(subs{:}, 2) = max(seg.y);
    end
end

function both = findBothBranches(xGrid, yRaw, params)
    both.rise = findMonotonicRun(xGrid, yRaw, params, "rising");
    both.desc = findMonotonicRun(xGrid, yRaw, params, "descending");
end

function seg = findMonotonicRun(xGrid, yRaw, params, direction)
    seg = struct("x", [], "y", [], "n", 0, "meanSlope", NaN, "invertible", false);

    good = ~isnan(yRaw);
    if nnz(good) < params.minSegLength, return, end
    x = xGrid(good);
    y = yRaw(good);

    if params.smoothWindow > 1 && numel(y) >= params.smoothWindow
        ys = movmean(y, params.smoothWindow);
    else
        ys = y;
    end

    slope = diff(ys) ./ diff(x);
    switch direction
        case "rising",     isValid = slope >  params.slopeTol;
        case "descending", isValid = slope < -params.slopeTol;
        otherwise
            error("findMonotonicRun:BadDirection", "%s", ...
                "direction must be 'rising' or 'descending', got: " + direction);
    end

    [runStart, runEnd, runLen] = longestTrueRun(isValid);
    if isempty(runStart) || runLen < params.minSegLength - 1, return, end

    sx = x(runStart:runEnd+1);
    sy = y(runStart:runEnd+1);

    % Trim to strictly-monotonic prefix on raw y
    if direction == "rising"
        bad = find(diff(sy) <= 0, 1, "first");
    else
        bad = find(diff(sy) >= 0, 1, "first");
    end
    if ~isempty(bad)
        sx = sx(1:bad);
        sy = sy(1:bad);
    end
    if numel(sy) < params.minSegLength, return, end

    % Flip descending branch so segY is monotonically increasing for interp1
    if direction == "descending"
        sx = flipud(sx(:));
        sy = flipud(sy(:));
    else
        sx = sx(:);
        sy = sy(:);
    end

    seg.x          = sx;
    seg.y          = sy;
    seg.n          = numel(sy);
    seg.meanSlope  = mean(diff(sy) ./ diff(sx));
    seg.invertible = seg.n >= params.minSegLength && abs(seg.meanSlope) > params.slopeTol;
end

function [runStart, runEnd, runLen] = longestTrueRun(v)
    runStart = []; runEnd = []; runLen = 0;
    d      = diff([false; v(:); false]);
    starts = find(d ==  1);
    ends   = find(d == -1) - 1;
    if isempty(starts), return, end
    lens = ends - starts + 1;
    [runLen, best] = max(lens);
    runStart = starts(best);
    runEnd   = ends(best);
end
