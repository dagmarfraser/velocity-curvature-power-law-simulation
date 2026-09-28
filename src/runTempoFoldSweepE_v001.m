%% runTempoFoldSweepE_v001.m
% Block E: do the registered compression floors (manuscript §3a; attractorLocation_v001.m)
% survive when tempo is held fixed? Floor statistics exactly as attractorLocation_v001
% L295-302: location = mean beta_rec over the grid's 22 beta_gen values, spread = max - min,
% slope = polyfit(beta_gen, beta_rec, 1), at fs 120 Hz.
%   REPLAY  v058 generator at VGF node 8 (181.3), sigma in canvas px (x pixelScale 4.8,
%           Toolchain_caller_v058.m L213/L545), v058 settings: edge clip 50, truth seeds.
%           Gate: reproduces perCoordinateSEM_v2_001 meanBetaRec within replicate SE.
%   FIXED   exact ellipse at the grid shape's least-squares fit (235 x 109 px / 4.8 = 49.0 x 22.7 mm; #236),
%           same noise-to-size ratio, fixed f0 in F0_E, same settings and beta_gen values.
% Each has sigma = 0 companions. Analysis: analyseTempoFoldSweepE_v001.m (local).
% Writes results/tempoFoldSweepE_v001.mat. Intended for BlueBEAR (all cores).
% After patching tempoFoldEngine_v001: delete(gcp('nocreate')) before re-running.

%% CONFIG
ROOT        = fileparts(fileparts(mfilename("fullpath")));
SEM_MAT     = fullfile(ROOT, "src", "perCoordinateSEM_v2_001.mat");
OUT_MAT     = fullfile(ROOT, "results", "tempoFoldSweepE_v001.mat");
N_WORKERS   = feature("numcores");
RNG_SEED    = 1729;
N_REPS      = 100;
FS          = 120;
VGF_IDX     = 8;                     % attractorLocation_v001 REF_VGF_I
ALPHAS      = [0 1 2 3 5];
SIGMAS_MM   = [10 20];
F0_E        = [0.25 0.5 1 2.5 4];
PIXEL_SCALE = 480 / 100;             % Toolchain_caller_v058.m L213
GRID_A_PX   = 235;  GRID_B_PX = 109; % grid ellipse fit, canvas px (extents 240 x 109.85; #236)
CLIP        = 50;  NCYC = 10;        % v058 settings (Toolchain_caller_v058.m L203, L207)

%% Setup
addpath(fullfile(ROOT, "src"));  addpath(genpath(fullfile(ROOT, "src", "functions")));  addpath(genpath(fullfile(ROOT, "src", "req")));
if ~isfile(SEM_MAT), error("sweepE:input", "%s", "Missing: " + SEM_MAT); end
if ~isfolder(fileparts(OUT_MAT)), error("sweepE:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
S = load(SEM_MAT, "coordTable");  T = S.coordTable;
for v = ["betaGen" "VGF" "fs"], T.(v) = double(T.(v)); end
vg = sort(unique(T.VGF));  VGF0 = vg(VGF_IDX);
BETA = sort(unique(T.betaGen));
aMM = GRID_A_PX / PIXEL_SCALE;  bMM = GRID_B_PX / PIXEL_SCALE;
fprintf("Workers %d, reps %d, VGF node %d = %.3f, %d beta_gen values, ellipse %.2f x %.2f mm\n", ...
    N_WORKERS, N_REPS, VGF_IDX, VGF0, numel(BETA), aMM, bMM);

%% Jobs
J = table();
[al, sg, bg] = ndgrid(ALPHAS, SIGMAS_MM, BETA);                          % replay, noisy
J = [J; jobs_local("REPLAY", "v058", NaN, NaN, NaN, VGF0, FS, bg(:), sg(:) * PIXEL_SCALE, sg(:), al(:), N_REPS, CLIP, NCYC)];
J = [J; jobs_local("REPLAY", "v058", NaN, NaN, NaN, VGF0, FS, BETA, 0, 0, 0, 1, CLIP, NCYC)];   % sigma = 0
[al, sg, f0, bg] = ndgrid(ALPHAS, SIGMAS_MM, F0_E, BETA);                % fixed tempo, noisy
J = [J; jobs_local("FIXED", "ellipse", aMM, bMM, f0(:), NaN, FS, bg(:), sg(:), sg(:), al(:), N_REPS, CLIP, NCYC)];
[f0, bg] = ndgrid(F0_E, BETA);
J = [J; jobs_local("FIXED", "ellipse", aMM, bMM, f0(:), NaN, FS, bg(:), 0, 0, 0, 1, CLIP, NCYC)];  % sigma = 0
J.jobID = (1:height(J))';
fprintf("Jobs: %d (trajectories: %d)\n", height(J), sum(J.nReps));

%% Run
if N_WORKERS > 0 && isempty(gcp("nocreate")), parpool(N_WORKERS); end
tic;  R = tempoFoldEngine_v001(J, struct("nWorkers", N_WORKERS, "rngSeed", RNG_SEED));
fprintf("Engine: %.1f min\n", toc/60);
fprintf("Length-excluded jobs: %d;  replicate-level failures: %s\n", sum(R.lengthExcluded), ...
    strjoin(compose("%s %d", R.pipelines(1, :)', sum(R.nFail, 1)'), ", "));

save(OUT_MAT, "R", "J", "VGF0", "BETA", "ALPHAS", "SIGMAS_MM", "F0_E", "PIXEL_SCALE", "aMM", "bMM", ...
     "CLIP", "NCYC", "FS", "N_REPS", "RNG_SEED", "-v7.3");
fprintf("Saved: %s\n", OUT_MAT);

%% =========================================================================
function J = jobs_local(block, gen, a, b, f0, vgf, fs, betaGen, sigma, sigmaMM, alpha, nReps, clip, nCyc)
    n = max([numel(f0) numel(betaGen) numel(sigma) numel(alpha)]);
    c = @(v) repmat(v(:), n / numel(v), 1);
    J = table(repmat(string(block), n, 1), repmat(string(gen), n, 1), c(a / b), c(a), c(b), c(f0), ...
        c(vgf), c(fs), c(betaGen), c(sigma), c(sigmaMM), c(alpha), c(nReps), c(clip), c(nCyc), ...
        repmat("truth", n, 1), zeros(n, 1), 'VariableNames', ["block" "generator" "ab" "a" "b" "f0" ...
        "vgf" "FS" "betaGen" "sigma" "sigmaMM" "alpha" "nReps" "edgeClip" "nCycles" "seedMode" "jobID"]);
end
