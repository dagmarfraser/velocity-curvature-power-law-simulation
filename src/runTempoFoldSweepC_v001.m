%% runTempoFoldSweepC_v001.m
% SPEC_TempoFoldSweep_v001 Block C (+ optional Block D). Gates V1-V5 passed in
% runTempoFoldSweep_v001 (2026-09-26).
%   C  each dataset's own geometry, fs and v004 noise centroid (Xu fGn, as v058), at
%      fixed tempo: does noise create a high-beta fold where sigma = 0 has none (P3)?
%      Each noisy condition has a sigma = 0 companion built identically.
%   D  bridge to Finding #2: a/b 2.16, fs 120, f0 1 Hz, alpha {0, 3}, sigma {0.1 0.5 2} mm.
% Fold classification is the production one: findBothBranches_v008 on the median curve
% with SD across replicates, runner default monoParams (runLoopClosureFftnoise_v012 L72-75).
% Writes results/tempoFoldSweepC_v001.mat. Intended for BlueBEAR (all cores).
% After patching tempoFoldEngine_v001: delete(gcp('nocreate')) before re-running.

%% CONFIG
ROOT      = fileparts(fileparts(mfilename("fullpath")));
OUT_MAT   = fullfile(ROOT, "results", "tempoFoldSweepC_v001.mat");
N_WORKERS = feature("numcores");
RNG_SEED  = 1729;
N_REPS    = 200;                       % matches v015 per-trial forward maps
% name, a (mm), b (mm), alpha, sigma (mm), fs: v015 median PCA axes; v004 centroids (EMPIRICAL_DATASETS)
DS = table(["Fraser"; "Cook_CTRL"; "Hickman_PLAC"], [28.768; 57.468; 56.749], [13.000; 27.126; 26.257], ...
           [4.289; 4.771; 5.34], [2.009; 8.15; 7.17], [240; 133; 133], ...
           'VariableNames', ["name" "a" "b" "alpha" "sigma" "fs"]);
F0_C   = [0.25 0.5 1 1.5 2.5 4];
BETA   = 0:1/60:0.75;
CLIP   = 20;  NCYC = 10;
RUN_D  = true;
D_AB = 2.16;  D_A = 50;  D_FS = 120;  D_F0 = 1;  D_ALPHA = [0 3];  D_SIGMA = [0.1 0.5 2];
MONO   = struct('slopeTol', 0.05, 'smoothWidth', 0.0625, 'minSegWidth', 0.0625, 'trimToleranceK', 1.0);
LABELS = ["BWFD-OLS" "BWFD-LMLS" "BWFD-IRLS" "SG-OLS" "SG-LMLS" "SG-IRLS"];   % engine order

%% Setup
addpath(fullfile(ROOT, "src"));  addpath(genpath(fullfile(ROOT, "src", "functions")));  addpath(genpath(fullfile(ROOT, "src", "req")));
if ~isfolder(fileparts(OUT_MAT)), error("sweepC:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
fprintf("Workers: %d   reps: %d\n", N_WORKERS, N_REPS);

%% Jobs
J = table();
for d = 1:height(DS)
    [f0, bg] = ndgrid(F0_C, BETA);
    for s = [DS.sigma(d) 0]                                     % noisy + sigma = 0 companion
        J = [J; jobs_local("C", DS.name(d), DS.a(d), DS.b(d), f0(:), DS.fs(d), bg(:), ...
                           s, DS.alpha(d), ifelse_local(s > 0, N_REPS, 1), CLIP, NCYC)]; %#ok<AGROW>
    end
end
if RUN_D
    [al, sg, bg] = ndgrid(D_ALPHA, [D_SIGMA 0], BETA);
    J = [J; jobs_local("D", "ab2.16", D_A, D_A/D_AB, D_F0, D_FS, bg(:), sg(:), al(:), ...
                       N_REPS * (sg(:) > 0) + (sg(:) == 0), CLIP, NCYC)];
end
J.jobID = (1:height(J))';
fprintf("Jobs: %d (trajectories: %d)\n", height(J), sum(J.nReps));

%% Run
if N_WORKERS > 0 && isempty(gcp("nocreate")), parpool(N_WORKERS); end
tic;  R = tempoFoldEngine_v001(J, struct("nWorkers", N_WORKERS, "rngSeed", RNG_SEED));
fprintf("Engine: %.1f min\n", toc/60);

%% Failures (fail loud)
fprintf("\nLength-excluded jobs: %d\n", sum(R.lengthExcluded));
nf = sum(R.nFail, 1);
fprintf("Pipeline failures (replicate-level): %s\n", strjoin(compose("%s %d", LABELS', nf'), ", "));
for p = find(nf > 0)
    m = unique(R.firstErr(R.firstErr(:, p) ~= "", p));
    fprintf("  %s first messages: %s\n", LABELS(p), strjoin(m(1:min(3, end)), " | "));
end

%% Classify every curve with the production monotonicity function
R.cond = R.dataset + "_a" + compose("%.2f", R.alpha) + "_s" + compose("%.3g", R.sigma);
S = classify_local(R(~R.lengthExcluded, :), ["block" "dataset" "cond" "sigma" "alpha" "f0"], BETA, MONO, LABELS);

%% P3: does noise create a fold where sigma = 0 has none? (Block C)
C = S(S.block == "C", :);
Z = C(C.sigma == 0, ["dataset" "f0" "pipeline" "hasDesc" "riseCov"]);
N = C(C.sigma > 0, :);
P3 = innerjoin(N, renamevars(Z, ["hasDesc" "riseCov"], ["hasDesc0" "riseCov0"]), "Keys", ["dataset" "f0" "pipeline"]);
fprintf("\nP3 (Block C): curves with a descending run, noise vs sigma = 0, by dataset x f0 (of 6 pipelines)\n");
P3.nWithNoise = double(P3.hasDesc);  P3.nAtZero = double(P3.hasDesc0);
disp(groupsummary(P3, ["dataset" "f0"], "sum", ["nWithNoise" "nAtZero"]))
fprintf("Noise-created folds (descending run with noise, none at sigma = 0): %d of %d curves\n", ...
    sum(P3.hasDesc & ~P3.hasDesc0), height(P3));
fprintf("\nSG-IRLS and BWFD-OLS with noise: floor = beta_rec(0), top = beta_rec(0.75), rising coverage, still rising at 2/3\n");
disp(N(ismember(N.pipeline, ["SG-IRLS" "BWFD-OLS"]), ["dataset" "f0" "pipeline" "floor" "top" "riseCov" "nRise" "nDesc" "risingAt23"]))

%% Block D: low-end dip on a fixed-tempo footing (Finding #2 bridge)
if RUN_D
    fprintf("\nBlock D (a/b 2.16, fs 120, f0 1 Hz): descending runs at the low end (beta_gen <= 0.2) and overall\n");
    disp(S(S.block == "D", ["alpha" "sigma" "pipeline" "floor" "top" "nRise" "nDesc" "lowDesc" "riseCov"]))
end

save(OUT_MAT, "R", "J", "S", "P3", "DS", "MONO", "LABELS", "RNG_SEED", "N_REPS", "-v7.3");
fprintf("Saved: %s\n", OUT_MAT);

%% =========================================================================
function J = jobs_local(block, dataset, a, b, f0, fs, betaGen, sigma, alpha, nReps, clip, nCyc)
    n = max([numel(f0) numel(betaGen) numel(sigma) numel(alpha) numel(nReps)]);
    c = @(v) repmat(v(:), n / numel(v), 1);
    J = table(repmat(string(block), n, 1), repmat(string(dataset), n, 1), repmat("ellipse", n, 1), ...
        c(a / b), c(a), c(b), c(f0), NaN(n, 1), c(fs), c(betaGen), c(sigma), c(alpha), c(nReps), ...
        c(clip), c(nCyc), repmat("runner", n, 1), zeros(n, 1), 'VariableNames', ["block" "dataset" ...
        "generator" "ab" "a" "b" "f0" "vgf" "FS" "betaGen" "sigma" "alpha" "nReps" "edgeClip" ...
        "nCycles" "seedMode" "jobID"]);
end

function v = ifelse_local(c, a, b)
    if c, v = a; else, v = b; end
end

function S = classify_local(Rb, keys, BETA, MONO, LABELS)
    % Per group x pipeline: production branch classification of the median curve.
    G = unique(Rb(:, keys), "rows");  S = table();
    for gi = 1:height(G)
        m = true(height(Rb), 1);
        for k = keys, m = m & Rb.(k) == G.(k)(gi); end
        c = sortrows(Rb(m, :), "betaGen");
        if height(c) ~= numel(BETA)
            error("sweepC:grid", "%s", sprintf("group %d has %d beta nodes, expected %d", gi, height(c), numel(BETA)));
        end
        br = c.betaRec;  if ~iscell(br), br = num2cell(br, 2); end
        for p = 1:numel(LABELS)
            reps = cellfun(@(x) x(:, p), br, "UniformOutput", false);
            med = cellfun(@(x) median(x, "omitnan"), reps)';
            sd  = cellfun(@(x) std(x, 0, "omitnan"), reps)';
            sd(~isfinite(sd)) = 0;                                  % single-rep (sigma = 0) rows
            ok = isfinite(med);
            if sum(ok) < 4, continue; end
            both = findBothBranches_v008(BETA(ok), med(ok), sd(ok), MONO);
            if ~isempty(both.rise) && ~isfield(both.rise, "x")
                error("sweepC:segFields", "findBothBranches_v008 segments lack field x");
            end
            span = max(BETA(ok)) - min(BETA(ok));
            riseCov = 0; lowDesc = false;
            for k = 1:numel(both.rise), riseCov = riseCov + (max(both.rise(k).x) - min(both.rise(k).x)) / span; end
            for k = 1:numel(both.desc), lowDesc = lowDesc || min(both.desc(k).x) <= 0.2; end
            [~, i23] = min(abs(BETA - 2/3));
            S = [S; [G(gi, :), table(LABELS(p), med(1), med(end), numel(both.rise), numel(both.desc), ...
                 ~isempty(both.desc), lowDesc, riseCov, med(i23 + 1) > med(i23), ...
                 'VariableNames', ["pipeline" "floor" "top" "nRise" "nDesc" "hasDesc" "lowDesc" ...
                 "riseCov" "risingAt23"])]]; %#ok<AGROW>
        end
    end
end
