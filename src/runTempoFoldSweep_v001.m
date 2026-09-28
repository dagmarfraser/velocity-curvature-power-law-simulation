%% runTempoFoldSweep_v001.m
% SPEC_TempoFoldSweep_v001, Blocks A and B with gates. Block C is a separate caller,
% run only after these gates pass.
%   A   sigma = 0, fixed tempo: eccentricity x f0 x beta_gen (+ fs, edge clip, nCycles sensitivity)
%   B1  v058's own generator at fixed VGF, v058 settings -> must reproduce the grid (V1a:
%       engine equivalence)
%   B2  analytic ellipse at each grid cell's REALISED tempo (from V0) -> does tempo alone
%       reproduce the grid's fold? (V1b)
%   V3  scale invariance at sigma = 0;  V4 analytic identity;  DET determinism
% Reads gridTempoAtFold_v003.mat (checkGridTempoAtFold_v003 output): results/ on RDS first,
% then the iMac mount; prints which. Writes results/tempoFoldSweep_v001.mat. No DB access;
% intended for BlueBEAR (uses every core of the node).
% After patching tempoFoldEngine_v001: delete(gcp('nocreate')) before re-running.

%% CONFIG
ROOT      = fileparts(fileparts(mfilename("fullpath")));
GRID_MAT_CANDS = [fullfile(ROOT, "results", "gridTempoAtFold_v003.mat")        % BlueBEAR (RDS)
                  "/Volumes/rdsprojects/f/fraserds-mpo-evaluation/2026_prereg/velocity-curvature-power-law-simulation-main/velocity-curvature-power-law-simulation-main/results/gridTempoAtFold_v003.mat"];  % iMac via mount
OUT_MAT   = fullfile(ROOT, "results", "tempoFoldSweep_v001.mat");
N_WORKERS = feature("numcores");   % all cores of the node (72 on BlueBEAR); 0 = serial
RNG_SEED  = 1729;
AB_LIST   = [1.25 1.7 2.16 3.0 4.1];      % a/b; 1.7 Dhieb, 2.16 grid ellipse fit (#236) and representative of Fraser/Cook/Hickman (2.12-2.21), 4.1 Zarandi
A_SIZE    = 50;                            % mm; sigma = 0 is scale-invariant (V3)
F0_LIST   = [0.25 0.5 0.75 1 1.5 2 2.5 3 4 6];
BETA_A    = 0:1/60:0.75;
FS_A = 120;  FS_EXTRA = [60 240];  CLIP_A = 20;  CLIP_SENS = 50;  NCYC_A = 10;  NCYC_SENS = [5 20];
VGF_B     = [90.017 164.02 330.3];         % nearest stored grid nodes are used
FS_B = 120;  CLIP_B = 50;  NCYC_B = 10;  A_GRID = 235;  B_GRID = 109;   % grid ellipse fit, px (#236)
V3_A      = 235;  V3_F0 = [0.5 2.5 4];  N_DET = 10;
FOLD_TOL  = 0.005;      % drop from peak to the last node that counts as a fold at sigma = 0
LABELS    = ["BWFD-OLS" "BWFD-LMLS" "BWFD-IRLS" "SG-OLS" "SG-LMLS" "SG-IRLS"];  % engine order

%% Setup
addpath(fullfile(ROOT, "src"));  addpath(genpath(fullfile(ROOT, "src", "functions")));  addpath(genpath(fullfile(ROOT, "src", "req")));
iG = find(isfile(GRID_MAT_CANDS), 1);
if isempty(iG), error("tempoSweep:gridMat", "%s", "V0 output not found at: " + strjoin(GRID_MAT_CANDS, " | ")); end
GRID_MAT = GRID_MAT_CANDS(iG);
fprintf("V0 grid file: %s\nWorkers: %d\n", GRID_MAT, N_WORKERS);
if ~isfolder(fileparts(OUT_MAT)), error("tempoSweep:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
G = load(GRID_MAT);  T = G.out.perCoordinate;
F = groupsummary(T(T.sampling_rate == FS_B, :), ["generated_beta" "vgf_value"], "mean", "f0");
vgAll = unique(F.vgf_value);
vgB = zeros(size(VGF_B));
for i = 1:numel(VGF_B), [~, ix] = min(abs(vgAll - VGF_B(i))); vgB(i) = vgAll(ix); end
betaGrid = unique(T.generated_beta);

%% Jobs
J = [makeJobs_local("A",        AB_LIST, A_SIZE, F0_LIST, FS_A,     BETA_A, CLIP_A,    NCYC_A)
     makeJobs_local("A_fs",     2.16,    A_SIZE, F0_LIST, FS_EXTRA, BETA_A, CLIP_A,    NCYC_A)
     makeJobs_local("A_clip50", 2.16,    A_SIZE, F0_LIST, FS_A,     BETA_A, CLIP_SENS, NCYC_A)
     makeJobs_local("A_ncyc",   2.16,    A_SIZE, F0_LIST, FS_A,     BETA_A, CLIP_A,    NCYC_SENS)
     makeJobs_local("V3",       2.16,    V3_A,   V3_F0,   FS_A,     BETA_A, CLIP_A,    NCYC_A)];
[vv, bb] = ndgrid(vgB, betaGrid);                                        % B1: v058 generator
J = [J; jobRows_local("B1", "v058", NaN, NaN, NaN, NaN, vv(:), FS_B, bb(:), CLIP_B, NCYC_B, "truth")];
F2 = F(ismember(F.vgf_value, vgB), :);                                   % B2: realised tempo
J = [J; jobRows_local("B2", "ellipse", A_GRID/B_GRID, A_GRID, B_GRID, F2.mean_f0, F2.vgf_value, ...
                      FS_B, F2.generated_beta, CLIP_B, NCYC_B, "truth")];
rng(RNG_SEED);  src = find(J.block == "A");  src = src(randperm(numel(src), N_DET));
D = J(src, :);  D.block(:) = "DET";  D.srcRow = src;  J.srcRow = NaN(height(J), 1);
J = [J; D];  J.jobID = (1:height(J))';
fprintf("Jobs: %d  (%s)\n", height(J), strjoin(compose("%s %d", string(unique(J.block)), ...
    groupcounts(J.block)), ", "));

%% Run
if N_WORKERS > 0 && isempty(gcp("nocreate")), parpool(N_WORKERS); end
tic;  R = tempoFoldEngine_v001(J, struct("nWorkers", N_WORKERS, "rngSeed", RNG_SEED));
fprintf("Engine: %.1f min\n", toc/60);

%% Failures (fail loud: reported, never absorbed)
fprintf("\nLength-excluded jobs: %d\n", sum(R.lengthExcluded));
nf = sum(R.nFail, 1);
fprintf("Pipeline failures: %s\n", strjoin(compose("%s %d", LABELS', nf'), ", "));
for p = find(nf > 0)
    m = unique(R.firstErr(R.firstErr(:, p) ~= "", p));
    fprintf("  %s first messages: %s\n", LABELS(p), strjoin(m(1:min(3, end)), " | "));
end

%% Gates
E = R(R.generator == "ellipse" & ~R.lengthExcluded & R.betaGen > 0, :);
fprintf("\nV4 analytic identity: max |beta_analytic - beta_gen| = %.3g (should be ~1e-15; LMLS/IRLS NaN at beta 0 excluded)\n", ...
    max(abs(E.betaAnalytic - E.betaGen), [], "all"));
Dr = R(R.block == "DET", :);
fprintf("DET determinism: max |diff| over %d repeated jobs = %.3g (should be 0)\n", height(Dr), ...
    max(abs(Dr.betaMean - R.betaMean(Dr.srcRow, :)), [], "all"));
K = ["f0" "betaGen"];
V3 = innerjoin(R(R.block == "V3", [K "betaMean"]), R(R.block == "A" & R.ab == 2.16, [K "betaMean"]), "Keys", K);
fprintf("V3 scale invariance (a = %g vs %g): max |diff| = %.3g over %d matched jobs\n", V3_A, A_SIZE, ...
    max(abs(V3.betaMean_left - V3.betaMean_right), [], "all"), height(V3));
cmpB1 = compareGrid_local(R(R.block == "B1", :), T, FS_B, LABELS);
cmpB2 = compareGrid_local(R(R.block == "B2", :), T, FS_B, LABELS);
fprintf("\nV1a  B1 (v058 generator, v058 settings) vs grid: should match to ~1e-6\n");  disp(cmpB1)
fprintf("V1b  B2 (analytic ellipse at realised grid tempo) vs grid: P4 predicts |diff| < 0.01\n");  disp(cmpB2)
foldB2 = foldSummary_local(R(R.block == "B2", :), "vgf", LABELS, FOLD_TOL);
fprintf("B2 fold by VGF node (compare grid: peak 0.43-0.63 and drops 0.05-0.22 for SG-IRLS):\n");
disp(foldB2(foldB2.pipeline == "SG-IRLS", :))

%% Block A readout: first f0 at which each (a/b, pipeline) folds, and where
blocks = ["A" "A_fs" "A_clip50" "A_ncyc"];
foldA = table();
for bk = blocks
    Rb = R(R.block == bk & ~R.lengthExcluded, :);
    Rb.cond = compose("ab%.2f_fs%d_clip%d_ncyc%d", Rb.ab, Rb.FS, Rb.edgeClip, Rb.nCycles);
    S = foldSummary_local(Rb, ["cond" "f0"], LABELS, FOLD_TOL);  S.block(:) = bk;
    foldA = [foldA; S]; %#ok<AGROW>
end
onset = groupsummary(foldA(foldA.isFold, :), ["block" "cond" "pipeline"], "min", "f0");
fprintf("\nBlock A: lowest f0 (Hz) at which the sigma = 0 map folds below beta_gen = 0.75 (absent = no fold up to %g Hz)\n", max(F0_LIST));
disp(unstack(onset(:, ["block" "cond" "pipeline" "min_f0"]), "min_f0", "pipeline"))

save(OUT_MAT, "R", "J", "cmpB1", "cmpB2", "foldB2", "foldA", "onset", "LABELS", "GRID_MAT", ...
     "FOLD_TOL", "RNG_SEED", "-v7.3");
fprintf("Saved: %s\n", OUT_MAT);

%% =========================================================================
function J = makeJobs_local(block, abList, aSize, f0List, fsList, betaList, clip, nCycList)
    [ab, f0, fs, bg, nc] = ndgrid(abList, f0List, fsList, betaList, nCycList);
    J = jobRows_local(block, "ellipse", ab(:), aSize, aSize ./ ab(:), f0(:), NaN, fs(:), bg(:), ...
                      clip, nc(:), "runner");
end

function J = jobRows_local(block, gen, ab, a, b, f0, vgf, FS, betaGen, clip, nCyc, seedMode)
    n = max([numel(ab) numel(f0) numel(vgf) numel(betaGen) numel(FS) numel(nCyc)]);
    c = @(v) repmat(v(:), n / numel(v), 1);          % expand scalars to n rows
    J = table(repmat(string(block), n, 1), repmat(string(gen), n, 1), c(ab), c(a), c(b), c(f0), ...
        c(vgf), c(FS), c(betaGen), zeros(n, 1), zeros(n, 1), ones(n, 1), c(clip), c(nCyc), ...
        repmat(string(seedMode), n, 1), zeros(n, 1), 'VariableNames', ["block" "generator" "ab" ...
        "a" "b" "f0" "vgf" "FS" "betaGen" "sigma" "alpha" "nReps" "edgeClip" "nCycles" "seedMode" "jobID"]);
end

function C = compareGrid_local(Rb, T, fs, LABELS)
    C = table();
    for p = 1:numel(LABELS)
        g = T(T.sampling_rate == fs & T.pipeline == LABELS(p), ["generated_beta" "vgf_value" "betaRec"]);
        r = table(Rb.betaGen, Rb.vgf, Rb.betaMean(:, p), 'VariableNames', ["generated_beta" "vgf_value" "eng"]);
        m = innerjoin(r(isfinite(r.eng), :), g, "Keys", ["generated_beta" "vgf_value"]);
        d = abs(m.eng - m.betaRec);
        C = [C; table(LABELS(p), height(m), max(d), median(d), 'VariableNames', ...
             ["pipeline" "nMatched" "maxAbsDiff" "medianAbsDiff"])]; %#ok<AGROW>
    end
end

function S = foldSummary_local(Rb, keys, LABELS, tol)
    % Per group x pipeline: peak location, drop from peak to the last node, fold flag,
    % and whether beta_rec is still rising at the node nearest beta_gen = 2/3.
    G = unique(Rb(:, keys), "rows");  S = table();
    for gi = 1:height(G)
        m = true(height(Rb), 1);
        for k = keys, m = m & Rb.(k) == G.(k)(gi); end
        c = sortrows(Rb(m, :), "betaGen");
        [~, i23] = min(abs(c.betaGen - 2/3));
        for p = 1:numel(LABELS)
            y = c.betaMean(:, p);  ok = isfinite(y);
            if sum(ok) < 3, continue; end
            b = c.betaGen(ok);  y = y(ok);
            [pk, ip] = max(y);
            rising = NaN;
            if i23 < height(c) && all(isfinite(c.betaMean([i23 i23+1], p)))
                rising = c.betaMean(i23+1, p) > c.betaMean(i23, p);
            end
            S = [S; [G(gi, :), table(LABELS(p), b(ip), pk, pk - y(end), pk - y(end) > tol, rising, ...
                 b(end), 'VariableNames', ["pipeline" "peakBeta" "peakRec" "drop" "isFold" ...
                 "risingAt23" "lastBeta"])]]; %#ok<AGROW>
        end
    end
end
