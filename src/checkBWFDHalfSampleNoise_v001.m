%% checkBWFDHalfSampleNoise_v001.m
% Does the BWFD pipeline's half-sample offset between velocity and acceleration matter under
% noise? differentiateKinematicsEBR case 2 stores a backward-difference velocity (L80, centred at
% i - 1/2) beside a centred second-difference acceleration (L84, centred at i), and curvature
% multiplies the two. At zero noise the offset shifts the gain onset (0.7% at 240 Hz, f_c 10 Hz;
% predictHarmonicGain_v002 P1) and leaves the fixed point at 1/3 unmoved. Under noise it also
% correlates the two estimates' errors on each axis (raw white-noise stencils: about -0.87; zero
% for a centred velocity), which could move the pull to the noise exponent.
% Three alignments, applied to identical noisy trajectories (paired):
%   REG  registered: differentiateKinematicsEBR case 2 as published (velocity at i - 1/2)
%   CEN  velocity averaged over i - 1/2 and i + 1/2 (= centred difference), acceleration as REG:
%        both at i
%   AVG  velocity as REG, acceleration averaged over i - 1 and i: both at i - 1/2
% CEN and AVG align the pair by averaging different quantities; if they agree, the effect is the
% alignment, not either stencil's noise gain. All three share one filtfilt stage (2nd-order
% Butterworth, f_c 10 Hz, zero lag, as the engine and the runner).
% Conditions: Block C of the tempo-fold sweep (runTempoFoldSweepC_v001 L20-25: each dataset's
% geometry, fs and v004 noise centroid, Xu fGn, fixed tempo) and a white-noise stress (Block D's
% geometry, a/b 2.16, a = 50 mm; alpha 0, sigma 0.5 mm) at 120 and 240 Hz. Each noisy condition
% has a sigma = 0 companion. Regressions OLS, LMLS, IRLS, runner seeding, limitBreak 0.
% Gain and target per median map: extractGainTarget_v002's definitions and reliability rules
% (its L33-39, L99-100, L206-230), copied here.
% Predictions, fixed before the run (2026-10-01), reported HELD / NOT HELD:
%  N1 null at empirical noise: in every Block C cell (dataset x f0 x beta_gen x regression),
%     |median paired delta| < MDC/2.77 for CEN and AVG, and |delta target| < N1_TARGET wherever
%     both targets are reliable.
%  N2 stress: the largest |median paired delta| under white noise exceeds the largest under
%     Block C noise.
%  N3 alignment: where |median delta_CEN| >= N3_MIN, delta_AVG has the same sign and
%     |delta_AVG - delta_CEN| < 0.5 |delta_CEN|. Not evaluable if no cell reaches N3_MIN.
% Anchors (error if broken): (1) the shared-filtfilt replica of case 2 equals
% differentiateKinematicsEBR exactly at every sampling rate used; (2) REG's sigma = 0 maps
% reproduce tempoFoldSweepC_v001's sigma = 0 companions (BWFD columns), within TOL_REG;
% (0, on resume only) two noisy jobs recomputed from scratch match the checkpoint within TOL_REG.
% TOL_REG is tiered because OLS is closed form while LMLS and IRLS stop at fitnlm's iteration
% tolerance: first run (2026-10-01), REG against Block C (BlueBEAR) on the lab PC, all 414 maps:
% max OLS 3.0e-15, LMLS 4.1e-11, IRLS 3.9e-9 (2 maps above 1e-9). A single 1e-9 tolerance
% stopped that run after the compute; the tiers sit two orders above those deviations and
% seven below MDC/2.77.
% Checkpoint: the raw run is saved to CKPT straight after the parfor, before any anchor. If CKPT
% exists and its job table matches this configuration exactly, the run is not repeated (anchor 0
% then re-derives two jobs). The first run's checkpoint was saved by hand from the workspace
% after anchor 2 stopped it (J only, no cfg); anchor 0 is what verifies it.
% About 22,000 noisy trajectories x 3 alignments x 3 regressions: roughly 30 min on 8 workers.
% Run from src/. Writes results/checkBWFDHalfSampleNoise_v001.mat,
% figures/checkBWFDHalfSampleNoise_v001.png

%% CONFIG
ROOT      = fileparts(fileparts(mfilename("fullpath")));
IN_C      = fullfile(ROOT, "results", "tempoFoldSweepC_v001.mat");
OUT_MAT   = fullfile(ROOT, "results", "checkBWFDHalfSampleNoise_v001.mat");
OUT_PNG   = fullfile(ROOT, "figures", "checkBWFDHalfSampleNoise_v001.png");
RNG_SEED  = 1729;
N_REPS    = 40;                                    % paired design: deltas need few reps
DS = table(["Fraser"; "Cook_CTRL"; "Hickman_PLAC"], [28.768; 57.468; 56.749], [13.000; 27.126; 26.257], ...
           [4.289; 4.771; 5.34], [2.009; 8.15; 7.17], [240; 133; 133], ...
           'VariableNames', ["name" "a" "b" "alpha" "sigma" "fs"]);   % as runTempoFoldSweepC_v001 L20-22
F0_C      = [0.25 0.5 1 1.5 2.5 4];                % as Block C
ST_A      = 50;  ST_AB = 2.16;  ST_FS = [120 240];  ST_F0 = [0.5 1 2];  ST_ALPHA = 0;  ST_SIGMA = 0.5;
BETA      = 0:1/30:0.75;
CLIP      = 20;  NCYC = 10;  FC = 10;  MIN_PTS = 10;
VARIANTS  = ["REG" "CEN" "AVG"];
REG       = [3 4 5];  REG_NAMES = ["OLS" "LMLS" "IRLS"];
FIT_WIN   = [0 0.75];  LOC_WIN = [0.2 0.5];        % as extractGainTarget_v002 L33-39
GAIN_MAX  = 0.9;  RMSE_MAX = 0.02;  AGREE_TOL = 0.02;  EXT = 0.05;  FAIL_MAX = 10;
N1_TARGET = 0.01;  N3_MIN = 0.005;
TOL_REG   = [1e-12 1e-9 1e-7];                      % OLS, LMLS, IRLS: anchors 0 and 2 (see header)
CKPT      = fullfile(ROOT, "results", "checkBWFDHalfSampleNoise_v001_runCheckpoint.mat");

%% Setup
addpath(fullfile(ROOT, "src"));  addpath(genpath(fullfile(ROOT, "src", "functions")));  addpath(genpath(fullfile(ROOT, "src", "req")));
if ~isfile(IN_C), error("halfSample:in", "%s", "Missing input: " + IN_C); end
for f = [OUT_MAT OUT_PNG], if ~isfolder(fileparts(f)), error("halfSample:outDir", "%s", "Missing folder: " + fileparts(f)); end, end
semAdequate = semAdequacyThreshold_v001;

%% Anchor 1: the shared-filtfilt replica equals differentiateKinematicsEBR case 2
for fs = unique([DS.fs; ST_FS(:)])'
    rs = RandStream("threefry4x64_20", "Seed", RNG_SEED);
    [x, y] = generatePowerLawEllipse_v001(ST_A, ST_A / ST_AB, 1, fs, 0.4, 'nCycles', NCYC);
    x = x(:) + 0.5 * randn(rs, numel(x), 1);  y = y(:) + 0.5 * randn(rs, numel(y), 1);
    [dx, dy] = differentiateKinematicsEBR(x, y, 2, [2 FC 1], fs);
    [vR, aR] = replica_local(x, y, fs, FC);
    N = numel(x);
    if ~isequal(vR(2:N, :), [dx(2:N, 2) dy(2:N, 2)]) || ~isequal(aR(2:N-1, :), [dx(2:N-1, 3) dy(2:N-1, 3)])
        error("halfSample:anchor1", "%s", sprintf("Replica differs from differentiateKinematicsEBR at fs %d Hz", fs));
    end
end
fprintf("Anchor 1 passed: replica equals differentiateKinematicsEBR case 2 exactly at fs %s Hz\n", mat2str(unique([DS.fs; ST_FS(:)])'));

%% Jobs
J = table();
for d = 1:height(DS)
    for s = [DS.sigma(d) 0]
        [f0, bg] = ndgrid(F0_C, BETA);
        J = [J; jobs_local("C", DS.name(d), DS.a(d), DS.b(d), DS.fs(d), f0(:), bg(:), s, DS.alpha(d), (s > 0) * N_REPS + (s == 0))]; %#ok<AGROW>
    end
end
for fs = ST_FS
    for s = [ST_SIGMA 0]
        [f0, bg] = ndgrid(ST_F0, BETA);
        J = [J; jobs_local("S", "white_fs" + fs, ST_A, ST_A / ST_AB, fs, f0(:), bg(:), s, ST_ALPHA, (s > 0) * N_REPS + (s == 0))]; %#ok<AGROW>
    end
end
J.jobID = (1:height(J))';
fprintf("Jobs: %d (noisy trajectories %d, sigma = 0 trajectories %d)\n", height(J), sum(J.nReps(J.sigma > 0)), sum(J.sigma == 0));

%% Run (paired: each noisy replicate goes through all three alignments), or resume
KEYS = ["block" "dataset" "a" "b" "FS" "f0" "betaGen" "sigma" "alpha" "nReps" "jobID"];
cfg = struct("seed", RNG_SEED, "clip", CLIP, "nCyc", NCYC, "fc", FC, "minPts", MIN_PTS, "reg", REG, "variants", VARIANTS);
resumed = isfile(CKPT);
if resumed
    S = load(CKPT, "J");
    if ~all(ismember(KEYS, string(S.J.Properties.VariableNames))) || ~isequaln(S.J(:, KEYS), J(:, KEYS))
        error("halfSample:ckpt", "%s", "Checkpoint jobs differ from this configuration; move it aside to rerun: " + CKPT);
    end
    J = S.J;
    fprintf("Resumed from checkpoint (run not repeated): %s\n", CKPT);
else
    Jc = table2struct(J);  out = cell(height(J), 1);
    nW = 0;  if ~isempty(gcp("nocreate")), nW = gcp("nocreate").NumWorkers; end
    tic;
    parfor (k = 1:height(J), nW)
        out{k} = runJob_local(Jc(k), cfg);
    end
    fprintf("Run: %.1f min\n", toc / 60);
    O = vertcat(out{:});
    J.lengthExcluded = [O.lengthExcluded]';  J.betaRec = {O.betaRec}';  J.nFail = vertcat(O.nFail);  J.firstErr = vertcat(O.firstErr);
    save(CKPT, "J", "cfg", "-v7.3");
    fprintf("Checkpoint written: %s\n", CKPT);
end

%% Anchor 0 (resume only): two noisy jobs recomputed from scratch match the checkpoint
if resumed
    chk = zeros(1, 2);  bl = ["C" "S"];
    for i = 1:2
        n = J(J.block == bl(i) & J.sigma > 0 & ~J.lengthExcluded, :);
        [~, w] = min(abs(n.f0 - max(n.f0)) + abs(n.betaGen - 1/3));  chk(i) = n.jobID(w);
    end
    for id = chk
        if J.jobID(id) ~= id, error("halfSample:jobID", "%s", "jobID is not the row index"); end
        o = runJob_local(table2struct(J(id, KEYS)), cfg);
        checkTiered_local(o.betaRec, J.betaRec{id}, TOL_REG, numel(VARIANTS), sprintf("anchor 0, job %d", id));
    end
    fprintf("Anchor 0 passed: jobs %s recomputed and match the checkpoint\n", mat2str(chk));
end
fprintf("Length-excluded jobs: %d\n", sum(J.lengthExcluded));
lab = reshape((VARIANTS' + "-" + REG_NAMES)', 1, []);           % betaRec column order: variant-major
nf = sum(J.nFail, 1);
fprintf("Replicate-level failures: %s\n", strjoin(arrayfun(@(i) sprintf("%s %d", lab(i), nf(i)), 1:numel(lab)), ", "));
for p = find(nf > 0)
    m = unique(J.firstErr(J.firstErr(:, p) ~= "", p));
    fprintf("  %s first messages: %s\n", lab(p), strjoin(m(1:min(3, end)), " | "));
end

%% Anchor 2: REG sigma = 0 maps reproduce Block C's sigma = 0 companions
Sc = load(IN_C, "R");  Rc = Sc.R(Sc.R.block == "C" & Sc.R.sigma == 0, :);
Z = J(J.block == "C" & J.sigma == 0 & ~J.lengthExcluded, :);
dev2 = NaN(height(Z), numel(REG));
for k = 1:height(Z)
    m = Rc.dataset == Z.dataset(k) & Rc.f0 == Z.f0(k) & round(Rc.betaGen * 60) == round(Z.betaGen(k) * 60);
    if nnz(m) ~= 1, error("halfSample:anchor2", "%s", sprintf("Block C has %d matches for job %d", nnz(m), Z.jobID(k))); end
    ref = Rc.betaRec{m};  ref = ref(1, 1:3);                       % engine order: BWFD-OLS, BWFD-LMLS, BWFD-IRLS
    got = Z.betaRec{k};  got = got(1, 1:3);                        % REG columns
    if ~isequal(isnan(ref), isnan(got))
        error("halfSample:anchor2", "%s", sprintf("Job %d: NaN pattern differs from Block C", Z.jobID(k)));
    end
    dev2(k, :) = abs(got - ref);
end
fprintf("Anchor 2: max |REG - Block C| over %d sigma = 0 maps: %s\n", height(Z), ...
    strjoin(arrayfun(@(q) sprintf("%s %.3g (tol %.0e)", REG_NAMES(q), max(dev2(:, q)), TOL_REG(q)), 1:numel(REG)), ", "));
if any(dev2 > TOL_REG, "all")
    [k, q] = find(dev2 > TOL_REG, 1);
    error("halfSample:anchor2", "%s", sprintf("Job %d %s deviates by %.3g (tol %.0e)", Z.jobID(k), REG_NAMES(q), dev2(k, q), TOL_REG(q)));
end
fprintf("Anchor 2 passed\n");

%% Paired deltas per cell (noisy jobs)
Nz = J(J.sigma > 0 & ~J.lengthExcluded, :);
D = table();
for k = 1:height(Nz)
    br = Nz.betaRec{k};
    for v = 2:numel(VARIANTS)
        for q = 1:numel(REG)
            d = br(:, col_local(v, q, numel(REG))) - br(:, col_local(1, q, numel(REG)));
            D = [D; table(Nz.block(k), Nz.dataset(k), Nz.FS(k), Nz.f0(k), Nz.betaGen(k), VARIANTS(v), REG_NAMES(q), ...
                median(d, "omitnan"), sum(isfinite(d)), 'VariableNames', ["block" "dataset" "FS" "f0" "betaGen" ...
                "variant" "regression" "medDelta" "nPaired"])]; %#ok<AGROW>
        end
    end
end
SD = groupsummary(D, ["block" "dataset" "f0" "variant" "regression"], {@(x) max(abs(x)), @(x) median(x, "omitnan")}, "medDelta");
SD = renamevars(SD, ["fun1_medDelta" "fun2_medDelta"], ["maxAbsMedDelta" "medianMedDelta"]);
fprintf("\nLargest |median paired delta| over beta_gen, per dataset x f0 x regression (MDC/2.77 = %.5f):\n", semAdequate);
disp(unstack(SD(:, ["block" "dataset" "f0" "regression" "variant" "maxAbsMedDelta"]), "maxAbsMedDelta", "variant"))

%% Gain and target of every median map (extractGainTarget_v002 definitions)
M = table();
G = unique(J(~J.lengthExcluded, ["block" "dataset" "FS" "sigma" "f0"]), "rows");
for g = 1:height(G)
    c = sortrows(J(J.block == G.block(g) & J.dataset == G.dataset(g) & J.sigma == G.sigma(g) & J.f0 == G.f0(g) & ...
        ~J.lengthExcluded, :), "betaGen");
    for v = 1:numel(VARIANTS)
        for q = 1:numel(REG)
            cl = col_local(v, q, numel(REG));
            med = cellfun(@(x) median(x(:, cl), "omitnan"), c.betaRec);
            nF = sum(cellfun(@(x) sum(~isfinite(x(:, cl))), c.betaRec));  failPct = 100 * nF / sum(c.nReps);
            r = gainTarget_local(c.betaGen, med, FIT_WIN, LOC_WIN);
            r.reliable = r.gain <= GAIN_MAX && r.rmse <= RMSE_MAX && r.nCross >= 1 && abs(r.targetLine - r.targetCross) <= AGREE_TOL && ...
                r.targetLine >= FIT_WIN(1) - EXT && r.targetLine <= FIT_WIN(2) + EXT && failPct <= FAIL_MAX;
            M = [M; [G(g, :), table(VARIANTS(v), REG_NAMES(q), failPct, 'VariableNames', ["variant" "regression" "failPct"]), ...
                struct2table(r)]]; %#ok<AGROW>
        end
    end
end
Mr = M(M.variant == "REG", ["block" "dataset" "sigma" "f0" "regression" "gain" "gainLocal" "targetLine" "reliable"]);
Mr = renamevars(Mr, ["gain" "gainLocal" "targetLine" "reliable"], ["gainREG" "gainLocalREG" "targetREG" "reliableREG"]);
MD = innerjoin(M(M.variant ~= "REG", :), Mr, "Keys", ["block" "dataset" "sigma" "f0" "regression"]);
MD.dGain = MD.gain - MD.gainREG;  MD.dGainLocal = MD.gainLocal - MD.gainLocalREG;
MD.dTarget = MD.targetLine - MD.targetREG;  MD.bothReliable = MD.reliable & MD.reliableREG;
MD.dTarget(~MD.bothReliable) = NaN;
fprintf("\nGain and target differences against REG (targets only where both are reliable):\n");
disp(MD(MD.sigma > 0, ["block" "dataset" "f0" "regression" "variant" "gainREG" "dGain" "dGainLocal" "targetREG" "dTarget" "failPct"]))
fprintf("Zero-noise companions, gain differences (the offset's effect without noise):\n");
disp(MD(MD.sigma == 0, ["block" "dataset" "f0" "regression" "variant" "gainREG" "dGain" "dGainLocal"]))

%% Predictions
P = struct();
dc = D(D.block == "C", :);  ds = D(D.block == "S", :);
tC = MD(MD.block == "C" & MD.sigma > 0, :);
P.N1 = all(abs(dc.medDelta) < semAdequate | isnan(dc.medDelta)) && all(abs(tC.dTarget) < N1_TARGET | isnan(tC.dTarget));
verdict_local("N1 null at empirical noise", P.N1, sprintf("Block C max |median delta| %.5f (MDC/2.77 %.5f), cells at or above: %d of %d; " + ...
    "max |delta target| %.4f over %d reliable pairs (threshold %.3f); NaN cells %d", max(abs(dc.medDelta)), semAdequate, ...
    sum(abs(dc.medDelta) >= semAdequate), height(dc), max(abs(tC.dTarget)), sum(isfinite(tC.dTarget)), N1_TARGET, sum(isnan(dc.medDelta))));
P.N2 = max(abs(ds.medDelta)) > max(abs(dc.medDelta));
verdict_local("N2 stress", P.N2, sprintf("max |median delta|: white noise %.5f, Block C %.5f", max(abs(ds.medDelta)), max(abs(dc.medDelta))));
W = unstack(D(:, ["block" "dataset" "FS" "f0" "betaGen" "regression" "variant" "medDelta"]), "medDelta", "variant");
big = abs(W.CEN) >= N3_MIN;
if any(big)
    agree = sign(W.AVG(big)) == sign(W.CEN(big)) & abs(W.AVG(big) - W.CEN(big)) < 0.5 * abs(W.CEN(big));
    P.N3 = all(agree);
    verdict_local("N3 alignment", P.N3, sprintf("%d of %d cells with |delta_CEN| >= %.3f agree", sum(agree), sum(big), N3_MIN));
else
    P.N3 = NaN;
    fprintf("N3 alignment NOT EVALUABLE: no cell has |median delta_CEN| >= %.3f\n", N3_MIN);
end

%% Figure
fg = figure("Color", "w", "Position", [80 80 1150 800]);  tl = tiledlayout(2, 2, "TileSpacing", "compact");
mk = ["o" "s" "^"];  ls = ["-" "--"];  cc = [0.05 0.27 0.49; 0.80 0.40 0.00; 0.20 0.60 0.30];
for blkName = ["C" "S"]
    nexttile; hold on;  sd = SD(SD.block == blkName, :);  dsn = unique(sd.dataset)';
    for i = 1:numel(dsn)
        for v = 2:numel(VARIANTS)
            for q = 1:numel(REG)
                t = sortrows(sd(sd.dataset == dsn(i) & sd.variant == VARIANTS(v) & sd.regression == REG_NAMES(q), :), "f0");
                plot(t.f0, t.maxAbsMedDelta, ls(v-1) + mk(q), "Color", cc(i, :), "DisplayName", ...
                    strrep(dsn(i), "_", " ") + " " + VARIANTS(v) + " " + REG_NAMES(q));
            end
        end
    end
    yline(semAdequate, ":", "MDC/2.77", "HandleVisibility", "off");  set(gca, "XScale", "log");  box on;
    xlabel("f_0 (Hz)");  ylabel("max |median paired \Delta\beta| over \beta_{gen}");
    if blkName == "C", title("Empirical noise (Block C centroids)"); else, title("White-noise stress (\sigma 0.5 mm)"); end
    legend("Location", "bestoutside", "FontSize", 6, "Interpreter", "none");
end
nexttile; hold on;  t = MD(MD.block == "C" & MD.sigma > 0, :);  dsn = unique(t.dataset)';
for i = 1:numel(dsn)
    for v = 2:numel(VARIANTS)
        u = sortrows(t(t.dataset == dsn(i) & t.variant == VARIANTS(v) & t.regression == "OLS", :), "f0");
        plot(u.f0, u.dTarget, ls(v-1) + "o", "Color", cc(i, :), "DisplayName", strrep(dsn(i), "_", " ") + " " + VARIANTS(v));
    end
end
yline([-N1_TARGET N1_TARGET], ":", "HandleVisibility", "off");  yline(0, "-", "Color", [0.7 0.7 0.7], "HandleVisibility", "off");
set(gca, "XScale", "log");  box on;  xlabel("f_0 (Hz)");  ylabel("\Delta target vs REG (OLS, reliable only)");
title("Target of the pull, empirical noise");  legend("Location", "best", "FontSize", 7, "Interpreter", "none");
nexttile; hold on;  ex = J(J.block == "C" & J.dataset == "Cook_CTRL" & J.f0 == 1 & J.sigma > 0, :);  ex = sortrows(ex, "betaGen");
for v = 1:numel(VARIANTS)
    med = cellfun(@(x) median(x(:, col_local(v, 1, numel(REG))), "omitnan"), ex.betaRec);
    plot(ex.betaGen, med, "-" + mk(v), "DisplayName", VARIANTS(v));
end
plot([0 0.75], [0 0.75], ":", "Color", [0.5 0.5 0.5], "HandleVisibility", "off");  box on;
xlabel("\beta_{gen}");  ylabel("median \beta_{rec} (BWFD-OLS)");  title("Example: Cook CTRL centroid, f_0 1 Hz");  legend("Location", "northwest");
title(tl, "BWFD half-sample offset between velocity and acceleration: paired effect under noise");
set(findall(fg, "Type", "axes"), "Toolbar", []);  exportgraphics(fg, OUT_PNG, "Resolution", 200);

save(OUT_MAT, "J", "D", "SD", "M", "MD", "P", "DS", "VARIANTS", "REG_NAMES", "BETA", "F0_C", "N_REPS", "RNG_SEED", ...
    "semAdequate", "N1_TARGET", "N3_MIN", "-v7.3");
fprintf("Saved: %s\nFigure: %s\n", OUT_MAT, OUT_PNG);

%% =========================================================================
function J = jobs_local(block, dataset, a, b, fs, f0, betaGen, sigma, alpha, nReps)
n = numel(f0);
J = table(repmat(string(block), n, 1), repmat(string(dataset), n, 1), repmat(a, n, 1), repmat(b, n, 1), repmat(fs, n, 1), ...
    f0(:), betaGen(:), repmat(sigma, n, 1), repmat(alpha, n, 1), repmat(nReps, n, 1), ...
    'VariableNames', ["block" "dataset" "a" "b" "FS" "f0" "betaGen" "sigma" "alpha" "nReps"]);
end

function o = runJob_local(j, cfg)
nV = numel(cfg.variants);  nR = numel(cfg.reg);
o = struct("lengthExcluded", false, "betaRec", NaN(j.nReps, nV * nR), "nFail", zeros(1, nV * nR), "firstErr", strings(1, nV * nR));
[x, y] = generatePowerLawEllipse_v001(j.a, j.b, j.f0, j.FS, j.betaGen, 'nCycles', cfg.nCyc);
x = x(:);  y = y(:);  M = numel(x);
if M - 2 * cfg.clip + 1 < cfg.minPts, o.lengthExcluded = true; return, end
if j.sigma > 0                                     % as tempoFoldEngine_v001 L85-89
    rs = RandStream("threefry4x64_20", "Seed", cfg.seed);  rs.Substream = j.jobID;  RandStream.setGlobalStream(rs);
end
for rep = 1:j.nReps
    xs = x;  ys = y;
    if j.sigma > 0
        xs = xs + reshape(generateCustomNoise_v003(M, j.alpha, j.sigma, j.FS), [], 1);
        ys = ys + reshape(generateCustomNoise_v003(M, j.alpha, j.sigma, j.FS), [], 1);
    end
    [dx, dy] = differentiateKinematicsEBR(xs, ys, 2, [2 cfg.fc 1], j.FS);          % REG, the registered path
    [vR, aR] = replica_local(xs, ys, j.FS, cfg.fc);
    N = M;  vC = NaN(N, 2);  vC(2:N-1, :) = (vR(2:N-1, :) + vR(3:N, :)) / 2;        % CEN: velocity to i
    aA = NaN(N, 2);  aA(3:N-1, :) = (aR(2:N-2, :) + aR(3:N-1, :)) / 2;              % AVG: acceleration to i - 1/2
    VA = {[dx(:, 2) dy(:, 2)], [dx(:, 3) dy(:, 3)]; vC, aR; vR, aA};
    c = cfg.clip;
    for v = 1:nV
        cols = (v - 1) * nR + (1:nR);
        vv = VA{v, 1}(c:end-c, :);  aa = VA{v, 2}(c:end-c, :);
        sp = hypot(vv(:, 1), vv(:, 2));  kp = curvatureKinematicEBR(vv(:, 1), vv(:, 2), aa(:, 1), aa(:, 2));
        ok = isfinite(sp) & isfinite(kp) & kp > 0;                                % as tempoFoldEngine_v001 L109
        if sum(ok) < cfg.minPts
            o = logFail_local(o, cols, sprintf("only %d valid points", sum(ok)));  continue
        end
        [b, errs] = regressThree_local(sp(ok), kp(ok), cfg.reg);
        o.betaRec(rep, cols) = b;
        for r = 1:nR, if errs(r) ~= "", o = logFail_local(o, cols(r), errs(r)); end, end
    end
end
end

function [vR, aR] = replica_local(x, y, fs, fc)
% Case 2 rebuilt from one filtfilt output, in differentiateKinematicsEBR's operation order
% (L80: diff after filtfilt; L84: diff of that velocity), so REG, CEN and AVG share it exactly.
[b, a] = butter(2, fc / (fs / 2));
xf = filtfilt(b, a, x);  yf = filtfilt(b, a, y);
N = numel(x);  vR = NaN(N, 2);  aR = NaN(N, 2);
vR(2:N, :) = [diff(xf) diff(yf)] * fs;
aR(2:N-1, :) = diff(vR(2:N, :)) * fs;
end

function [b, errs] = regressThree_local(sp, kp, REG)
% As tempoFoldEngine_v001 L126-140, seedMode "runner".
b = NaN(1, numel(REG));  errs = strings(1, numel(REG));  lm = [1, -1/3];
for r = 1:numel(REG)
    try
        [bb, vv] = regressDataEBR(sp, kp, REG(r), lm, 0, 0);
        b(r) = bb;
        if ~isfinite(bb), errs(r) = "non-finite beta"; end
        if r == 1 && isfinite(bb) && isfinite(vv), lm = [vv, bb]; end
    catch ME
        errs(r) = "regress: " + ME.message;
    end
end
end

function o = logFail_local(o, cols, msg)
o.nFail(cols) = o.nFail(cols) + 1;
for c = cols, if o.firstErr(c) == "", o.firstErr(c) = string(msg); end, end
end

function c = col_local(v, q, nR)
c = (v - 1) * nR + q;
end

function checkTiered_local(a, b, tol, nV, tag)
% Compare two nReps x (nV * nReg) replicate matrices (variant-major) with per-regression tolerances.
tc = repmat(tol(:)', 1, nV);
if ~isequal(size(a), size(b)) || ~isequal(isnan(a), isnan(b)) || any(abs(a - b) > tc, "all")
    error("halfSample:anchor0", "%s", sprintf("%s: replicates differ (max %.3g)", tag, max(abs(a - b), [], "all", "omitnan")));
end
end

function r = gainTarget_local(b, y, fitWin, locWin)
% extractGainTarget_v002 row_local (L206-230), measures only.
b = b(:);  y = y(:);
ok = isfinite(y) & b >= fitWin(1) & b <= fitWin(2);  b = b(ok);  y = y(ok);
[g, c, rm, gl, Tl, Tc] = deal(NaN);  nC = 0;
if numel(b) >= 4
    p = polyfit(b, y, 1);  g = p(1);  c = p(2);  res = y - polyval(p, b);  rm = sqrt(mean(res.^2));
    w = b >= locWin(1) & b <= locWin(2);
    if sum(w) >= 3, q = polyfit(b(w), y(w), 1); gl = q(1); end
    if g < 1, Tl = c / (1 - g); end
    d = y - b;  s = find(d(1:end-1) .* d(2:end) < 0 | d(1:end-1) == 0);  nC = numel(s);
    if nC >= 1
        xs = zeros(nC, 1);
        for k = 1:nC
            j = s(k);
            if d(j) == 0, xs(k) = b(j); else, xs(k) = b(j) - d(j) * (b(j+1) - b(j)) / (d(j+1) - d(j)); end
        end
        if isfinite(Tl), [~, kk] = min(abs(xs - Tl)); else, kk = 1; end
        Tc = xs(kk);
    end
end
r = struct("gain", g, "intercept", c, "targetLine", Tl, "targetCross", Tc, "nCross", nC, "rmse", rm, "gainLocal", gl);
end

function verdict_local(tag, held, detail)
s = "NOT HELD";  if held, s = "HELD"; end
fprintf("%s %s: %s\n", tag, s, detail);
end
