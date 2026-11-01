%% obsBetaRateCutoff_v002.m
% D-D16 stage S1 (docs/SPEC_D16_SamplingRateCutoff_v001.md, amendment of 2026-10-24 13:04). v001 changed the rate
% with decimate(), whose 8th-order Chebyshev low-pass is not what a slower digitiser does and cannot show aliasing.
% Dagmar, 2026-10-24 13:04: "decimation should not add filtering - we subsample cleanly". v002 subsamples cleanly:
% x(1:D:end), y(1:D:end), D = 2 (120 Hz) and 4 (60 Hz), no filter, starting at the first sample. Everything else is
% v001's. The 240 Hz settings (A240p and the cutoff sweep) involve no rate change, so their beta_obs are taken from
% v001's results file (results/obsBetaRateCutoff_v001.mat, which passed anchor A) after checking that the setting
% definitions and trial order match; only the subsampled arms, the noise at 120 and 60 Hz, R1, R3, Q1 and Q4 are new.
% On the Fraser trials of the v015 corpus, without forward maps:
%   arms     A120p/A60p (edge clip 20 samples, production), A120d/A60d (clip 10/5, 0.083 s), clean subsampling;
%   noise    characteriseBiologicalNoise_v004(Method="irasa") on the subsampled trials (as rerunNoiseCharacterisation_v005);
%   beta_obs exactly as runLoopClosureFftnoise_v015's processOneTrial_local (first block): differentiate, clip,
%            regress OLS, LMLS, IRLS with OLS start values; SG curvature residue (below RESID) counted per trial (RS).
% Rules R1 (alpha per octave, person bootstrap, negligible if 95% CI within +-0.15) and R3 (clip arms to the HPC
% if the median |p - d| >= 0.0108 in any pipeline at either rate) as in the spec; Q1/Q4 descriptive; Q3 (240 Hz
% sweep and its zero-noise comparator) recomputed from v001's values, so identical to v001's.
% Anchor (fail loud): B (noise at 240 Hz on N_ANCHOR random trials equals noiseCharacterisation_fraser.mat).
% Limit: one subsampling phase (the first sample); the other D - 1 phases are not run.
% Reads results/obsBetaRateCutoff_v001.mat, src/loopClosureResults_Fraser_all_shaped_xu_v015.mat,
% src/noiseCharacterisation_fraser.mat, results/filterCutoffCollapse_v002.mat. Writes results/obsBetaRateCutoff_v002.mat
% (checkpoint _ckpt.mat while running) and a dated log in logs/. Run from src/. Starts a 3-worker pool if none is open.

%% CONFIG
ROOT    = fileparts(fileparts(mfilename("fullpath")));  SRC = fullfile(ROOT, "src");
IN_V015 = fullfile(SRC, "loopClosureResults_Fraser_all_shaped_xu_v015.mat");
IN_BIO  = fullfile(SRC, "noiseCharacterisation_fraser.mat");
IN_COLL = fullfile(ROOT, "results", "filterCutoffCollapse_v002.mat");
OUT_MAT = fullfile(ROOT, "results", "obsBetaRateCutoff_v002.mat");  CKPT = fullfile(ROOT, "results", "obsBetaRateCutoff_v002_ckpt.mat");
IN_V1   = fullfile(ROOT, "results", "obsBetaRateCutoff_v001.mat");
FS0 = 240;  SIG_TO_MM = 1 / 10.41793;
SET = struct("name", {"A240p" "A120p" "A120d" "A60p" "A60d" "S240fc05" "S240fc07" "S240fc10" "S240fc14" "S240fc20"}, ...
             "dec",  {1 2 2 4 4 1 1 1 1 1}, "clip", {20 20 10 20 5 41 41 41 41 41}, ...
             "fc",   {10 10 10 10 10 5 7 10 14 20}, "ref", {17 17 17 17 17 34 24 17 12 8.5});
PIPE = ["BWFD-OLS" "SG-OLS" "BWFD-LMLS" "SG-LMLS" "BWFD-IRLS" "SG-IRLS"];   % column-major order of betaObs(d, r), as v015
REG  = [3 4 5];  N_ANCHOR = 50;  TOL_NOISE = 1e-8;  N_BOOT = 2000;  R1_BOUND = 0.15;  SEM = 0.0108;  SEED = 20261024;
RESID = 1e-10;          % SG curvature below this is floating-point residue (see header)
ZARANDI_ALPHA = 3.10;   % PATCH_DesignRecs_TempoMatching_v002, "Correction" (v015 IRASA median), for Q4 only

%% Setup
logFile = fullfile(ROOT, "logs", "obsBetaRateCutoff_v002_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log");
LG = fopen(logFile, "w");  if LG < 0, error("d16:log", "%s", "Cannot open log: " + logFile); end
say_local(LG, "obsBetaRateCutoff_v002 run %s (clean subsampling)\n", string(datetime("now")));
addpath(SRC);  addpath(genpath(fullfile(SRC, "functions")));  addpath(genpath(fullfile(SRC, "req")));
for f = [IN_V015 IN_BIO IN_COLL IN_V1], if ~isfile(f), error("d16:input", "%s", "Missing: " + f); end, end
if isempty(gcp("nocreate")), parpool("Processes", 3); end
V = load(IN_V015, "results").results;  bio = load(IN_BIO, "bioResults").bioResults;
ids = string({V.trialID})';  nT = numel(ids);
trials = importDB_fraser_v001(Shapes=3, Visits=[1 2], Verbose=false);
[ok, loc] = ismember(ids, string({trials.trialID}));
if ~all(ok), error("d16:import", "%s", sprintf("%d v015 trials not imported", sum(~ok))); end
trials = trials(loc);
if any([trials.fs] ~= FS0), error("d16:fs", "%s", "Imported trials are not all at 240 Hz"); end
say_local(LG, "Trials: %d (v015 Fraser set), %d persons.\n", nT, numel(unique(string({V.subjectID}))));

%% Subsampled copies: every D-th sample from the first, no filter
X = cell(nT, 3);  Y = cell(nT, 3);  DECS = [1 2 4];
for k = 1:nT
    for q = 1:3
        x = double(trials(k).x);  y = double(trials(k).y);  X{k, q} = x(1:DECS(q):end);  Y{k, q} = y(1:DECS(q):end);
    end
end

%% The 240 Hz settings from v001 (no rate change; v001 passed anchor A)
P1 = load(IN_V1, "B", "RS", "SET", "T", "frames");   % v001 saves the trial IDs in T
if ~isequal(P1.T.trialID, ids), error("d16:v1ids", "%s", "v001 trial order differs"); end
if ~isequal(struct2table(P1.SET), struct2table(SET)), error("d16:v1set", "%s", "v001 setting definitions differ from v002's"); end

%% beta_obs in every arm and sweep setting
B = NaN(numel(SET), nT, 6);  RS = NaN(numel(SET), nT);  frames = NaN(numel(SET), 1);
for s = 1:numel(SET)
    if SET(s).dec == 1
        B(s, :, :) = P1.B(s, :, :);  RS(s, :) = P1.RS(s, :);  frames(s) = P1.frames(s);
        say_local(LG, "%-9s fs 240 Hz: taken from v001 (SG curvature residue in %d trials)\n", SET(s).name, sum(RS(s, :) > 0));
        continue
    end
    t0 = tic;  q = find(DECS == SET(s).dec);  fs = FS0 / SET(s).dec;
    xs = X(:, q);  ys = Y(:, q);  cfg = SET(s);  out = NaN(nT, 7);
    parfor k = 1:nT
        out(k, :) = betaObs_local(xs{k}, ys{k}, fs, cfg, REG, RESID); %#ok<PFBNS>
    end
    B(s, :, :) = out(:, 1:6);  RS(s, :) = out(:, 7);  frames(s) = sgFrame_local(SET(s).ref, fs);
    say_local(LG, "%-9s fs %3d Hz, clip %2d, BW %4.1f Hz, SG frame %2d (%.3f s): %d trials in %.0f s; finite SG-IRLS %d; SG curvature residue in %d trials\n", ...
        SET(s).name, fs, SET(s).clip, SET(s).fc, frames(s), frames(s) / fs, nT, toc(t0), sum(isfinite(out(:, 6))), sum(out(:, 7) > 0));
    save(CKPT, "B", "RS", "frames", "SET", "ids", "-v7.3");
end

%% Noise: anchor at 240 Hz, then the subsampled rates
rng(SEED);  ia = sort(randperm(nT, N_ANCHOR));
nb = noise_local(trials(ia), X(ia, 1), Y(ia, 1), FS0);
[~, jb] = ismember(string(nb.trialID), string(bio.trialID));
d = [abs(nb.ira_alphaMean - bio.ira_alphaMean(jb)); abs(nb.sigmaMean - bio.sigmaMean(jb))];
if any(jb == 0) || max(d(isfinite(d))) > TOL_NOISE || any(isnan(nb.ira_alphaMean) ~= isnan(bio.ira_alphaMean(jb)))
    fail_local(LG, "d16:anchorB", sprintf("Noise at 240 Hz does not reproduce the noise file on %d trials (max diff %.3g)", N_ANCHOR, max(d(isfinite(d)))));
end
say_local(LG, "Anchor B passed: noise at 240 Hz reproduces noiseCharacterisation_fraser.mat on %d random trials (max diff %.3g).\n", N_ANCHOR, max(d(isfinite(d))));
[~, j0] = ismember(ids, string(bio.trialID));
A = NaN(nT, 3);  Sg = NaN(nT, 3);  A(:, 1) = bio.ira_alphaMean(j0);  Sg(:, 1) = bio.sigmaMean(j0) * SIG_TO_MM;
for q = 2:3
    t0 = tic;  nq = noise_local(trials, X(:, q), Y(:, q), FS0 / DECS(q));
    [~, jq] = ismember(ids, string(nq.trialID));  A(:, q) = nq.ira_alphaMean(jq);  Sg(:, q) = nq.sigmaMean(jq) * SIG_TO_MM;
    say_local(LG, "Noise at %d Hz: %d trials in %.0f s; finite alpha %d\n", FS0 / DECS(q), nT, toc(t0), sum(isfinite(A(:, q))));
    save(CKPT, "B", "RS", "frames", "SET", "ids", "A", "Sg", "-v7.3");
end

%% R1: alpha per octave of rate (person-level bootstrap of the mean per-trial slope)
lr = log2(FS0 ./ DECS);  slope = NaN(nT, 1);
for k = 1:nT, g = isfinite(A(k, :)); if sum(g) >= 2, p = polyfit(lr(g), A(k, g), 1); slope(k) = p(1); end, end
person = string({V.subjectID})';  [up, ~, pIdx] = unique(person);  rng(SEED);  bs = NaN(N_BOOT, 1);
for b = 1:N_BOOT
    pick = randi(numel(up), numel(up), 1);  w = accumarray(pick, 1, [numel(up) 1]);  wt = w(pIdx);
    g = isfinite(slope) & wt > 0;  bs(b) = sum(slope(g) .* wt(g)) / sum(wt(g));
end
r1 = struct("mean", mean(slope, "omitnan"), "lo", prctile(bs, 2.5), "hi", prctile(bs, 97.5), "nTrials", sum(isfinite(slope)), "nPersons", numel(up));
r1.negligible = r1.lo > -R1_BOUND && r1.hi < R1_BOUND;
say_local(LG, "\nR1: alpha changes by %+.3f per octave of rate (95%% person-bootstrap CI %+.3f to %+.3f; %d trials, %d persons): %s (bound +-%.2f)\n", ...
    r1.mean, r1.lo, r1.hi, r1.nTrials, r1.nPersons, ifelse_local(r1.negligible, "NEGLIGIBLE", "NOT negligible"), R1_BOUND);

%% R3: does the edge clip matter?
R3 = table();
for rate = [120 60]
    sp = find([SET.name] == "A" + rate + "p");  sd = find([SET.name] == "A" + rate + "d");
    md = squeeze(median(abs(B(sp, :, :) - B(sd, :, :)), 2, "omitnan"))';
    R3 = [R3; array2table(md, "VariableNames", PIPE, "RowNames", "fs" + rate)]; %#ok<AGROW>
end
r3trigger = any(R3{:, :} >= SEM, "all");
say_local(LG, "\nR3: median |beta_obs(clip in samples) - beta_obs(clip in duration)| per pipeline:\n%s-> %s\n", evalc("disp(R3)"), ...
    ifelse_local(r3trigger, "TRIGGERED: the duration-matched arms also go to the HPC", "not triggered: the clip is immaterial"));

%% Q1 and Q4: the matched-rate table (production arms)
Q1 = table();  PROD = ["A240p" "A120p" "A60p"];
for q = 1:3
    s = find([SET.name] == PROD(q));
    Q1 = [Q1; [table(FS0 / DECS(q), median(A(:, q), "omitnan"), median(Sg(:, q), "omitnan"), 'VariableNames', ["fs" "alphaMedian" "sigmaMMMedian"]), ...
        array2table(squeeze(median(B(s, :, :), 2, "omitnan"))', "VariableNames", "betaObs_" + replace(PIPE, "-", "_"))]]; %#ok<AGROW>
end
say_local(LG, "\nQ1: medians over trials, production pipelines:\n%s", evalc("disp(Q1)"));
say_local(LG, "Q4: Fraser median alpha %.2f (240), %.2f (120), %.2f (60) Hz against Zarandi %.2f (100 Hz, processed export).\n", Q1.alphaMedian, ZARANDI_ALPHA);

%% Q3, observed side: the cutoff sweep, and the zero-noise harmonic-gain comparator
C = load(IN_COLL, "G").G;  f0 = [V.f0]';  bgs = [V.betaGenStarMed]';
gBW = curve_local(C(C.family == "BW" & C.fc == 10, :), "gainOLS");
gSG = curve_local(C(C.family == "SG" & abs(C.Tw - 41 / 240) < 1e-9, :), "gainIRLS");
s10 = find([SET.name] == "S240fc10");  Q3 = table();
for s = find(startsWith([SET.name], "S240") & [SET.name] ~= "S240fc10")
    obsBW = B(s, :, 1)' - B(s10, :, 1)';  obsSG = B(s, :, 6)' - B(s10, :, 6)';
    prBW = (gBW(f0 / SET(s).fc) - gBW(f0 / 10)) .* (bgs - 1/3);
    prSG = (gSG(f0 * frames(s) / FS0) - gSG(f0 * frames(s10) / FS0)) .* (bgs - 1/3);
    for pz = ["BWFD-OLS" "SG-IRLS"]
        if pz == "BWFD-OLS", o = obsBW; pr = prBW; else, o = obsSG; pr = prSG; end
        g = isfinite(o) & isfinite(pr);
        Q3 = [Q3; table(string(SET(s).name), pz, sum(g), median(o(g)), median(pr(g)), ccc_local(o(g), pr(g)), median(abs(o(g) - pr(g))), ...
            'VariableNames', ["setting" "pipeline" "nTrials" "medianObsChange" "medianPredChange" "CCC" "medianAbsErr"])]; %#ok<AGROW>
    end
end
say_local(LG, "\nQ3 (observed side; comparator = zero-noise harmonic gain x (beta_gen* - 1/3), no rule):\n%s", evalc("disp(Q3)"));

%% Save
T = table(ids, person, f0, bgs, A, Sg, slope, 'VariableNames', ["trialID" "person" "f0" "betaGenStarMed" "alpha" "sigmaMM" "alphaSlopePerOct"]);
save(OUT_MAT, "B", "RS", "RESID", "SET", "PIPE", "frames", "ids", "T", "r1", "bs", "R3", "r3trigger", "Q1", "Q3", "SEED", "-v7.3");
if isfile(CKPT), delete(CKPT); end
say_local(LG, "\nSaved: %s\nLog: %s\n", OUT_MAT, logFile);  fclose(LG);

%% =========================================================================
function b = betaObs_local(x, y, fs, cfg, REG, RESID)
% runLoopClosureFftnoise_v015 processOneTrial_local, first block, with this setting's filters and clip.
% Element 7: the number of SG curvature samples below RESID (floating-point residue; see header).
    dc = struct("filterType", {2, 6}, "filterParams", {[2 cfg.fc 1], [4 cfg.ref]});  bo = NaN(2, 3);  C = cfg.clip;
    for dI = 1:2
        [dx, dy] = differentiateKinematicsEBR(x, y, dc(dI).filterType, dc(dI).filterParams, fs);
        vx = dx(C:end-C, 2);  vy = dy(C:end-C, 2);  ax = dx(C:end-C, 3);  ay = dy(C:end-C, 3);
        sp = hypot(vx, vy);  kp = curvatureKinematicEBR(vx, vy, ax, ay);  lm = [1, -1/3];  if dI == 2, nres = sum(kp < RESID); end
        for rI = 1:3
            try
                [bb, v] = regressDataEBR(sp, kp, REG(rI), lm, 0, 0);  bo(dI, rI) = bb;
                if rI == 1, lm = [v, bb]; end
            catch
            end
        end
    end
    b = [bo(:)' nres];
end

function nb = noise_local(trials, xs, ys, fs)
% characteriseBiologicalNoise_v004(Method="irasa") on the given coordinates, in parallel chunks.
    n = numel(trials);  nc = 12;  edges = round(linspace(0, n, nc + 1));  parts = cell(nc, 1);
    for k = 1:n, trials(k).x = xs{k};  trials(k).y = ys{k};  trials(k).fs = fs;  trials(k).t = (0:numel(xs{k}) - 1)' / fs; end
    parfor c = 1:nc
        parts{c} = characteriseBiologicalNoise_v004(trials(edges(c) + 1:edges(c + 1)), Method="irasa", Verbose=false); %#ok<PFBNS>
    end
    nb = vertcat(parts{:});
end

function n = sgFrame_local(ref, fs)
% The frame differentiateKinematicsEBR case 6 uses (reference fs 100 Hz, form 4k+1), for logging and the comparator.
    n = round(ref / 100 * fs);  n = max(n, 5);  r = mod(n - 1, 4);
    if r == 0 && mod(n, 2) == 0, n = n + 4; elseif r == 1, n = n - 1; elseif r == 2, n = n + 2; elseif r == 3, n = n + 1; end
    if n < 5, n = 5; end
end

function g = curve_local(Cs, col)
% Gain against x (f0/f_c or f0*T_w) on the collapse nodes; log-x interpolation, clamped at the ends.
    Cs = sortrows(Cs(isfinite(Cs.(col)) & isfinite(Cs.x), :), "x");
    lx = log(Cs.x);  gv = Cs.(col);
    g = @(x) interp1(lx, gv, min(max(log(x), lx(1)), lx(end)), "linear");
end

function c = ccc_local(a, b)
    m = mean(a);  n = mean(b);  c = 2 * mean((a - m) .* (b - n)) / (var(a, 1) + var(b, 1) + (m - n) ^ 2);
end

function v = ifelse_local(c, a, b)
    if c, v = a; else, v = b; end
end

function fail_local(fid, id, msg)
% Write the failure to the log, then stop (an error alone does not reach the log).
    fprintf(fid, "FAILED %s: %s\n", id, msg);  fclose(fid);  error(id, "%s", msg);
end

function say_local(fid, fmt, varargin)
    fprintf(fmt, varargin{:});  fprintf(fid, fmt, varargin{:});
end
