%% obsAlphaRateLadder_v001.m
% Sampling-rate ladder: noise colour (alpha, sigma) and beta_obs on cleanly subsampled real trials, every integer divisor of
% 240 Hz. Rules were fixed in docs/SPEC_AlphaRateLadder_v001.md before the first run; read that first.
% Why: Finding #98 (Pilot, 744 integer-stride schemes, ellipses) gave +0.061 alpha per octave and Finding #260 (Fraser, 240/120/60 Hz,
% ellipses) gave +0.187. They differ in the template (global 12 harmonics fitted at 240 Hz then subsampled, against v004's sliding
% 4 harmonics refitted at each rate), in phases (all against the first) and in the rate range fitted.
% Noise arms:  A4 = production (subsample, then characteriseBiologicalNoise_v004, 4 harmonics); A12 = the same with 12 harmonics;
%              B = the Pilot method (edge clip 20, global 12-harmonic template at 240 Hz, subsample the residuals, iraAlphaSigma_v001).
% beta_obs arms (OLS, LMLS, IRLS each): BW = [2 10 1]; SGp = differentiateKinematicsEBR case 6 [4 17] (production: fs-scaled, 4k+1,
%              floor 5 samples); SGw = case 4 with the odd frame nearest 41/240 s (window-matched, floor 5 samples).
%              Edge clip = max(round(fs/12), SG half-frame + 1) samples, so it is 20 samples at 240 Hz as in production.
% Stages (STAGE below): 0 smoke test; 1 every trial, first phase of every rung (iMac); 2 every phase of every rung (HPC; restartable).
% Reads importDB_fraser_v001 or importDB_pilot_v001 and, for the Fraser anchor, src/noiseCharacterisation_fraser.mat.
% Writes results/obsAlphaRateLadder_v001_<corpus>_stage<k>.mat (checkpoint beside it while running) and a dated log in logs/.
% Run from src/. Starts a pool if none is open. Fail loud: a failed anchor stops the run; IRASA failures are counted and logged.

%% CONFIG
STAGE     = 0;                      % 0 smoke test | 1 all trials, first phase | 2 all phases (see header)
CORPUS    = "fraser";               % "fraser" | "pilot"
SHAPES    = 3;                      % imageNumber(s); the project copy of Fraser holds 3 (ellipse) only
DATA_ROOT = "";                     % "" = importer default; set to an all-shapes root to use other SHAPES (spec section 2)
FS0 = 240;  REF = 17;  BW_FC = 10;  WM_SEC = 41 / FS0;  CLIP_SEC = 20 / FS0;  MIN_USABLE = 30;
NH_LIST = [4 12];  PILOT_NH = 12;  PILOT_CLIP = 20;  PILOT_MIN_BINS = 5;
REG = [3 4 5];  RESID = 1e-10;  SEED = 20261102;  N_WORKERS = 3;  MAX_PER_PERSON = Inf;
ANCHOR_FILE = "noiseCharacterisation_fraser.mat";  TOL_NOISE = 1e-8;   % "" disables the anchor (the log says so)
PIPE = ["BW-OLS" "BW-LMLS" "BW-IRLS" "SGp-OLS" "SGp-LMLS" "SGp-IRLS" "SGw-OLS" "SGw-LMLS" "SGw-IRLS"];
ALL_FS = FS0 ./ find(mod(FS0, 1:FS0) == 0);
switch STAGE
    case 0, MAX_TRIALS = 24;   RUN_FS = [240 120 60 24 10];  PHASES = 1;
    case 1, MAX_TRIALS = Inf;  RUN_FS = ALL_FS;              PHASES = 1;
    case 2, MAX_TRIALS = Inf;  RUN_FS = ALL_FS;              PHASES = Inf;
    otherwise, error("ladder:stage", "%s", "STAGE must be 0, 1 or 2");
end

%% Setup
ROOT = fileparts(fileparts(mfilename("fullpath")));  SRC = fullfile(ROOT, "src");
TAG  = sprintf("%s_stage%d", CORPUS, STAGE);
OUT_MAT = fullfile(ROOT, "results", "obsAlphaRateLadder_v001_" + TAG + ".mat");  CKPT = replace(OUT_MAT, ".mat", "_ckpt.mat");
logFile = fullfile(ROOT, "logs", "obsAlphaRateLadder_v001_" + TAG + "_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log");
LG = fopen(logFile, "w");  if LG < 0, error("ladder:log", "%s", "Cannot open log: " + logFile); end
say_local(LG, "obsAlphaRateLadder_v001 %s, stage %d, %s\n", CORPUS, STAGE, string(datetime("now")));
addpath(SRC);  addpath(genpath(fullfile(SRC, "functions")));  addpath(genpath(fullfile(SRC, "req")));
if isempty(gcp("nocreate")), parpool("Processes", N_WORKERS); end
chk = [240 41; 120 21; 80 13; 60 9; 48 9; 40 9; 30 5; 24 5; 20 5; 10 5];   % case 6 frames, checked against differentiateKinematicsEBR 2026-11-02
for i = 1:size(chk, 1)
    if sgFrame_local(REF, chk(i, 1)) ~= chk(i, 2), fail_local(LG, "ladder:sgframe", sprintf("SG frame at %d Hz is %d, expected %d", chk(i, 1), sgFrame_local(REF, chk(i, 1)), chk(i, 2))); end
end

%% Trials
args = {};  if DATA_ROOT ~= "", args = {char(DATA_ROOT)}; end
if CORPUS == "fraser", trials = importDB_fraser_v001(args{:}, Shapes=SHAPES, Visits=[1 2], Verbose=false);
elseif CORPUS == "pilot", trials = importDB_pilot_v001(args{:}, Shapes=SHAPES, Visits=[1 2], Verbose=false);
else, fail_local(LG, "ladder:corpus", "CORPUS must be fraser or pilot"); end
trials = trials(:);  if isempty(trials), fail_local(LG, "ladder:empty", "No trials imported (shape not in this data root?)"); end
if any([trials.fs] ~= FS0), fail_local(LG, "ladder:fs", "Imported trials are not all at 240 Hz"); end
rng(SEED);  person = string({trials.subjectID})';
if isfinite(MAX_PER_PERSON)
    keep = false(numel(trials), 1);
    for u = unique(person)', j = find(person == u);  j = j(randperm(numel(j), min(MAX_PER_PERSON, numel(j))));  keep(j) = true; end
    trials = trials(keep);  person = person(keep);
end
if isfinite(MAX_TRIALS) && numel(trials) > MAX_TRIALS, j = sort(randperm(numel(trials), MAX_TRIALS));  trials = trials(j);  person = person(j); end
nT = numel(trials);  ids = string({trials.trialID})';
shapeCode = NaN(nT, 1);  if isfield(trials, "shapeCode"), shapeCode = [trials.shapeCode]'; end
if numel(unique(SHAPES)) > 1 && ~isfield(trials, "shapeCode"), fail_local(LG, "ladder:shape", "Several shapes requested but trials carry no shapeCode"); end
s2mm = NaN(nT, 1);  if isfield(trials, "SigmaToMM"), s2mm = [trials.SigmaToMM]'; end
say_local(LG, "Trials: %d, persons: %d, shapes: %s\n", nT, numel(unique(person)), mat2str(unique(shapeCode(isfinite(shapeCode)))'));
X0 = cellfun(@double, {trials.x}, "UniformOutput", false)';  Y0 = cellfun(@double, {trials.y}, "UniformOutput", false)';

%% Schemes: every rung, phases 1..min(D, PHASES)
Dv = FS0 ./ RUN_FS(:);  if any(Dv ~= round(Dv)), fail_local(LG, "ladder:rungs", "RUN_FS must be integer divisors of 240"); end
rows = zeros(0, 2);
for d = Dv', rows = [rows; repmat(d, min(d, PHASES), 1) (1:min(d, PHASES))']; end %#ok<AGROW>
nS = size(rows, 1);  S = table(rows(:, 1), rows(:, 2), FS0 ./ rows(:, 1), 'VariableNames', ["D" "phase" "fs"]);
for i = 1:nS, RG(i) = rung_local(S.D(i), FS0, REF, BW_FC, WM_SEC, CLIP_SEC); end %#ok<SAGROW>
S.frameProd = [RG.frameProd]';  S.frameWM = [RG.frameWM]';  S.durProd = S.frameProd ./ S.fs;  S.durWM = S.frameWM ./ S.fs;
S.floorBound = [RG.floorBound]';  S.bwValid = [RG.bwValid]';  S.clip = [RG.clip]';
say_local(LG, "Schemes: %d (rungs %d, phases per rung up to %s)\n%s", nS, numel(Dv), string(PHASES), evalc("disp(S(S.phase == 1, :))"));

%% Pilot-method residuals at 240 Hz (arm B), computed once
xn = X0;  yn = Y0;  rXc = cell(nT, 1);  rYc = cell(nT, 1);  f0n = NaN(nT, 1);
parfor k = 1:nT
    [rXc{k}, rYc{k}, f0n(k)] = pilotResid_local(xn{k}, yn{k}, FS0, PILOT_NH, PILOT_CLIP);
end
say_local(LG, "Arm B residuals: %d of %d trials usable; f0 median %.3f Hz\n", sum(~cellfun(@isempty, rXc)), nT, median(f0n, "omitnan"));

%% Result arrays and checkpoint
nNH = numel(NH_LIST);
aA = NaN(nT, nS, nNH);  sA = aA;  aB = NaN(nT, nS);  sB = aB;  bet = NaN(nT, nS, 9);  rs = aB;  nUse = aB;
nBadB = zeros(nS, 1);  tSec = NaN(nS, 1);  done = false(nS, 1);
if isfile(CKPT)
    K = load(CKPT);
    if ~isequal(K.ids, ids) || ~isequal(K.S.D, S.D) || ~isequal(K.S.phase, S.phase)
        fail_local(LG, "ladder:ckpt", "Checkpoint " + CKPT + " does not match this configuration; delete it or change STAGE/CORPUS");
    end
    aA = K.aA;  sA = K.sA;  aB = K.aB;  sB = K.sB;  bet = K.bet;  rs = K.rs;  nUse = K.nUse;  nBadB = K.nBadB;  tSec = K.tSec;  done = K.done;
    say_local(LG, "Resumed from checkpoint: %d of %d schemes done.\n", sum(done), nS);
end
anchorOn = CORPUS == "fraser" && ANCHOR_FILE ~= "";
if ~anchorOn, say_local(LG, "NOTE: anchor against the stored noise file is OFF (CORPUS %s, ANCHOR_FILE '%s').\n", CORPUS, ANCHOR_FILE); end

%% Ladder
for s = 1:nS
    if done(s), continue, end
    D = S.D(s);  p = S.phase(s);  fs = S.fs(s);  Rg = RG(s);  t0 = tic;
    xs = cell(nT, 1);  ys = cell(nT, 1);
    for k = 1:nT, xs{k} = X0{k}(p:D:end);  ys{k} = Y0{k}(p:D:end); end
    bk = NaN(nT, 11);
    parfor k = 1:nT
        bk(k, :) = betaObs_local(xs{k}, ys{k}, fs, Rg, REG, RESID, MIN_USABLE); %#ok<PFBNS>
    end
    bet(:, s, :) = bk(:, 1:9);  rs(:, s) = bk(:, 10);  nUse(:, s) = bk(:, 11);
    for h = 1:nNH
        nb = noise_local(trials, xs, ys, fs, NH_LIST(h));
        [okn, jn] = ismember(ids, string(nb.trialID));
        if ~all(okn), fail_local(LG, "ladder:noiseids", "characteriseBiologicalNoise_v004 dropped trials"); end
        aA(:, s, h) = nb.ira_alphaMean(jn);  sA(:, s, h) = nb.sigmaMean(jn);
    end
    ab = NaN(nT, 1);  sb = NaN(nT, 1);  bad = false(nT, 1);
    parfor k = 1:nT
        [ab(k), sb(k), bad(k)] = armB_local(rXc{k}, rYc{k}, D, p, fs, PILOT_MIN_BINS);
    end
    aB(:, s) = ab;  sB(:, s) = sb;  nBadB(s) = sum(bad);
    if anchorOn && D == 1 && p == 1 && NH_LIST(1) == 4, anchor_local(LG, ANCHOR_FILE, ids, aA(:, s, 1), sA(:, s, 1), TOL_NOISE); end
    tSec(s) = toc(t0);  done(s) = true;
    say_local(LG, "%3d Hz phase %d: SG frame %2d/%2d (%.3f/%.3f s)%s, clip %2d, usable median %4.0f; alpha A4 %.2f A12 %.2f B %.2f; SGp-IRLS %.3f; IRASA failures B %d; %.0f s\n", ...
        fs, p, S.frameProd(s), S.frameWM(s), S.durProd(s), S.durWM(s), ifelse_local(S.floorBound(s), " floor", ""), S.clip(s), median(nUse(:, s), "omitnan"), ...
        median(aA(:, s, 1), "omitnan"), median(aA(:, s, nNH), "omitnan"), median(aB(:, s), "omitnan"), median(bet(:, s, 6), "omitnan"), nBadB(s), tSec(s));
    save(CKPT, "ids", "S", "aA", "sA", "aB", "sB", "bet", "rs", "nUse", "nBadB", "tSec", "done", "-v7.3");
end
if anchorOn && ~any(S.D == 1 & S.phase == 1), say_local(LG, "NOTE: anchor skipped because the 240 Hz rung was not run.\n"); end

%% Save
CFG = struct("STAGE", STAGE, "CORPUS", CORPUS, "SHAPES", SHAPES, "DATA_ROOT", DATA_ROOT, "FS0", FS0, "REF", REF, "BW_FC", BW_FC, "WM_SEC", WM_SEC, ...
    "CLIP_SEC", CLIP_SEC, "MIN_USABLE", MIN_USABLE, "NH_LIST", NH_LIST, "PILOT_NH", PILOT_NH, "PILOT_CLIP", PILOT_CLIP, "REG", REG, "SEED", SEED, "PHASES", PHASES, "anchorOn", anchorOn);
f0 = f0n;  s2mmTrial = s2mm; %#ok<NASGU>
save(OUT_MAT, "CFG", "ids", "person", "shapeCode", "s2mmTrial", "f0", "S", "PIPE", "aA", "sA", "aB", "sB", "bet", "rs", "nUse", "nBadB", "tSec", "-v7.3");
if isfile(CKPT), delete(CKPT); end
say_local(LG, "Saved: %s\nLog: %s\nTotal %.0f min\n", OUT_MAT, logFile, sum(tSec, "omitnan") / 60);  fclose(LG);

%% =========================================================================
function R = rung_local(D, FS0, REF, BW_FC, WM_SEC, CLIP_SEC)
% Per-rung filter settings and the edge clip. The SG frames are logged, never silently altered.
    fs = FS0 / D;  fp = sgFrame_local(REF, fs);  fw = round(WM_SEC * fs);  fw = max(fw + mod(fw + 1, 2), 5);
    R = struct("D", D, "fs", fs, "refProd", REF, "bwFc", BW_FC, "bwValid", BW_FC < fs / 2, "frameProd", fp, "frameWM", fw, ...
        "floorBound", round(REF / 100 * fs) < 5, "clip", max(round(CLIP_SEC * fs), (max(fp, fw) - 1) / 2 + 1));
end

function n = sgFrame_local(ref, fs)
% The frame differentiateKinematicsEBR case 6 uses (reference fs 100 Hz, form 4k+1, floor order + 1 = 5).
    n = round(ref / 100 * fs);  n = max(n, 5);  r = mod(n - 1, 4);
    if r == 0 && mod(n, 2) == 0, n = n + 4; elseif r == 1, n = n - 1; elseif r == 2, n = n + 2; elseif r == 3, n = n + 1; end
    if n < 5, n = 5; end
end

function b = betaObs_local(x, y, fs, R, REG, RESID, MIN_USABLE)
% runLoopClosureFftnoise_v015 processOneTrial_local, first block, for BW, SGp, SGw. Elements 1-9 beta_obs (BW, SGp, SGw; OLS, LMLS, IRLS),
% 10 the number of SGp curvature samples below RESID (floating-point residue), 11 the usable samples after the clip.
    b = NaN(1, 11);  C = R.clip;  b(11) = numel(x) - 2 * C + 1;
    if b(11) < MIN_USABLE, return, end
    dc = struct("filterType", {2, 6, 4}, "filterParams", {[2 R.bwFc 1], [4 R.refProd], [4 R.frameWM]});
    for dI = 1:3
        if dI == 1 && ~R.bwValid, continue, end
        if dI == 3 && R.frameWM == R.frameProd, b(7:9) = b(4:6);  continue, end
        [dx, dy] = differentiateKinematicsEBR(x, y, dc(dI).filterType, dc(dI).filterParams, fs);
        vx = dx(C:end-C, 2);  vy = dy(C:end-C, 2);  ax = dx(C:end-C, 3);  ay = dy(C:end-C, 3);
        sp = hypot(vx, vy);  kp = curvatureKinematicEBR(vx, vy, ax, ay);  lm = [1, -1/3];  if dI == 2, b(10) = sum(kp < RESID); end
        for rI = 1:3
            try
                [bb, v] = regressDataEBR(sp, kp, REG(rI), lm, 0, 0);  b(3 * (dI - 1) + rI) = bb;
                if rI == 1, lm = [v, bb]; end
            catch    % non-convergence (limitBreak 0) stays NaN, as in the production runner
            end
        end
    end
end

function nb = noise_local(trials, xs, ys, fs, nh)
% characteriseBiologicalNoise_v004(Method="irasa") on the given coordinates, in parallel chunks.
    n = numel(trials);  nc = min(12, n);  edges = round(linspace(0, n, nc + 1));  parts = cell(nc, 1);
    for k = 1:n, trials(k).x = xs{k};  trials(k).y = ys{k};  trials(k).fs = fs;  trials(k).t = (0:numel(xs{k}) - 1)' / fs; end
    parfor c = 1:nc
        parts{c} = characteriseBiologicalNoise_v004(trials(edges(c) + 1:edges(c + 1)), Method="irasa", NHarmonics=nh, Verbose=false); %#ok<PFBNS>
    end
    nb = vertcat(parts{:});
end

function [rX, rY, f0] = pilotResid_local(x, y, fs, nHarm, edgeClip)
% subsamplingAnalysisPilot_v002 processTrial_local: f0 from the whole trial, clip, global harmonic template, residuals at the native rate.
    rX = [];  rY = [];  f0 = estimateF0_local(x(:), y(:), fs);  N = numel(x);
    if N - 2 * edgeClip < 40 || ~isfinite(f0) || f0 <= 0, return, end
    xC = x(edgeClip + 1:N - edgeClip);  yC = y(edgeClip + 1:N - edgeClip);  xC = xC(:);  yC = yC(:);
    t = (0:numel(xC) - 1)' / fs;  A = [ones(numel(t), 1), cos(2 * pi * (1:nHarm) .* f0 .* t), sin(2 * pi * (1:nHarm) .* f0 .* t)];
    rX = xC - A * (A \ xC);  rY = yC - A * (A \ yC);
end

function [a, s, bad] = armB_local(rX, rY, D, p, fs, minBins)
% Subsample the native-rate residuals and run IRASA; sigma is the residual SD (the Pilot script used the IRASA sigma).
    a = NaN;  s = NaN;  bad = false;
    if isempty(rX), return, end
    idx = p:D:numel(rX);  if numel(idx) < 8, return, end
    [aX, b1] = irasa_local(rX(idx), fs, minBins);  [aY, b2] = irasa_local(rY(idx), fs, minBins);  bad = b1 || b2;
    if isfinite(aX) && isfinite(aY), a = (aX + aY) / 2;  s = (std(rX(idx)) + std(rY(idx))) / 2; end
end

function [alpha, bad] = irasa_local(x, fs, minBins)
% iraAlphaSigma_v001 as in the Pilot script: band [1, min(20, fs/2 - 1)] Hz, N >= 30, at least minBins frequency bins.
    alpha = NaN;  bad = false;  N = numel(x);  fHigh = min(20, fs / 2 - 1);
    if N < 30 || fHigh <= 1 || (fHigh - 1) / (fs / N) < minBins, return, end
    try
        alpha = iraAlphaSigma_v001(x(:), fs, 1.0, fHigh, 1.1:0.05:1.9);
    catch
        bad = true;    % counted and logged per scheme
    end
end

function f0 = estimateF0_local(x, y, fs)
% Dominant FFT peak in [0.1, 5] Hz, as runLoopClosureFftnoise_v007 and the Pilot script.
    N = numel(x);  nfft = 2 ^ nextpow2(4 * N);  Xf = abs(fft(detrend(x, "linear"), nfft));  Yf = abs(fft(detrend(y, "linear"), nfft));
    fAx = (0:nfft - 1)' * fs / nfft;  idx = find(fAx > 0.1 & fAx < 5);  [mx, pkX] = max(Xf(idx));  [my, pkY] = max(Yf(idx));
    f0x = fAx(idx(pkX));  f0y = fAx(idx(pkY));
    if abs(f0x - f0y) / max(f0x, 0.01) < 0.2, f0 = mean([f0x, f0y]); elseif mx > my, f0 = f0x; else, f0 = f0y; end
end

function anchor_local(LG, file, ids, a, sg, tol)
% Fail loud unless the 240 Hz production arm reproduces the stored noise file on every trial.
    if ~isfile(file), file = fullfile(fileparts(mfilename("fullpath")), file); end
    if ~isfile(file), fail_local(LG, "ladder:anchorfile", "Anchor file missing: " + file); end
    bio = load(file, "bioResults").bioResults;  [ok, j] = ismember(ids, string(bio.trialID));
    if ~all(ok), fail_local(LG, "ladder:anchorids", sprintf("%d trials are not in the anchor file (set ANCHOR_FILE to '' to run without it)", sum(~ok))); end
    d = [abs(a - bio.ira_alphaMean(j)); abs(sg - bio.sigmaMean(j))];
    if any(isnan(a) ~= isnan(bio.ira_alphaMean(j))) || max(d(isfinite(d))) > tol
        fail_local(LG, "ladder:anchor", sprintf("Noise at 240 Hz does not reproduce the noise file (max diff %.3g)", max(d(isfinite(d)))));
    end
    say_local(LG, "Anchor passed: noise at 240 Hz reproduces %s on %d trials (max diff %.3g).\n", file, numel(ids), max(d(isfinite(d))));
end

function say_local(fid, varargin)
    fprintf(varargin{:});  fprintf(fid, varargin{:});
end

function fail_local(fid, id, msg)
    say_local(fid, "FAIL %s: %s\n", id, msg);  fclose(fid);  error(id, "%s", msg);
end

function v = ifelse_local(c, a, b)
    if c, v = a; else, v = b; end
end
