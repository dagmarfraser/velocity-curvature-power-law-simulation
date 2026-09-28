%% positiveControlAdditiveVsRetimed_v001
% Positive control and robustness for the generative-model check (Findings
% #222/#223, D20) that places Zarandi outside the forward model's validity
% domain. Session 125, 2026-09-28: the exclusion must rest on a validated
% estimator, not on a pattern.
%
%   0. SELF-CHECK  checkAdditiveVsRetimed_v003 reproduces v002 exactly on
%                  real data (Zarandi, per_subject). Fail loud otherwise.
%   1. CALIBRATE   synthetic trials of KNOWN construction, Zarandi-like
%                  (100 Hz, f0 0.65 Hz, 4:1 ellipse at 45 deg, alpha 3.18,
%                  sigma 0.477 cm, 20 s). Expect lambda ~ 0 for additive
%                  deviations on power-law timing (the forward map's own
%                  model) and lambda ~ 1 for ENACTED geometry (the path
%                  deviates, then is traversed at power-law speed along the
%                  path actually drawn), built independently of the check's
%                  own re-timer (dense arc-length traversal here, not
%                  retime_local).
%   2. PROCESSING  additive trials generated at 200 Hz, then exported to
%                  100 Hz four ways: anti-aliased resample; decimation with no
%                  filter; the authors' zero-phase Butterworth at normalised
%                  0.07 (7 Hz at 200 Hz) then decimation; event-driven
%                  sampling (timing jitter, 0.001 cm quantisation) then pchip
%                  onto a uniform 100 Hz grid. Plus enacted geometry resampled.
%                  Question: can export-style processing alone produce
%                  lambda > 1 (beyond the re-timed ceiling), as Zarandi shows?
%   3. ZARANDI     all trials, band cut at 2, 3 and 4 x f0; per-trial share
%                  beyond the ceiling; the 3 x f0 run must reproduce #223
%                  (subject-median lambda 1.94, lambdaHigh 2.13).
%
% Synthetic trials carry a small white sensor noise (0.003 cm) in every
% condition. Unit of analysis for synthetic conditions: trial (one trial per
% synthetic "subject"). Everything is re-runnable; nothing is hand-entered
% except the published values the self-checks compare against.
%
% OUTPUT: results/positiveControlAdditiveVsRetimed_v001_<stamp>.mat (summary
%   tables) plus the v002/v003 trial tables the check itself saves.
% USAGE: run from src/.  Cost: a few minutes on the M5 (parfor).
%
% Fraser, D.S. (2026)

%% CONFIG ---------------------------------------------------------------------
CFG.NSurr      = 10;          % as the full-corpus runner
CFG.NPerCond   = 20;          % synthetic trials per condition
CFG.FS         = 100;         % Zarandi export rate
CFG.FSHI       = 200;         % Zarandi reported rate
CFG.DUR        = 20;          % s
CFG.F0         = 0.65;        % Hz, Zarandi median f0 (checkEmpiricalTempoVsGridVGF_v001)
CFG.A          = 8;  CFG.B = 2;  CFG.TH = pi/4;   % cm; Zarandi template 16 x 4 cm at 45 deg
CFG.BETA       = 1/3;
CFG.ALPHA      = 3.18;        % Zarandi residual colour (v004)
CFG.SIGMA      = 0.477;       % cm, Zarandi residual sigma (4.77 mm)
CFG.SENSOR     = 0.003;       % cm, white sensor noise
CFG.JITTER     = 0.0015;      % s, event-timing jitter SD
CFG.QUANT      = 0.001;       % cm, position quantisation for the event-driven condition
CFG.BUTTER_WN  = 0.07;        % authors' stated (normalised) cut-off
CFG.SEED       = 20260928;
CFG.CAL_ADD    = [-0.30 0.30];   % calibration window, additive (median lambda)
CFG.CAL_RET    = [ 0.70 1.30];   % calibration window, enacted geometry
CFG.Z_LAMBDA   = 1.94;  CFG.Z_LAMBDAHI = 2.13;  CFG.Z_TOL = 0.005;   % Finding #223
CFG.CBM        = [2 3 4];

srcDir = fileparts(mfilename("fullpath"));
cd(srcDir);
addpath(genpath(fullfile(srcDir, "functions")));
stamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));

%% 0. SELF-CHECK: v003 == v002 on real data ------------------------------------
fprintf("\n=== 0. Self-check: v003 reproduces v002 (Zarandi, per_subject) ===\n");
T2 = checkAdditiveVsRetimed_v002(Datasets="Zarandi", TrialSelection="per_subject", NSurr=5);
T3 = checkAdditiveVsRetimed_v003(Datasets="Zarandi", TrialSelection="per_subject", NSurr=5, SaveTag="selfcheck_");
cols = ["lambda","lambdaLow","lambdaHigh","rhoReal","rhoAddMed","rhoRetMed","f0"];
for c = cols
    if ~isequaln(T2.(c), T3.(c))
        error("posCtrl:SelfCheck", "%s", "v003 does not reproduce v002 on column " + c + ".");
    end
end
if ~isequal(T2.status, T3.status), error("posCtrl:SelfCheck", "%s", "v003 status differs from v002."); end
fprintf("Self-check passed: %d trials, columns %s identical.\n", height(T3), strjoin(cols, ", "));

%% 1-2. SYNTHETIC TRIALS ------------------------------------------------------------
rng(CFG.SEED, "twister");
conds = ["A_additive", "R_enacted", "A_resample200", "A_decimate200", ...
         "A_butter007_200", "A_eventInterp200", "R_enacted_resample200"];
S = struct("x", {}, "y", {}, "trialID", {}, "subjectID", {}, "dataset", {});
N  = round(CFG.DUR * CFG.FS);  NH = round(CFG.DUR * CFG.FSHI);
for c = conds
    for k = 1:CFG.NPerCond
        switch c
            case "A_additive"
                [x, y] = additive_local(CFG, CFG.FS, N);
            case "R_enacted"
                [x, y] = enacted_local(CFG, CFG.FS, N);
            case "A_resample200"
                [x, y] = additive_local(CFG, CFG.FSHI, NH);
                x = resample(x, 1, 2); y = resample(y, 1, 2);
            case "A_decimate200"
                [x, y] = additive_local(CFG, CFG.FSHI, NH);
                x = x(1:2:end); y = y(1:2:end);
            case "A_butter007_200"
                [x, y] = additive_local(CFG, CFG.FSHI, NH);
                [bb, aa] = butter(2, CFG.BUTTER_WN);
                x = filtfilt(bb, aa, x); y = filtfilt(bb, aa, y);
                x = x(1:2:end); y = y(1:2:end);
            case "A_eventInterp200"
                [x, y] = eventInterp_local(CFG, N);
            case "R_enacted_resample200"
                [x, y] = enacted_local(CFG, CFG.FSHI, NH);
                x = resample(x, 1, 2); y = resample(y, 1, 2);
        end
        if numel(x) ~= N || any(~isfinite([x; y]))
            error("posCtrl:BadTrial", "%s", sprintf("%s #%d: %d samples (expected %d) or non-finite.", c, k, numel(x), N));
        end
        S(end+1) = struct("x", x(:), "y", y(:), "trialID", sprintf("%s_%02d", c, k), ...
            "subjectID", sprintf("%s_%02d", c, k), "dataset", c); %#ok<SAGROW>
    end
end
fprintf("\n=== 1-2. %d synthetic trials, %d conditions ===\n", numel(S), numel(conds));
TS = checkAdditiveVsRetimed_v003(Trials=S, TrialsFS=CFG.FS, NSurr=CFG.NSurr, UseParfor=true, SaveTag="posctrl_");
Csum = summarise_local(TS, conds);
disp(Csum);

mAdd = Csum.lambda_med(Csum.condition == "A_additive");
mRet = Csum.lambda_med(Csum.condition == "R_enacted");
calOK = mAdd >= CFG.CAL_ADD(1) && mAdd <= CFG.CAL_ADD(2) && mRet >= CFG.CAL_RET(1) && mRet <= CFG.CAL_RET(2);
fprintf("\nCALIBRATION: additive median lambda %.3f (window %s); enacted %.3f (window %s) -> %s\n", ...
    mAdd, mat2str(CFG.CAL_ADD), mRet, mat2str(CFG.CAL_RET), string(ternary_local(calOK, "PASS", "FAIL")));

%% 3. ZARANDI ROBUSTNESS ----------------------------------------------------------
fprintf("\n=== 3. Zarandi, all trials, band cut x f0 in %s ===\n", mat2str(CFG.CBM));
Zsum = table();
for cbm = CFG.CBM
    TZ = checkAdditiveVsRetimed_v003(Datasets="Zarandi", TrialSelection="all", NSurr=CFG.NSurr, ...
        CycleBandMult=cbm, UseParfor=true, SaveTag=sprintf("zarandi_cbm%d_", cbm));
    ok = isfinite(TZ.lambda); okH = isfinite(TZ.lambdaHigh);
    sL = splitapply(@median, TZ.lambda(ok),  findgroups(TZ.subjectID(ok)));
    sH = splitapply(@median, TZ.lambdaHigh(okH), findgroups(TZ.subjectID(okH)));
    Zsum = [Zsum; table(cbm, height(TZ), sum(ok), median(sL), sum(sL > 1), numel(sL), mean(TZ.lambda(ok) > 1), ...
        median(sH), sum(sH > 1), numel(sH), mean(TZ.lambdaHigh(okH) > 1), ...
        'VariableNames', ["cbm","nTrials","nFinite","subjMed_lambda","nSubjAbove1","nSubj", ...
        "shareTrialsAbove1","subjMed_lambdaHigh","nSubjHighAbove1","nSubjHigh","shareTrialsHighAbove1"])]; %#ok<AGROW>
end
disp(Zsum);
z3 = Zsum(Zsum.cbm == 3, :);
if abs(z3.subjMed_lambda - CFG.Z_LAMBDA) > CFG.Z_TOL || abs(z3.subjMed_lambdaHigh - CFG.Z_LAMBDAHI) > CFG.Z_TOL
    error("posCtrl:Z223", "%s", sprintf("3 x f0 run gives lambda %.3f / lambdaHigh %.3f, not #223's %.2f / %.2f.", ...
        z3.subjMed_lambda, z3.subjMed_lambdaHigh, CFG.Z_LAMBDA, CFG.Z_LAMBDAHI));
end
fprintf("Self-check passed: 3 x f0 reproduces #223 (%.3f, %.3f).\n", z3.subjMed_lambda, z3.subjMed_lambdaHigh);

%% SAVE, then fail loud on calibration --------------------------------------------
resDir = fullfile(fileparts(srcDir), "results");
if ~isfolder(resDir), error("posCtrl:NoResults", "%s", "FAILED PATH: " + resDir); end
outFile = fullfile(resDir, "positiveControlAdditiveVsRetimed_v001_" + stamp + ".mat");
save(outFile, "CFG", "Csum", "Zsum", "calOK");
fprintf("\nSaved: %s\n", outFile);
if ~calOK
    error("posCtrl:Calibration", "%s", "Calibration FAILED: the estimator does not recover known constructions. Do not use lambda for exclusion until resolved.");
end

%% ============================== local functions ==================================
function [x, y] = additive_local(CFG, fs, n)
% Power-law ellipse template + additive coloured deviations + white sensor noise.
[x, y] = template_local(CFG, fs, n);
x = x + coloured_local(n, CFG.ALPHA, CFG.SIGMA) + CFG.SENSOR*randn(n, 1);
y = y + coloured_local(n, CFG.ALPHA, CFG.SIGMA) + CFG.SENSOR*randn(n, 1);
end

function [x, y] = template_local(CFG, fs, n)
% Arc-length power-law ellipse (same construction as the runner's template).
phi  = linspace(0, 2*pi, 10000)';
dsdp = sqrt((CFG.A*sin(phi)).^2 + (CFG.B*cos(phi)).^2);
kp   = (CFG.A*CFG.B) ./ max((CFG.A^2*sin(phi).^2 + CFG.B^2*cos(phi).^2).^1.5, eps);
cumT = cumsum(dsdp .* kp.^CFG.BETA); cumT = cumT / cumT(end);
tN   = min(mod((0:n-1)'/fs*CFG.F0, 1), 1 - 1e-9);
pA   = interp1(cumT, phi, tN, "linear", "extrap");
xE = CFG.A*cos(pA); yE = CFG.B*sin(pA);
x = xE*cos(CFG.TH) - yE*sin(CFG.TH);  y = xE*sin(CFG.TH) + yE*cos(CFG.TH);
end

function [x, y] = enacted_local(CFG, fs, n)
% Enacted geometry: the PATH deviates (coloured deviations indexed by path
% parameter), then is traversed at v = K*kappa_path^-beta along the path
% actually drawn (dense arc-length integration), sampled at fs; + sensor noise.
nCyc = CFG.F0 * CFG.DUR;
nP   = round(nCyc * 4000);
phi  = linspace(0, 2*pi*nCyc, nP)';
xE = CFG.A*cos(phi); yE = CFG.B*sin(phi);
px = xE*cos(CFG.TH) - yE*sin(CFG.TH) + coloured_local(nP, CFG.ALPHA, CFG.SIGMA);
py = xE*sin(CFG.TH) + yE*cos(CFG.TH) + coloured_local(nP, CFG.ALPHA, CFG.SIGMA);
xp = gradient(px); yp = gradient(py); xpp = gradient(xp); ypp = gradient(yp);
k  = abs(xp.*ypp - yp.*xpp) ./ max((xp.^2 + yp.^2).^1.5, eps);
k  = max(k, prctile(k(k > 0), 1));
w  = k.^(-CFG.BETA);
ds = max(hypot(diff(px), diff(py)), eps);
t  = [0; cumsum(ds ./ (0.5*(w(1:end-1) + w(2:end))))];
t  = t * CFG.DUR / t(end);
tq = (0:n-1)'/fs;
x  = interp1(t, px, tq, "pchip") + CFG.SENSOR*randn(n, 1);
y  = interp1(t, py, tq, "pchip") + CFG.SENSOR*randn(n, 1);
end

function [x, y] = eventInterp_local(CFG, n)
% Additive trial evaluated at jittered 200 Hz event times, quantised, then
% interpolated (pchip) onto a uniform 100 Hz grid.
fsD = 1000; nD = round(CFG.DUR * fsD);
[xD, yD] = template_local(CFG, fsD, nD);
xD = xD + coloured_local(nD, CFG.ALPHA, CFG.SIGMA);
yD = yD + coloured_local(nD, CFG.ALPHA, CFG.SIGMA);
tD = (0:nD-1)'/fsD;
tE = (0:1/CFG.FSHI:CFG.DUR - 1/fsD)' + CFG.JITTER*randn(round(CFG.DUR*CFG.FSHI), 1);
tE = sort(min(max(tE, 0), tD(end)));
[tE, iu] = unique(tE);
xE = round(interp1(tD, xD, tE, "pchip") / CFG.QUANT) * CFG.QUANT + CFG.SENSOR*randn(numel(tE), 1);
yE = round(interp1(tD, yD, tE, "pchip") / CFG.QUANT) * CFG.QUANT + CFG.SENSOR*randn(numel(tE), 1);
if numel(iu) < 0.99*round(CFG.DUR*CFG.FSHI)
    error("posCtrl:EventTimes", "%s", "Too many duplicate event times after jitter.");
end
tq = (0:n-1)'/CFG.FS;
x = interp1(tE, xE, tq, "pchip", "extrap"); y = interp1(tE, yE, tq, "pchip", "extrap");
end

function z = coloured_local(n, alpha, sigma)
% 1/f^alpha Gaussian noise by spectral shaping, zero mean, SD = sigma.
m  = 2^nextpow2(2*n);
f  = [0, 1:m/2, -(m/2-1):-1]' / m;
amp = zeros(m, 1); nz = f ~= 0; amp(nz) = abs(f(nz)).^(-alpha/2);
Z  = (randn(m, 1) + 1i*randn(m, 1)) .* amp;
z  = real(ifft(Z)); z = z(1:n); z = z - mean(z);
s  = std(z);
if ~(s > 0), error("posCtrl:Noise", "%s", "Degenerate coloured noise."); end
z  = sigma * z / s;
end

function C = summarise_local(T, conds)
C = table();
for c = conds
    m = T.dataset == c; L = T.lambda(m); LL = T.lambdaLow(m); LH = T.lambdaHigh(m);
    st = T.status(m & ~isfinite(T.lambda));
    reasons = "";
    if ~isempty(st), [u,~,g] = unique(st); reasons = strjoin(u + "=" + string(accumarray(g,1)), "; "); end
    C = [C; table(c, sum(m), sum(isfinite(L)), median(L,"omitnan"), prctile(L,25), prctile(L,75), mean(L(isfinite(L)) > 1), ...
        median(LL,"omitnan"), median(LH,"omitnan"), mean(LH(isfinite(LH)) > 1), reasons, ...
        'VariableNames', ["condition","n","nFinite","lambda_med","lambda_q25","lambda_q75","shareAbove1", ...
        "lambdaLow_med","lambdaHigh_med","shareHighAbove1","nanReasons"])]; %#ok<AGROW>
end
end

function out = ternary_local(c, a, b)
if c, out = a; else, out = b; end
end
