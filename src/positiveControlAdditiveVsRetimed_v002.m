%% positiveControlAdditiveVsRetimed_v002
% Calibration of the generative-model check (Findings #222/#223, D20), v002.
% v001 (2026-09-28) FAILED its pre-set windows: additive synthetic trials read
% lambda 0.36 (not 0), enacted geometry 1.36 (not 1), a common offset with the
% intended ~1.0 separation; no processed-additive trial exceeded 1; Zarandi
% 14/14 subjects and 100% of evaluable trials above 1 at 2, 3, 4 x f0.
%
% v002 isolates the offset and calibrates:
%   - NOISE CONSTRUCTION: random-phase 1/f^alpha (v001) versus shaped_xu
%     (empirical magnitude on Xu's fractional phase), which is how the check's
%     own anchors are built (generateLoopClosureNoise_v003). Hypothesis: the
%     offset is phase, so shaped_xu additive trials should read ~0.
%   - SENSOR NOISE: each additive and enacted world with and without 0.003 cm
%     white noise.
%   - POINTER ACCELERATION: additive (shaped_xu + sensor) trials passed through
%     a speed-dependent gain on displacement (as OS/driver pointer
%     acceleration, e.g. Dhieb's EPP). Exploratory: can this processing
%     produce lambda > 1?
%   - CALIBRATED lambda*: lambda* = (lambda - medA) / (medR - medA), with medA,
%     medR the matched (shaped_xu + sensor; enacted + sensor) condition medians.
%     Applied to Zarandi (geometry-matched, so valid) and, for orientation
%     only, to every dataset in the saved full-corpus table (their geometry,
%     fs and colour differ from the calibration's, so their lambda* is
%     approximate and flagged as such).
%
% Pre-set calibration criterion (matched conditions): |medA| <= 0.30 and
% 0.70 <= medR <= 1.30. Written before the run; FAIL is reported and saved,
% then raised as an error (Fail Loud).
%
% Synthetic geometry as v001 (Zarandi-like: 100 Hz, f0 0.65 Hz, 16 x 4 cm
% ellipse at 45 deg, alpha 3.18, sigma 0.477 cm, 20 s). Unit: trial.
% OUTPUT: results/positiveControlAdditiveVsRetimed_v002_<stamp>.mat
% USAGE: run from src/.
%
% Fraser, D.S. (2026)

%% CONFIG ---------------------------------------------------------------------
CFG.NSurr      = 10;
CFG.NPerCond   = 20;
CFG.FS         = 100;
CFG.DUR        = 20;
CFG.F0         = 0.65;
CFG.A          = 8;  CFG.B = 2;  CFG.TH = pi/4;
CFG.BETA       = 1/3;
CFG.ALPHA      = 3.18;
CFG.SIGMA      = 0.477;
CFG.SENSOR     = 0.003;
CFG.ACCEL_G    = 1.0;         % pointer acceleration: gain = 1 + G*tanh(speed/median speed)
CFG.SEED       = 20260929;
CFG.CAL_ADD    = [-0.30 0.30];
CFG.CAL_RET    = [ 0.70 1.30];
CFG.ZARANDI_TABLE = "checkAdditiveVsRetimed_v003_zarandi_cbm3_20260928_173727.mat";   % v001 step 3
CFG.CORPUS_TABLE  = "checkAdditiveVsRetimed_v002_20260922_205130.mat";                % full corpus (#222)
CFG.CORPUS_N      = 3945;

srcDir = fileparts(mfilename("fullpath"));
cd(srcDir);
addpath(genpath(fullfile(srcDir, "functions")));
stamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));

%% SYNTHETIC TRIALS -------------------------------------------------------------
rng(CFG.SEED, "twister");
conds = ["A_fft_clean", "A_fft_sensor", "A_sxu_clean", "A_sxu_sensor", ...
         "R_enacted_clean", "R_enacted_sensor", "A_sxu_pointerAccel"];
S = struct("x", {}, "y", {}, "trialID", {}, "subjectID", {}, "dataset", {});
N = round(CFG.DUR * CFG.FS);
for c = conds
    for k = 1:CFG.NPerCond
        switch c
            case "A_fft_clean",        [x, y] = additive_local(CFG, N, "fft", 0);
            case "A_fft_sensor",       [x, y] = additive_local(CFG, N, "fft", CFG.SENSOR);
            case "A_sxu_clean",        [x, y] = additive_local(CFG, N, "sxu", 0);
            case "A_sxu_sensor",       [x, y] = additive_local(CFG, N, "sxu", CFG.SENSOR);
            case "R_enacted_clean",    [x, y] = enacted_local(CFG, N, 0);
            case "R_enacted_sensor",   [x, y] = enacted_local(CFG, N, CFG.SENSOR);
            case "A_sxu_pointerAccel"
                [x, y] = additive_local(CFG, N, "sxu", CFG.SENSOR);
                [x, y] = pointerAccel_local(x, y, CFG.FS, CFG.ACCEL_G);
        end
        if numel(x) ~= N || any(~isfinite([x; y]))
            error("posCtrl2:BadTrial", "%s", sprintf("%s #%d bad (%d samples).", c, k, numel(x)));
        end
        S(end+1) = struct("x", x(:), "y", y(:), "trialID", sprintf("%s_%02d", c, k), ...
            "subjectID", sprintf("%s_%02d", c, k), "dataset", c); %#ok<SAGROW>
    end
end
fprintf("\n=== %d synthetic trials, %d conditions ===\n", numel(S), numel(conds));
TS = checkAdditiveVsRetimed_v003(Trials=S, TrialsFS=CFG.FS, NSurr=CFG.NSurr, UseParfor=true, SaveTag="posctrl2_");
Csum = summarise_local(TS, conds);
disp(Csum(:, 1:end-1));
for i = 1:height(Csum)
    if strlength(Csum.nanReasons(i)) > 0, fprintf("  NaN %-20s %s\n", Csum.condition(i), Csum.nanReasons(i)); end
end

medA = Csum.lambda_med(Csum.condition == "A_sxu_sensor");
medR = Csum.lambda_med(Csum.condition == "R_enacted_sensor");
calOK = abs(medA) <= CFG.CAL_ADD(2) && medR >= CFG.CAL_RET(1) && medR <= CFG.CAL_RET(2);
fprintf("\nCALIBRATION (matched): additive %.3f, enacted %.3f, separation %.3f -> %s\n", ...
    medA, medR, medR - medA, string(ternary_local(calOK, "PASS", "FAIL")));
if ~(medR - medA > 0)
    error("posCtrl2:NoSeparation", "%s", "Enacted does not exceed additive: lambda* undefined.");
end
toStar = @(L) (L - medA) / (medR - medA);

%% CALIBRATED lambda* -----------------------------------------------------------
resSrc = fullfile(srcDir, "results");
fz = fullfile(resSrc, CFG.ZARANDI_TABLE);
if ~isfile(fz), error("posCtrl2:NoZ", "%s", "FAILED PATH: " + fz); end
TZ = load(fz, "T").T;
okZ = isfinite(TZ.lambda);
sZ  = splitapply(@median, TZ.lambda(okZ), findgroups(TZ.subjectID(okZ)));
Zstar = table(numel(sZ), median(sZ), median(toStar(sZ)), prctile(toStar(sZ), 25), prctile(toStar(sZ), 75), ...
    sum(toStar(sZ) > 1), mean(toStar(TZ.lambda(okZ)) > 1), ...
    'VariableNames', ["nSubj","subjMed_lambda","subjMed_lambdaStar","q25_star","q75_star","nSubjStarAbove1","shareTrialsStarAbove1"]);
fprintf("\nZarandi (geometry-matched calibration):\n"); disp(Zstar);

fc = fullfile(resSrc, CFG.CORPUS_TABLE);
if ~isfile(fc), error("posCtrl2:NoCorpus", "%s", "FAILED PATH: " + fc); end
TC = load(fc, "T").T;
if height(TC) ~= CFG.CORPUS_N
    error("posCtrl2:CorpusN", "%s", sprintf("Corpus table has %d rows, expected %d (#222).", height(TC), CFG.CORPUS_N));
end
dsets = unique(TC.dataset, "stable");
Dstar = table();
for d = dsets'
    m = TC.dataset == d & isfinite(TC.lambda);
    sL = splitapply(@median, TC.lambda(m), findgroups(TC.subjectID(m)));
    Dstar = [Dstar; table(d, numel(sL), median(sL), median(toStar(sL)), ...
        'VariableNames', ["dataset","nSubj","subjMed_lambda","subjMed_lambdaStar_APPROX"])]; %#ok<AGROW>
end
fprintf("All datasets (lambda* APPROXIMATE outside Zarandi: calibration is Zarandi-geometry):\n"); disp(Dstar);

%% SAVE, then fail loud ---------------------------------------------------------
resDir = fullfile(fileparts(srcDir), "results");
if ~isfolder(resDir), error("posCtrl2:NoResults", "%s", "FAILED PATH: " + resDir); end
outFile = fullfile(resDir, "positiveControlAdditiveVsRetimed_v002_" + stamp + ".mat");
save(outFile, "CFG", "Csum", "medA", "medR", "calOK", "Zstar", "Dstar");
fprintf("\nSaved: %s\n", outFile);
if ~calOK
    error("posCtrl2:Calibration", "%s", sprintf("Calibration FAILED (additive %.3f, enacted %.3f). lambda* above is offset-corrected but the pre-set criterion was not met.", medA, medR));
end

%% ============================== local functions ==================================
function [x, y] = additive_local(CFG, n, kind, sensor)
[x, y] = template_local(CFG, CFG.FS, n);
rMaj = coloured_local(n, CFG.ALPHA, CFG.SIGMA); rMin = coloured_local(n, CFG.ALPHA, CFG.SIGMA);
switch kind
    case "fft"
        nMaj = rMaj; nMin = rMin;
    case "sxu"   % the check's own anchor construction, reference = rMaj/rMin
        [nMaj, nMin] = generateLoopClosureNoise_v003("shaped_xu", rMaj, rMin, CFG.FS, CFG.ALPHA, CFG.ALPHA);
        if any(~isfinite([nMaj; nMin])), error("posCtrl2:SXU", "%s", "shaped_xu returned non-finite noise."); end
end
x = x + nMaj*cos(CFG.TH) - nMin*sin(CFG.TH) + sensor*randn(n, 1);
y = y + nMaj*sin(CFG.TH) + nMin*cos(CFG.TH) + sensor*randn(n, 1);
end

function [x, y] = template_local(CFG, fs, n)
phi  = linspace(0, 2*pi, 10000)';
dsdp = sqrt((CFG.A*sin(phi)).^2 + (CFG.B*cos(phi)).^2);
kp   = (CFG.A*CFG.B) ./ max((CFG.A^2*sin(phi).^2 + CFG.B^2*cos(phi).^2).^1.5, eps);
cumT = cumsum(dsdp .* kp.^CFG.BETA); cumT = cumT / cumT(end);
tN   = min(mod((0:n-1)'/fs*CFG.F0, 1), 1 - 1e-9);
pA   = interp1(cumT, phi, tN, "linear", "extrap");
xE = CFG.A*cos(pA); yE = CFG.B*sin(pA);
x = xE*cos(CFG.TH) - yE*sin(CFG.TH);  y = xE*sin(CFG.TH) + yE*cos(CFG.TH);
end

function [x, y] = enacted_local(CFG, n, sensor)
nCyc = CFG.F0 * CFG.DUR;  nP = round(nCyc * 4000);
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
tq = (0:n-1)'/CFG.FS;
x  = interp1(t, px, tq, "pchip") + sensor*randn(n, 1);
y  = interp1(t, py, tq, "pchip") + sensor*randn(n, 1);
end

function [x, y] = pointerAccel_local(x, y, fs, G)
% Screen displacement = gain(speed) x device displacement, gain rising with speed.
d  = [diff(x), diff(y)];
s  = vecnorm(d, 2, 2) * fs;
g  = 1 + G * tanh(s / median(s));
d  = d .* g;
x  = x(1) + [0; cumsum(d(:,1))];  y = y(1) + [0; cumsum(d(:,2))];
end

function z = coloured_local(n, alpha, sigma)
m  = 2^nextpow2(2*n);
f  = [0, 1:m/2, -(m/2-1):-1]' / m;
amp = zeros(m, 1); nz = f ~= 0; amp(nz) = abs(f(nz)).^(-alpha/2);
Z  = (randn(m, 1) + 1i*randn(m, 1)) .* amp;
z  = real(ifft(Z)); z = z(1:n); z = z - mean(z);
s  = std(z);
if ~(s > 0), error("posCtrl2:Noise", "%s", "Degenerate coloured noise."); end
z  = sigma * z / s;
end

function C = summarise_local(T, conds)
C = table();
for c = conds
    m = T.dataset == c; L = T.lambda(m); LH = T.lambdaHigh(m);
    st = T.status(m & ~isfinite(T.lambda)); reasons = "";
    if ~isempty(st), [u,~,g] = unique(st); reasons = strjoin(u + "=" + string(accumarray(g,1)), "; "); end
    C = [C; table(c, sum(m), sum(isfinite(L)), median(L,"omitnan"), prctile(L,25), prctile(L,75), ...
        mean(L(isfinite(L)) > 1), median(LH,"omitnan"), reasons, ...
        'VariableNames', ["condition","n","nFinite","lambda_med","lambda_q25","lambda_q75", ...
        "shareAbove1","lambdaHigh_med","nanReasons"])]; %#ok<AGROW>
end
end

function out = ternary_local(c, a, b)
if c, out = a; else, out = b; end
end
