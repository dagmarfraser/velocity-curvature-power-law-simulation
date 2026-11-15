%% positiveControlMatchedByDataset_v005  v002 with ONE SEED PER DATASET (2026-11-07, Dagmar: option A). v002 drew every dataset from one rng(20261031) stream, Fraser first, so a change in Fraser's trial count moved every other dataset's controls by Monte Carlo order alone. v005 seeds each dataset as rng(SEED + k, 'twister'), k = its position in CFG.Datasets, before its trial draw; nothing else changes. Writes results/positiveControlMatchedByDataset_v005_<stamp>.mat.
%% (v002 header follows)  Deep Review GD17/lambda*: per-dataset matched controls for the generative-model check (D20).
% WHY: lambda* = (lambda - medA)/(medR - medA) was normalised by ONE control built on Zarandi geometry, rate and noise colour. v001 built matched controls per
% dataset but ASSUMED a sensor floor and missed its windows outside Fraser. v002 replaces the assumed floor with MEASURED bounds (checkSensorFloor_v001):
% lower bound = quantisation SD (grid step/sqrt(12)); upper bound = high-frequency PSD level above 25 Hz (a plateau only for Fraser).
% DECIDED BEFORE THE RUN: each real trial's IRASA sigma and the shaped_xu reference both come from the whole real residual, which already contains any sensor
% noise, so adding more would double-count. The PRIMARY matched pair is therefore A_sxu_clean / R_enacted_clean (no added sensor noise), tested against the
% unchanged windows |medA| <= 0.30 and 0.70 <= medR <= 1.30. SECONDARY: the same pair with the quantisation floor added (A_sxu_quant / R_enacted_quant).
% STRESS: R_enacted_hf adds the high-frequency level (an upper bound). Everything else as v001: each real trial supplies ellipse, f0, rate, length (the
% trimmed residual is wrapped to N), betaTpl and per-axis IRASA (alpha, sigma); analysed by the unchanged checkAdditiveVsRetimed_v003.
% Run from src/. Writes results/positiveControlMatchedByDataset_v005_<stamp>.mat.   Fraser, D.S. (2026)  v002
CFG.Datasets = ["Fraser","Cook CTRL","Cook ASD","Hickman PLAC","Hickman HALO","Zarandi","Dhieb"];
CFG.NTrials = 60;  CFG.NSurr = 10;  CFG.SEED = 20261031; 
CFG.CAL_ADD = [-0.30 0.30];  CFG.CAL_RET = [0.70 1.30];  CFG.RULE = 0.5;
CFG.FLOOR_FILE = "checkSensorFloor_v001.mat";
CFG.CORPUS_TABLE = "checkAdditiveVsRetimed_v002_20260922_205130.mat";
srcDir = fileparts(mfilename("fullpath"));  cd(srcDir);  addpath(genpath(fullfile(srcDir, "functions")));
stamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
FL = load(fullfile(fileparts(srcDir), "results", CFG.FLOOR_FILE), "out").out;
TC = load(fullfile(srcDir, "results", CFG.CORPUS_TABLE), "T").T;
out = table();  perTrial = table();
for dName = CFG.Datasets
    rng(CFG.SEED + find(CFG.Datasets == dName), "twister");   % v005: one seed per dataset
    cfg = datasetConfig_local(dName, srcDir);  FS = cfg.FS;
    trials = cfg.importFn();  tids = string({trials.trialID}');
    cand = find(TC.dataset == dName & isfinite(TC.lambda));
    if numel(cand) > CFG.NTrials, cand = cand(randperm(numel(cand), CFG.NTrials)); end
    S = struct("x", {}, "y", {}, "trialID", {}, "subjectID", {}, "dataset", {});  chk = [];
    for q = cand'
        j = find(tids == string(TC.trialID(q)), 1);
        if isempty(j), error("ctrlM:NoTrial", "%s", "Trial not in import: " + string(TC.trialID(q))); end
        x = double(trials(j).x(:));  y = double(trials(j).y(:));  N = numel(x);
        f0 = estimateF0_local(x, y, FS);
        if abs(f0 - TC.f0(q)) > 1e-9, error("ctrlM:F0", "%s", sprintf("f0 mismatch for %s: %.6f vs saved %.6f", string(TC.trialID(q)), f0, TC.f0(q))); end
        [xRes, yRes] = templateSubtract_local(x, y, FS, f0, 4, 4);
        xc = x - mean(x); yc = y - mean(y);  [V, D] = eig(cov(xc, yc));
        [lams, ord] = sort(diag(D), "descend"); V = V(:, ord);
        a = sqrt(2*max(lams(1), 0)); b = sqrt(2*max(lams(2), 0)); th = atan2(V(2,1), V(1,1));
        rMaj =  xRes*cos(th) + yRes*sin(th);  rMin = -xRes*sin(th) + yRes*cos(th);
        fHi = min(20, FS/2 - 1);  hset = 1.1:0.05:1.9;
        [aMaj, sMaj] = iraAlphaSigma_v001(rMaj, FS, 1.0, fHi, hset);  [aMin, sMin] = iraAlphaSigma_v001(rMin, FS, 1.0, fHi, hset);
        if ~all(isfinite([aMaj aMin sMaj sMin])), continue, end
        chk(end+1, :) = [sMaj/std(rMaj), sMin/std(rMin)]; %#ok<AGROW>
        beta = TC.betaTpl(q);  fr = FL(FL.dataset == dName, :);  sQ = fr.quantSD;  sH = fr.hfFloorSD;
        if isempty(fr) || ~isfinite(sQ) || ~isfinite(sH), error("ctrlM:NoFloor", "%s", "No sensor floor for " + dName); end
        [xT, yT] = ellipse_local(a, b, th, f0, FS, beta, N);
        wrap = mod((0:N-1)', numel(rMaj)) + 1;   % templateSubtract_local returns the residual trimmed to length M < N; wrap it to N so the additive trial has the real trial's length
        [nMaj, nMin] = generateLoopClosureNoise_v003("shaped_xu", rMaj(wrap), rMin(wrap), FS, aMaj, aMin);
        if any(~isfinite([nMaj; nMin])), continue, end
        if ~isequal(size(xT), size(nMaj), size(nMin)), error("ctrlM:Size", "%s", sprintf("size mismatch N=%d xT %s nMaj %s nMin %s", N, mat2str(size(xT)), mat2str(size(nMaj)), mat2str(size(nMin)))); end
        xA0 = xT + nMaj*cos(th) - nMin*sin(th);  yA0 = yT + nMaj*sin(th) + nMin*cos(th);
        [xE0, yE0] = enacted_local(a, b, th, f0, FS, N, beta, aMaj, aMin, sMaj, sMin, 0);
        [xEQ, yEQ] = enacted_local(a, b, th, f0, FS, N, beta, aMaj, aMin, sMaj, sMin, sQ);
        [xEH, yEH] = enacted_local(a, b, th, f0, FS, N, beta, aMaj, aMin, sMaj, sMin, sH);
        id = string(TC.trialID(q));  sb = string(TC.subjectID(q));
        add = @(x, y, c) struct("x", x, "y", y, "trialID", id, "subjectID", sb, "dataset", c);
        S(end+1) = add(xA0, yA0, "A_sxu_clean"); %#ok<AGROW>
        S(end+1) = add(xA0 + sQ*randn(N, 1), yA0 + sQ*randn(N, 1), "A_sxu_quant"); %#ok<AGROW>
        S(end+1) = add(xE0, yE0, "R_enacted_clean"); %#ok<AGROW>
        S(end+1) = add(xEQ, yEQ, "R_enacted_quant"); %#ok<AGROW>
        S(end+1) = add(xEH, yEH, "R_enacted_hf"); %#ok<AGROW>
    end
    fprintf("\n=== %s: %d source trials -> %d synthetic; IRASA sigma / std(residual), median %.2f ===\n", dName, numel(cand), numel(S), median(chk(:)));
    TS = checkAdditiveVsRetimed_v003(Trials=S, TrialsFS=FS, NSurr=CFG.NSurr, UseParfor=true, SaveTag="ctrlM_");
    g = @(c) median(TS.lambda(TS.dataset == c & isfinite(TS.lambda)), "omitnan");
    n = @(c) sum(TS.dataset == c & isfinite(TS.lambda));
    medA = g("A_sxu_clean");  medR = g("R_enacted_clean");  medAq = g("A_sxu_quant");  medRq = g("R_enacted_quant");  medRh = g("R_enacted_hf");
    calOK = abs(medA) <= CFG.CAL_ADD(2) && medR >= CFG.CAL_RET(1) && medR <= CFG.CAL_RET(2);
    calOKq = abs(medAq) <= CFG.CAL_ADD(2) && medRq >= CFG.CAL_RET(1) && medRq <= CFG.CAL_RET(2);
    mm = TC.dataset == dName & isfinite(TC.lambda);
    sL = splitapply(@median, TC.lambda(mm), findgroups(TC.subjectID(mm)));
    lamSubj = median(sL);  lamStar = (lamSubj - medA) / (medR - medA);
    out = [out; table(dName, numel(cand), n("A_sxu_clean"), n("R_enacted_clean"), medA, medR, calOK, medAq, medRq, calOKq, medRh, lamSubj, lamStar, lamStar > CFG.RULE, ...
        'VariableNames', ["dataset","nSource","nFiniteA","nFiniteR","medA","medR","calOK","medA_quant","medR_quant","calOK_quant","medR_hf","subjMed_lambda","lambdaStar_own","outside"])]; %#ok<AGROW>
    TS.srcDataset = repmat(dName, height(TS), 1);
    perTrial = [perTrial; TS(:, ["srcDataset","dataset","trialID","subjectID","lambda","lambdaHigh","status"])]; %#ok<AGROW>
    disp(out(end, :));
end
disp(out);
save(fullfile(fileparts(srcDir), "results", "positiveControlMatchedByDataset_v005_" + stamp + ".mat"), "CFG", "out", "perTrial");
if any(~out.calOK), warning("ctrlM:Calibration", "%s", "Pre-set windows missed for: " + strjoin(out.dataset(~out.calOK), ", ") + " (reported as run)."); end

function [x, y] = ellipse_local(a, b, th, f0, FS, beta, M)
phi  = linspace(0, 2*pi, 10000)';
dsdp = sqrt((a*sin(phi)).^2 + (b*cos(phi)).^2);
kp   = (a*b) ./ max((a^2*sin(phi).^2 + b^2*cos(phi).^2).^1.5, eps);
cumT = cumsum(dsdp .* kp.^beta); cumT = cumT / cumT(end);
tN   = min(mod((0:M-1)'/FS*f0, 1), 1 - 1e-9);
pA   = interp1(cumT, phi, tN, "linear", "extrap");
xE = a*cos(pA); yE = b*sin(pA);
x = xE*cos(th) - yE*sin(th);  y = xE*sin(th) + yE*cos(th);
end

function [x, y] = enacted_local(a, b, th, f0, FS, n, beta, aMaj, aMin, sMaj, sMin, sensor)
% Noise enters the drawn path; speed follows the drawn curvature^-beta (positiveControlAdditiveVsRetimed_v002 enacted_local, matched per trial).
dur = n / FS;  nCyc = f0 * dur;  nP = round(nCyc * 4000);
phi = linspace(0, 2*pi*nCyc, nP)';
xE = a*cos(phi); yE = b*sin(phi);
nMaj = coloured_local(nP, aMaj, sMaj);  nMin = coloured_local(nP, aMin, sMin);
px = xE*cos(th) - yE*sin(th) + nMaj*cos(th) - nMin*sin(th);
py = xE*sin(th) + yE*cos(th) + nMaj*sin(th) + nMin*cos(th);
xp = gradient(px); yp = gradient(py); xpp = gradient(xp); ypp = gradient(yp);
k  = abs(xp.*ypp - yp.*xpp) ./ max((xp.^2 + yp.^2).^1.5, eps);
k  = max(k, prctile(k(k > 0), 1));
w  = k.^(-beta);
ds = max(hypot(diff(px), diff(py)), eps);
t  = [0; cumsum(ds ./ (0.5*(w(1:end-1) + w(2:end))))];
t  = t * dur / t(end);
tq = (0:n-1)'/FS;
x  = interp1(t, px, tq, "pchip") + sensor*randn(n, 1);
y  = interp1(t, py, tq, "pchip") + sensor*randn(n, 1);
end

function z = coloured_local(n, alpha, sigma)
m  = 2^nextpow2(2*n);
f  = [0, 1:m/2, -(m/2-1):-1]' / m;
amp = zeros(m, 1); nz = f ~= 0; amp(nz) = abs(f(nz)).^(-alpha/2);
Z  = (randn(m, 1) + 1i*randn(m, 1)) .* amp;
z  = real(ifft(Z)); z = z(1:n); z = z - mean(z);
s  = std(z);
if ~(s > 0), error("ctrlM:Noise", "%s", "Degenerate coloured noise."); end
z  = sigma * z / s;
end

% ---- dataset configuration, mirrors checkAdditiveVsRetimed_v003 datasetConfig_local ----
function c = datasetConfig_local(name, srcDir)
% Mirrors runLoopClosureFftnoise_v012 L118-162.
    switch name
        case "Cook CTRL",    f = "cook";         s = 0.248;  fs = 133; im = @() importDB_cook_v002(Group="CTRL", Tasks=7, Verbose=false);
        case "Cook ASD",     f = "cookASD";      s = 0.248;  fs = 133; im = @() importDB_cook_v002(Group="ASD",  Tasks=7, Verbose=false);
        case "Hickman PLAC", f = "hickmanPLAC";  s = 0.248;  fs = 133; im = @() importDB_hickman_v003(Study=2, Group="PLAC", Verbose=false);
        case "Hickman HALO", f = "hickmanHALO";  s = 0.248;  fs = 133; im = @() importDB_hickman_v003(Study=2, Group="HALO", Verbose=false);
        case "Zarandi",      f = "zarandi";      s = 10.0;   fs = 100; im = @() importDB_zarandi_v001(Verbose=false);
        case "Dhieb",        f = "dhieb";        s = 0.1478; fs = 100; im = @() importDB_dhieb_v001(Verbose=false);
        case "Fraser",       f = "fraser";       s = 1/10.41793; fs = 240; im = @() importDB_fraser_v001(Verbose=false);
        otherwise
            error("ctrlM:UnknownDataset", "%s", "Unknown dataset: " + name);
    end
    c = struct("matFile", fullfile(srcDir, "noiseCharacterisation_" + f + ".mat"), ...
        "sigToMM", s, "FS", fs, "importFn", im);
    if ~isfile(c.matFile)
        error("ctrlM:NoNoiseMat", "%s", "Missing " + c.matFile);
    end
end

% ---- Copied verbatim from checkAdditiveVsRetimed_v003 (itself copied from runLoopClosureFftnoise_v012) ----
function f0 = estimateF0_local(x, y, fs)
    N    = numel(x); nfft = 2^nextpow2(4*N);
    Xf   = abs(fft(detrend(x(:), 'linear'), nfft));
    Yf   = abs(fft(detrend(y(:), 'linear'), nfft));
    fAx  = (0:nfft-1)' * fs / nfft;
    band = fAx > 0.1 & fAx < 5; idx = find(band);
    [~, pkX] = max(Xf(band)); [~, pkY] = max(Yf(band));
    f0x = fAx(idx(pkX)); f0y = fAx(idx(pkY));
    if abs(f0x - f0y) / max(f0x, 0.01) < 0.2
        f0 = mean([f0x, f0y]);
    elseif max(Xf(band)) > max(Yf(band))
        f0 = f0x;
    else
        f0 = f0y;
    end
end

function [resX, resY, fitX, fitY] = templateSubtract_local(x, y, fs, f0, nH, nCW)
    N   = numel(x);
    win = round(nCW/f0*fs);
    hop = round(win/2);
    win = min(win, N);
    win = max(win, round(2/f0*fs));
    resX = zeros(N,1); resY = zeros(N,1); ww = zeros(N,1);
    s = 1;
    while s + win - 1 <= N
        e   = s + win - 1;
        tw  = (0:win-1)'/fs;
        han = 0.5*(1 - cos(2*pi*(0:win-1)'/(win-1)));
        D   = [ones(win,1), tw];
        for h = 1:nH
            D(:,end+1) = cos(2*pi*h*f0*tw); %#ok<AGROW>
            D(:,end+1) = sin(2*pi*h*f0*tw); %#ok<AGROW>
        end
        Dw        = D .* han;
        bX        = Dw \ (x(s:e) .* han);
        bY        = Dw \ (y(s:e) .* han);
        resX(s:e) = resX(s:e) + (x(s:e) - D*bX) .* han;
        resY(s:e) = resY(s:e) + (y(s:e) - D*bY) .* han;
        ww(s:e)   = ww(s:e) + han;
        s         = s + hop;
    end
    ok = ww > 0;
    resX(ok) = resX(ok) ./ ww(ok);
    resY(ok) = resY(ok) ./ ww(ok);
    cl   = max(round(win/4), 1);
    cs   = cl;
    ce   = min(max(N - cl, cs + 100), N);
    fitX = x(cs:ce) - resX(cs:ce);
    fitY = y(cs:ce) - resY(cs:ce);
    resX = resX(cs:ce);
    resY = resY(cs:ce);
end
