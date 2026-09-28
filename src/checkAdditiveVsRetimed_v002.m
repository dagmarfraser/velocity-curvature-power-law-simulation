function T = checkAdditiveVsRetimed_v002(opts)
% checkAdditiveVsRetimed_v002  Does speed follow the EXECUTED path's curvature
% or the TEMPLATE's? Tests the additive-entry assumption of the forward map.
%
% The loop-closure forward map (runLoopClosureFftnoise_v012 L806-846) builds
% each synthetic trial as  template(beta_gen) + additive shaped_xu noise
% (Maoz et al., 2005, generalised). R7 checks the surrogate's alpha/sigma, a
% per-axis property of the POSITION residual, which is identical whether
% motor deviations are (A) additive on top of template timing or (R) enacted
% as path geometry with power-law timing along the path actually drawn.
% This script measures the property that separates them.
%
% STATISTIC (per trace, SG differentiation, EdgeClip samples trimmed):
%   log v = c + bT*log(kTpl) + bD*delta,   delta = log(kObs) - log(kTpl)
%   kTpl from the OLA harmonic fit (templateSubtract_local), kObs from the
%   trace itself.  rho = bD/bT.  (A) predicts rho ~ 0; (R) predicts rho ~ 1.
%   Both are attenuated by differentiation noise, so rho is NOT read raw:
%   each real trial is placed between two surrogate anchors built from its
%   own geometry, f0 and residual spectrum, analysed identically:
%     rhoAdd - template + shaped_xu noise (the forward map's own model)
%     rhoRet - the same curve re-timed so v = K*kappa_path^-beta, duration kept
%   lambda = (rhoReal - med(rhoAdd)) / (med(rhoRet) - med(rhoAdd))
%   lambda ~ 0: data behave as the forward map assumes.  ~ 1: re-timed world.
%
% LIMITATIONS (stated, not hidden):
%   - rhoRet re-times ALL surrogate noise into geometry, including any
%     sensor-like component; it is an idealised upper anchor.
%   - The OLA fit absorbs within-window path deviation, so delta sees only
%     deviation beyond 4 harmonics; identical for real and both anchors.
%   - The two-predictor fit uses plain OLS (mldivide): a diagnostic ratio, not
%     a beta estimate, so regressDataEBR (single predictor) does not apply.
%   - Anchors that fail to separate give lambda = NaN, status reported.
%   - Surrogates are built at the real residual's trimmed length M (~N/2
%     for slow trials). If the OLA window exceeds M (fewer than ~2 cycles
%     left), coupling cannot be computed: status surrogate_olaWindowExceedsTrace.
%     This excludes SLOW drawers specifically; report the count per dataset.
%
% v002 ADDITIONS (v001 statistic retained unchanged as rho/lambda):
%   - f0 of the real trace recorded per trial (tempo; within-dataset lambda~f0).
%   - Band-split delta: delta = dLow + dHigh, dLow = zero-phase 4th-order
%     Butterworth low-pass of delta at fc = CycleBandMult*f0 (default 3),
%     dHigh = delta - dLow (exactly complementary). Fit
%       log v = c + bT*log(kTpl) + bL*dLow + bH*dHigh
%     rhoLow = bL/bT, rhoHigh = bH/bT, each placed between its own anchors
%     (lambdaLow, lambdaHigh; own separation gates, statusLow/statusHigh).
%     Enacted cycle-scale geometry predicts lambdaLow >> lambdaHigh; upstream
%     processing predicts coupling carried in the high band too.
%     The filter acts on delta (a diagnostic decomposition), not on
%     positions, so the Fraser et al. (2025) filtering caveat on beta does
%     not apply; delta gaps are linearly filled before filtering, fit on
%     valid samples only.
%
% USAGE:
%   T = checkAdditiveVsRetimed_v002();
%   T = checkAdditiveVsRetimed_v002(Datasets="Zarandi", NSurr=10);
%
% Fraser, D.S. (2026) v002 (from v001)

arguments
    opts.Datasets (1,:) string = ["Fraser","Cook CTRL","Cook ASD", ...
        "Hickman PLAC","Hickman HALO","Zarandi","Dhieb"]
    opts.TrialSelection (1,1) string {mustBeMember(opts.TrialSelection, ...
        ["per_subject","all"])} = "per_subject"
    opts.NSurr    (1,1) double {mustBeInteger, mustBePositive} = 5
    opts.EdgeClip (1,1) double {mustBeInteger, mustBePositive} = 20
    opts.RngSeed  (1,1) double = 1729
    opts.SepK     (1,1) double {mustBePositive} = 2   % anchor separation, in pooled SDs
    opts.CycleBandMult (1,1) double {mustBePositive} = 3   % fc = CycleBandMult*f0
    opts.UseParfor (1,1) logical = false
end

srcDir = fileparts(mfilename("fullpath"));
addpath(srcDir); addpath(genpath(fullfile(srcDir, "functions")));
addpath(genpath(fullfile(srcDir, "req")));
SG = struct("filterType", 6, "filterParams", [4 17]);   % runner's SG config

if opts.UseParfor
    nW = max(1, feature("numcores") - 1);                % leave one core free
    p  = gcp("nocreate");
    if isempty(p) || p.NumWorkers > nW
        delete(p); parpool("Processes", nW);
    end
    fprintf("Pool: %d workers (%d cores)\n", gcp("nocreate").NumWorkers, feature("numcores"));
end

rows = {};
for dName = opts.Datasets
    cfg = datasetConfig_local(dName, srcDir);
    bio = load(cfg.matFile, "bioResults").bioResults;
    selIdx = selectTrials_local(bio, cfg.sigToMM, opts.TrialSelection);
    trials = cfg.importFn();
    tids   = string({trials.trialID}');
    [ok, loc] = ismember(string(bio.trialID(selIdx)), tids);
    if any(~ok)
        warning("checkAdditiveVsRetimed_v002:NoMatch", "%s", sprintf( ...
            "%s: %d selected trials not found in import, excluded.", dName, sum(~ok)));
    end
    trialsSel = trials(loc(ok)); selIdx = selIdx(ok);
    nT = numel(trialsSel);
    fprintf("=== %s: %d trials (%s), NSurr=%d ===\n", dName, nT, opts.TrialSelection, opts.NSurr);

    out = cell(nT, 1);
    FS = cfg.FS; NS = opts.NSurr; EC = opts.EdgeClip; seed0 = opts.RngSeed; sepK = opts.SepK; cbm = opts.CycleBandMult;
    if opts.UseParfor
        parfor k = 1:nT
            out{k} = oneTrial_local(trialsSel(k), FS, SG, NS, EC, seed0 + k, sepK, cbm);
        end
    else
        for k = 1:nT
            out{k} = oneTrial_local(trialsSel(k), FS, SG, NS, EC, seed0 + k, sepK, cbm);
        end
    end
    R = struct2table([out{:}]');
    R.dataset   = repmat(dName, nT, 1);
    R.subjectID = string(bio.subjectID(selIdx));
    rows{end+1} = R; %#ok<AGROW>
end
T = vertcat(rows{:});

%% ---- Summary (subject is the unit when several trials per subject) -------
fprintf("\n%-13s %5s %5s %8s %8s %8s %8s %16s\n", "dataset", "nOK", "nNaN", ...
    "rhoReal", "rhoAdd", "rhoRet", "lambda", "lambda IQR");
for dName = opts.Datasets
    m = T.dataset == dName; okm = m & isfinite(T.lambda);
    if ~any(okm)
        fprintf("%-13s %5d %5d   NO FINITE LAMBDA -- see NaN reasons\n", dName, 0, sum(m));
        subjLam = NaN;
    else
        subjLam = splitapply(@median, T.lambda(okm), findgroups(T.subjectID(okm)));
    end
    fprintf("%-13s %5d %5d %8.3f %8.3f %8.3f %8.3f   [%6.3f %6.3f]\n", dName, ...
        sum(okm), sum(m & ~okm), median(T.rhoReal(okm)), median(T.rhoAddMed(okm)), ...
        median(T.rhoRetMed(okm)), median(subjLam), prctile(subjLam, 25), prctile(subjLam, 75));
    st = T.status(m & ~okm);
    if ~isempty(st)
        [u, ~, g] = unique(st);
        fprintf("   NaN reasons: %s\n", strjoin(u + "=" + string(accumarray(g, 1)), ", "));
    end
end

resDir = fullfile(srcDir, "results");
if ~isfolder(resDir), mkdir(resDir); end
outFile = fullfile(resDir, "checkAdditiveVsRetimed_v002_" + ...
    string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".mat");
save(outFile, "T", "opts");
fprintf("\nSaved: %s\n", outFile);
end

% =========================================================================
function r = oneTrial_local(tr, FS, SG, NS, EC, seed, sepK, cbm)
    rng(seed, "twister");   % per-trial stream: parfor-deterministic
    r = struct("trialID", string(tr.trialID), "f0", NaN, "betaTpl", NaN, "bTReal", NaN, ...
        "rhoReal", NaN, "rhoAddMed", NaN, "rhoRetMed", NaN, "rhoAddSD", NaN, "rhoRetSD", NaN, "lambda", NaN, ...
        "rhoRealLow", NaN, "rhoAddLowMed", NaN, "rhoRetLowMed", NaN, "lambdaLow", NaN, ...
        "rhoRealHigh", NaN, "rhoAddHighMed", NaN, "rhoRetHighMed", NaN, "lambdaHigh", NaN, ...
        "status", "ok", "statusLow", "ok", "statusHigh", "ok");
    x = double(tr.x(:)); y = double(tr.y(:));

    cR = coupling_local(x, y, FS, SG, EC, cbm);
    r.f0 = cR.f0;
    if ~isfinite(cR.rho(1)), r.status = "real_" + cR.why; return, end
    r.rhoReal = cR.rho(1); r.rhoRealLow = cR.rho(2); r.rhoRealHigh = cR.rho(3);
    r.bTReal = cR.bT;
    betaT = min(max(-cR.bT, 0.05), 0.70);  r.betaTpl = betaT;

    % Geometry and per-axis residual, as runner L735-750 (native units)
    xc = x - mean(x); yc = y - mean(y);
    [V, D] = eig(cov(xc, yc));
    [lams, ord] = sort(diag(D), "descend"); V = V(:, ord);
    a = sqrt(2*max(lams(1), 0)); b = sqrt(2*max(lams(2), 0));
    th = atan2(V(2,1), V(1,1));
    rMaj =  cR.xRes*cos(th) + cR.yRes*sin(th);
    rMin = -cR.xRes*sin(th) + cR.yRes*cos(th);
    fHi  = min(20, FS/2 - 1);
    aMaj = iraAlphaSigma_v001(rMaj, FS, 1.0, fHi, 1.1:0.05:1.9);
    aMin = iraAlphaSigma_v001(rMin, FS, 1.0, fHi, 1.1:0.05:1.9);
    if ~isfinite(aMaj) || ~isfinite(aMin), r.status = "iraFail"; return, end

    [xT, yT] = powerLawEllipse_local(a, b, th, cR.f0, FS, betaT, numel(rMaj));
    rhoA = NaN(NS, 3); rhoRt = NaN(NS, 3); whyS = strings(0, 1);
    for s = 1:NS
        [nMaj, nMin] = generateLoopClosureNoise_v003("shaped_xu", rMaj, rMin, FS, aMaj, aMin);
        if any(~isfinite([nMaj; nMin])), continue, end
        xA = xT + nMaj*cos(th) - nMin*sin(th);
        yA = yT + nMaj*sin(th) + nMin*cos(th);
        cA = coupling_local(xA, yA, FS, SG, EC, cbm);
        [xR, yR] = retime_local(xA, yA, FS, SG, betaT);
        cRt = coupling_local(xR, yR, FS, SG, EC, cbm);
        rhoA(s, :) = cA.rho; rhoRt(s, :) = cRt.rho;
        whyS = [whyS; cA.why; cRt.why]; %#ok<AGROW>
    end
    if sum(isfinite(rhoA(:,1))) < 2 || sum(isfinite(rhoRt(:,1))) < 2
        whyS = whyS(whyS ~= "");
        if isempty(whyS), r.status = "surrogate_noiseFail";
        else, r.status = "surrogate_" + string(mode(categorical(whyS))); end
        return
    end
    r.rhoAddSD = std(rhoA(:,1), "omitnan"); r.rhoRetSD = std(rhoRt(:,1), "omitnan");
    [r.lambda, r.rhoAddMed, r.rhoRetMed, r.status] = place_local(r.rhoReal, rhoA(:,1), rhoRt(:,1), sepK);
    [r.lambdaLow, r.rhoAddLowMed, r.rhoRetLowMed, r.statusLow] = place_local(r.rhoRealLow, rhoA(:,2), rhoRt(:,2), sepK);
    [r.lambdaHigh, r.rhoAddHighMed, r.rhoRetHighMed, r.statusHigh] = place_local(r.rhoRealHigh, rhoA(:,3), rhoRt(:,3), sepK);
end

function [lam, mA, mR, status] = place_local(rhoReal, rA, rR, sepK)
% Place a real rho between its additive and re-timed anchors; NaN + reason if not separable.
    lam = NaN; status = "ok";
    mA = median(rA, "omitnan"); mR = median(rR, "omitnan");
    if ~isfinite(rhoReal), status = "realBandFail"; return, end
    if sum(isfinite(rA)) < 2 || sum(isfinite(rR)) < 2, status = "anchorBandFail"; return, end
    sep = mR - mA;
    if ~(sep > 0) || sep < sepK * sqrt((std(rA, "omitnan")^2 + std(rR, "omitnan")^2)/2)
        status = "anchorsUnseparated"; return
    end
    lam = (rhoReal - mA) / sep;
end

% -------------------------------------------------------------------------
function c = coupling_local(x, y, FS, SG, EC, cbm)
% rho(1) = bD/bT from log v ~ log kTpl + delta   (v001 statistic, unchanged)
% rho(2:3) = bL/bT, bH/bT from log v ~ log kTpl + dLow + dHigh (band split)
% NaN + reason on failure.
    c = struct("rho", NaN(1,3), "bT", NaN, "f0", NaN, "xRes", [], "yRes", [], "why", "");
    c.f0 = estimateF0_local(x, y, FS);
    [c.xRes, c.yRes, xFit, yFit] = templateSubtract_local(x, y, FS, c.f0, 4, 4);
    % OLA window can exceed the trace for slow, short traces: residual all zeros.
    if all(c.xRes == 0) && all(c.yRes == 0), c.why = "olaWindowExceedsTrace"; return, end
    [vO, kO] = kin_local(xFit + c.xRes, yFit + c.yRes, FS, SG, EC);
    [~,  kT] = kin_local(xFit, yFit, FS, SG, EC);
    m = isfinite(vO) & isfinite(kO) & isfinite(kT) & vO > 0 & kO > 0 & kT > 0;
    if sum(m) < 50, c.why = "tooFewValidSamples"; return, end
    lkT = log(kT(m)); d = log(kO(m)) - lkT; lv = log(vO(m));
    X = [ones(sum(m),1), lkT, d];
    if rank(X) < 3, c.why = "rankDeficient"; return, end
    b = X \ lv;
    if abs(b(2)) <= eps, c.why = "zeroTemplateSlope"; return, end
    c.bT = b(2); c.rho(1) = b(3) / b(2);

    % Band split on the contiguous delta series (gaps filled, fit on valid only)
    fc = cbm * c.f0;
    if ~(fc > 0) || fc >= 0.9 * FS/2, return, end        % rho(2:3) stay NaN
    dAll = NaN(numel(m), 1); dAll(m) = d;
    dAll = fillmissing(dAll, "linear", "EndValues", "nearest");
    [bb, aa] = butter(4, fc / (FS/2));
    if numel(dAll) <= 3 * 8, return, end                  % filtfilt length guard
    dLow = filtfilt(bb, aa, dAll); dHigh = dAll - dLow;
    X4 = [ones(sum(m),1), lkT, dLow(m), dHigh(m)];
    if rank(X4) < 4, return, end
    b4 = X4 \ lv;
    if abs(b4(2)) > eps, c.rho(2:3) = b4(3:4)' / b4(2); end
end

function [sp, kp] = kin_local(x, y, FS, SG, EC)
    [dx, dy] = differentiateKinematicsEBR(x, y, SG.filterType, SG.filterParams, FS);
    vx = dx(EC:end-EC, 2); vy = dy(EC:end-EC, 2);
    ax = dx(EC:end-EC, 3); ay = dy(EC:end-EC, 3);
    sp = hypot(vx, vy);
    kp = curvatureKinematicEBR(vx, vy, ax, ay);
end

% -------------------------------------------------------------------------
function [xR, yR] = retime_local(xC, yC, FS, SG, beta)
% Traverse the curve (xC,yC) with v = K*kappa_path^-beta, keeping duration.
% Curvature is parametrisation-invariant, so it is taken from the curve as
% sampled; a 1st-percentile floor guards near-inflection blow-up.
    [dx, dy] = differentiateKinematicsEBR(xC, yC, SG.filterType, SG.filterParams, FS);
    k = curvatureKinematicEBR(dx(:,2), dy(:,2), dx(:,3), dy(:,3));
    k = fillmissing(k, "nearest");
    kPos = k(k > 0 & isfinite(k));
    if isempty(kPos)
        error("checkAdditiveVsRetimed_v002:NoCurvature", "%s", "Curve has no positive curvature.");
    end
    k  = max(k, prctile(kPos, 1));
    w  = k.^(-beta);                                   % relative speed
    ds = max(hypot(diff(xC), diff(yC)), eps);
    tC = [0; cumsum(ds ./ (0.5*(w(1:end-1) + w(2:end))))];
    N  = numel(xC);
    tC = tC * ((N-1)/FS) / tC(end);                    % preserve duration
    tq = (0:N-1)' / FS;
    xR = interp1(tC, xC, tq, "pchip");
    yR = interp1(tC, yC, tq, "pchip");
end

function [x, y] = powerLawEllipse_local(a, b, th, f0, FS, beta, M)
% Arc-length power-law ellipse, as runner L806-821.
    phi  = linspace(0, 2*pi, 10000)';
    dsdp = sqrt((a*sin(phi)).^2 + (b*cos(phi)).^2);
    kp   = (a*b) ./ max((a^2*sin(phi).^2 + b^2*cos(phi).^2).^1.5, eps);
    cumT = cumsum(dsdp .* kp.^beta); cumT = cumT / cumT(end);
    tN   = min(mod((0:M-1)'/FS*f0, 1), 1 - 1e-9);
    pA   = interp1(cumT, phi, tN, "linear", "extrap");
    xE = a*cos(pA); yE = b*sin(pA);
    x = xE*cos(th) - yE*sin(th);  y = xE*sin(th) + yE*cos(th);
end

% -------------------------------------------------------------------------
function sel = selectTrials_local(bio, sigToMM, mode)
% Runner v012 L203-245 per_subject logic: trial nearest subject (alpha,sigma) median.
    if mode == "all", sel = (1:height(bio))'; return, end
    aV = bio.ira_alphaMean; sV = bio.sigmaMean * sigToMM;
    subj = unique(bio.subjectID); sel = zeros(numel(subj), 1);
    for s = 1:numel(subj)
        rowsS = find(bio.subjectID == subj(s));
        a = aV(rowsS); g = sV(rowsS); v = isfinite(a) & isfinite(g);
        if ~any(v), sel(s) = rowsS(1); continue, end
        d = (a - median(a(v))).^2 / max(var(a(v)), eps) + ...
            (g - median(g(v))).^2 / max(var(g(v)), eps);
        d(~v) = Inf; [~, j] = min(d); sel(s) = rowsS(j);
    end
end

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
            error("checkAdditiveVsRetimed_v002:UnknownDataset", "%s", "Unknown dataset: " + name);
    end
    c = struct("matFile", fullfile(srcDir, "noiseCharacterisation_" + f + ".mat"), ...
        "sigToMM", s, "FS", fs, "importFn", im);
    if ~isfile(c.matFile)
        error("checkAdditiveVsRetimed_v002:NoNoiseMat", "%s", "Missing " + c.matFile);
    end
end

% ---- Copied verbatim from runLoopClosureFftnoise_v012 (L903-963) ----------
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
