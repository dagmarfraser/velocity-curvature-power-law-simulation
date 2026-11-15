function R = positiveControlMatchedByDataset_v006(opts)
% POSITIVECONTROLMATCHEDBYDATASET_V006  Per-dataset matched controls for the generative-model check (D20), REPLICATED over seeds (2026-11-07).
% WHY: v005 (one seed per dataset) showed the single-draw result is seed-sensitive at NTrials = 60 (which datasets meet the windows, Dhieb's
% control separating or not). v006 follows the seed practice settled for the fractal-noise validation (threefry master stream, one substream
% per cell, a seed manifest, summaries over many replicates, a master-seed stability check; never a single draw):
%   - master stream RandStream("threefry4x64_20", Seed = Master); cell (dataset d, replicate r) uses Substream = d*100000 + r, set as the
%     global stream before the trial draw, the quantisation noise and the enacted noise; the inner checkAdditiveVsRetimed_v003 surrogates take
%     RngSeed = a 31-bit integer drawn from the same substream, so every random draw of a cell is fixed by (Master, d, r) and independent of
%     run order, other datasets and parallel scheduling;
%   - a manifest row per cell (master, dataset, replicate, substream, innerSeed); per-replicate rows as v005's `out`;
%   - summaries per dataset over replicates: median and 5-95% range of medA, medR, lambdaStar_own, the share of replicates meeting the pre-set
%     windows (|medA| <= 0.30, 0.70 <= medR <= 1.30) and the share with lambdaStar_own > 0.5 (outside).
% The per-trial construction is v005's (v002's) unchanged. Stability across masters is checked by summariseMatchedControlSeeds (separate).
%   R = positiveControlMatchedByDataset_v006(Master=1, Reps=1:20, Datasets=[...], UseParfor=true, Tag="")
%   Writes results/positiveControlMatchedByDataset_v006_m<Master>_r<first>-<last><Tag>_<stamp>.mat (CFG, rows, manifest, summary).   Fraser, D.S. (2026)  v006
arguments
    opts.Master (1,1) double {mustBeInteger, mustBePositive} = 1
    opts.Reps (1,:) double {mustBeInteger, mustBePositive} = 1:20
    opts.Datasets (1,:) string = ["Fraser","Cook CTRL","Cook ASD","Hickman PLAC","Hickman HALO","Zarandi","Dhieb"]
    opts.UseParfor (1,1) logical = true
    opts.Tag (1,1) string = ""
end
CFG.Datasets = ["Fraser","Cook CTRL","Cook ASD","Hickman PLAC","Hickman HALO","Zarandi","Dhieb"];   % fixed order: d indexes substreams
CFG.NTrials = 60;  CFG.NSurr = 10;  CFG.Master = opts.Master;  CFG.Reps = opts.Reps;
CFG.CAL_ADD = [-0.30 0.30];  CFG.CAL_RET = [0.70 1.30];  CFG.RULE = 0.5;
CFG.FLOOR_FILE = "checkSensorFloor_v001.mat";
CFG.CORPUS_TABLE = "checkAdditiveVsRetimed_v002_20260922_205130.mat";
if any(~ismember(opts.Datasets, CFG.Datasets)), error("ctrlM6:dataset", "%s", "Unknown dataset in Datasets"); end
srcDir = fileparts(mfilename("fullpath"));  cd(srcDir);  addpath(genpath(fullfile(srcDir, "functions")));
stamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
FL = load(fullfile(fileparts(srcDir), "results", CFG.FLOOR_FILE), "out").out;
TC = load(fullfile(srcDir, "results", CFG.CORPUS_TABLE), "T").T;
rows = table();  manifest = table();  perTrial = table();
for dName = opts.Datasets
    d = find(CFG.Datasets == dName);
    cfg = datasetConfig_local(dName, srcDir);  trials = cfg.importFn();   % import once per dataset
    for r = opts.Reps
        st = RandStream("threefry4x64_20", "Seed", CFG.Master);  st.Substream = d*100000 + r;  RandStream.setGlobalStream(st);
        innerSeed = randi(2^31 - 1);
        t0 = tic;
        cellTag = sprintf("ctrlM_v006_m%d_%s_r%d_", CFG.Master, strrep(dName, " ", ""), r);
        [row, TS] = oneCell_local(CFG, opts, dName, cfg, trials, TC, FL, srcDir, innerSeed, cellTag);
        cf = dir(fullfile(srcDir, "results", "checkAdditiveVsRetimed_v003_" + cellTag + "*.mat"));   % the inner check saves one file per call: file it per cell
        cellDir = fullfile(srcDir, "results", "ctrlM_v006_cells");  if ~isfolder(cellDir), mkdir(cellDir); end
        for k = 1:numel(cf), movefile(fullfile(cf(k).folder, cf(k).name), fullfile(cellDir, cf(k).name)); end
        row.master = CFG.Master;  row.rep = r;  row.substream = d*100000 + r;  row.innerSeed = innerSeed;  row.seconds = toc(t0);
        rows = [rows; row]; %#ok<AGROW>
        manifest = [manifest; table(CFG.Master, dName, r, d*100000 + r, innerSeed, 'VariableNames', ["master","dataset","rep","substream","innerSeed"])]; %#ok<AGROW>
        TS.srcDataset = repmat(dName, height(TS), 1);  TS.rep = repmat(r, height(TS), 1);  perTrial = [perTrial; TS]; %#ok<AGROW>
        fprintf("rep %d %-13s medA %+.3f medR %.3f calOK %d lambdaStar %+.3f (%.0f s)\n", r, dName, row.medA, row.medR, row.calOK, row.lambdaStar_own, row.seconds);
    end
end
summary = summarise_local(rows, CFG);  disp(summary);
f = fullfile(fileparts(srcDir), "results", sprintf("positiveControlMatchedByDataset_v006_m%d_r%d-%d%s_%s.mat", CFG.Master, min(opts.Reps), max(opts.Reps), opts.Tag, stamp));
save(f, "CFG", "rows", "manifest", "summary", "perTrial");  fprintf("Saved %s\n", f);
R = struct("rows", rows, "manifest", manifest, "summary", summary, "file", f);
end

function S = summarise_local(rows, CFG)
S = table();
for dName = unique(rows.dataset, "stable")'
    x = rows(rows.dataset == dName, :);  q = @(v) prctile(v, [5 50 95]);
    S = [S; table(dName, height(x), q(x.medA), q(x.medR), q(x.lambdaStar_own), mean(x.calOK), mean(x.calOK_quant), mean(x.outside), ...
        'VariableNames', ["dataset","nReps","medA_p5_50_95","medR_p5_50_95","lambdaStar_p5_50_95","shareCalOK","shareCalOK_quant","shareOutside"])]; %#ok<AGROW>
end
end

function [row, TS] = oneCell_local(CFG, opts, dName, cfg, trials, TC, FL, srcDir, innerSeed, cellTag) %#ok<INUSD>
FS = cfg.FS;  tids = string({trials.trialID}');
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
    TS = checkAdditiveVsRetimed_v003(Trials=S, TrialsFS=FS, NSurr=CFG.NSurr, UseParfor=opts.UseParfor, RngSeed=innerSeed, SaveTag=cellTag);
    g = @(c) median(TS.lambda(TS.dataset == c & isfinite(TS.lambda)), "omitnan");
    n = @(c) sum(TS.dataset == c & isfinite(TS.lambda));
    medA = g("A_sxu_clean");  medR = g("R_enacted_clean");  medAq = g("A_sxu_quant");  medRq = g("R_enacted_quant");  medRh = g("R_enacted_hf");
    calOK = abs(medA) <= CFG.CAL_ADD(2) && medR >= CFG.CAL_RET(1) && medR <= CFG.CAL_RET(2);
    calOKq = abs(medAq) <= CFG.CAL_ADD(2) && medRq >= CFG.CAL_RET(1) && medRq <= CFG.CAL_RET(2);
    mm = TC.dataset == dName & isfinite(TC.lambda);
    sL = splitapply(@median, TC.lambda(mm), findgroups(TC.subjectID(mm)));
    lamSubj = median(sL);  lamStar = (lamSubj - medA) / (medR - medA);
    row = table(dName, numel(cand), n("A_sxu_clean"), n("R_enacted_clean"), medA, medR, calOK, medAq, medRq, calOKq, medRh, lamSubj, lamStar, lamStar > CFG.RULE, ...
        'VariableNames', ["dataset","nSource","nFiniteA","nFiniteR","medA","medR","calOK","medA_quant","medR_quant","calOK_quant","medR_hf","subjMed_lambda","lambdaStar_own","outside"]);
    TS = TS(:, ["dataset","trialID","subjectID","lambda","lambdaHigh","status"]);
end

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

