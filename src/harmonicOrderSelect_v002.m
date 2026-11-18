function out = harmonicOrderSelect_v002(opts)
%HARMONICORDERSELECT_V002 Information-criterion harmonic order, global-minimum row and first-minimum row (Amendments 01 and 02), acceptance test and real-data selection (2026-11-17).
%   v002 = v001 + the Amendment 02 first-minimum row + Mode "real" on all seven datasets. v001 is not edited. Position data only: no beta_obs, beta_gen* or forward map.
%   Per trial and axis, for K = 1 .. Kmax (Kmax = min(20, floor(10 Hz / f0)), at least 1): production end-anchored template with K harmonics (constant, trend and K harmonic pairs per Hann-weighted
%   window of 4 cycles, overlap-add; as runLoopClosureFftnoise_v028 templateSubtract_local) gives the residual r_K of length N; sigma2_K(h) is the one-step prediction variance of an AR(h) fitted to
%   r_K by least squares, H = min(20, round((ln N)^2));  C(K, h) = N ln sigma2_K(h) + (2K + h) ln N  (sigma2_K(0) = var(r_K)).
%   GLOBAL row (Amendment 01, unchanged from v001): K = argmin over K and h = 1 .. H of C(K, h); the trial's K is the larger of the two axes' values.
%   FIRST-MINIMUM row (Amendment 02, Kavalieris & Hannan 1994 Theorem 2): Cmin(K) = min over h = 0 .. H of C(K, h) (the frozen text of Amendment 02 says h = 0 to H; the global row follows Amendment 01 as
%       written, h = 1 to H; the two ranges differ, so Kfirst <= Kglobal is not strictly guaranteed); the axis's K is the first K with Cmin(K) <= Cmin(K + 1), or Kmax if the profile falls throughout;
%       the trial's Kfirst is the larger of the two axes' values.
%   BIC row (sensitivity, unchanged): argmin over K of N ln var(r_K) + 2 K ln N.
%   Dataset K (Amendment 01 section 3 item 3) = median of the per-trial K, rounded, per dataset.
%   Modes: out = harmonicOrderSelect_v002(Mode="acceptance");   % semi-synthetic trials as harmonicRecoveryTest (real Fraser residual as noise; sine with 0 or 10% 2nd and 3rd harmonics); descriptive for the first-minimum row
%          out = harmonicOrderSelect_v002(Mode="real");         % the three IC rows per dataset on real positions (trials of the v028 primary file: trial IDs only are read from it)
%   Regression: with Mode="acceptance" the global K must equal harmonicOrderSelect_v001's Kic on results/harmonicOrderSelect_v001_acceptance.mat trial for trial (checked inside, error if not, when that file exists).
%   Writes results/harmonicOrderSelect_v002_<Mode><Tag>.mat
arguments
    opts.Mode (1,1) string {mustBeMember(opts.Mode, ["acceptance", "real"])} = "acceptance"
    opts.MaxTrials (1,1) double = 300
    opts.Distortions (1,:) double = [0 0.1]
    opts.Datasets (1,:) string = ["FraserN48" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Zarandi" "Dhieb"]
    opts.Tag (1,1) string = ""
    opts.UseParfor (1,1) logical = true
    opts.Seed (1,1) double = 1729
end
here = fileparts(mfilename('fullpath'));  addpath(genpath(fullfile(here, 'functions')));
resDir = fullfile(fileparts(here), 'results');  if ~isfolder(resDir), mkdir(resDir); end
ws = warning('off', 'all');  wCleanup = onCleanup(@() warning(ws)); %#ok<NASGU>
if opts.UseParfor, try, ensurePool_v001(); catch, end, end %#ok<TRYNC>
if opts.Mode == "real"
    S = table();  det = struct();
    for dsn = opts.Datasets
        switch dsn
            case "FraserN48",      tr = importDB_fraser_v002(Verbose=false);
            case "Cook_CTRL",      tr = importDB_cook_v002(Group="CTRL", Tasks=7, Verbose=false);
            case "Cook_ASD",       tr = importDB_cook_v002(Group="ASD",  Tasks=7, Verbose=false);
            case "Hickman_PLAC",   tr = importDB_hickman_v003(Study=2, Group="PLAC", Verbose=false);
            case "Hickman_HALO",   tr = importDB_hickman_v003(Study=2, Group="HALO", Verbose=false);
            case "Zarandi",        tr = importDB_zarandi_v001(Verbose=false);
            case "Dhieb",          tr = importDB_dhieb_v001(Verbose=false);
        end
        R = loadV028_v001(strrep(dsn, "FraserN48", "Fraser"));  keep = string({R.trialID});   % trial IDs only
        tr = tr(ismember(string({tr.trialID}), keep));  n = numel(tr);
        X = {tr.x};  Y = {tr.y};  FS = [tr.fs];  Kg = NaN(n, 1);  Kb = NaN(n, 1);  Kf = NaN(n, 1);  f0 = NaN(n, 1);
        if opts.UseParfor
            parfor i = 1:n, [Kg(i), Kb(i), Kf(i), f0(i)] = selectOne_local(double(X{i}), double(Y{i}), FS(i)); end
        else
            for i = 1:n, [Kg(i), Kb(i), Kf(i), f0(i)] = selectOne_local(double(X{i}), double(Y{i}), FS(i)); end
        end
        ok = isfinite(Kg);
        fprintf('%-13s n %5d usable %5d  median K: global %.1f, first minimum %.1f, plain BIC %.1f  -> dataset K (rounded median): global %d, first minimum %d\n', dsn, n, nnz(ok), ...
            median(Kg, 'omitnan'), median(Kf, 'omitnan'), median(Kb, 'omitnan'), round(median(Kg, 'omitnan')), round(median(Kf, 'omitnan')));
        S = [S; table(dsn, n, nnz(ok), median(f0, 'omitnan'), median(Kg, 'omitnan'), round(median(Kg, 'omitnan')), median(Kf, 'omitnan'), round(median(Kf, 'omitnan')), median(Kb, 'omitnan'), round(median(Kb, 'omitnan')), ...
            'VariableNames', {'dataset', 'nTrials', 'nUsable', 'f0Median', 'KglobalMedian', 'KglobalDataset', 'KfirstMedian', 'KfirstDataset', 'KbicMedian', 'KbicDataset'})]; %#ok<AGROW>
        det.(matlab.lang.makeValidName(dsn)) = struct('Kglobal', Kg, 'Kfirst', Kf, 'Kbic', Kb, 'f0', f0, 'trialID', string({tr.trialID})');
    end
    disp(S);
    out = struct('summary', S, 'detail', det, 'opts', opts, 'note', "global row h = 1..H (Amendment 01); first-minimum row h = 0..H (Amendment 02)", 'generated', string(datetime('now')));
else
    trials = importDB_fraser_v002(Verbose=false);  n0 = numel(trials);
    sel = unique(round(linspace(1, n0, min(opts.MaxTrials, n0))));  trials = trials(sel);  n = numel(trials);
    X = {trials.x};  Y = {trials.y};  FS = [trials.fs];
    nD = numel(opts.Distortions);  truth = 1 + 2 * (opts.Distortions > 0);
    Kg = NaN(n, nD);  Kb = NaN(n, nD);  Kf = NaN(n, nD);
    rs = RandStream('threefry4x64_20', 'Seed', opts.Seed);  PH = 2*pi*rand(rs, n, nD, 2);   % same phases as v001
    for d = 1:nD
        tg = NaN(n, 1);  tb = NaN(n, 1);  tf = NaN(n, 1);  dist = opts.Distortions(d);  phd = squeeze(PH(:, d, :));
        if opts.UseParfor
            parfor i = 1:n, [tg(i), tb(i), tf(i)] = acceptOne_local(double(X{i}), double(Y{i}), FS(i), dist, phd(i, :)); end
        else
            for i = 1:n, [tg(i), tb(i), tf(i)] = acceptOne_local(double(X{i}), double(Y{i}), FS(i), dist, phd(i, :)); end
        end
        Kg(:, d) = tg;  Kb(:, d) = tb;  Kf(:, d) = tf;
    end
    % regression against v001
    f1 = fullfile(resDir, 'harmonicOrderSelect_v001_acceptance.mat');
    if isfile(f1)
        A = load(f1, 'Kic');  same = isequaln(A.Kic, Kg);
        fprintf('REGRESSION vs harmonicOrderSelect_v001 acceptance Kic (global minimum): identical = %d\n', same);
        if ~same, error('hos2:regression', 'v002 global-minimum K differs from v001 on the acceptance trials'); end
    end
    rows = table();  usable = isfinite(Kg(:, 1));
    for d = 1:nD
        for m = ["IC global", "IC first minimum", "plain BIC"]
            k = Kg(:, d);  if m == "IC first minimum", k = Kf(:, d); elseif m == "plain BIC", k = Kb(:, d); end
            ok = isfinite(k) & usable;   % same 166-trial set for every row
            rows = [rows; table(opts.Distortions(d), truth(d), m, median(k(ok)), mean(k(ok) == truth(d)), mean(k(ok) > truth(d)), mean(k(ok) < truth(d)), nnz(ok), ...
                'VariableNames', {'distortion', 'trueK', 'criterion', 'medianK', 'shareExact', 'shareOver', 'shareUnder', 'n'})]; %#ok<AGROW>
        end
    end
    disp(rows);
    out = struct('summary', rows, 'Kglobal', Kg, 'Kfirst', Kf, 'Kbic', Kb, 'truth', truth, 'opts', opts, 'note', "global row h = 1..H (Amendment 01); first-minimum row h = 0..H (Amendment 02); descriptive only", 'generated', string(datetime('now')));
end
f = fullfile(resDir, "harmonicOrderSelect_v002_" + opts.Mode + opts.Tag + ".mat");  save(f, '-struct', 'out', '-v7.3');  fprintf('Saved: %s\n', f);
end

function [Kg, Kb, Kf] = acceptOne_local(x, y, fs, dist, ph)
% semi-synthetic trial: real 4-harmonic residual as the noise, known movement added (as harmonicRecoveryTest, no drift)
Kg = NaN;  Kb = NaN;  Kf = NaN;  f0 = estimateF0_local(x, y, fs);
if ~isfinite(f0) || f0 <= 0, return, end
[nx, ny] = templateSubtract_local(x, y, fs, f0, 4, 4);  N = numel(nx);  t = (0:N-1)' / fs;
D = [ones(numel(x), 1), (0:numel(x)-1)'/fs, cos(2*pi*f0*(0:numel(x)-1)'/fs), sin(2*pi*f0*(0:numel(x)-1)'/fs)];
bx = D \ detrend(x(:), 'linear');  by = D \ detrend(y(:), 'linear');  Ax = hypot(bx(3), bx(4));  Ay = hypot(by(3), by(4));
phs = 2*pi*f0*t;  mx = Ax * cos(phs);  my = Ay * sin(phs);
for h = 2:3, mx = mx + dist * Ax * cos(h*phs + ph(1)) / h;  my = my + dist * Ay * sin(h*phs + ph(2)) / h;  end
[Kg, Kb, Kf] = selectCore_local(mx + nx(:), my + ny(:), fs, f0);
end

function [Kg, Kb, Kf, f0] = selectOne_local(x, y, fs)
Kg = NaN;  Kb = NaN;  Kf = NaN;  f0 = estimateF0_local(x, y, fs);
if ~isfinite(f0) || f0 <= 0, return, end
[Kg, Kb, Kf] = selectCore_local(x, y, fs, f0);
end

function [Kg, Kb, Kf] = selectCore_local(x, y, fs, f0)
Kmax = max(1, min(20, floor(10 / f0)));  kg = NaN(1, 2);  kb = NaN(1, 2);  kf = NaN(1, 2);
for a = 1:2
    s = ifelse_local(a == 1, x, y);  Cg = Inf(1, Kmax);  Cf = Inf(1, Kmax);  CB = Inf(1, Kmax);
    for K = 1:Kmax
        [rx, ~] = templateSubtract_local(s, s, fs, f0, K, 4);   % same fit on the chosen axis (the second output is a duplicate)
        r = rx(:) - mean(rx);  N = numel(r);  lnN = log(N);  H = min(20, round(lnN^2));
        v0 = var(r);
        if v0 > 0, CB(K) = N * log(v0) + 2 * K * lnN;  Cf(K) = N * log(v0) + (2*K + 0) * lnN;  end   % h = 0 enters the first-minimum profile only
        for h = 1:H
            s2 = arPredVar_local(r, h);  if ~isfinite(s2) || s2 <= 0, continue, end
            C = N * log(s2) + (2*K + h) * lnN;
            if C < Cg(K), Cg(K) = C; end
            if C < Cf(K), Cf(K) = C; end
        end
    end
    if all(~isfinite(Cg)), continue, end   % no fitted template (series shorter than the window): NaN, excluded
    [~, kg(a)] = min(Cg);  [~, kb(a)] = min(CB);
    k1 = Kmax;  for K = 1:Kmax-1, if Cf(K) <= Cf(K + 1), k1 = K;  break, end, end
    kf(a) = k1;
end
Kg = max(kg);  Kb = max(kb);  Kf = max(kf);
end

function s2 = arPredVar_local(r, h)
% least-squares AR(h) one-step prediction variance
N = numel(r);  if N <= 3*h, s2 = NaN; return, end
Z = zeros(N - h, h);  for j = 1:h, Z(:, j) = r(h-j+1:N-j); end
y = r(h+1:N);  b = Z \ y;  e = y - Z*b;  s2 = mean(e .^ 2);
end

function v = ifelse_local(c, a, b)
if c, v = a; else, v = b; end
end

function f0 = estimateF0_local(x, y, fs)   % as runLoopClosureFftnoise_v028
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

function [resX, resY] = templateSubtract_local(x, y, fs, f0, nH, nCW)   % as runLoopClosureFftnoise_v028, END-ANCHORED
    N   = numel(x);
    win = round(nCW/f0*fs);
    hop = round(win/2);
    win = min(win, N);
    win = max(win, round(2/f0*fs));
    resX = zeros(N,1); resY = zeros(N,1); ww = zeros(N,1);
    starts = 1:hop:(N - win + 1);
    if win > N, starts = []; end
    if ~isempty(starts) && starts(end) + win - 1 < N, starts(end+1) = N - win + 1; end
    for s = starts
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
    end
    ok = ww > 0;
    resX(ok) = resX(ok) ./ ww(ok);
    resY(ok) = resY(ok) ./ ww(ok);
    cl   = max(round(win/4), 1);
    cs   = cl;
    ce   = min(max(N - cl, cs + 100), N);
    resX = resX(cs:ce);
    resY = resY(cs:ce);
end
