function out = harmonicOrderSelect_v001(opts)
%HARMONICORDERSELECT_V001 Information-criterion harmonic order for the template, and its acceptance test (2026-11-17).
%   Implements section 3 of docs/V028_PreCommitment_AMENDMENT01_20261117_DRAFT_v002.md (DRAFT, not frozen when written). Position data only: no beta_obs, beta_gen* or forward map.
%   Per trial and per axis, for K = 1 .. Kmax (Kmax = min(20, floor(10 Hz / f0)), at least 1): the production end-anchored template with K harmonics is subtracted (constant, trend and K harmonic
%   sine/cosine pairs per Hann-weighted window of 4 cycles, overlap-add; the same fit as runLoopClosureFftnoise_v028 templateSubtract_local), giving the order-K residual r_K of length N.
%   An AR(h) model is fitted to r_K by least squares for h = 1 .. H, H = min(20, round((ln N)^2)); sigma2_K(h) is its one-step prediction variance. The criterion is
%       C(K, h) = N ln sigma2_K(h) + (2 K + h) ln N          (Christensen et al. 2011, penalty 2 per harmonic because f0 is fixed; Kavalieris & Hannan form for the AR order)
%   minimised jointly over K and h; the axis's K is the argmin and the trial's K the larger of the two axes' selections. Sensitivity row: the plain white-noise criterion
%       C0(K) = N ln var(r_K) + 2 K ln N.
%   Modes:  out = harmonicOrderSelect_v001(Mode="acceptance");   % semi-synthetic trials as in harmonicRecoveryTest_v001 (real Fraser residual as noise; sine with 0 or 10% 2nd and 3rd harmonics);
%                                                                 % pass: median K equals the true content and at least 70% of trials select it exactly; over- and under-selection rates reported
%           out = harmonicOrderSelect_v001(Mode="real");         % the selection on real Fraser N=48 positions (NOT to be run before the amendment is frozen)
%   Writes results/harmonicOrderSelect_v001_<Mode><Tag>.mat
arguments
    opts.Mode (1,1) string {mustBeMember(opts.Mode, ["acceptance", "real"])} = "acceptance"
    opts.MaxTrials (1,1) double = 300
    opts.Distortions (1,:) double = [0 0.1]
    opts.Tag (1,1) string = ""
    opts.UseParfor (1,1) logical = true
    opts.Seed (1,1) double = 1729
end
here = fileparts(mfilename('fullpath'));  addpath(genpath(fullfile(here, 'functions')));
resDir = fullfile(fileparts(here), 'results');  if ~isfolder(resDir), mkdir(resDir); end
ws = warning('off', 'all');  wCleanup = onCleanup(@() warning(ws)); %#ok<NASGU>
trials = importDB_fraser_v002(Verbose=false);  n0 = numel(trials);
sel = unique(round(linspace(1, n0, min(opts.MaxTrials, n0))));  trials = trials(sel);  n = numel(trials);
X = {trials.x};  Y = {trials.y};  FS = [trials.fs];
if opts.Mode == "real"
    Kic = NaN(n, 1);  Kbic = NaN(n, 1);  f0 = NaN(n, 1);
    if opts.UseParfor, try, ensurePool_v001(); catch, end, end %#ok<TRYNC>
    if opts.UseParfor
        parfor i = 1:n, [Kic(i), Kbic(i), f0(i)] = selectOne_local(double(X{i}), double(Y{i}), FS(i)); end
    else
        for i = 1:n, [Kic(i), Kbic(i), f0(i)] = selectOne_local(double(X{i}), double(Y{i}), FS(i)); end
    end
    out = struct('Kic', Kic, 'Kbic', Kbic, 'f0', f0, 'trialID', string({trials.trialID})', 'subjectID', string({trials.subjectID})', 'opts', opts, 'generated', string(datetime('now')));
    fprintf('real Fraser N=48: median K (IC) %.0f, median K (plain BIC) %.0f over %d trials\n', median(Kic, 'omitnan'), median(Kbic, 'omitnan'), n);
else
    nD = numel(opts.Distortions);  truth = 1 + 2 * (opts.Distortions > 0);   % true harmonic content: 1 for the pure sine, 3 with weak 2nd and 3rd harmonics
    Kic = NaN(n, nD);  Kbic = NaN(n, nD);
    rs = RandStream('threefry4x64_20', 'Seed', opts.Seed);  PH = 2*pi*rand(rs, n, nD, 2);   % fixed random phases for the weak harmonics
    for d = 1:nD
        tmpIC = NaN(n, 1);  tmpB = NaN(n, 1);
        if opts.UseParfor, try, ensurePool_v001(); catch, end, end %#ok<TRYNC>
        dist = opts.Distortions(d);  phd = squeeze(PH(:, d, :));
        if opts.UseParfor
            parfor i = 1:n, [tmpIC(i), tmpB(i)] = acceptOne_local(double(X{i}), double(Y{i}), FS(i), dist, phd(i, :)); end
        else
            for i = 1:n, [tmpIC(i), tmpB(i)] = acceptOne_local(double(X{i}), double(Y{i}), FS(i), dist, phd(i, :)); end
        end
        Kic(:, d) = tmpIC;  Kbic(:, d) = tmpB;
    end
    rows = table();
    for d = 1:nD
        for m = ["IC", "plainBIC"]
            k = Kic(:, d);  if m == "plainBIC", k = Kbic(:, d); end
            ok = isfinite(k);
            rows = [rows; table(opts.Distortions(d), truth(d), m, median(k(ok)), mean(k(ok) == truth(d)), mean(k(ok) > truth(d)), mean(k(ok) < truth(d)), nnz(ok), ...
                'VariableNames', {'distortion', 'trueK', 'criterion', 'medianK', 'shareExact', 'shareOver', 'shareUnder', 'n'})]; %#ok<AGROW>
        end
    end
    disp(rows);
    pass = all(rows.medianK(rows.criterion == "IC") == truth(:)) && all(rows.shareExact(rows.criterion == "IC") >= 0.70);
    fprintf('ACCEPTANCE (IC rule: median K equals the truth in each scenario and >= 70%% exact): %s\n', ifelse_local(pass, 'PASS', 'FAIL'));
    out = struct('summary', rows, 'Kic', Kic, 'Kbic', Kbic, 'truth', truth, 'pass', pass, 'opts', opts, 'generated', string(datetime('now')));
end
f = fullfile(resDir, "harmonicOrderSelect_v001_" + opts.Mode + opts.Tag + ".mat");  save(f, '-struct', 'out', '-v7.3');  fprintf('Saved: %s\n', f);
end

function [Kic, Kbic, f0] = acceptOne_local(x, y, fs, dist, ph)
% semi-synthetic trial: real 4-harmonic residual as the noise, known movement added (as harmonicRecoveryTest_v001, no drift)
Kic = NaN;  Kbic = NaN;  f0 = estimateF0_local(x, y, fs);
if ~isfinite(f0) || f0 <= 0, return, end
[nx, ny] = templateSubtract_local(x, y, fs, f0, 4, 4);  N = numel(nx);  t = (0:N-1)' / fs;
D = [ones(numel(x), 1), (0:numel(x)-1)'/fs, cos(2*pi*f0*(0:numel(x)-1)'/fs), sin(2*pi*f0*(0:numel(x)-1)'/fs)];
bx = D \ detrend(x(:), 'linear');  by = D \ detrend(y(:), 'linear');  Ax = hypot(bx(3), bx(4));  Ay = hypot(by(3), by(4));
phs = 2*pi*f0*t;  mx = Ax * cos(phs);  my = Ay * sin(phs);
for h = 2:3, mx = mx + dist * Ax * cos(h*phs + ph(1)) / h;  my = my + dist * Ay * sin(h*phs + ph(2)) / h;  end
[Kic, Kbic] = selectCore_local(mx + nx(:), my + ny(:), fs, f0);
end

function [Kic, Kbic, f0] = selectOne_local(x, y, fs)
Kic = NaN;  Kbic = NaN;  f0 = estimateF0_local(x, y, fs);
if ~isfinite(f0) || f0 <= 0, return, end
[Kic, Kbic] = selectCore_local(x, y, fs, f0);
end

function [Kic, Kbic] = selectCore_local(x, y, fs, f0)
Kmax = max(1, min(20, floor(10 / f0)));  ka = NaN(1, 2);  kb = NaN(1, 2);
for a = 1:2
    s = ifelse_local(a == 1, x, y);  bestC = Inf;  bestK = NaN;  bestB = Inf;  bestKb = NaN;
    for K = 1:Kmax
        [rx, ry] = templateSubtract_local(s, s, fs, f0, K, 4);   % same fit on the chosen axis (the second output is a duplicate)
        r = rx(:) - mean(rx);  N = numel(r);  lnN = log(N);  H = min(20, round(lnN^2));
        cB = N * log(var(r)) + 2 * K * lnN;  if cB < bestB, bestB = cB;  bestKb = K; end
        for h = 1:H
            s2 = arPredVar_local(r, h);  if ~isfinite(s2) || s2 <= 0, continue, end
            C = N * log(s2) + (2*K + h) * lnN;
            if C < bestC, bestC = C;  bestK = K; end
        end
    end
    ka(a) = bestK;  kb(a) = bestKb;
end
Kic = max(ka);  Kbic = max(kb);
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
