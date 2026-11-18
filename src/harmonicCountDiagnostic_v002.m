function out = harmonicCountDiagnostic_v002(opts)
%HARMONICCOUNTDIAGNOSTIC_V002 Does the template's harmonic count matter? Pooled spectra in units of f0, per-harmonic amplitude against the noise between harmonics, f0 precision, residual colour (2026-11-17).
%   Replaces v001's per-trial line test, which cannot work on a record of about 4 cycles (Hann line width about +/-f0/2). Fraser N=48 trials (functions/importDB_fraser_v002).
%   Per trial, no beta_obs, beta_gen* or forward map is used:
%     spec     log10 Hann periodogram of the raw detrended record, mean over the two axes, interpolated onto a grid in units of f0 (g = 0.5 : 0.05 : 20.5), minus its median, so trials weigh equally;
%              pooled over trials as the mean of the log (a geometric mean). Harmonics appear as peaks at whole numbers, the floor between them at half integers. specRes{nH}: the same for the
%              residual after the production end-anchored template with nH harmonics (nH in opts.ResCounts): leftover lines (peaks) and over-removal (notches at whole numbers).
%     amp, ampNull   one least-squares fit over the whole record (constant, linear trend, then a cosine and sine at each of the frequencies below), per axis, amplitude = hypot(cos coef, sin coef):
%              amp(k) at k f0 and ampNull(k) the geometric mean of the amplitudes at (k - 0.5) f0 and (k + 0.5) f0 (the noise between harmonics). snr(k) = amp(k) / ampNull(k), mean of the two axes.
%              snr near 1 means the harmonic is no larger than the noise next to it. Pooled: the median over trials, and the share of trials with snr > opts.SnrCut (default 2).
%     f0, f0fit, relErr   f0 from the runner's FFT peak estimate; f0fit by maximising the power of a 3-harmonic fit over f0 * (0.9 : 0.002 : 1.1); relErr = (f0fit - f0) / f0.
%                         k * |relErr| against 0.25 (about one resolution cell of a 4-cycle record, in f0 units) shows from which harmonic on the template's frequency is mis-specified.
%     alpha(nH)  PWLW exponent of the template residual over [1, min(20, fs/2 - 1)] Hz, nH in opts.Counts (as the runner's alpha selection), mean of the axes.
%   out = harmonicCountDiagnostic_v002();  out = harmonicCountDiagnostic_v002(MaxTrials=40, UseParfor=false, Tag="_test");   Writes results/harmonicCountDiagnostic_v002<Tag>.mat
arguments
    opts.Counts (1,:) double {mustBeInteger, mustBePositive} = [1 2 3 4 5 6 8 10 15 20]
    opts.ResCounts (1,:) double {mustBeInteger, mustBeNonnegative} = [0 4 10 20]
    opts.MaxTrials (1,1) double = Inf
    opts.SnrCut (1,1) double = 2
    opts.Tag (1,1) string = ""
    opts.UseParfor (1,1) logical = true
end
here = fileparts(mfilename('fullpath'));  addpath(genpath(fullfile(here, 'functions')));
resDir = fullfile(fileparts(here), 'results');  if ~isfolder(resDir), mkdir(resDir); end
trials = importDB_fraser_v002(Verbose=false);
n0 = numel(trials);  sel = 1:n0;
if isfinite(opts.MaxTrials), sel = unique(round(linspace(1, n0, min(opts.MaxTrials, n0)))); end
trials = trials(sel);  n = numel(trials);  counts = opts.Counts;  nC = numel(counts);  rc = opts.ResCounts;  nR = numel(rc);
G = 0.5:0.05:20.5;  nG = numel(G);  K = 20;
fprintf('harmonicCountDiagnostic_v002: %d of %d trials\n', n, n0);
f0 = NaN(n, 1);  f0fit = NaN(n, 1);  alphaN = NaN(n, nC);  spec = NaN(n, nG);  specRes = NaN(n, nR, nG);  amp = NaN(n, K);  ampNull = NaN(n, K);
X = {trials.x};  Y = {trials.y};  FS = [trials.fs];  ids = strings(n, 1);  subj = strings(n, 1);
for i = 1:n, ids(i) = string(trials(i).trialID);  subj(i) = string(trials(i).subjectID); end
tStart = tic;
if opts.UseParfor
    try, ensurePool_v001(); catch, end   %#ok<TRYNC>
    parfor i = 1:n
        [f0(i), f0fit(i), alphaN(i, :), spec(i, :), specRes(i, :, :), amp(i, :), ampNull(i, :)] = oneTrial_local(double(X{i}), double(Y{i}), FS(i), counts, rc, G, K);
    end
else
    for i = 1:n
        [f0(i), f0fit(i), alphaN(i, :), spec(i, :), specRes(i, :, :), amp(i, :), ampNull(i, :)] = oneTrial_local(double(X{i}), double(Y{i}), FS(i), counts, rc, G, K);
    end
end
fprintf('done in %.0f s\n', toc(tStart));
relErr = (f0fit - f0) ./ f0;  snr = amp ./ ampNull;
T = table(ids, subj, f0, f0fit, relErr, 'VariableNames', {'trialID', 'subjectID', 'f0', 'f0fit', 'relErr'});
T = [T array2table(alphaN, 'VariableNames', "alpha_nH" + counts)];
S = table(counts', median(alphaN, 1, 'omitnan')', prctile(alphaN, 25, 1)', prctile(alphaN, 75, 1)', 'VariableNames', {'nH', 'alphaMedian', 'alphaP25', 'alphaP75'});
H = table((1:K)', median(snr, 1, 'omitnan')', mean(snr > opts.SnrCut, 1, 'omitnan')', median(amp, 1, 'omitnan')', median(ampNull, 1, 'omitnan')', ...
    median(abs(relErr), 'omitnan') * (1:K)', 'VariableNames', {'k', 'snrMedian', 'shareSnrAboveCut', 'ampMedian', 'ampNullMedian', 'kTimesMedianAbsRelErr'});
fprintf('\nresidual alpha by harmonic count:\n');  disp(S);
fprintf('per-harmonic amplitude over the noise between harmonics (median over trials) and f0 precision:\n');  disp(H);
fprintf('median |relErr| of f0 = %.4f; harmonic at which k * |relErr| first exceeds 0.25: k = %s\n', median(abs(relErr), 'omitnan'), string(find(H.kTimesMedianAbsRelErr > 0.25, 1)));
specMean = mean(spec - median(spec, 2, 'omitnan'), 1, 'omitnan');
resMean = squeeze(mean(specRes - median(specRes, 3, 'omitnan'), 1, 'omitnan'));
out = struct('T', T, 'alphaSummary', S, 'harmonicSummary', H, 'G', G, 'specMean', specMean, 'resMean', resMean, 'resCounts', rc, 'amp', amp, 'ampNull', ampNull, 'snr', snr, 'counts', counts, 'opts', opts, 'nTrials', n, 'generated', string(datetime('now')));
f = fullfile(resDir, "harmonicCountDiagnostic_v002" + opts.Tag + ".mat");  save(f, '-struct', 'out', '-v7.3');  fprintf('Saved: %s\n', f);
end

function [f0, f0fit, alphaN, spec, specRes, amp, ampNull] = oneTrial_local(x, y, fs, counts, rc, G, K)
nC = numel(counts);  nR = numel(rc);  alphaN = NaN(1, nC);  spec = NaN(1, numel(G));  specRes = NaN(1, nR, numel(G));  amp = NaN(1, K);  ampNull = NaN(1, K);  f0fit = NaN;
f0 = estimateF0_local(x, y, fs);
if ~isfinite(f0) || f0 <= 0, return, end
spec = specOnGrid_local(detrend(x(:), 'linear'), detrend(y(:), 'linear'), fs, f0, G);
f0fit = fitF0_local(x, y, fs, f0);
[amp, ampNull] = harmAmp_local(x, y, fs, f0, K);
band = [1, min(20, fs/2 - 1)];
for c = 1:nC
    [rx, ry] = templateSubtract_local(x, y, fs, f0, counts(c), 4);
    ax = NaN(1, 2);
    for a = 1:2, s = ifelse_local(a == 1, rx, ry);  s = s(:) - mean(s);  try, ax(a) = alphaPWLW_v001(s, fs, band); catch, ax(a) = NaN; end, end
    alphaN(c) = mean(ax, 'omitnan');
end
for c = 1:nR
    if rc(c) == 0, rx = detrend(x(:), 'linear');  ry = detrend(y(:), 'linear'); else, [rx, ry] = templateSubtract_local(x, y, fs, f0, rc(c), 4); end
    specRes(1, c, :) = specOnGrid_local(rx, ry, fs, f0, G);
end
end

function s = specOnGrid_local(a, b, fs, f0, G)
% mean over the two axes of log10 Hann periodogram, interpolated at G * f0 (NaN above Nyquist)
N = numel(a);  nfft = 2^nextpow2(8 * N);  w = 0.5 * (1 - cos(2*pi*(0:N-1)' / max(N - 1, 1)));  fa = (0:nfft/2-1)' * fs / nfft;
Pa = abs(fft((a(:) - mean(a)) .* w, nfft)) .^ 2;  Pb = abs(fft((b(:) - mean(b)) .* w, nfft)) .^ 2;
L = 0.5 * (log10(Pa(1:nfft/2) + eps) + log10(Pb(1:nfft/2) + eps));
fq = G * f0;  s = NaN(1, numel(G));  ok = fq < fs / 2 - 0.5;
s(ok) = interp1(fa, L, fq(ok), 'linear');
end

function f0fit = fitF0_local(x, y, fs, f0)
% maximise the explained power of a 3-harmonic fit over f0 * (0.9 : 0.002 : 1.1)
x = detrend(x(:), 'linear');  y = detrend(y(:), 'linear');  t = (0:numel(x)-1)' / fs;  grid = f0 * (0.9:0.002:1.1);  e = zeros(size(grid));
for g = 1:numel(grid)
    D = [];  for h = 1:3, D = [D, cos(2*pi*h*grid(g)*t), sin(2*pi*h*grid(g)*t)]; end %#ok<AGROW>
    bx = D \ x;  by = D \ y;  e(g) = sum((D*bx).^2) + sum((D*by).^2);
end
[~, i] = max(e);  f0fit = grid(i);
end

function [amp, ampNull] = harmAmp_local(x, y, fs, f0, K)
x = detrend(x(:), 'linear');  y = detrend(y(:), 'linear');  t = (0:numel(x)-1)' / fs;
fr = (0.5:0.5:K + 0.5) * f0;  fr = fr(fr < fs / 2 - 0.5);  nf = numel(fr);  A = NaN(2, nf);
D0 = [ones(numel(t), 1), t];
for j = 1:nf
    D = [D0, cos(2*pi*fr(j)*t), sin(2*pi*fr(j)*t)];
    bx = D \ x;  by = D \ y;  A(1, j) = hypot(bx(3), bx(4));  A(2, j) = hypot(by(3), by(4));
end
Am = mean(A, 1);  amp = NaN(1, K);  ampNull = NaN(1, K);
for k = 1:K
    ik = 2 * k;  if ik + 1 > nf, break, end        % index of k f0 is 2k (frequencies are 0.5, 1, 1.5, ...); neighbours 2k-1 and 2k+1
    amp(k) = Am(ik);  ampNull(k) = sqrt(Am(ik - 1) * Am(ik + 1));
end
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

function [resX, resY, tailN] = templateSubtract_local(x, y, fs, f0, nH, nCW)   % as runLoopClosureFftnoise_v028, END-ANCHORED
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
    lastEnd = 0;  if ~isempty(starts), lastEnd = starts(end) + win - 1; end
    tailN = max(0, ce - lastEnd);
    resX = resX(cs:ce);
    resY = resY(cs:ce);
end
