function out = harmonicCountDiagnostic_v001(opts)
%HARMONICCOUNTDIAGNOSTIC_V001 Does the template's harmonic count matter, and what would a data-driven stopping rule choose? (2026-11-17)
%   Dagmar 2026-11-17: rule (b) of docs/TemplateJustification_LitCheck_20261115_v001.md section 3 item 2, run on the Fraser N=48 corpus, using criteria that
%   do NOT involve beta_gen*: (1) leftover periodic power at k*f0 after subtraction, (2) notches at k*f0 (over-removal), (3) residual noise colour (alpha) against the
%   harmonic count. Per trial, for each count nH in opts.Counts, the production end-anchored template (same fit as runLoopClosureFftnoise_v028 templateSubtract_local:
%   constant, linear trend, nH harmonics, Hann-weighted windows of 4 cycles, hop half a window, final end-anchored window) is subtracted, and on the residual of each axis:
%     peakRatio(nH, k)  periodogram power at k f0 over the median power in the flanking bands (0.25-0.75 f0 either side of k f0), averaged over the two axes;
%                       a value well above 1 at k <= nH means the template left a periodic line behind; below 1 at k <= nH means a notch (power removed from the noise floor);
%     alpha(nH)         PWLW exponent of the residual over the band [1, min(20, fs/2-1)] Hz (as the runner's alpha selection), averaged over the axes.
%   Rule (b), on the RAW detrended record: the number of consecutive harmonics k = 1, 2, ... whose power at k f0 exceeds the flanking median by the factor opts.Factor
%   (default 4.6 = -log(0.01) for an exponentially distributed periodogram ordinate, p < 0.01 each, using the median-based floor mean = median / log 2), stopping at the
%   first harmonic that fails; nHb is that count (at least 1). Counted on each axis, then the larger of the two.
%   Nothing here uses beta_obs, beta_gen* or any forward map. Fraser N=48 trials come from functions/importDB_fraser_v002 (2,877 trials).
%   out = harmonicCountDiagnostic_v001();                       all trials, parfor
%   out = harmonicCountDiagnostic_v001(MaxTrials=40);           quick local test (evenly spaced trials)
%   Writes results/harmonicCountDiagnostic_v001.mat (or _<Tag>) with out.T (per trial) and out.summary.
arguments
    opts.Counts (1,:) double {mustBeInteger, mustBePositive} = [1 2 3 4 5 6 8 10 15 20]
    opts.MaxTrials (1,1) double = Inf
    opts.Factor (1,1) double = 4.6
    opts.MaxK (1,1) double {mustBeInteger, mustBePositive} = 20
    opts.Tag (1,1) string = ""
    opts.UseParfor (1,1) logical = true
end
here = fileparts(mfilename('fullpath'));  addpath(genpath(fullfile(here, 'functions')));
resDir = fullfile(fileparts(here), 'results');  if ~isfolder(resDir), mkdir(resDir); end
trials = importDB_fraser_v002(Verbose=false);
n0 = numel(trials);  sel = 1:n0;
if isfinite(opts.MaxTrials), sel = unique(round(linspace(1, n0, min(opts.MaxTrials, n0)))); end
trials = trials(sel);  n = numel(trials);  counts = opts.Counts;  nC = numel(counts);
fprintf('harmonicCountDiagnostic_v001: %d of %d trials, counts %s, factor %.2f\n', n, n0, mat2str(counts), opts.Factor);
f0 = NaN(n, 1);  nHb = NaN(n, 1);  rawRatio = NaN(n, opts.MaxK);  alphaN = NaN(n, nC);  peak = NaN(n, nC, 6);  nTail = NaN(n, nC);
ids = strings(n, 1);  subj = strings(n, 1);
X = {trials.x};  Y = {trials.y};  FS = [trials.fs];
for i = 1:n, ids(i) = string(trials(i).trialID);  subj(i) = string(trials(i).subjectID); end
tStart = tic;
if opts.UseParfor
    try, ensurePool_v001(); catch, end   %#ok<TRYNC>
    parfor i = 1:n
        [f0(i), nHb(i), rawRatio(i, :), alphaN(i, :), peak(i, :, :), nTail(i, :)] = oneTrial_local(double(X{i}), double(Y{i}), FS(i), counts, opts.Factor, opts.MaxK);
    end
else
    for i = 1:n
        [f0(i), nHb(i), rawRatio(i, :), alphaN(i, :), peak(i, :, :), nTail(i, :)] = oneTrial_local(double(X{i}), double(Y{i}), FS(i), counts, opts.Factor, opts.MaxK);
    end
end
fprintf('done in %.0f s\n', toc(tStart));
T = table(ids, subj, f0, nHb, 'VariableNames', {'trialID', 'subjectID', 'f0', 'nHb'});
cap = floor(10 ./ f0);   % harmonics up to the 10 Hz cutoff (the 'fc' rule)
T.nFc = cap;
alphaT = array2table(alphaN, 'VariableNames', "alpha_nH" + counts);  T = [T alphaT];
% summary: median alpha by count, median leftover peak ratio at k = 1..3 after subtracting >= k harmonics, and notch share at k <= nH
S = table(counts', median(alphaN, 1, 'omitnan')', 'VariableNames', {'nH', 'alphaMedian'});
S.alphaP25 = prctile(alphaN, 25, 1)';  S.alphaP75 = prctile(alphaN, 75, 1)';
for k = 1:3
    v = NaN(nC, 1);
    for c = 1:nC, if counts(c) >= k, v(c) = median(peak(:, c, k), 'omitnan'); end, end
    S.("peakRatio_k" + k) = v;
end
v = NaN(nC, 1);  for c = 1:nC, kk = 1:min(counts(c), 6);  pk = peak(:, c, kk);  v(c) = mean(pk(isfinite(pk)) < 1); end
S.shareBelowOne_kLEnH = v;   % share of (trial, k <= min(nH, 6)) with the line at k f0 below the flanking floor: notch
fprintf('\nresidual alpha and leftover lines by harmonic count (median over trials):\n');  disp(S);
fprintf('rule (b) nHb: median %.1f, IQR %.1f-%.1f, share nHb <= 4: %.1f%%, share nHb >= 10: %.1f%%; fc rule (10 Hz / f0): median %.1f\n', ...
    median(nHb, 'omitnan'), prctile(nHb, 25), prctile(nHb, 75), 100 * mean(nHb <= 4, 'omitnan'), 100 * mean(nHb >= 10, 'omitnan'), median(cap, 'omitnan'));
out = struct('T', T, 'summary', S, 'peak', peak, 'rawRatio', rawRatio, 'counts', counts, 'opts', opts, 'nTrials', n, 'generated', string(datetime('now')));
f = fullfile(resDir, "harmonicCountDiagnostic_v001" + opts.Tag + ".mat");  save(f, '-struct', 'out', '-v7.3');  fprintf('Saved: %s\n', f);
end

function [f0, nHb, rawRatio, alphaN, peak, nTail] = oneTrial_local(x, y, fs, counts, factor, maxK)
nC = numel(counts);  alphaN = NaN(1, nC);  peak = NaN(1, nC, 6);  nTail = NaN(1, nC);  rawRatio = NaN(1, maxK);  nHb = NaN;
f0 = estimateF0_local(x, y, fs);
if ~isfinite(f0) || f0 <= 0, return, end
% rule (b) on the raw detrended record (both axes; the larger count)
nb = zeros(1, 2);
for a = 1:2
    s = ifelse_local(a == 1, x, y);  r = lineRatio_local(detrend(s(:), 'linear'), fs, f0, maxK);
    if a == 1, rawRatio = r; else, rawRatio = max(rawRatio, r); end
    fail = find(~(r > factor), 1);  if isempty(fail), nb(a) = maxK; else, nb(a) = fail - 1; end
end
nHb = max(1, max(nb));
band = [1, min(20, fs/2 - 1)];
for c = 1:nC
    [rx, ry, tailN] = templateSubtract_local(x, y, fs, f0, counts(c), 4);
    nTail(c) = tailN;
    ax = NaN(1, 2);  pk = NaN(2, 6);
    for a = 1:2
        s = ifelse_local(a == 1, rx, ry);  s = s(:) - mean(s);
        try, ax(a) = alphaPWLW_v001(s, fs, band); catch, ax(a) = NaN; end
        r = lineRatio_local(s, fs, f0, 6);  pk(a, 1:numel(r)) = r;
    end
    alphaN(c) = mean(ax, 'omitnan');
    peak(1, c, :) = mean(pk, 1, 'omitnan');
end
end

function r = lineRatio_local(s, fs, f0, K)
% power at k f0 (nearest bin of a Hann periodogram, zero-padded x8) over the flanking floor (median power at 0.25-0.75 f0 either side of k f0, scaled by 1/log 2)
N = numel(s);  nfft = 2^nextpow2(8 * N);  w = 0.5 * (1 - cos(2*pi*(0:N-1)' / max(N - 1, 1)));
P = abs(fft((s(:) - mean(s)) .* w, nfft)) .^ 2;  P = P(1:nfft/2);  fa = (0:nfft/2-1)' * fs / nfft;
r = NaN(1, K);
for k = 1:K
    fk = k * f0;  if fk >= fs / 2 - 0.5, break, end
    [~, i0] = min(abs(fa - fk));
    fl = (fa >= fk - 0.75*f0 & fa <= fk - 0.25*f0) | (fa >= fk + 0.25*f0 & fa <= fk + 0.75*f0);
    fl = fl & fa > 0;
    if nnz(fl) < 5, continue, end
    floorMean = median(P(fl)) / log(2);
    r(k) = max(P(max(i0-2, 1):min(i0+2, numel(P)))) / floorMean;
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
