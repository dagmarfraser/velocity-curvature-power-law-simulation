function out = harmonicRecoveryTest_v001(opts)
%HARMONICRECOVERYTEST_V001 Truth test for the template's harmonic count (2026-11-17). Which count recovers a KNOWN noise best?
%   Semi-synthetic trials: the noise is the REAL residual of a Fraser N=48 trial (its end-anchored 4-harmonic template residual, per axis), so its colour and size are the data's; the movement is
%   synthetic and known: a sine at the trial's own f0 with the trial's own fundamental amplitude (one cycle-phase per axis, quarter period apart, an ellipse), optionally with weak 2nd and 3rd harmonics
%   (opts.Distortion, as a fraction of the fundamental, random phases) and optionally a tempo drift (opts.Drift, relative linear change of frequency across the record).
%   Each template count nH in opts.Counts is then run on x_syn = movement + noise and the recovered residual is compared with the injected noise on the samples the template returns:
%     alphaBias   PWLW exponent of the recovered residual minus that of the injected noise (band [1, 20] Hz), mean over the two axes
%     sdRatio     std(recovered) / std(injected)
%     errRel      std(recovered - injected) / std(injected)    (0 = perfect recovery of the noise sample by sample)
%   Scenarios (Distortion, Drift): none, distorted, drifting, both. Nothing uses beta_obs, beta_gen* or any forward map.
%   out = harmonicRecoveryTest_v001();  out = harmonicRecoveryTest_v001(MaxTrials=60, Tag="_test");   Writes results/harmonicRecoveryTest_v001<Tag>.mat
arguments
    opts.Counts (1,:) double {mustBeInteger, mustBePositive} = [1 2 3 4 6 10 20]
    opts.MaxTrials (1,1) double = 300
    opts.Scenarios (:,2) double = [0 0; 0.1 0; 0 0.03; 0.1 0.03]
    opts.Tag (1,1) string = ""
    opts.Seed (1,1) double = 1729
end
here = fileparts(mfilename('fullpath'));  addpath(genpath(fullfile(here, 'functions')));
ws = warning('off', 'all');  wCleanup = onCleanup(@() warning(ws)); %#ok<NASGU>   % alphaPWLW's DFA pre-estimate warns about short segments; the estimate is still returned
resDir = fullfile(fileparts(here), 'results');  if ~isfolder(resDir), mkdir(resDir); end
trials = importDB_fraser_v002(Verbose=false);  n0 = numel(trials);
sel = unique(round(linspace(1, n0, min(opts.MaxTrials, n0))));  trials = trials(sel);  n = numel(trials);
counts = opts.Counts;  nC = numel(counts);  nS = size(opts.Scenarios, 1);
rs = RandStream('threefry4x64_20', 'Seed', opts.Seed);
alphaBias = NaN(n, nC, nS);  sdRatio = NaN(n, nC, nS);  errRel = NaN(n, nC, nS);  alphaTrue = NaN(n, 1);  f0v = NaN(n, 1);
fprintf('harmonicRecoveryTest_v001: %d trials, %d scenarios, counts %s\n', n, nS, mat2str(counts));
tic;
for i = 1:n
    fs = trials(i).fs;  x = double(trials(i).x);  y = double(trials(i).y);  f0 = estimateF0_local(x, y, fs);  f0v(i) = f0;
    if ~isfinite(f0) || f0 <= 0, continue, end
    [nx, ny, ~, cs0] = templateSubtract_local(x, y, fs, f0, 4, 4);   % real residual, the noise
    N = numel(nx);  t = (0:N-1)' / fs;
    % fundamental amplitude per axis from a whole-record fit at f0
    D = [ones(numel(x), 1), (0:numel(x)-1)'/fs, cos(2*pi*f0*(0:numel(x)-1)'/fs), sin(2*pi*f0*(0:numel(x)-1)'/fs)];
    bx = D \ detrend(x(:), 'linear');  by = D \ detrend(y(:), 'linear');  Ax = hypot(bx(3), bx(4));  Ay = hypot(by(3), by(4));
    band = [1, min(20, fs/2 - 1)];
    try, aT = mean([alphaPWLW_v001(nx(:) - mean(nx), fs, band), alphaPWLW_v001(ny(:) - mean(ny), fs, band)]); catch, aT = NaN; end
    alphaTrue(i) = aT;
    for sc = 1:nS
        dist = opts.Scenarios(sc, 1);  drift = opts.Scenarios(sc, 2);
        ph = 2*pi*(t * f0 + 0.5 * drift * f0 * t.^2 / max(t(end), eps));      % instantaneous phase with a linear frequency change of drift * f0 across the record
        mx = Ax * cos(ph);  my = Ay * sin(ph);
        for h = 2:3
            p = 2*pi*rand(rs, 1, 2);  mx = mx + dist * Ax * cos(h*ph + p(1)) / h;  my = my + dist * Ay * sin(h*ph + p(2)) / h;
        end
        xs = mx + nx(:);  ys = my + ny(:);
        for c = 1:nC
            [rx, ry, ~, cs] = templateSubtract_local(xs, ys, fs, f0, counts(c), 4);
            m = numel(rx);  tx = nx(cs:cs+m-1);  ty = ny(cs:cs+m-1);
            try, aR = mean([alphaPWLW_v001(rx(:) - mean(rx), fs, band), alphaPWLW_v001(ry(:) - mean(ry), fs, band)]);
                 aN = mean([alphaPWLW_v001(tx(:) - mean(tx), fs, band), alphaPWLW_v001(ty(:) - mean(ty), fs, band)]);
                 alphaBias(i, c, sc) = aR - aN; catch, end
            sdRatio(i, c, sc) = mean([std(rx) / std(tx), std(ry) / std(ty)]);
            errRel(i, c, sc) = mean([std(rx(:) - tx(:)) / std(tx), std(ry(:) - ty(:)) / std(ty)]);
        end
    end
end
fprintf('done in %.0f s\n', toc);
S = table();
for sc = 1:nS
    for c = 1:nC
        S = [S; table(opts.Scenarios(sc, 1), opts.Scenarios(sc, 2), counts(c), median(alphaBias(:, c, sc), 'omitnan'), prctile(alphaBias(:, c, sc), 25), prctile(alphaBias(:, c, sc), 75), ...
            median(sdRatio(:, c, sc), 'omitnan'), median(errRel(:, c, sc), 'omitnan'), 'VariableNames', {'distortion', 'drift', 'nH', 'alphaBiasMedian', 'alphaBiasP25', 'alphaBiasP75', 'sdRatioMedian', 'errRelMedian'})]; %#ok<AGROW>
    end
end
disp(S);
out = struct('S', S, 'alphaBias', alphaBias, 'sdRatio', sdRatio, 'errRel', errRel, 'alphaTrue', alphaTrue, 'f0', f0v, 'counts', counts, 'opts', opts, 'generated', string(datetime('now')));
f = fullfile(resDir, "harmonicRecoveryTest_v001" + opts.Tag + ".mat");  save(f, '-struct', 'out', '-v7.3');  fprintf('Saved: %s\n', f);
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

function [resX, resY, tailN, cs] = templateSubtract_local(x, y, fs, f0, nH, nCW)   % as runLoopClosureFftnoise_v028, END-ANCHORED; also returns cs (first kept sample)
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
