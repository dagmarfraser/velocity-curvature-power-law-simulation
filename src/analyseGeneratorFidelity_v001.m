function T = analyseGeneratorFidelity_v001(opts)
%ANALYSEGENERATORFIDELITY_V001  Which noise generator departs from the real residuals in Zarandi? (SESSION_LOG Addendum 38)
%   Per trial (runner-identical residuals, PCA frame) and per model (shaped_xu, bootstrap), K draws are made with generateLoopClosureNoise_v003 and four properties compared with the REAL residuals, reported as
%   generated / real (1 = faithful; medians over trials of the median over draws):  sdRatio = sd(maj)/sd(min);  rhoMM = |zero-lag corr(maj, min)| (difference, not ratio);  cyc = autocorrelation of maj at one cycle lag (1/f0),
%   i.e. cycle-locked residual structure left after template subtraction;  alpha = IRASA alpha of maj (difference).  shaped_xu randomises phase (Xu), bootstrap keeps the real temporal structure, so a cyc gap points at
%   cycle-locked residual as the carrier of the bootstrap offset.  Descriptive; a random trial subset per dataset (SubsetN, seed 1729).  Saves results/analyseGeneratorFidelity_v001.mat.   Fraser, D.S. (2026)  v001
arguments
    opts.Datasets (1,:) string = ["Zarandi" "Hickman PLAC" "Fraser"]
    opts.SubsetN (1,1) double = 60
    opts.K (1,1) double = 20
end
srcDir = fileparts(mfilename('fullpath'));  addpath(srcDir);  addpath(genpath(fullfile(srcDir, 'functions')));  addpath(genpath(fullfile(srcDir, 'req')));
MODELS = ["shaped_xu" "bootstrap"];  T = struct('dataset', {}, 'model', {}, 'n', {}, 'sdRatioReal', {}, 'sdRatioGen', {}, 'dRhoMM', {}, 'cycReal', {}, 'cycGen', {}, 'dAlpha', {});
fprintf('%-13s %-10s %3s | %9s %9s | %8s | %8s %8s | %8s\n', 'dataset', 'model', 'n', 'anisReal', 'anisGen', 'd|rhoMM|', 'cycReal', 'cycGen', 'dAlphaMaj');
for d = opts.Datasets
    [sigToMM, FS, importFn] = dataset_local(d);  trials = importFn();  rs = RandStream('twister', 'Seed', 1729);  pick = sort(randperm(rs, numel(trials), min(opts.SubsetN, numel(trials))));
    R = nan(numel(pick), 2, 6);   % per trial, model: anisReal, anisGen, dRhoMM, cycReal, cycGen, dAlpha
    for i = 1:numel(pick)
        x = double(trials(pick(i)).x(:));  y = double(trials(pick(i)).y(:));  f0 = estimateF0_local(x, y, FS);  [rx, ry] = templateSubtract_local(x, y, FS, f0, 4, 4);
        xm = x * sigToMM;  xm = xm - mean(xm);  ym = y * sigToMM;  ym = ym - mean(ym);  [V, D] = eig(cov(xm, ym));  [~, o] = sort(diag(D), 'descend');  V = V(:, o);  th = atan2(V(2, 1), V(1, 1));
        rmaj = rx * cos(th) + ry * sin(th);  rmin = -rx * sin(th) + ry * cos(th);  lag = round(FS / f0);
        if numel(rmaj) < 4 * lag || lag < 1, continue, end
        aj = iraAlphaSigma_v001(rmaj, FS, 1.0, min(20, FS / 2 - 1), 1.1:0.05:1.9);  an = iraAlphaSigma_v001(rmin, FS, 1.0, min(20, FS / 2 - 1), 1.1:0.05:1.9);
        if ~isfinite(aj) || ~isfinite(an), continue, end
        for m = 1:2
            g = nan(opts.K, 4);
            for k = 1:opts.K
                if MODELS(m) == "shaped_xu", [gj, gn] = generateLoopClosureNoise_v003('shaped_xu', rmaj, rmin, FS, aj, an);
                else, [gj, gn] = generateLoopClosureNoise_v003('bootstrap', rmaj, rmin, FS); end
                g(k, :) = [std(gj) / std(gn), abs(corr(gj, gn)), acf_local(gj, lag), iraAlphaSigma_v001(gj, FS, 1.0, min(20, FS / 2 - 1), 1.1:0.05:1.9)];
            end
            R(i, m, :) = [std(rmaj) / std(rmin), median(g(:, 1), 'omitnan'), median(g(:, 2), 'omitnan') - abs(corr(rmaj, rmin)), acf_local(rmaj, lag), median(g(:, 3), 'omitnan'), median(g(:, 4), 'omitnan') - aj];
        end
    end
    for m = 1:2
        Q = squeeze(R(:, m, :));  ok = all(isfinite(Q), 2);  q = median(Q(ok, :), 1);
        T(end + 1) = struct('dataset', d, 'model', MODELS(m), 'n', nnz(ok), 'sdRatioReal', q(1), 'sdRatioGen', q(2), 'dRhoMM', q(3), 'cycReal', q(4), 'cycGen', q(5), 'dAlpha', q(6)); %#ok<AGROW>
        fprintf('%-13s %-10s %3d | %9.2f %9.2f | %+8.3f | %8.3f %8.3f | %+8.3f\n', d, MODELS(m), nnz(ok), q);
    end
end
save(fullfile(srcDir, '..', 'results', 'analyseGeneratorFidelity_v001.mat'), 'T', 'opts');  fprintf('\nSaved results/analyseGeneratorFidelity_v001.mat\n=== GENERATOR FIDELITY DONE ===\n');
end

function r = acf_local(v, lag)
    v = v(:) - mean(v);  r = sum(v(1:end - lag) .* v(1 + lag:end)) / sum(v .^ 2);
end

function [sigToMM, FS, importFn] = dataset_local(name)
    switch name
        case "Fraser",       sigToMM = 1/10.41793; FS = 240; importFn = @() importDB_fraser_v001(Verbose=false);
        case "Hickman PLAC", sigToMM = 0.248;      FS = 133; importFn = @() importDB_hickman_v003(Study=2, Group="PLAC", Verbose=false);
        case "Zarandi",      sigToMM = 10.0;       FS = 100; importFn = @() importDB_zarandi_v001(Verbose=false);
        otherwise, error('xcoh:dataset', 'Dataset "%s" not supported', name);
    end
end

% estimateF0_local and templateSubtract_local: copied from runLoopClosureFftnoise_v017.m (the runner's own, so the residuals match the runs).
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

function [resX, resY] = templateSubtract_local(x, y, fs, f0, nH, nCW)
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
    resX = resX(cs:ce);
    resY = resY(cs:ce);
end
