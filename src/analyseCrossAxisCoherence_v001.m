function T = analyseCrossAxisCoherence_v001(opts)
%ANALYSECROSSAXISCOHERENCE_V001  Why is bootstrap lower than shaped_xu only in Zarandi? Tests the cross-axis route (SESSION_LOG Addendum 37).
%   shaped_xu draws the PCA major and minor residual axes independently; bootstrap keeps their cross-axis coherence exactly (generateLoopClosureNoise_v002 header). Zarandi's template is rotated 45 degrees (theta median 42.5
%   degrees in the runs), so independent tablet-frame x/y noise of unequal size becomes CORRELATED in the PCA frame, a loss that shaped_xu cannot reproduce and bootstrap does. Axis-aligned datasets lose nothing.
%   Per trial (same chain as the runner: estimateF0, template subtraction 4 harmonics / 4 cycle window, PCA theta from the mm trajectory): zero-lag rho of the tablet-frame residuals (x, y) and of the PCA-frame residuals
%   (maj, min); mean magnitude-squared coherence of maj/min over 1 Hz to min(20, fs/2 - 1); sigma_x/sigma_y. Joined by trialID to the per-trial H1 difference d = betaGenStarMed(bootstrap v017) - betaGenStarMed(shaped_xu v015),
%   and reported per dataset: medians, Spearman(d, |rho_majmin|) and Spearman(d, coherence). Descriptive. Datasets: Zarandi, Hickman PLAC, Fraser (importers as verifyTemplateColumns_v001).
%   T = analyseCrossAxisCoherence_v001(Dir=<folder with the loopClosureResults H1 files>). Saves results/analyseCrossAxisCoherence_v001.mat.   Fraser, D.S. (2026)  v001
arguments
    opts.Dir (1,1) string = ""
    opts.Datasets (1,:) string = ["Zarandi" "Hickman PLAC" "Fraser"]
end
srcDir = fileparts(mfilename('fullpath'));  addpath(srcDir);  addpath(genpath(fullfile(srcDir, 'functions')));  addpath(genpath(fullfile(srcDir, 'req')));
if opts.Dir == "", opts.Dir = "/Volumes/rdsprojects/f/fraserds-mpo-evaluation/2026_prereg/velocity-curvature-power-law-simulation-main/velocity-curvature-power-law-simulation-main/src"; end   % bootstrap v017 files live on the RDS only
T = struct('dataset', {}, 'n', {}, 'rhoXY', {}, 'rhoMajMin', {}, 'coh', {}, 'sdRatio', {}, 'dMed', {}, 'spRho', {}, 'spCoh', {}, 'perTrial', {});
fprintf('%-13s %4s | %9s %11s %7s %8s | %8s | %12s %12s\n', 'dataset', 'n', 'med rhoXY', 'med|rhoMaMi|', 'cohMM', 'sdx/sdy', 'med d', 'Sp(d,|rho|)', 'Sp(d,coh)');
for d = opts.Datasets
    tag = strrep(d, ' ', '_');  xv = "v015";  if d == "Zarandi", xv = "v017"; end   % shaped_xu v017 exists on the RDS for Zarandi only; others compare v015 (disclosed in the printout)
    fprintf('%s: shaped_xu %s vs bootstrap v017\n', d, xv);  A = load_local(fullfile(opts.Dir, sprintf('loopClosureResults_%s_all_shaped_xu_%s.mat', tag, xv)));  B = load_local(fullfile(opts.Dir, sprintf('loopClosureResults_%s_all_bootstrap_v017.mat', tag)));
    [~, ia, ib] = intersect(A.ids, B.ids, 'stable');  dd = B.bgs(ib) - A.bgs(ia);  ids = A.ids(ia);
    [sigToMM, FS, importFn] = dataset_local(d);  trials = importFn();  tids = string({trials.trialID})';  [tf, loc] = ismember(ids, tids);
    if ~all(tf), error('xcoh:match', '%s: %d result trials not found in the importer', d, nnz(~tf)), end
    n = numel(ids);  M = nan(n, 4);                                   % rhoXY, rhoMajMin, coherence, sdx/sdy
    for i = 1:n
        x = double(trials(loc(i)).x(:));  y = double(trials(loc(i)).y(:));  f0 = estimateF0_local(x, y, FS);  [rx, ry] = templateSubtract_local(x, y, FS, f0, 4, 4);
        xm = x * sigToMM;  xm = xm - mean(xm);  ym = y * sigToMM;  ym = ym - mean(ym);  [V, D] = eig(cov(xm, ym));  [~, o] = sort(diag(D), 'descend');  V = V(:, o);  th = atan2(V(2, 1), V(1, 1));
        rmaj = rx * cos(th) + ry * sin(th);  rmin = -rx * sin(th) + ry * cos(th);  ns = min(256, floor(numel(rmaj) / 4));
        if ns < 32 || ~all(isfinite([rx; ry])), continue, end
        [Cxy, fq] = mscohere(rmaj, rmin, hamming(ns), [], [], FS);  band = fq >= 1 & fq <= min(20, FS / 2 - 1);
        M(i, :) = [corr(rx, ry), corr(rmaj, rmin), mean(Cxy(band)), std(rx) / std(ry)];
    end
    ok = isfinite(dd) & all(isfinite(M), 2);  sp1 = corr(dd(ok), abs(M(ok, 2)), 'Type', 'Spearman');  sp2 = corr(dd(ok), M(ok, 3), 'Type', 'Spearman');
    T(end + 1) = struct('dataset', d, 'n', nnz(ok), 'rhoXY', median(M(ok, 1)), 'rhoMajMin', median(abs(M(ok, 2))), 'coh', median(M(ok, 3)), 'sdRatio', median(M(ok, 4)), 'dMed', median(dd(ok)), ...
        'spRho', sp1, 'spCoh', sp2, 'perTrial', table(ids(ok), dd(ok), M(ok, 1), M(ok, 2), M(ok, 3), M(ok, 4), 'VariableNames', {'id', 'd', 'rhoXY', 'rhoMajMin', 'coh', 'sdRatio'})); %#ok<AGROW>
    fprintf('%-13s %4d | %9.3f %11.3f %7.3f %8.2f | %+8.3f | %12.3f %12.3f\n', d, nnz(ok), T(end).rhoXY, T(end).rhoMajMin, T(end).coh, T(end).sdRatio, T(end).dMed, sp1, sp2);
end
save(fullfile(srcDir, '..', 'results', 'analyseCrossAxisCoherence_v001.mat'), 'T', 'opts');  fprintf('\nSaved results/analyseCrossAxisCoherence_v001.mat\n=== CROSS-AXIS DONE ===\n');
end

function M = load_local(f)
    if ~isfile(f), error('xcoh:missing', '%s not found', f), end
    r = load(f, 'results').results(:);  M.ids = string({r.trialID})';  M.bgs = [r.betaGenStarMed]';
    if numel(unique(M.ids)) ~= numel(M.ids), error('xcoh:dup', 'duplicate trialIDs in %s', f), end
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
