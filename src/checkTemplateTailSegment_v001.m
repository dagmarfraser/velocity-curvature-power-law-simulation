function S = checkTemplateTailSegment_v001(opts)
%CHECKTEMPLATETAILSEGMENT_V001 Do the end-matched segments of the v023 runs include the production template's uncovered tail?
%   Coherence check, 2026-11-15. The production template (endAnchored = false; the v015 basis and the v023 runs) lays its
%   overlap-add windows from the start of the trial, so the residual between the last window's end and the end clip is never
%   reached: it stays exactly zero, after a step (K-MF5, docs/TODO_PreSubmission_v007: "the residual ends in a single large
%   step followed by an exactly flat stretch"). The end-matched noise takes its spectrum from a segment of that residual; a
%   segment that contains the step and the zero stretch would carry the step's broadband energy into the noise.
%   Per trial (all trials of each arm): uncovered tail length inside the kept residual, the commonLength segments
%   (endToEndSegment_v004, Tol 1, as v023), and whether either axis's segment overlaps the tail. Same template and rotation
%   code as runLoopClosureFftnoise_v023 (production path). Deterministic.
%   Writes results/endMismatchSensitivity/templateTailSegment_v001.csv (per trial) and prints the summary per arm.
arguments
    opts.Arms (1,:) string = ["FraserN48", "Hickman PLAC", "Hickman HALO", "Cook CTRL", "Cook ASD", "Zarandi", "Dhieb"]
end
here = fileparts(mfilename('fullpath'));  addpath(here, fullfile(here, 'functions'));
ensurePool_v001();  wait(parfevalOnAll(gcp, @addpath, 0, here, fullfile(here, 'functions')));
T = table();
for arm = opts.Arms
    switch arm
        case "FraserN48",    FS = 240; tr = importDB_fraser_v002(Verbose=false);
        case "Hickman PLAC", FS = 133; tr = importDB_hickman_v003(Study=2, Group="PLAC", Verbose=false);
        case "Hickman HALO", FS = 133; tr = importDB_hickman_v003(Study=2, Group="HALO", Verbose=false);
        case "Cook CTRL",    FS = 133; tr = importDB_cook_v002(Group="CTRL", Tasks=7, Verbose=false);
        case "Cook ASD",     FS = 133; tr = importDB_cook_v002(Group="ASD", Tasks=7, Verbose=false);
        case "Zarandi",      FS = 100; tr = importDB_zarandi_v001(Verbose=false);
        case "Dhieb",        FS = 100; tr = importDB_dhieb_v001(Verbose=false);
    end
    n = numel(tr);  rows = cell(n, 1);
    parfor t = 1:n
        rows{t} = oneTrial(tr(t), FS);
    end
    Ta = vertcat(rows{:});  Ta.arm = repmat(arm, height(Ta), 1);  T = [T; Ta]; %#ok<AGROW>
    ok = Ta.status == "ok";  hasTail = ok & Ta.tail > 0;
    fprintf('%-13s %4d trials (%d with a segment): tail in %d (%.0f%%), median %.0f samples (%.1f%% of residual); segment overlaps the tail in %d trials (%.0f%% of tailed), median overlap %.0f samples\n', ...
        arm, n, nnz(ok), nnz(hasTail), 100 * nnz(hasTail) / max(nnz(ok), 1), median(Ta.tail(hasTail)), 100 * median(Ta.tail(hasTail) ./ Ta.M(hasTail)), ...
        nnz(hasTail & Ta.overlap > 0), 100 * nnz(hasTail & Ta.overlap > 0) / max(nnz(hasTail), 1), median(Ta.overlap(Ta.overlap > 0)));
end
od = fullfile(here, '..', 'results', 'endMismatchSensitivity');  writetable(T, fullfile(od, 'templateTailSegment_v001.csv'));
S = T;
end

function R = oneTrial(t, FS)
x = double(t.x(:));  y = double(t.y(:));  N = numel(x);  tail = NaN;  ov = NaN;  M = NaN;  L = NaN;  st = "ok";
try
    f0 = estimateF0_b(x, y, FS);
    [xR, yR, cs, lastEnd] = templateSubtract_b(x, y, FS, f0, 4, 4);
    M = numel(xR);  ce = cs + M - 1;  tail = max(0, ce - lastEnd);  t0 = M - tail + 1;     % tail = residual indices t0..M
    if M < 280, error('tail:short', 'residual below the runner minimum'); end
    xc = x - mean(x);  yc = y - mean(y);  [V, D] = eig(cov(xc, yc));  [~, o] = sort(diag(D), 'descend');  V = V(:, o);  th = atan2(V(2,1), V(1,1));
    rMaj = xR*cos(th) + yR*sin(th);  rMin = -xR*sin(th) + yR*cos(th);
    [i0, i1, info] = endToEndSegment_v004([rMaj, rMin], Tol=1, Mode="commonLength");  L = info.L;
    ov = 0;  if tail > 0, ov = max(max(0, i1 - max(i0, t0) + 1)); end
catch ME
    st = "FAIL: " + string(ME.message);
end
R = table(string(t.trialID), N, M, tail, L, ov, st, 'VariableNames', {'trialID', 'N', 'M', 'tail', 'L', 'overlap', 'status'});
end

% ---- runLoopClosureFftnoise_v023 production template path, with the clip and last window end returned ----
function f0 = estimateF0_b(x, y, fs)
    N    = numel(x); nfft = 2^nextpow2(4*N);
    Xf   = abs(fft(detrend(x(:), 'linear'), nfft));
    Yf   = abs(fft(detrend(y(:), 'linear'), nfft));
    fAx  = (0:nfft-1)' * fs / nfft;
    band = fAx > 0.1 & fAx < 5; idx = find(band);
    [~, pkX] = max(Xf(band)); [~, pkY] = max(Yf(band));
    f0x = fAx(idx(pkX)); f0y = fAx(idx(pkY));
    if abs(f0x - f0y) / max(f0x, 0.01) < 0.2, f0 = mean([f0x, f0y]);
    elseif max(Xf(band)) > max(Yf(band)), f0 = f0x;
    else, f0 = f0y; end
end

function [resX, resY, cs, lastEnd] = templateSubtract_b(x, y, fs, f0, nH, nCW)
    N   = numel(x);  win = round(nCW/f0*fs);  hop = round(win/2);  win = min(win, N);  win = max(win, round(2/f0*fs));
    resX = zeros(N,1); resY = zeros(N,1); ww = zeros(N,1);
    starts = 1:hop:(N - win + 1);  if win > N, starts = []; end
    lastEnd = 0;
    for s = starts
        e = s + win - 1;  lastEnd = e;  tw = (0:win-1)'/fs;  han = 0.5*(1 - cos(2*pi*(0:win-1)'/(win-1)));
        D = [ones(win,1), tw];
        for h = 1:nH, D(:,end+1) = cos(2*pi*h*f0*tw); D(:,end+1) = sin(2*pi*h*f0*tw); end %#ok<AGROW>
        Dw = D .* han;  bX = Dw \ (x(s:e) .* han);  bY = Dw \ (y(s:e) .* han);
        resX(s:e) = resX(s:e) + (x(s:e) - D*bX) .* han;  resY(s:e) = resY(s:e) + (y(s:e) - D*bY) .* han;  ww(s:e) = ww(s:e) + han;
    end
    ok = ww > 0;  resX(ok) = resX(ok) ./ ww(ok);  resY(ok) = resY(ok) ./ ww(ok);
    cl = max(round(win/4), 1);  cs = cl;  ce = min(max(N - cl, cs + 100), N);
    resX = resX(cs:ce);  resY = resY(cs:ce);
end
