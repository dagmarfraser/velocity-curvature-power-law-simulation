function T = templateComparison_v002(opts)
%TEMPLATECOMPARISON_V002 Which harmonic template should define the noise? Stage 1: residual measures (2026-11-15).
%   Plan and decision rule: docs/TemplateJustification_LitCheck_20261115_v001.md section 4 (stated before this ran):
%   prefer the template with full coverage that leaves no significant movement power at k*f0 below the filter cutoff while
%   removing the least noise.
%   Templates (all: constant + linear trend + nH harmonics of f0 in Hann-weighted windows of 4 cycles, at least 2, hop half a
%   window, overlap-added; f0 and the end clip exactly as runLoopClosureFftnoise_v024):
%     prod4     production windows (from the start), nH 4, uncovered tail cut (= the running v024 runs)
%     end4      end-anchored windows (one final window ending at the last sample), nH 4
%     endFc     end-anchored, nH = floor(fc / f0), fc = 10 Hz (the BWFD cutoff): every harmonic the pipeline passes
%     endSel    end-anchored, nH chosen per trial: the smallest nH >= 1 at which no remaining harmonic k (nH < k <= fc/f0) of
%               the residual stands out by more than ExcessDb above its neighbourhood (rule 2b)
%   Per trial and template (major and minor axis, reported as their mean unless stated):
%     nH, residual SD, PWLW alpha (1 to min(20, fs/2 - 1) Hz), fairShare common length (endToEndSegment_v004) as a share of
%     the residual, and the harmonic profile of the residual: for each k = 1..floor(fc/f0), the power in a band of +/-0.1 f0
%     around k*f0 against the median power in the flanks (0.25 to 0.45 f0 either side), in dB, from a Hann-tapered periodogram
%     zero-padded to 4x (resolution well below f0). excessMedDb = median over k; nOver3dB = harmonics standing out by > 3 dB
%     (leftover movement); nUnderMinus3dB = harmonics more than 3 dB BELOW their flanks (notches: noise removed with movement).
%   v002 (same day, after a Cook CTRL smoke run of v001): the dB counts flag narrow features but not how much variance they
%   carry, so v002 adds variance shares of the residual spectrum, for position, velocity ((2 pi f)^2) and acceleration
%   ((2 pi f)^4) weighting, below fc: leftover* = power standing above the flank level in the +/-0.1 f0 bands around k*f0
%   (movement left in), notch* = power missing below the flank level in the same bands (noise removed with the movement).
%   Each as a share of the residual's total power below fc under the same weighting.
%   Writes results/templateComparison/templateComparison_v002.csv (one row per trial x template) and _summary_v002.csv.
%   Deterministic (no random draws).
arguments
    opts.Arms (1,:) string = ["FraserN48", "Hickman PLAC", "Hickman HALO", "Cook CTRL", "Cook ASD"]
    opts.ExcessDb (1,1) double = 3
    opts.Fc (1,1) double = 10
end
here = fileparts(mfilename('fullpath'));  addpath(here, fullfile(here, 'functions'));
od = fullfile(here, '..', 'results', 'templateComparison');  if ~isfolder(od), mkdir(od); end
ensurePool_v001();  wait(parfevalOnAll(gcp, @addpath, 0, here, fullfile(here, 'functions')));
T = table();  exDb = opts.ExcessDb;  fc = opts.Fc;
for arm = opts.Arms
    switch arm
        case "FraserN48",    FS = 240; tr = importDB_fraser_v002(Verbose=false);
        case "Hickman PLAC", FS = 133; tr = importDB_hickman_v003(Study=2, Group="PLAC", Verbose=false);
        case "Hickman HALO", FS = 133; tr = importDB_hickman_v003(Study=2, Group="HALO", Verbose=false);
        case "Cook CTRL",    FS = 133; tr = importDB_cook_v002(Group="CTRL", Tasks=7, Verbose=false);
        case "Cook ASD",     FS = 133; tr = importDB_cook_v002(Group="ASD", Tasks=7, Verbose=false);
    end
    n = numel(tr);  rows = cell(n, 1);  t0 = tic;
    parfor t = 1:n
        rows{t} = oneTrial(tr(t), FS, fc, exDb);
    end
    Ta = vertcat(rows{:});  Ta.arm = repmat(arm, height(Ta), 1);  T = [T; Ta]; %#ok<AGROW>
    fprintf('%s: %d trials in %.1f min\n', arm, n, toc(t0) / 60);
end
writetable(T, fullfile(od, 'templateComparison_v002.csv'));
ok = T.status == "ok";
S = groupsummary(T(ok, :), ["arm", "template"], {@(v) median(v, 'omitnan'), @(v) prctile(v, 10), @(v) prctile(v, 90)}, ...
    ["nH", "sdRes", "alpha", "Lfrac", "excessMedDb", "nOver3dB", "nUnderMinus3dB", "nHarm", "leftPos", "leftVel", "leftAcc", "notchPos", "notchVel", "notchAcc"]);
S.Properties.VariableNames = replace(replace(replace(S.Properties.VariableNames, "fun1_", "med_"), "fun2_", "p10_"), "fun3_", "p90_");
m = groupsummary(T(ok, :), ["arm", "template"], "mean", ["anyOver3dB", "anyUnderMinus3dB"]);
S = join(S, removevars(m, "GroupCount"), 'Keys', ["arm", "template"]);
writetable(S, fullfile(od, 'templateComparison_summary_v002.csv'));
disp(S(:, ["arm", "template", "GroupCount", "med_nH", "med_sdRes", "med_alpha", "med_Lfrac", "med_excessMedDb", "med_nOver3dB", "med_nUnderMinus3dB", "med_leftPos", "med_leftVel", "med_leftAcc", "med_notchPos", "med_notchVel", "med_notchAcc"]));
fprintf('FAIL rows: %d of %d\n', nnz(~ok), height(T));
end

function R = oneTrial(t, FS, fc, exDb)
x = double(t.x(:));  y = double(t.y(:));
names = ["prod4", "end4", "endFc", "endSel"];  R = table();
try
    f0 = estimateF0_b(x, y, FS);
    xc = x - mean(x);  yc = y - mean(y);  [V, D] = eig(cov(xc, yc));  [~, o] = sort(diag(D), 'descend');  V = V(:, o);  th = atan2(V(2,1), V(1,1));
    kMax = max(1, floor(fc / f0));
    for tpl = names
        switch tpl
            case "prod4",  nH = 4;     [xR, yR, tailN] = tpl_b(x, y, FS, f0, nH, false);  xR = xR(1:end - tailN);  yR = yR(1:end - tailN);
            case "end4",   nH = 4;     [xR, yR] = tpl_b(x, y, FS, f0, nH, true);
            case "endFc",  nH = kMax;  [xR, yR] = tpl_b(x, y, FS, f0, nH, true);
            case "endSel"
                for nH = 1:kMax
                    [xR, yR] = tpl_b(x, y, FS, f0, nH, true);
                    [exA] = harmProfile([xR*cos(th) + yR*sin(th), -xR*sin(th) + yR*cos(th)], FS, f0, kMax);
                    if nH == kMax || all(max(exA(nH + 1:end, :), [], 2) <= exDb), break, end
                end
        end
        rMaj = xR*cos(th) + yR*sin(th);  rMin = -xR*sin(th) + yR*cos(th);  M = numel(rMaj);
        [ex, sh] = harmProfile([rMaj, rMin], FS, f0, kMax, fc);   % kMax x 2, dB; variance shares below fc
        a = NaN(1, 2);  try, a = [alphaPWLW_v001(rMaj, FS, [1 min(20, FS/2 - 1)]), alphaPWLW_v001(rMin, FS, [1 min(20, FS/2 - 1)])]; catch, end
        Lf = NaN;  try, [~, ~, inf4] = endToEndSegment_v004([rMaj, rMin], Tol=1, Mode="commonLength");  Lf = inf4.Lfrac; catch, end
        exM = mean(ex, 2);
        R = [R; table(string(t.trialID), tpl, f0, nH, kMax, M, mean([std(rMaj), std(rMin)]), mean(a, 'omitnan'), Lf, ...
            median(exM), nnz(exM > 3), nnz(exM < -3), kMax, any(exM > 3), any(exM < -3), sh(1), sh(2), sh(3), sh(4), sh(5), sh(6), "ok", ...
            'VariableNames', {'trialID','template','f0','nH','kMax','M','sdRes','alpha','Lfrac','excessMedDb','nOver3dB','nUnderMinus3dB','nHarm','anyOver3dB','anyUnderMinus3dB','leftPos','leftVel','leftAcc','notchPos','notchVel','notchAcc','status'})]; %#ok<AGROW>
    end
catch ME
    R = table(string(t.trialID), "all", NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, false, false, NaN, NaN, NaN, NaN, NaN, NaN, "FAIL: " + string(ME.message), ...
        'VariableNames', {'trialID','template','f0','nH','kMax','M','sdRes','alpha','Lfrac','excessMedDb','nOver3dB','nUnderMinus3dB','nHarm','anyOver3dB','anyUnderMinus3dB','leftPos','leftVel','leftAcc','notchPos','notchVel','notchAcc','status'});
end
end

function [ex, sh] = harmProfile(Rm, FS, f0, kMax, fc)
% power in +/-0.1 f0 around k*f0 against the median of the flanks (0.25-0.45 f0 either side), dB, per column
[M, C] = size(Rm);  nfft = 2^nextpow2(4 * M);  w = 0.5 * (1 - cos(2*pi*(0:M-1)'/(M-1)));
P = abs(fft((Rm - mean(Rm)) .* w, nfft)).^2;  f = (0:nfft-1)' * FS / nfft;  ex = NaN(kMax, C);
if nargin < 5, fc = FS / 2; end
inFc = f > 0 & f <= fc;  wt = [ones(size(f)), (2*pi*f).^2, (2*pi*f).^4];   % position, velocity, acceleration weighting
tot = zeros(3, C);  lo = zeros(3, C);  no = zeros(3, C);
for q = 1:3, tot(q, :) = sum(P(inFc, :) .* wt(inFc, q), 1); end
for k = 1:kMax
    c = k * f0;  inB = abs(f - c) <= 0.1 * f0 & f < FS/2;  fl = abs(f - c) >= 0.25 * f0 & abs(f - c) <= 0.45 * f0 & f < FS/2 & f > 0;
    if ~any(inB) || nnz(fl) < 3, continue, end
    lvl = median(P(fl, :), 1);
    ex(k, :) = 10 * log10(mean(P(inB, :), 1) ./ lvl);
    if c <= fc
        d = P(inB, :) - lvl;
        for q = 1:3, lo(q, :) = lo(q, :) + sum(max(d, 0) .* wt(inB, q), 1);  no(q, :) = no(q, :) + sum(max(-d, 0) .* wt(inB, q), 1); end
    end
end
sh = [mean(lo ./ tot, 2)', mean(no ./ tot, 2)'];   % [leftPos leftVel leftAcc notchPos notchVel notchAcc], mean over axes
end

% ---- runLoopClosureFftnoise_v024 template code, generalised to nH and returning the uncovered tail ----
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

function [resX, resY, tailN] = tpl_b(x, y, fs, f0, nH, endAnchored)
    N   = numel(x);  nCW = 4;  win = round(nCW/f0*fs);  hop = round(win/2);  win = min(win, N);  win = max(win, round(2/f0*fs));
    resX = zeros(N,1); resY = zeros(N,1); ww = zeros(N,1);
    starts = 1:hop:(N - win + 1);  if win > N, starts = []; end
    if endAnchored && ~isempty(starts) && starts(end) + win - 1 < N, starts(end+1) = N - win + 1; end
    tw = (0:win-1)'/fs;  han = 0.5*(1 - cos(2*pi*(0:win-1)'/(win-1)));
    D = [ones(win,1), tw];  for h = 1:nH, D = [D, cos(2*pi*h*f0*tw), sin(2*pi*h*f0*tw)]; end %#ok<AGROW>
    Dw = D .* han;
    for s = starts
        e = s + win - 1;  bX = Dw \ (x(s:e) .* han);  bY = Dw \ (y(s:e) .* han);
        resX(s:e) = resX(s:e) + (x(s:e) - D*bX) .* han;  resY(s:e) = resY(s:e) + (y(s:e) - D*bY) .* han;  ww(s:e) = ww(s:e) + han;
    end
    ok = ww > 0;  resX(ok) = resX(ok) ./ ww(ok);  resY(ok) = resY(ok) ./ ww(ok);
    cl = max(round(win/4), 1);  cs = cl;  ce = min(max(N - cl, cs + 100), N);
    lastEnd = 0;  if ~isempty(starts), lastEnd = starts(end) + win - 1; end
    tailN = max(0, ce - lastEnd);
    resX = resX(cs:ce);  resY = resY(cs:ce);
end
