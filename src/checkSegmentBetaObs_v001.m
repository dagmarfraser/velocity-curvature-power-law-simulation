function [R, S] = checkSegmentBetaObs_v001(opts)
%CHECKSEGMENTBETAOBS_V001 Does cutting a real trial to its end-matched noise segment change beta_obs?
%   Dagmar (2026-11-14) chose to shorten the synthetic trajectory to the end-matched noise segment (fairShare rule, one
%   segment for both axes: endToEndSegment_v003 on [resid_maj, resid_min], Tol 1). The forward map is then built at length
%   L while beta_obs comes from the whole trial (length N). This script measures what length alone does to beta_obs: per
%   trial and pipeline, beta_obs on the whole trial against beta_obs on the trial cut to the samples its noise segment
%   covers (residual index i is trial index i + cs - 1, cs the template-subtraction clip). If the difference is small
%   against the MDC (0.03), the shorter synthetic is comparable; if not, beta_obs should come from the same segment.
%   Residuals: production template (endAnchored = false; the v015 reporting basis), nH = 4, nCW = 4, functions copied
%   verbatim from runLoopClosureFftnoise_v021 as in endMismatchSensitivity_v006. Pipelines as the runner: BWFD [2 10 1] and
%   SG [4 17] x OLS, LMLS, IRLS (regressDataEBR, limitBreak 0, OLS fit seeds LMLS/IRLS), EDGE_CLIP 20. Trials whose residual
%   is shorter than 4*(EDGE_CLIP + 50) are skipped as the runner skips them. Deterministic: no random draws.
%   Writes results/endMismatchSensitivity/segmentBetaObs_v001.csv (one row per trial x pipeline) and _summary_v001.csv.
arguments
    opts.Arms (1,:) string = ["FraserN48", "Hickman PLAC", "Hickman HALO", "Cook CTRL", "Cook ASD"]
    opts.MaxTrials (1,1) double {mustBePositive} = Inf
end
here = fileparts(mfilename('fullpath'));  addpath(here, fullfile(here, 'functions'));
od = fullfile(here, '..', 'results', 'endMismatchSensitivity');  if ~isfolder(od), mkdir(od); end
pool = ensurePool_v001();  wait(parfevalOnAll(pool, @addpath, 0, here, fullfile(here, 'functions')));
pl = ["BWFD-OLS", "SG-OLS", "BWFD-LMLS", "SG-LMLS", "BWFD-IRLS", "SG-IRLS"];   % runner order: (rIdx-1)*2 + dIdx
R = table();
for arm = opts.Arms
    switch arm
        case "FraserN48",    FS = 240; tr = importDB_fraser_v002(Verbose=false);
        case "Hickman PLAC", FS = 133; tr = importDB_hickman_v003(Study=2, Group="PLAC", Verbose=false);
        case "Hickman HALO", FS = 133; tr = importDB_hickman_v003(Study=2, Group="HALO", Verbose=false);
        case "Cook CTRL",    FS = 133; tr = importDB_cook_v002(Group="CTRL", Tasks=7, Verbose=false);
        case "Cook ASD",     FS = 133; tr = importDB_cook_v002(Group="ASD", Tasks=7, Verbose=false);
        otherwise, error('checkSegmentBetaObs:arm', 'Unknown arm "%s".', arm);
    end
    tr = tr(1:min(numel(tr), opts.MaxTrials));  nT = numel(tr);  out = cell(nT, 1);
    parfor t = 1:nT
        out{t} = oneTrial(tr(t), FS, pl);
    end
    Ra = vertcat(out{:});  Ra.arm = repmat(arm, height(Ra), 1);  R = [R; Ra]; %#ok<AGROW>
end
writetable(R, fullfile(od, "segmentBetaObs_v001.csv"));
ok = R.status == "ok" & isfinite(R.dBeta);
S = groupsummary(R(ok, :), ["arm", "pipeline"], {@(v) median(v), @(v) prctile(v, 25), @(v) prctile(v, 75), @(v) mean(abs(v) > 0.03), @(v) median(abs(v))}, "dBeta");
S.Properties.VariableNames = replace(S.Properties.VariableNames, ["fun1_dBeta", "fun2_dBeta", "fun3_dBeta", "fun4_dBeta", "fun5_dBeta"], ["medianDBeta", "q25", "q75", "shareAbsAboveMDC", "medianAbsDBeta"]);
L = groupsummary(R(ok & R.pipeline == pl(1), :), "arm", {@(v) median(v), @(v) prctile(v, 10), @(v) prctile(v, 90)}, "segFracOfTrial");
S = join(S, L(:, ["arm", "fun1_segFracOfTrial", "fun2_segFracOfTrial", "fun3_segFracOfTrial"]), 'Keys', "arm");
writetable(S, fullfile(od, "segmentBetaObs_summary_v001.csv"));
fprintf('%d trial x pipeline rows, %d not ok (status) or non-finite\n', height(R), nnz(~ok));
disp(groupcounts(R(R.pipeline == pl(1), :), ["arm", "status"]));  disp(S);
end

function T = oneTrial(tr, FS, pl)
x = double(tr.x(:));  y = double(tr.y(:));  N = numel(x);  EC = 20;
bFull = NaN(1, 6);  bSeg = NaN(1, 6);  st = "ok";  segFrac = NaN;  i0 = NaN;  i1 = NaN;
try
    f0 = estimateF0_local(x, y, FS);
    [xR, yR, ~, ~, cs] = templateSubtract_local(x, y, FS, f0, 4, 4, false);
    M = numel(xR);
    if M < 4 * (EC + 50), error('checkSegmentBetaObs:short', 'residual of %d samples is below the runner minimum %d', M, 4 * (EC + 50)); end
    xc = x - mean(x);  yc = y - mean(y);
    [V, D] = eig(cov(xc, yc));  [~, ord] = sort(diag(D), 'descend');  V = V(:, ord);  th = atan2(V(2,1), V(1,1));
    rMaj = xR*cos(th) + yR*sin(th);  rMin = -xR*sin(th) + yR*cos(th);
    [i0, i1] = endToEndSegment_v003([rMaj, rMin], Tol=1);
    a = i0 + cs - 1;  b = i1 + cs - 1;  segFrac = (b - a + 1) / N;
    bFull = betaObs_local(x, y, FS, EC);  bSeg = betaObs_local(x(a:b), y(a:b), FS, EC);
catch ME
    st = "FAIL: " + ME.message;
end
T = table(repmat(string(tr.trialID), 6, 1), pl(:), repmat(N, 6, 1), repmat(i0, 6, 1), repmat(i1, 6, 1), repmat(segFrac, 6, 1), ...
    bFull(:), bSeg(:), bSeg(:) - bFull(:), repmat(st, 6, 1), ...
    'VariableNames', {'trialID', 'pipeline', 'N', 'segI0', 'segI1', 'segFracOfTrial', 'betaObsFull', 'betaObsSeg', 'dBeta', 'status'});
end

function b6 = betaObs_local(x, y, FS, EC)
% the runner's beta_obs block (runLoopClosureFftnoise_v021 L944-966): 2 filters x 3 regressions, ppCM = (rIdx-1)*2 + dIdx
dc = {2, [2 10 1]; 6, [4 17]};  b6 = NaN(1, 6);
for dIdx = 1:2
    [dx, dy] = differentiateKinematicsEBR(x, y, dc{dIdx, 1}, dc{dIdx, 2}, FS);
    vx = dx(EC:end-EC, 2);  vy = dy(EC:end-EC, 2);  ax = dx(EC:end-EC, 3);  ay = dy(EC:end-EC, 3);
    sp = hypot(vx, vy);  kp = curvatureKinematicEBR(vx, vy, ax, ay);  lm = [1, -1/3];
    for rIdx = 1:3
        try
            [b, v] = regressDataEBR(sp, kp, rIdx + 2, lm, 0, 0);  b6((rIdx - 1) * 2 + dIdx) = b;
            if rIdx == 1, lm = [v, b]; end
        catch
        end
    end
end
end

% ---- copied verbatim from runLoopClosureFftnoise_v021.m, except that templateSubtract_local also returns cs ----
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

function [resX, resY, fitX, fitY, cs] = templateSubtract_local(x, y, fs, f0, nH, nCW, endAnchored)
    N   = numel(x);
    win = round(nCW/f0*fs);
    hop = round(win/2);
    win = min(win, N);
    win = max(win, round(2/f0*fs));
    resX = zeros(N,1); resY = zeros(N,1); ww = zeros(N,1);
    starts = 1:hop:(N - win + 1);   % the production windows (identical to: s = 1; while s + win - 1 <= N, s = s + hop)
    if win > N, starts = []; end
    if endAnchored && ~isempty(starts) && starts(end) + win - 1 < N, starts(end+1) = N - win + 1; end   % v019: END-ANCHORED final window; v020: only when Template = 'endAnchored'
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
    fitX = x(cs:ce) - resX(cs:ce);
    fitY = y(cs:ce) - resY(cs:ce);
    resX = resX(cs:ce);
    resY = resY(cs:ce);
end
