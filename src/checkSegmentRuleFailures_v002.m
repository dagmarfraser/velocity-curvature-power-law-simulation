function T = checkSegmentRuleFailures_v002(opts)
%CHECKSEGMENTRULEFAILURES_V002 checkSegmentRuleFailures_v001 for the end-anchored template, adding the detectability rule.
%   Per trial (all trials of seven arms), end-anchored template with 4 harmonics, three rules via endToEndSegment_v006:
%   Budget 1%, Budget 2% and "detect" (parameter-free). Also checks v006 equals v005 for the budget rule on every trial.
%   Reports per arm: share ok / below 280 / no segment, median L/M, and for "detect" the implied energy share of the wrap
%   (cJump/(L-1), max over axes; median and 90th percentile).
%   Writes results/templateComparison/segmentRuleFailures_v002.csv and _summary_v002.csv.
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
end
od = fullfile(here, '..', 'results', 'templateComparison');  writetable(T, fullfile(od, 'segmentRuleFailures_v002.csv'));
if any(T.eqV005 == 0), error('segRule2:equivalence', 'v006 budget differs from v005 on %d rows', nnz(T.eqV005 == 0)); end
G = groupsummary(T(T.residualOk, :), ["arm", "rule"], {@(v) mean(v == "ok"), @(v) mean(v == "short"), @(v) mean(v == "none")}, "status");
G.Properties.VariableNames = replace(G.Properties.VariableNames, ["fun1_status", "fun2_status", "fun3_status"], ["shareOk", "shareShort", "shareNone"]);
L = groupsummary(T(T.residualOk & T.status == "ok", :), ["arm", "rule"], {@(v) median(v), @(v) median(v), @(v) prctile(v, 90)}, ["Lfrac", "wrapShare", "wrapShare"]);
G = join(G, L(:, [1 2 4 6 8]), 'Keys', ["arm", "rule"]);  G.Properties.VariableNames(end-2:end) = {'medLfrac', 'medWrapShare', 'p90WrapShare'};
writetable(G, fullfile(od, 'segmentRuleFailures_summary_v002.csv'));  disp(G);
fprintf('v006 (budget) = v005 on all %d rows\n', nnz(T.rule == "budget1"));
end

function R = oneTrial(t, FS)
x = double(t.x(:));  y = double(t.y(:));  R = table();
f0 = estimateF0_b(x, y, FS);
xc = x - mean(x);  yc = y - mean(y);  [V, D] = eig(cov(xc, yc));  [~, o] = sort(diag(D), 'descend');  V = V(:, o);  th = atan2(V(2,1), V(1,1));
[xR, yR] = tpl_b(x, y, FS, f0, 4, true);
Rm = [xR*cos(th) + yR*sin(th), -xR*sin(th) + yR*cos(th)];  M = size(Rm, 1);
for rule = ["budget1", "budget2", "detect"]
    st = "ok";  Lx = NaN;  ws = NaN;  eq = NaN;
    try
        switch rule
            case "budget1"
                [a6, b6, inf6] = endToEndSegment_v006(Rm, Budget=0.01);  [a5, b5] = endToEndSegment_v005(Rm, Budget=0.01);  eq = isequal([a6 b6], [a5 b5]);
            case "budget2", [~, ~, inf6] = endToEndSegment_v006(Rm, Budget=0.02);
            case "detect",  [~, ~, inf6] = endToEndSegment_v006(Rm, Rule="detect");
        end
        Lx = inf6.L;  ws = max(inf6.cJump) / (Lx - 1);
        if Lx < 280, st = "short"; end
    catch
        st = "none";
        if rule == "budget1"
            e5 = "";  e6 = "";
            try, endToEndSegment_v005(Rm, Budget=0.01); catch ME5, e5 = ME5.identifier; end
            try, endToEndSegment_v006(Rm, Budget=0.01); catch ME6, e6 = ME6.identifier; end
            eq = strcmp(e5, e6);
        end
    end
    R = [R; table(string(t.trialID), rule, M, M >= 280, Lx, Lx / M, ws, st, double(eq), ...
        'VariableNames', {'trialID', 'rule', 'M', 'residualOk', 'L', 'Lfrac', 'wrapShare', 'status', 'eqV005'})]; %#ok<AGROW>
end
end

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

function [resX, resY] = tpl_b(x, y, fs, f0, nH, endAnchored)
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
    resX = resX(cs:ce);  resY = resY(cs:ce);
end
