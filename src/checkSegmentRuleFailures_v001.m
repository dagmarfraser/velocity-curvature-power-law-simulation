function T = checkSegmentRuleFailures_v001(opts)
%CHECKSEGMENTRULEFAILURES_V001 How many trials lose their end-matched noise under each segment rule and template? (2026-11-15)
%   After the template tail was excluded (v024) about a quarter of Hickman trials had no common length of at least 280 samples
%   (the runner minimum, 4*(EDGE_CLIP + 50)). Per trial, all trials of seven arms, for two templates (production 4 harmonics
%   with the uncovered tail cut, as v024/v025; end-anchored 4 harmonics) and three rules (Tol 1 = endToEndSegment_v004;
%   Budget 1% and 2% = endToEndSegment_v005): the common length L, L as a share of the residual, and the status
%   (ok, below 280, or no segment). v005 with no budget is first checked to equal v004 on every trial.
%   Writes results/templateComparison/segmentRuleFailures_v001.csv and prints the failure table.
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
od = fullfile(here, '..', 'results', 'templateComparison');  if ~isfolder(od), mkdir(od); end
writetable(T, fullfile(od, 'segmentRuleFailures_v001.csv'));
if any(T.v005eqV004 == 0), error('segRule:equivalence', 'endToEndSegment_v005 without a budget differs from v004 on %d trials', nnz(T.v005eqV004 == 0)); end
G = groupsummary(T(T.residualOk, :), ["arm", "template", "rule"], {@(v) mean(v == "ok"), @(v) mean(v == "short"), @(v) mean(v == "none")}, "status");
G.Properties.VariableNames = replace(G.Properties.VariableNames, ["fun1_status", "fun2_status", "fun3_status"], ["shareOk", "shareShort", "shareNone"]);
L = groupsummary(T(T.residualOk & T.status == "ok", :), ["arm", "template", "rule"], @(v) median(v), "Lfrac");
G = join(G, L(:, ["arm", "template", "rule", "fun1_Lfrac"]), 'Keys', ["arm", "template", "rule"]);
writetable(G, fullfile(od, 'segmentRuleFailures_summary_v001.csv'));  disp(G);
fprintf('v005 (no budget) = v004 on all %d trial x template rows\n', nnz(T.rule == "tol1"));
end

function R = oneTrial(t, FS)
x = double(t.x(:));  y = double(t.y(:));  R = table();
f0 = estimateF0_b(x, y, FS);
xc = x - mean(x);  yc = y - mean(y);  [V, D] = eig(cov(xc, yc));  [~, o] = sort(diag(D), 'descend');  V = V(:, o);  th = atan2(V(2,1), V(1,1));
for tpl = ["prod4", "end4"]
    [xR, yR, tailN] = tpl_b(x, y, FS, f0, 4, tpl == "end4");
    if tpl == "prod4", xR = xR(1:end - tailN);  yR = yR(1:end - tailN); end
    Rm = [xR*cos(th) + yR*sin(th), -xR*sin(th) + yR*cos(th)];  M = size(Rm, 1);  resOk = M >= 280;
    eq = NaN;
    for rule = ["tol1", "budget1", "budget2"]
        st = "ok";  Lx = NaN;
        try
            switch rule
                case "tol1"
                    [a4, b4, i4] = endToEndSegment_v004(Rm, Tol=1, Mode="commonLength");  [a5, b5] = endToEndSegment_v005(Rm, Tol=1);
                    eq = isequal([a4 b4], [a5 b5]);  Lx = i4.L;
                case "budget1", [~, ~, i5] = endToEndSegment_v005(Rm, Budget=0.01);  Lx = i5.L;
                case "budget2", [~, ~, i5] = endToEndSegment_v005(Rm, Budget=0.02);  Lx = i5.L;
            end
            if Lx < 280, st = "short"; end
        catch
            st = "none";
            if rule == "tol1"
                e4 = "";  e5 = "";
                try, endToEndSegment_v004(Rm, Tol=1, Mode="commonLength"); catch ME4, e4 = ME4.identifier; end
                try, endToEndSegment_v005(Rm, Tol=1); catch ME5, e5 = ME5.identifier; end
                eq = strcmp(e4, e5);
            end
        end
        R = [R; table(string(t.trialID), tpl, rule, M, resOk, Lx, Lx / M, st, double(eq), ...
            'VariableNames', {'trialID', 'template', 'rule', 'M', 'residualOk', 'L', 'Lfrac', 'status', 'v005eqV004'})]; %#ok<AGROW>
    end
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
