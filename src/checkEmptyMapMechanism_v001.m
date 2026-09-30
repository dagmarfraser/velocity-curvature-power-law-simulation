%% checkEmptyMapMechanism_v001.m
% Follow-up to checkEmptyMaps_v001, which found that the 84 in-domain trials without a forward map
% (15.7% of the failing trial-cells in the in-domain FAIL cells) are 29 SHORT trials (template-
% subtracted length M < 280 samples) and 55 IRA-FAIL trials (per-axis IRASA exponent not finite),
% the latter all with an analysed stretch of at most one orbit. This script:
%   A. shows the code path. templateSubtract_local (runner v013 L1000-1035, copied here verbatim
%      and checked against the runner's text) fits windows of at least two cycles
%      (win >= round(2/f0*fs)). A trial shorter than that fits no window, its residual is
%      identically zero, and iraAlphaSigma_v001 finds no positive PSD in 1-20 Hz and returns NaN
%      (L60-65) without an error, so the runner's warning (L840-842, on a throw only) never fires.
%      Demonstrated on a synthetic ellipse at N = W2 - 1, W2 and W4 samples.
%   B. tests that account on every real trial with M >= 280: from its saved M, f0 and fs, which
%      branch of templateSubtract_local can produce M (M >= 280 fixes M = N - 2*cl + 1).
%      Errors unless every IRA-FAIL trial fits the zero-residual branch and every mapped trial
%      a fitting branch.
%   C. checks which downstream analyses already exclude these trials: the tempo decomposition's
%      filter (analyseTempoDecomposition_v002 L69) and anything that needs beta_gen*.
%   D. recomputes the Part 2 failure composition (checkPart2VerdictBasis_v001, section 1b) with
%      the trials without a map removed. Anchor: V0 reproduces that script's saved per-cell table.
% Reads:  src/loopClosureResults_<dataset>_all_shaped_xu_v015.mat; results/checkEmptyMaps_v001.mat,
%         results/checkPart2VerdictBasis_v001.mat, results/tempoDecomposition_v002.mat;
%         src/runLoopClosureFftnoise_v013.m (text only, for the verbatim check)
% Writes: results/checkEmptyMapMechanism_v001.mat
% USAGE:  from the project root: checkEmptyMapMechanism_v001
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT     = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(ROOT, "src", "functions")));
DATASETS = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Dhieb" "Zarandi"];
IN_DOM   = [true true true true true false false];                    % #242
FS       = [240 133 133 133 133 100 100];                               % FS_DS, runner v013 L185-218
PL       = ["BWFD-OLS" "SG-OLS" "BWFD-LMLS" "SG-LMLS" "BWFD-IRLS" "SG-IRLS"];
M_FLOOR  = 280;  SLOW_HZ = 0.35;  PUB_NOMAP = 15.7;                     % runner L806; checkPart2VerdictBasis_v001
IRA      = {1.0, 1.1:0.05:1.9};                                         % runner L834-836 (fHigh = min(20, fs/2 - 1))
RUNNER   = fullfile(ROOT, "src", "runLoopClosureFftnoise_v013.m");
EMP_MAT  = fullfile(ROOT, "results", "checkEmptyMaps_v001.mat");
P2_MAT   = fullfile(ROOT, "results", "checkPart2VerdictBasis_v001.mat");
TD_MAT   = fullfile(ROOT, "results", "tempoDecomposition_v002.mat");
OUT_MAT  = fullfile(ROOT, "results", "checkEmptyMapMechanism_v001.mat");
for f = [RUNNER EMP_MAT P2_MAT TD_MAT], if ~isfile(f), error("emptyMech:input", "%s", "FAILED PATH: " + f); end, end

%% A. The code path
SIG = "function [resX, resY, fitX, fitY] = templateSubtract_local(x, y, fs, f0, nH, nCW)";
if ~isequal(fnBlock_local(RUNNER, SIG), fnBlock_local(string(mfilename("fullpath")) + ".m", SIG))
    error("emptyMech:copy", "%s", "templateSubtract_local here differs from runLoopClosureFftnoise_v013's");
end
fprintf("VERBATIM CHECK passed: templateSubtract_local is runLoopClosureFftnoise_v013 L1000-1035.\n");
fsD = 240;  f0D = 0.25;  W2 = round(2/f0D*fsD);  W4 = round(4/f0D*fsD);  rng(1729);
demo = table();
for N = [W2 - 1, W2, W4]
    t = (0:N - 1)' / fsD;
    [rx, ry] = templateSubtract_local(100*cos(2*pi*f0D*t) + randn(N, 1), 45*sin(2*pi*f0D*t) + randn(N, 1), fsD, f0D, 4, 4);
    demo = [demo; table(N, N*f0D/fsD, numel(rx), all(rx == 0 & ry == 0), iraAlphaSigma_v001(rx, fsD, IRA{1}, 20, IRA{2}), ...
        'VariableNames', ["N" "orbits" "M" "residualAllZero" "alphaIRASA"])]; %#ok<AGROW>
end
disp(demo)
if ~(demo.residualAllZero(1) && isnan(demo.alphaIRASA(1)) && ~any(demo.residualAllZero(2:3)) && all(isfinite(demo.alphaIRASA(2:3))))
    error("emptyMech:demo", "%s", "synthetic demonstration does not show zero residual and NaN alpha below two cycles only");
end
fprintf("DEMO passed: below %d samples (two cycles at %.2f Hz) the residual is identically zero and IRASA returns NaN; at %d it is not.\n\n", W2, f0D, W2);

%% Per trial (v015)
E = load(EMP_MAT, "T").T;  P = load(P2_MAT, "comp").comp;  D = load(TD_MAT, "T", "PIPES");
B = table();  C = table();  TC = table();
for k = 1:numel(DATASETS)
    f = fullfile(ROOT, "src", "loopClosureResults_" + DATASETS(k) + "_all_shaped_xu_v015.mat");
    if ~isfile(f), error("emptyMech:input", "%s", "FAILED PATH: " + f); end
    R = load(f, "results").results;  n = numel(R);  e = E(E.dataset == DATASETS(k), :);
    M = arrayfun(@(r) double(r.M), R(:));  f0 = arrayfun(@(r) double(r.f0), R(:));
    if height(e) ~= n || ~isequal(e.M, M), error("emptyMech:join", "%s", DATASETS(k) + ": checkEmptyMaps_v001 rows do not match the corpus"); end
    noMap = e.cause ~= "mapped";
    [zOK, lOK] = branch_local(M, f0, FS(k));
    B = [B; table(repmat(DATASETS(k), n, 1), e.cause, M, f0, M .* f0 / FS(k), zOK, lOK, (M + 2*max(round(round(2 ./ f0 * FS(k)) / 4), 1) - 1) .* f0 / FS(k), ...
        'VariableNames', ["dataset" "cause" "M" "f0" "orbitsM" "zeroBranch" "loopBranch" "rawOrbitsIfZero"])]; %#ok<AGROW>
    a = arrayfun(@(r) double(r.a_mm), R(:));  s = arrayfun(@(r) double(r.sigmaMM), R(:));  al = arrayfun(@(r) double(r.alphaIRA), R(:));
    keep69 = isfinite(f0) & f0 > 0 & isfinite(s) & isfinite(a) & a > 0 & isfinite(al);   % analyseTempoDecomposition_v002 L69
    gs = cell2mat(arrayfun(@(r) double(r.betaGenStar(:))', R(:), "UniformOutput", false));
    gm = arrayfun(@(r) double(r.betaGenStarMed), R(:));
    nTD = NaN;  if any(D.T.dataset == DATASETS(k)), nTD = nnz(D.T.dataset == DATASETS(k) & D.T.pipeline == D.PIPES(1)); end
    C = [C; table(DATASETS(k), n, nnz(noMap), nnz(~keep69), isequal(~keep69, noMap), nTD, nnz(isfinite(gs(noMap, :))) + nnz(isfinite(gm(noMap))), ...
        'VariableNames', ["dataset" "n" "nNoMap" "nDroppedByL69" "sameTrials" "nInSavedDecomposition" "nFiniteGenStarNoMap"])]; %#ok<AGROW>
    St = arrayfun(@(r) string(r.invertStatus(:))', R(:), "UniformOutput", false);  st = vertcat(St{:});
    pos = repmat("nomap", n, 6);                                           % as checkPart2VerdictBasis_v001 L54-63
    for i = 1:n
        Mc = R(i).betaRecCurveMed;
        if isempty(Mc), continue; end
        for q = 1:6
            c = Mc(q, :);  b = R(i).betaObs(q);
            if all(isnan(c)) || ~isfinite(b), pos(i, q) = "nan";
            elseif b < min(c), pos(i, q) = "below";  elseif b > max(c), pos(i, q) = "above";  else, pos(i, q) = "within"; end
        end
    end
    TC = [TC; table(repelem(DATASETS(k), 6*n, 1), repelem(IN_DOM(k), 6*n, 1), repelem(f0, 6), repmat(PL(:), n, 1), ...
        reshape(st', [], 1), reshape(pos', [], 1), repelem(noMap, 6), ...
        'VariableNames', ["dataset" "inDomain" "f0" "pipeline" "status" "obsPos" "noMap"])]; %#ok<AGROW>
end

%% B. Which branch produced each real trial's M
b = B(B.M >= M_FLOOR & isfinite(B.f0) & B.f0 > 0, :);
if any(~b.zeroBranch & ~b.loopBranch), error("emptyMech:branch", "%s", sprintf("%d trials fit no branch", nnz(~b.zeroBranch & ~b.loopBranch))); end
badI = b.cause == "iraFail" & ~b.zeroBranch;  badM = b.cause == "mapped" & ~b.loopBranch;
if any(badI) || any(badM)
    error("emptyMech:branchCause", "%s", sprintf("%d IRA-FAIL trials cannot be zero-residual; %d mapped trials cannot be windowed", nnz(badI), nnz(badM)));
end
iF = b(b.cause == "iraFail", :);
fprintf("BRANCH CHECK passed on %d trials with M >= %d: all %d IRA-FAIL trials fit only the zero-residual branch (recording shorter than two cycles);\n", ...
    height(b), M_FLOOR, height(iF));
fprintf("  every mapped trial fits a windowed branch; %d trials fit both (boundary rounding). IRA-FAIL recordings span %.2f-%.2f orbits (median %.2f).\n", ...
    nnz(b.zeroBranch & b.loopBranch), min(iF.rawOrbitsIfZero), max(iF.rawOrbitsIfZero), median(iF.rawOrbitsIfZero));
if any(iF.zeroBranch & iF.loopBranch), warning("emptyMech:ambiguous", "%s", sprintf("%d IRA-FAIL trials also fit a windowed branch", nnz(iF.zeroBranch & iF.loopBranch))); end

%% C. Downstream exclusion
disp(C)
chk = C(~isnan(C.nInSavedDecomposition), :);
if ~all(chk.sameTrials) || any(chk.n - chk.nNoMap ~= chk.nInSavedDecomposition) || any(C.nFiniteGenStarNoMap > 0)
    error("emptyMech:downstream", "%s", "the decomposition filter or beta_gen* does not exclude exactly the trials without a map");
end
fprintf("DOWNSTREAM CHECK passed: analyseTempoDecomposition_v002's L69 filter drops exactly the trials without a map (its saved table\n");
fprintf("  holds n - nNoMap trials per dataset), and no trial without a map has a beta_gen*. The tempo decomposition, the selection GLMM\n");
fprintf("  and every beta_gen* summary already exclude them; the coverage verdicts count them as failures.\n\n");

%% D. Part 2 failure composition without the trials that have no map (v015)
Vt = TC(TC.status ~= "no_beta_obs", :);  Vt.rise = Vt.status == "rise";
G = groupsummary(Vt, ["dataset" "pipeline"], "mean", "rise");
G.verdict = arrayfun(@verdict_local, G.mean_rise);
G.inDomain = ismember(G.dataset, DATASETS(IN_DOM));
F = G(G.inDomain & G.verdict == "FAIL", :);
comp = table();  pf = table();
for r = 1:height(F)
    x = TC(TC.dataset == F.dataset(r) & TC.pipeline == F.pipeline(r) & TC.status ~= "no_beta_obs" & isfinite(TC.f0), :);
    x1 = x(~x.noMap, :);
    [m0, m1] = deal(cellMetrics_local(x, SLOW_HZ), cellMetrics_local(x1, SLOW_HZ));
    comp = [comp; table(F.dataset(r), F.pipeline(r), m0, m1, nnz(x.noMap & x.status ~= "rise"), 'VariableNames', ["dataset" "pipeline" "V0" "V1" "nFailNoMap"])]; %#ok<AGROW>
    pf = [pf; x(x.status ~= "rise", :)]; %#ok<AGROW>
end
[ok, loc] = ismember(comp.dataset + "|" + comp.pipeline, P.dataset + "|" + P.pipeline);
if ~all(ok) || height(comp) ~= height(P), error("emptyMech:anchorCells", "%s", "FAIL cells differ from checkPart2VerdictBasis_v001's"); end
Pm = P(loc, :);
ref = [Pm.nValid Pm.nFail Pm.f0MedFail Pm.f0MedRise Pm.pctFailSlow Pm.pctFailOffMap Pm.pctFailFold];
if any(abs(comp.V0 - ref) > 1e-3, "all")
    error("emptyMech:anchorComp", "%s", "V0 does not reproduce checkPart2VerdictBasis_v001's per-cell table");
end
nb = pf.status == "neither";
sh0 = 100 * [mean(nb & pf.obsPos == "above"), mean(nb & pf.obsPos == "below"), mean(nb & ismember(pf.obsPos, ["nan" "nomap"])), ...
    mean(nb & pf.obsPos == "within"), mean(ismember(pf.status, ["ambiguous" "desc"]))];
if abs(sh0(3) - PUB_NOMAP) > 0.05, error("emptyMech:anchorPooled", "%s", sprintf("no-map share %.2f%%, expected %.1f%%", sh0(3), PUB_NOMAP)); end
fprintf("REGRESSION ANCHOR passed: V0 reproduces checkPart2VerdictBasis_v001's %d FAIL cells and its %.1f%% no-map share.\n", height(comp), sh0(3));
p1 = pf(~pf.noMap, :);  nb1 = p1.status == "neither";
sh1 = 100 * [mean(nb1 & p1.obsPos == "above"), mean(nb1 & p1.obsPos == "below"), 0, mean(nb1 & p1.obsPos == "within"), ...
    mean(ismember(p1.status, ["ambiguous" "desc"]))];
fprintf("Pooled failing trial-cells, %% above the map's top / below its floor / no map / within, no monotonic run / fold branch:\n");
fprintf("  V0 all %d:           %5.1f / %5.1f / %5.1f / %4.1f / %5.1f\n", height(pf), sh0);
fprintf("  V1 with a map, %d:   %5.1f / %5.1f / %5.1f / %4.1f / %5.1f  (off the map %.1f%%)\n", height(p1), sh1, sh1(1) + sh1(2));
cols = ["nValid" "nFail" "f0MedFail" "f0MedRise" "pctFailSlow" "pctFailOffMap" "pctFailFold"];
fprintf("\nPer FAIL cell (V0 -> V1, trials without a map removed): failing trials, their median f0, %% of them below %.2f Hz\n", SLOW_HZ);
for r = 1:height(comp)
    fprintf("  %-13s %-9s fail %4d -> %4d (no map %2d); f0 fail %.3f -> %.3f Hz (rise %.3f); slow %5.1f%% -> %5.1f%%\n", comp.dataset(r), comp.pipeline(r), ...
        comp.V0(r, 2), comp.V1(r, 2), comp.nFailNoMap(r), comp.V0(r, 3), comp.V1(r, 3), comp.V0(r, 4), comp.V0(r, 5), comp.V1(r, 5));
end

runDate = string(datetime("now", "Format", "yyyy-MM-dd"));
save(OUT_MAT, "demo", "B", "C", "comp", "cols", "sh0", "sh1", "runDate");
fprintf("\nSaved: %s\n", OUT_MAT);

%% =========================================================================
function [zOK, lOK] = branch_local(M, f0, fs)
% Which branch of templateSubtract_local (nCW = 4) can give length M when M >= 280 (ce = N - cl).
    W4 = round(4 ./ f0 * fs);  W2 = round(2 ./ f0 * fs);
    zOK = M + 2*max(round(W2 / 4), 1) - 1 < W2;                            % N < W2: win = W2, no window fits
    lOK = M + 2*max(round(W4 / 4), 1) - 1 >= W4;                           % N >= W4: win = W4
    for i = find(~lOK & isfinite(W2))'                                     % W2 <= N < W4: win = N
        Nc = W2(i):W4(i) - 1;
        lOK(i) = any(Nc - 2*max(round(Nc / 4), 1) + 1 == M(i));
    end
end

function m = cellMetrics_local(x, slowHz)
% nValid, nFail, median f0 failing / rising, % failing below slowHz, % failing off the map, % on a fold
    bad = x(x.status ~= "rise", :);  good = x(x.status == "rise", :);
    m = [height(x), height(bad), median(bad.f0), median(good.f0), 100 * mean(bad.f0 < slowHz), ...
        100 * mean(bad.status == "neither" & ismember(bad.obsPos, ["below" "above"])), 100 * mean(ismember(bad.status, ["ambiguous" "desc"]))];
end

function v = verdict_local(c)
    if c >= 0.95, v = "PASS"; elseif c >= 0.90, v = "CONDITIONAL"; else, v = "FAIL"; end
end

function blk = fnBlock_local(file, sig)
% The lines of local function sig in file, from its signature to the first unindented "end".
    L = splitlines(string(fileread(file)));
    i0 = find(strtrim(L) == sig, 1);
    if isempty(i0), error("emptyMech:sig", "%s", "signature not found in " + file); end
    i1 = i0 + find(deblank(L(i0 + 1:end)) == "end", 1);
    blk = strtrim(L(i0:i1));
end

function [resX, resY, fitX, fitY] = templateSubtract_local(x, y, fs, f0, nH, nCW)
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
    fitX = x(cs:ce) - resX(cs:ce);
    fitY = y(cs:ce) - resY(cs:ce);
    resX = resX(cs:ce);
    resY = resY(cs:ce);
end
