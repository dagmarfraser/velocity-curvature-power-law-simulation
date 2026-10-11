%% checkSelectionDefect_v002.m
% v001 (run 2026-10-11) measured how the invertibility gate selects on observables (Meng 2018,
% Ann. Appl. Stat. 12:685-726: mean(G | recorded) - mean(G) = rho * sqrt((1-f)/f) * sd(G)) and
% subsampled Fraser by SUBJECT to each smaller dataset's subject count. v002 keeps Part A unchanged
% (regression anchor against results/checkSelectionDefect_v001.mat) and adds:
%  Part B  Fraser subsampled to each smaller dataset's DESIGN: its subject count AND its own
%          trials-per-subject counts (mode "subjects+trials"), beside the subjects-only mode.
%          v001 matched subjects only, so Fraser kept about 1100 trials against Cook CTRL's 94.
%  Part C  per-pipeline pair differences in the selection shift (v001 tested the all-pipeline gate
%          only; FINDINGS_REFERENCE.md L825 says the bias is "symmetric within each pair").
%  Part D  the Cook ASD - CTRL tempo contrast on ALL trials (the basis of Finding #235's +0.50 log2,
%          analyseTempoDecomposition_v002 L144) against the same contrast on GATED trials (the sample
%          that carries beta_gen*): does the gate change the tempo gap that the beta_gen* contrast sees?
% R (all-pipeline) = isfinite(betaGenStarMed). R (pipeline p) = invertStatus(p) == "rise" over trials
% with a beta_obs (compare42CellVerdict_v002 verdictOneCell_local, L73-90). In-domain datasets only.
% Prior statements: FINDINGS_REFERENCE.md L825 (selection bias; symmetric within each pair) and
% Finding #142 (L3464; N_fin untraced). In this file's v015 input, Cook's N_fin is the finite
% betaGenStarMed count; Fraser and Hickman differ from the hub's N_fin by 2-6 trials (basis not traced).
% What this can and cannot show: selection on OBSERVABLES only; Meng's identity covers selection, not
% measurement bias; the permutation null ignores subject clustering; the constellation median is not a
% mean (Meng footnote 4). Part B models sampling of Fraser's own subjects and trials, not Fraser's
% different device, noise or paradigm, so "N alone does not explain it" is all it can say.
% Predictions, fixed 2026-10-11 before v002 was run (P1-P5 are v001's, already run and HELD):
%  P4  (both modes) the median of subsample medians stays within MDC/2.77 of the full-Fraser median.
%  P6  trial thinning widens the 95% interval: at the Cook CTRL design (19 subjects) the interval for
%      "subjects+trials" is at least 1.5 times as wide as for "subjects".
%  P7  Cook CTRL's and Cook ASD's own medians still fall below the 5th percentile of the thinned draws
%      (N at their design does not by itself produce them).
%  P8  not a prediction: v001's printed table suggested 3 of 12 per-pipeline pair differences exceed
%      MDC/2.77. The run shows 4 of 14 pair differences (12 single-pipeline plus 2 all-pipeline); all four
%      are single-pipeline: Cook SG-LMLS, BWFD-IRLS, SG-IRLS and Hickman BWFD-OLS. Reported with the flag per cell.
%  P9  among gated (all-pipeline) trials the Cook ASD - CTRL tempo contrast is larger than on all
%      trials (v001: gated trials are faster by 0.296 log2 units in ASD, 0.184 in CTRL).
% Anchors (error if broken): dataset medians span 0.272-0.328 (Finding #247); the registered rule gives
% 4/2/24 over the 30 cells; Meng's identity holds to 1e-9; Part A rho and shift reproduce v001's saved
% rows to 1e-10 (regression anchor).
% Fail loud: nothing is filled; shortfalls in trials per subject when thinning are counted and printed.
% Run from src/. Reads src/loopClosureResults_<dataset>_all_shaped_xu_v015.mat,
% results/checkSelectionDefect_v001.mat, results/tempoDecomposition_v002.mat (printed for reference).
% Writes results/checkSelectionDefect_v002.mat

%% CONFIG
ROOT      = fileparts(fileparts(mfilename("fullpath")));
DATASETS  = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO"];
FILE_FMT  = "loopClosureResults_%s_all_shaped_xu_v015.mat";
V1_MAT    = fullfile(ROOT, "results", "checkSelectionDefect_v001.mat");
TEMPO_MAT = fullfile(ROOT, "results", "tempoDecomposition_v002.mat");
OUT_MAT   = fullfile(ROOT, "results", "checkSelectionDefect_v002.mat");
B_PERM    = 5000;
B_SUB     = 2000;
RNG_SEED  = 1729;
ANCHOR    = [0.272 0.328];
VERDICTS  = [4 2 24];
PASS_COV  = 0.95;
COND_COV  = 0.90;
MDC       = 0.03;
G_NAMES   = ["betaObs" "log2 f0" "sigmaMM" "alphaIRA"];
PAIRS     = ["Cook_ASD" "Cook_CTRL"; "Hickman_HALO" "Hickman_PLAC"];

%% Setup
addpath(genpath(fullfile(ROOT, "src", "functions")));
if ~isfolder(fileparts(OUT_MAT)), error("selDef:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
if ~isfile(V1_MAT), error("selDef:v1", "%s", "Missing regression anchor: " + V1_MAT); end
if ~isfile(TEMPO_MAT), error("selDef:tempo", "%s", "Missing input: " + TEMPO_MAT); end
semAdequate = semAdequacyThreshold_v001;
rng(RNG_SEED, "twister");

%% Load
nD = numel(DATASETS);  data = repmat(struct(), 1, nD);
for di = 1:nD
    f = fullfile(ROOT, "src", sprintf(FILE_FMT, DATASETS(di)));
    if ~isfile(f), error("selDef:in", "%s", "Missing input: " + f); end
    S = load(f, "results", "pipelineLabels");  r = S.results(:);
    st = vertcat(r.invertStatus);  obs = vertcat(r.betaObs);
    if size(st, 2) ~= 6 || size(obs, 2) ~= 6 || size(st, 1) ~= numel(r), error("selDef:shape", "%s", "Unexpected shape in " + DATASETS(di)); end
    data(di).name = DATASETS(di);  data(di).labels = string(S.pipelineLabels(:))';
    data(di).subj = string({r.subjectID})';  data(di).nSub = numel(unique(data(di).subj));
    data(di).obs = obs;  data(di).obsMed = median(obs, 2, "omitnan");  data(di).status = st;
    data(di).genMed = [r.betaGenStarMed]';  data(di).gen = horzcat(r.betaGenStar)';
    data(di).f0 = [r.f0]';  data(di).sigma = [r.sigmaMM]';  data(di).alpha = [r.alphaIRA]';
end

%% Anchors
cMed = arrayfun(@(d) median(d.genMed, "omitnan"), data);
if round(min(cMed), 3) ~= ANCHOR(1) || round(max(cMed), 3) ~= ANCHOR(2)
    error("selDef:anchor", "Dataset medians of betaGenStarMed span %.4f-%.4f, not Finding #247's %.3f-%.3f.", min(cMed), max(cMed), ANCHOR(1), ANCHOR(2));
end
verd = strings(nD, 6);  cov = zeros(nD, 6);
for di = 1:nD, for p = 1:6, [verd(di, p), cov(di, p)] = covVerdict_local(data(di).status(:, p), PASS_COV, COND_COV); end, end
counts = [sum(verd(:) == "PASS") sum(verd(:) == "CONDITIONAL") sum(verd(:) == "FAIL")];
if ~isequal(counts, VERDICTS)
    error("selDef:verdicts", "Registered rule gives PASS/COND/FAIL = %s, not %s.", mat2str(counts), mat2str(VERDICTS));
end
fprintf("Anchors passed: medians %.4f-%.4f; verdicts PASS/COND/FAIL = %s over 30 in-domain cells\n", min(cMed), max(cMed), mat2str(counts));

%% Part A: Meng decomposition of the gate (unchanged from v001)
rowsA = table();
for di = 1:nD
    d = data(di);  R = isfinite(d.genMed);
    Gs = {d.obsMed, log2(d.f0), d.sigma, d.alpha};
    for k = 1:numel(G_NAMES)
        rowsA = [rowsA; [table(d.name, "all-pipeline", G_NAMES(k), 'VariableNames', ["dataset" "gate" "G"]), mengRow_local(R, Gs{k}, B_PERM)]]; %#ok<AGROW>
    end
    for p = 1:6
        ok = d.status(:, p) ~= "no_beta_obs";
        rowsA = [rowsA; [table(d.name, d.labels(p), "betaObs", 'VariableNames', ["dataset" "gate" "G"]), ...
            mengRow_local(d.status(ok, p) == "rise", d.obs(ok, p), B_PERM)]]; %#ok<AGROW>
    end
end
rowsA.shiftMDC = rowsA.shift / MDC;
V1 = load(V1_MAT, "rowsA");
if height(V1.rowsA) ~= height(rowsA) || any(V1.rowsA.dataset ~= rowsA.dataset) || any(V1.rowsA.gate ~= rowsA.gate) || any(V1.rowsA.G ~= rowsA.G)
    error("selDef:v1rows", "Part A rows do not line up with v001's saved rows.");
end
dRho = max(abs(V1.rowsA.rho - rowsA.rho));  dShift = max(abs(V1.rowsA.shift - rowsA.shift));
if dRho > 1e-10 || dShift > 1e-10, error("selDef:v1regress", "Part A differs from v001: max |d rho| %.2e, max |d shift| %.2e", dRho, dShift); end
fprintf("Regression anchor passed: Part A rho and shift reproduce v001 (max |d rho| %.1e, max |d shift| %.1e; max |d pPerm| %.1e)\n", dRho, dShift, max(abs(V1.rowsA.pPerm - rowsA.pPerm)));

%% Part B: Fraser subsampled to each smaller dataset's design
fr = data(1);  us = unique(fr.subj);  nF = numel(us);
frIdx = arrayfun(@(u) find(fr.subj == u), us, 'UniformOutput', false);
fprintf("\nFraser: %d subjects, trials per subject min %d, median %d, max %d\n", nF, min(cellfun(@numel, frIdx)), median(cellfun(@numel, frIdx)), max(cellfun(@numel, frIdx)));
MODES = ["subjects" "subjects+trials"];
rowsB = table();  covB = table();
for m = 2:nD
    t = data(m);  ns = t.nSub;
    [~, ~, g] = unique(t.subj);  tc = accumarray(g, 1);                  % the target's own trials per subject
    for mode = MODES
        med = zeros(B_SUB, 1);  nTr = med;  cv = zeros(B_SUB, 6);  nShort = 0;
        for b = 1:B_SUB
            pick = randperm(nF, ns);
            if mode == "subjects"
                idx = vertcat(frIdx{pick});
            else
                k = tc(randperm(ns));  parts = cell(ns, 1);
                for j = 1:ns
                    a = frIdx{pick(j)};  kk = k(j);
                    if kk > numel(a), nShort = nShort + 1; kk = numel(a); end
                    parts{j} = a(randperm(numel(a), kk));
                end
                idx = vertcat(parts{:});
            end
            med(b) = median(fr.genMed(idx), "omitnan");  nTr(b) = numel(idx);
            for p = 1:6, [~, cv(b, p)] = covVerdict_local(fr.status(idx, p), PASS_COV, COND_COV); end
        end
        if any(~isfinite(med)), error("selDef:nanMedian", "%s", sprintf("Non-finite subsample median, target %s, mode %s", t.name, mode)); end
        q = prctile(med, [2.5 5 97.5]);
        rowsB = [rowsB; table(t.name, mode, ns, sum(tc), mean(nTr), median(med), q(1), q(3), q(3) - q(1), median(med) - cMed(1), cMed(m), 100 * mean(med <= cMed(m)), q(2), nShort, ...
            'VariableNames', ["target" "mode" "nSubjects" "targetTrials" "meanTrials" "medOfMedians" "lo95" "hi95" "width95" "shiftFromFull" "targetMedian" "targetPctile" "p5" "nShort"])]; %#ok<AGROW>
        if mode == "subjects+trials"
            for p = 1:6
                c = prctile(cv(:, p), [2.5 97.5]);
                covB = [covB; table(t.name, fr.labels(p), cov(1, p), median(cv(:, p)), c(1), c(2), mean(cv(:, p) >= PASS_COV), mean(cv(:, p) >= COND_COV), ...
                    'VariableNames', ["target" "pipeline" "fullCoverage" "medCoverage" "covLo" "covHi" "pPass" "pAtLeastCond"])]; %#ok<AGROW>
            end
        end
    end
end
fprintf("\nPart B: Fraser full-sample median %.4f; %d subsamples per cell, without replacement. targetPctile = percentile of the target's own median among the draws; nShort = subject-draws where Fraser had fewer trials than the target's count\n", cMed(1), B_SUB);
disp(rowsB)
fprintf("Coverage under the registered rule for Fraser thinned to each target's design (mode subjects+trials)\n");
disp(covB)

%% Part C: selection shift within each clinical pair, per gate
Bp = rowsA(rowsA.G == "betaObs", :);  gates = unique(Bp.gate, "stable");
rowsC = table();
for k = 1:size(PAIRS, 1)
    for gi = 1:numel(gates)
        a = Bp(Bp.dataset == PAIRS(k, 1) & Bp.gate == gates(gi), :);  c = Bp(Bp.dataset == PAIRS(k, 2) & Bp.gate == gates(gi), :);
        rowsC = [rowsC; table(PAIRS(k, 1) + " - " + PAIRS(k, 2), gates(gi), a.f, c.f, a.rho, c.rho, a.shift, c.shift, a.shift - c.shift, abs(a.shift - c.shift) > semAdequate, ...
            'VariableNames', ["pair" "gate" "fA" "fC" "rhoA" "rhoC" "shiftA" "shiftC" "diff" "exceedsMDC277"])]; %#ok<AGROW>
    end
end
fprintf("\nPart C: selection shift in beta_obs, arm A minus arm C, per gate (exceedsMDC277: |diff| > %.5f)\n", semAdequate);
disp(rowsC)

%% Part D: Cook tempo contrast, all trials against gated trials
ck = find(ismember(DATASETS, ["Cook_CTRL" "Cook_ASD"]));
W = table();
for di = ck
    d = data(di);
    W = [W; table(log2(d.f0), repmat(double(d.name == "Cook_ASD"), numel(d.f0), 1), d.name + "_" + d.subj, isfinite(d.genMed), isfinite(d.gen(:, 6)), ...
        'VariableNames', ["lf" "group" "person" "gAll" "gSG"])]; %#ok<AGROW>
end
if any(~isfinite(W.lf)), error("selDef:lf", "%d non-finite log2 f0 among Cook trials", sum(~isfinite(W.lf))); end
W.person = categorical(W.person);
sets = {"all trials", true(height(W), 1); "gated, all-pipeline", W.gAll; "gated, SG-IRLS", W.gSG};
rowsD = table();
for s = 1:size(sets, 1)
    V = W(sets{s, 2}, :);
    mt = fitlme(V, "lf ~ 1 + group + (1|person)");
    c = mt.Coefficients;  i = find(strcmp(string(c.Name), "group"), 1);
    if isempty(i), error("selDef:coef", "%s", "No group coefficient"); end
    rowsD = [rowsD; table(string(sets{s, 1}), height(V), numel(unique(V.person)), c.Estimate(i), c.Lower(i), c.Upper(i), c.pValue(i), ...
        'VariableNames', ["trialSet" "nTrials" "nSubjects" "b_log2" "ciLo" "ciHi" "p"])]; %#ok<AGROW>
end
fprintf("\nPart D: Cook ASD - CTRL log2 f0 contrast, lf ~ group + (1|subject), by trial set\n");
disp(rowsD)
fprintf("Reference: analyseTempoDecomposition_v002's saved tempo contrast rows (its trial filter also requires finite sigma, a, alpha)\n");
Tm = load(TEMPO_MAT, "C3");  disp(Tm.C3(Tm.C3.outcome == "log2 f0 (tempo)", :))

%% Predictions
verdict_local("P4 subsampling does not move Fraser's centre (both modes)", all(abs(rowsB.shiftFromFull) <= semAdequate), ...
    sprintf("max |shift| %.4f (limit MDC/2.77 = %.5f)", max(abs(rowsB.shiftFromFull)), semAdequate));
w1 = rowsB.width95(rowsB.target == "Cook_CTRL" & rowsB.mode == "subjects");  w2 = rowsB.width95(rowsB.target == "Cook_CTRL" & rowsB.mode == "subjects+trials");
verdict_local("P6 trial thinning widens the interval at the Cook CTRL design", w2 >= 1.5 * w1, sprintf("width %.4f (subjects) -> %.4f (subjects+trials), ratio %.2f", w1, w2, w2 / w1));
cb = rowsB(ismember(rowsB.target, ["Cook_CTRL" "Cook_ASD"]) & rowsB.mode == "subjects+trials", :);
verdict_local("P7 Cook's own medians sit below the 5th percentile of the thinned draws", all(cb.targetMedian < cb.p5), ...
    sprintf("%s", strjoin(cb.target + " median " + compose("%.4f", cb.targetMedian) + " vs p5 " + compose("%.4f", cb.p5) + " (percentile " + compose("%.1f", cb.targetPctile) + ")", "; ")));
fprintf("P8 (not a prediction) %d of %d per-pipeline-and-all-pipeline pair differences exceed MDC/2.77\n", sum(rowsC.exceedsMDC277), height(rowsC));
verdict_local("P9 gating widens the Cook tempo gap", rowsD.b_log2(2) > rowsD.b_log2(1), sprintf("b: all %.3f, gated all-pipeline %.3f, gated SG-IRLS %.3f log2 units", rowsD.b_log2(1), rowsD.b_log2(2), rowsD.b_log2(3)));

save(OUT_MAT, "rowsA", "rowsB", "covB", "rowsC", "rowsD", "cMed", "verd", "cov", "DATASETS", "B_PERM", "B_SUB", "RNG_SEED", "semAdequate", "-v7.3");
fprintf("Saved: %s\n", OUT_MAT);

%% =========================================================================
function t = mengRow_local(R, G, B)
% Meng's identity for one recorded-indicator R and one quantity G (finite G only). Errors if it breaks.
R = R(:);  G = G(:);  ok = isfinite(G);  R = logical(R(ok));  G = G(ok);
N = numel(G);  n = sum(R);
if n < 2 || n >= N - 1, error("selDef:degenerate", "%s", sprintf("n = %d of N = %d", n, N)); end
f = n / N;  sG = std(G, 1);
if sG == 0, error("selDef:constG", "%s", "G is constant"); end
Rz = (double(R) - f) / sqrt(f * (1 - f));  Gz = (G - mean(G)) / sG;
rho = mean(Rz .* Gz);
lhs = mean(G(R)) - mean(G);  rhs = rho * sqrt((1 - f) / f) * sG;
if abs(lhs - rhs) > 1e-9 * max(1, sG), error("selDef:identity", "Meng identity broken: %.3e vs %.3e", lhs, rhs); end
pr = zeros(B, 1);
for b = 1:B, pr(b) = mean(Rz(randperm(N)) .* Gz); end
pPerm = (1 + sum(abs(pr) >= abs(rho))) / (B + 1);
nStar = f / ((1 - f) * rho^2);
if rho == 0, nEff = N; else, nEff = nStar / (1 + (nStar - 1) / N); end
t = table(N, n, f, rho, sG, lhs, sqrt(N - 1) * rho, nEff, pPerm, 'VariableNames', ["N" "n" "f" "rho" "sdG" "shift" "Z" "nEff" "pPerm"]);
end

function [verdict, coverage] = covVerdict_local(statusCol, passCov, condCov)
% The registered coverage rule (compare42CellVerdict_v002 verdictOneCell_local).
valid = statusCol ~= "no_beta_obs";
if ~any(valid), verdict = "no-data"; coverage = NaN; return, end
coverage = sum(statusCol(valid) == "rise") / sum(valid);
if coverage >= passCov, verdict = "PASS"; elseif coverage >= condCov, verdict = "CONDITIONAL"; else, verdict = "FAIL"; end
end

function verdict_local(tag, held, detail)
s = "NOT HELD";  if held, s = "HELD"; end
fprintf("%s %s: %s\n", tag, s, detail);
end
