%% checkSelectionDefect_v001.m
% Meng (2018, Ann. Appl. Stat. 12:685-726) decomposes the error of a sample mean over a non-random
% subset as  mean(G | recorded) - mean(G) = rho(R,G) * sqrt((1-f)/f) * sd(G)   (data quality x data
% quantity x difficulty; exact, no assumptions). Here R = the trial entered beta_gen* (the
% invertibility gate), f = share of trials that did, G = any per-trial quantity observed for ALL trials.
% Truth for beta_gen is not observed on non-invertible trials, but beta_obs, tempo (f0) and the noise
% descriptors are, so the selection the gate imposes on them can be measured exactly (Part A).
% Part B asks the N question directly: subsample Fraser by subject to each smaller dataset's subject
% count (without replacement) and see whether the centre of the constellation median moves, and how
% often the registered coverage rule would still return PASS.
% Prior statements this quantifies (checked 2026-10-11, ledger-first): FINDINGS_REFERENCE.md L825 (under
% item 45) says inversion fails preferentially for the most compressed (low beta_obs) trials, biasing
% reported means upward, and that "group comparisons are unaffected, bias direction is symmetric within
% each pair"; Finding #142 (L3464) records N_fin as a smaller, untraced subset of the raw trial count.
% Part A measures the first and tests the second; for the Cook arms and the Hickman arms the count
% here is the finite betaGenStarMed count in the v015 file, which differs slightly from the hub's N_fin.
% Cell definitions: in-domain datasets only (Fraser, Cook CTRL, Cook ASD, Hickman PLAC, Hickman HALO).
% R (all-pipeline) = isfinite(betaGenStarMed). R (pipeline p) = invertStatus(p) == "rise", over trials
% with a beta_obs (the rule in compare42CellVerdict_v002 verdictOneCell_local, L73-90).
% What this can and cannot show: Part A measures selection on OBSERVABLES. It does not measure
% selection on the unobserved generator exponent, and Meng's identity covers selection, not
% measurement bias (the pull to the noise exponent is a separate channel). The permutation null
% shuffles R across trials and ignores subject clustering, so p-values are optimistic. Meng's identity
% is for means; the constellation median is an approximation (his footnote 4).
% Predictions, fixed before the run (2026-10-11):
%  P1 gate passers have a HIGHER all-pipeline median beta_obs than the dataset's trials overall
%     (shift > 0) in Cook CTRL, Cook ASD and Hickman HALO, permutation p < .05 in Cook CTRL and ASD.
%     (Non-invertible trials are mostly slow and off the map's floor, #247.) If P1 fails, selection
%     cannot be what separates Fraser from the smaller datasets in either direction.
%  P2 the gate selects on tempo: rho(R, log2 f0) > 0 in all five datasets (#235).
%  P3 Fraser's |shift| in beta_obs is the smallest of the five (Meng's quantity term).
%  P4 Fraser subsampled to 19, 21 and 35 subjects keeps its centre: the median of the subsample
%     medians stays within MDC/2.77 of the full-Fraser median at every subject count (N alone does not
%     move it). Where each small dataset's own median falls in that distribution is reported, no
%     prediction. rho with sigmaMM and alphaIRA is exploratory, no prediction.
%  P5 the selection shift in beta_obs is symmetric within each clinical pair: |shift(ASD) - shift(CTRL)|
%     and |shift(HALO) - shift(PLAC)| are each below MDC/2.77 (all-pipeline gate). If P5 fails, the
%     "bias is symmetric within each pair" sentence at L825 does not hold for that pair.
% Anchors (error if broken): dataset medians of betaGenStarMed span 0.272-0.328 (Finding #247); the
% registered rule gives PASS/CONDITIONAL/FAIL = 4/2/24 over these 30 cells (N_REPS = 200, Finding #242);
% Meng's identity holds to 1e-9 in every row.
% Fail loud: nothing is filled; the hub's N_fin is printed beside this run's count, not forced to match.
% Run from src/. Reads src/loopClosureResults_<dataset>_all_shaped_xu_v015.mat; writes
% results/checkSelectionDefect_v001.mat

%% CONFIG
ROOT      = fileparts(fileparts(mfilename("fullpath")));
DATASETS  = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO"];
FILE_FMT  = "loopClosureResults_%s_all_shaped_xu_v015.mat";
OUT_MAT   = fullfile(ROOT, "results", "checkSelectionDefect_v001.mat");
HUB_NFIN  = [2744 72 81 343 301];
B_PERM    = 5000;
B_SUB     = 2000;
RNG_SEED  = 1729;
ANCHOR    = [0.272 0.328];
VERDICTS  = [4 2 24];
PASS_COV  = 0.95;
COND_COV  = 0.90;
MDC       = 0.03;
G_NAMES   = ["betaObs" "log2 f0" "sigmaMM" "alphaIRA"];

%% Setup
addpath(genpath(fullfile(ROOT, "src", "functions")));
if ~isfolder(fileparts(OUT_MAT)), error("selDef:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
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

%% Anchors and counts
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
nFin = arrayfun(@(d) sum(isfinite(d.genMed)), data);
finNotRise = arrayfun(@(d) sum(isfinite(d.gen(:)) & d.status(:) ~= "rise"), data);
cnt = table(DATASETS', arrayfun(@(d) numel(d.genMed), data)', nFin', nFin' ./ arrayfun(@(d) numel(d.genMed), data)', HUB_NFIN', ...
    nFin' - HUB_NFIN', [data.nSub]', finNotRise', 'VariableNames', ["dataset" "nTrials" "nFinite" "f" "hubNfin" "diffFromHub" "nSubjects" "finiteNotRise"]);
disp(cnt)
fprintf("diffFromHub: v015 count minus the hub's N_fin (a different file/basis); finiteNotRise counts finite beta_gen* cells whose status is not 'rise'\n");

%% Part A: Meng decomposition of the gate
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
cols = ["dataset" "gate" "G" "N" "n" "f" "rho" "shift" "shiftMDC" "Z" "nEff" "pPerm"];
fprintf("\nPart A: selection by the gate. shift = mean(G | passed) - mean(G); shiftMDC is meaningful for betaObs only; nEff = SRS-equivalent trials (Meng 3.5)\n");
disp(rowsA(rowsA.gate == "all-pipeline", cols))
fprintf("Per pipeline, G = betaObs (R = status 'rise' among trials with a beta_obs)\n");
disp(rowsA(rowsA.gate ~= "all-pipeline", cols))

%% Part B: Fraser subsampled by subject
fr = data(1);  us = unique(fr.subj);  nsT = unique([data(2:end).nSub]);
rowsB = table();  sumB = table();
for ns = nsT
    med = zeros(B_SUB, 1);  nTr = med;  cv = zeros(B_SUB, 6);
    for b = 1:B_SUB
        idx = ismember(fr.subj, us(randperm(numel(us), ns)));
        med(b) = median(fr.genMed(idx), "omitnan");  nTr(b) = sum(idx);
        for p = 1:6, [~, cv(b, p)] = covVerdict_local(fr.status(idx, p), PASS_COV, COND_COV); end
    end
    q = prctile(med, [2.5 97.5]);
    sumB = [sumB; table(ns, mean(nTr), median(med), q(1), q(2), median(med) - cMed(1), 'VariableNames', ["nSubjects" "meanTrials" "medOfMedians" "lo95" "hi95" "shiftFromFull"])]; %#ok<AGROW>
    for p = 1:6
        c = prctile(cv(:, p), [2.5 97.5]);
        rowsB = [rowsB; table(ns, fr.labels(p), cov(1, p), median(cv(:, p)), c(1), c(2), mean(cv(:, p) >= PASS_COV), mean(cv(:, p) >= COND_COV), ...
            'VariableNames', ["nSubjects" "pipeline" "fullCoverage" "medCoverage" "covLo" "covHi" "pPass" "pAtLeastCond"])]; %#ok<AGROW>
    end
    for m = find([data(2:end).nSub] == ns) + 1
        fprintf("Fraser at %d subjects: %s's own median %.4f sits at percentile %.1f of the subsample medians\n", ns, DATASETS(m), cMed(m), 100 * mean(med <= cMed(m)));
    end
end
fprintf("\nPart B: Fraser full-sample median %.4f; %d subsamples per subject count, without replacement\n", cMed(1), B_SUB);
disp(sumB)
disp(rowsB)

%% Predictions
A = rowsA(rowsA.gate == "all-pipeline" & rowsA.G == "betaObs", :);
sm = ismember(A.dataset, ["Cook_CTRL" "Cook_ASD" "Hickman_HALO"]);  ck = ismember(A.dataset, ["Cook_CTRL" "Cook_ASD"]);
p1 = all(A.shift(sm) > 0) && all(A.pPerm(ck) < 0.05);
verdict_local("P1 gate passers have higher beta_obs in the smaller datasets", p1, sprintf("shift [p]: %s", strjoin(A.dataset + " " + compose("%+.4f", A.shift) + " [" + compose("%.3f", A.pPerm) + "]", "; ")));
T = rowsA(rowsA.gate == "all-pipeline" & rowsA.G == "log2 f0", :);
verdict_local("P2 the gate selects on tempo", all(T.rho > 0), sprintf("rho: %s", strjoin(T.dataset + " " + compose("%+.3f", T.rho), "; ")));
isF = A.dataset == "Fraser";
verdict_local("P3 Fraser's selection shift is the smallest", abs(A.shift(isF)) < min(abs(A.shift(~isF))), sprintf("|shift|: %s", strjoin(A.dataset + " " + compose("%.4f", abs(A.shift)), "; ")));
verdict_local("P4 subsampling Fraser does not move its centre", all(abs(sumB.shiftFromFull) <= semAdequate), sprintf("shift at %s subjects: %s (limit MDC/2.77 = %.5f)", strjoin(string(sumB.nSubjects), "/"), strjoin(compose("%+.4f", sumB.shiftFromFull), "/"), semAdequate));

PAIRS = ["Cook_ASD" "Cook_CTRL"; "Hickman_HALO" "Hickman_PLAC"];  dPair = zeros(size(PAIRS, 1), 1);  pairTxt = strings(size(PAIRS, 1), 1);
for k = 1:size(PAIRS, 1)
    a = A(A.dataset == PAIRS(k, 1), :);  c = A(A.dataset == PAIRS(k, 2), :);
    dPair(k) = a.shift - c.shift;
    pairTxt(k) = sprintf("%s - %s: f %.3f vs %.3f, rho %+.3f vs %+.3f, shift %+.4f vs %+.4f, difference %+.4f", PAIRS(k, 1), PAIRS(k, 2), a.f, c.f, a.rho, c.rho, a.shift, c.shift, dPair(k));
end
verdict_local("P5 the selection shift is symmetric within each clinical pair", all(abs(dPair) < semAdequate), strjoin(pairTxt, " | "));

save(OUT_MAT, "rowsA", "rowsB", "sumB", "cnt", "cMed", "verd", "cov", "DATASETS", "B_PERM", "B_SUB", "RNG_SEED", "semAdequate", "-v7.3");
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
