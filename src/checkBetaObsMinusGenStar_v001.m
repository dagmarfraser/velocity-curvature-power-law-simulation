%% checkBetaObsMinusGenStar_v001.m
% Sol Round 4, row 1: v004 L2749-2759 says the artefacts "would have hidden a one-third law, not
% created one". That sentence is a claim about the sign and size of beta_obs minus beta_gen*: it
% needs a script or it goes (docs/CritiqueResponse_ResearchTrajectory_v001.md, row 1).
% For each in-domain dataset (Fraser, Cook CTRL, Cook ASD, Hickman PLAC, Hickman HALO; Zarandi and
% Dhieb are outside the validity domain) and each pipeline, plus the all-pipeline median:
%   beta_obs      the pipeline's output on the empirical trial;
%   beta_gen*     the pointwise inversion's generator estimate for that trial and pipeline (NaN where
%                 the trial is not invertible for that pipeline);
%   difference    beta_obs - beta_gen*, paired by trial, finite pairs only. Median, subject-cluster
%                 bootstrap 95% CI (resampling subjects, B_BOOT draws, percentile), share of pairs
%                 below zero, and the distance of each median from 1/3.
% Unit of analysis: trial within subject; inference clusters on subject. For the all-pipeline row,
% beta_obs is the per-trial median over the six pipelines and beta_gen* is betaGenStarMed.
% What this can and cannot show: it measures how far inversion moves a measured exponent and in
% which direction. It cannot test model discrepancy: inversion assumes a power-law generator and
% the six pipelines share that assumption (Sol Round 4, rows 1 and 2), so a discrepancy common to
% all of them is invisible here. A HELD below is a statement about the calibrated pipelines, not
% about the generator.
% Predictions, fixed before the run (2026-10-01), on the all-pipeline rows:
%  P1 inversion does not move any dataset's median away from 1/3 by more than MDC/2.77:
%     |median beta_gen* - 1/3| <= |median beta_obs - 1/3| + MDC/2.77 in all five datasets.
%  P2 the artefacts pull the measured exponent down, in every dataset: the upper bound of the
%     cluster-bootstrap CI of the median difference is below zero.
% "Hidden, not created" needs both. If P2 fails in a dataset, the sentence cannot be applied to it.
% Anchor (error if broken): the dataset medians of betaGenStarMed span 0.272 to 0.328 at three
% decimals, the v015 range of Finding #247 (checkAbstractBetaGenRange_v001).
% Fail loud: non-finite counts are printed per cell, nothing is filled; fewer than three subjects
% in a cell is an error.
% Run from src/. Reads src/loopClosureResults_<dataset>_all_shaped_xu_v015.mat; writes
% results/checkBetaObsMinusGenStar_v001.mat

%% CONFIG
ROOT      = fileparts(fileparts(mfilename("fullpath")));
DATASETS  = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO"];
FILE_FMT  = "loopClosureResults_%s_all_shaped_xu_v015.mat";
OUT_MAT   = fullfile(ROOT, "results", "checkBetaObsMinusGenStar_v001.mat");
TARGET    = 1/3;
B_BOOT    = 2000;
RNG_SEED  = 1729;
ANCHOR    = [0.272 0.328];

%% Setup
addpath(genpath(fullfile(ROOT, "src", "functions")));
if ~isfolder(fileparts(OUT_MAT)), error("obsGen:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
semAdequate = semAdequacyThreshold_v001;
rng(RNG_SEED, "twister");

%% Per dataset and pipeline
rows = table();  cMed = zeros(1, numel(DATASETS));
for di = 1:numel(DATASETS)
    f = fullfile(ROOT, "src", sprintf(FILE_FMT, DATASETS(di)));
    if ~isfile(f), error("obsGen:in", "%s", "Missing input: " + f); end
    S = load(f, "results", "pipelineLabels");
    r = S.results;  labels = string(S.pipelineLabels(:))';
    obs = vertcat(r.betaObs);  gen = horzcat(r.betaGenStar)';  genMed = [r.betaGenStarMed]';
    subj = string({r.subjectID})';
    if size(obs, 2) ~= numel(labels) || size(gen, 2) ~= numel(labels) || size(obs, 1) ~= numel(r) || size(gen, 1) ~= numel(r)
        error("obsGen:shape", "%s", "Unexpected betaObs/betaGenStar shape in " + DATASETS(di));
    end
    cMed(di) = median(genMed, "omitnan");
    for p = 1:numel(labels) + 1
        if p <= numel(labels), o = obs(:, p);  g = gen(:, p);  lab = labels(p);
        else, o = median(obs, 2, "omitnan");  g = genMed;  lab = "all-pipeline median"; end
        ok = isfinite(o) & isfinite(g);
        [gi, ~] = findgroups(subj(ok));  ns = max([gi; 0]);
        if ns < 3, error("obsGen:subjects", "%s", sprintf("%s, %s: %d subjects with finite pairs", DATASETS(di), lab, ns)); end
        d = o(ok) - g(ok);
        dsub = accumarray(gi, d, [], @(x) {x});
        bm = zeros(B_BOOT, 1);
        for b = 1:B_BOOT, bm(b) = median(vertcat(dsub{randi(ns, ns, 1)})); end
        ci = prctile(bm, [2.5 97.5]);
        mo = median(o(ok));  mg = median(g(ok));
        rows = [rows; table(DATASETS(di), lab, numel(r), sum(isfinite(o)), sum(isfinite(g)), sum(ok), ns, mo, mg, median(d), ci(1), ci(2), ...
            mean(d < 0), abs(mo - TARGET), abs(mg - TARGET), 'VariableNames', ["dataset" "pipeline" "nTrials" "nObsFinite" "nGenFinite" ...
            "nPaired" "nSubjects" "medObs" "medGen" "medDiff" "ciLo" "ciHi" "shareBelowZero" "devObsFromThird" "devGenFromThird"])]; %#ok<AGROW>
    end
end

%% Anchor
if round(min(cMed), 3) ~= ANCHOR(1) || round(max(cMed), 3) ~= ANCHOR(2)
    error("obsGen:anchor", "Dataset medians of betaGenStarMed span %.4f-%.4f, not Finding #247's %.3f-%.3f.", min(cMed), max(cMed), ANCHOR(1), ANCHOR(2));
end
fprintf("Anchor passed: dataset medians of betaGenStarMed span %.4f-%.4f (Finding #247: %.3f-%.3f)\n", min(cMed), max(cMed), ANCHOR(1), ANCHOR(2));
fprintf("MDC/2.77 = %.5f; target %.4f; %d subject-cluster bootstrap draws, seed %d\n", semAdequate, TARGET, B_BOOT, RNG_SEED);

%% Report
for di = 1:numel(DATASETS)
    fprintf("\n%s: beta_obs - beta_gen*, paired by trial (finite pairs), subject-cluster bootstrap CI\n", DATASETS(di));
    disp(rows(rows.dataset == DATASETS(di), ["pipeline" "nTrials" "nObsFinite" "nGenFinite" "nPaired" "nSubjects" "medObs" "medGen" "medDiff" "ciLo" "ciHi" "shareBelowZero"]))
end

%% Predictions
A = rows(rows.pipeline == "all-pipeline median", :);
p1 = A.devGenFromThird <= A.devObsFromThird + semAdequate;
verdict_local("P1 inversion does not move the median away from 1/3", all(p1), sprintf("distance from 1/3, obs then gen*: %s%s", ...
    strjoin(A.dataset + " " + compose("%.3f", A.devObsFromThird) + " -> " + compose("%.3f", A.devGenFromThird), "; "), ...
    ternary_local(all(p1), "", "; fails in " + strjoin(A.dataset(~p1), ", "))));
p2 = A.ciHi < 0;
verdict_local("P2 the artefacts pull the measured exponent down", all(p2), sprintf("median difference [CI]: %s%s", ...
    strjoin(A.dataset + " " + compose("%.3f", A.medDiff) + " [" + compose("%.3f", A.ciLo) + ", " + compose("%.3f", A.ciHi) + "]", "; "), ...
    ternary_local(all(p2), "", "; CI not below zero in " + strjoin(A.dataset(~p2), ", "))));

save(OUT_MAT, "rows", "cMed", "DATASETS", "TARGET", "B_BOOT", "RNG_SEED", "semAdequate", "-v7.3");
fprintf("Saved: %s\n", OUT_MAT);

%% =========================================================================
function verdict_local(tag, held, detail)
s = "NOT HELD";  if held, s = "HELD"; end
fprintf("%s %s: %s\n", tag, s, detail);
end

function s = ternary_local(c, a, b)
if c, s = string(a); else, s = string(b); end
end
