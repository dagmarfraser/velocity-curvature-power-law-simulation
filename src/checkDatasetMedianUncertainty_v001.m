%% checkDatasetMedianUncertainty_v001.m
% checkSelectionDefect_v002 Part B showed that Fraser, subsampled to each smaller dataset's design,
% keeps its centre near 0.328 and the smaller datasets' own medians fall in the extreme lower tail.
% That comparison uses FRASER's between-subject spread. A clinical or smaller sample may be more
% heterogeneous, so the proper test uses each dataset's OWN uncertainty. This script:
%  (1) subject-cluster bootstrap (resample subjects with replacement, keep all their trials) of each
%      in-domain dataset's median beta_gen*, for two estimators: the constellation median
%      (betaGenStarMed, the basis of the 0.272-0.328 anchor) and SG-IRLS (betaGenStar column 6);
%  (2) Fraser minus each smaller dataset, from independent bootstraps (different participants), with a
%      percentile 95% CI and a two-sided bootstrap p;
%  (3) the clinical contrasts on the same footing: Cook ASD minus CTRL (independent groups) and
%      Hickman HALO minus PLAC (paired: the same subjects resampled in both arms).
% Unit of analysis: trial within subject; inference clusters on subject. Non-invertible trials (NaN
% beta_gen*) are excluded from the median, as in the anchor; selection by the gate is the subject of
% checkSelectionDefect_v001/v002 and is not corrected here.
% What this can and cannot show: it asks whether the observed differences exceed what each dataset's
% own sampling variability produces. It cannot separate population, device, noise and paradigm, which
% every dataset bundles. Percentile bootstraps with 19-21 subjects are rough (few distinct resamples of
% the tails); read the CIs as indicative.
% Predictions, fixed 2026-10-11 before the run:
%  Q1 for Cook CTRL and Cook ASD the 95% CI of (Fraser - dataset), constellation median, includes zero:
%     on their own variance the data cannot separate them from Fraser.
%  Q2 for at least one Hickman arm (35 subjects) the CI excludes zero.
%  Q3 the Cook ASD - CTRL constellation-median difference has a CI half-width between 0.03 and 0.08
%     (the "CIs about +/-0.05" in Finding #235 / the hub).
% Anchor (error if broken): the dataset constellation medians span 0.272-0.328 (Finding #247).
% Fail loud: a bootstrap resample with no finite beta_gen* is an error; counts are printed.
% Run from src/. Reads src/loopClosureResults_<dataset>_all_shaped_xu_v015.mat; writes
% results/checkDatasetMedianUncertainty_v001.mat

%% CONFIG
ROOT      = fileparts(fileparts(mfilename("fullpath")));
DATASETS  = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO"];
FILE_FMT  = "loopClosureResults_%s_all_shaped_xu_v015.mat";
OUT_MAT   = fullfile(ROOT, "results", "checkDatasetMedianUncertainty_v001.mat");
EST       = ["constellation" "SG-IRLS"];
B_BOOT    = 5000;
RNG_SEED  = 1729;
ANCHOR    = [0.272 0.328];

%% Setup
if ~isfolder(fileparts(OUT_MAT)), error("dsUnc:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
rng(RNG_SEED, "twister");

%% Load
nD = numel(DATASETS);  data = repmat(struct(), 1, nD);
for di = 1:nD
    f = fullfile(ROOT, "src", sprintf(FILE_FMT, DATASETS(di)));
    if ~isfile(f), error("dsUnc:in", "%s", "Missing input: " + f); end
    S = load(f, "results", "pipelineLabels");  r = S.results(:);
    if string(S.pipelineLabels(6)) ~= "SG-IRLS", error("dsUnc:labels", "%s", "Pipeline 6 is not SG-IRLS in " + DATASETS(di)); end
    g = horzcat(r.betaGenStar)';
    data(di).subj = string({r.subjectID})';
    data(di).vals = [[r.betaGenStarMed]', g(:, 6)];
    data(di).ids  = unique(data(di).subj);
end

%% Anchor
cMed = arrayfun(@(d) median(d.vals(:, 1), "omitnan"), data);
if round(min(cMed), 3) ~= ANCHOR(1) || round(max(cMed), 3) ~= ANCHOR(2)
    error("dsUnc:anchor", "Dataset medians span %.4f-%.4f, not Finding #247's %.3f-%.3f.", min(cMed), max(cMed), ANCHOR(1), ANCHOR(2));
end
fprintf("Anchor passed: constellation medians span %.4f-%.4f\n", min(cMed), max(cMed));

%% (1) Own-variance bootstrap of each dataset's median
boot = cell(nD, numel(EST));  rowsM = table();
for di = 1:nD
    for e = 1:numel(EST)
        d = data(di);  vs = bySubj_local(d.subj, d.vals(:, e), d.ids);
        boot{di, e} = bootMedian_local(vs, B_BOOT, DATASETS(di));
        q = prctile(boot{di, e}, [2.5 97.5]);  m = median(d.vals(:, e), "omitnan");
        rowsM = [rowsM; table(DATASETS(di), EST(e), numel(d.ids), sum(isfinite(d.vals(:, e))), m, q(1), q(2), (q(2) - q(1)) / 2, ...
            'VariableNames', ["dataset" "estimator" "nSubjects" "nFinite" "median" "ciLo" "ciHi" "halfWidth"])]; %#ok<AGROW>
    end
end
fprintf("\n(1) Subject-cluster bootstrap of each dataset's median beta_gen* (%d draws, seed %d)\n", B_BOOT, RNG_SEED);
disp(rowsM)

%% (2) Fraser minus each smaller dataset
rowsD = table();
for di = 2:nD
    for e = 1:numel(EST)
        rowsD = [rowsD; diffRow_local("Fraser - " + DATASETS(di), EST(e), boot{1, e} - boot{di, e}, ...
            median(data(1).vals(:, e), "omitnan") - median(data(di).vals(:, e), "omitnan"))]; %#ok<AGROW>
    end
end

%% (3) Clinical contrasts
iC = find(DATASETS == "Cook_CTRL");  iA = find(DATASETS == "Cook_ASD");
iP = find(DATASETS == "Hickman_PLAC");  iH = find(DATASETS == "Hickman_HALO");
ids = union(data(iP).ids, data(iH).ids);
fprintf("Hickman paired bootstrap: %d subjects in PLAC, %d in HALO, %d in the union, %d in both\n", numel(data(iP).ids), numel(data(iH).ids), numel(ids), numel(intersect(data(iP).ids, data(iH).ids)));
for e = 1:numel(EST)
    rowsD = [rowsD; diffRow_local("Cook_ASD - Cook_CTRL", EST(e), boot{iA, e} - boot{iC, e}, ...
        median(data(iA).vals(:, e), "omitnan") - median(data(iC).vals(:, e), "omitnan"))]; %#ok<AGROW>
    [mH, mP] = bootMedianPair_local(bySubj_local(data(iH).subj, data(iH).vals(:, e), ids), bySubj_local(data(iP).subj, data(iP).vals(:, e), ids), B_BOOT);
    rowsD = [rowsD; diffRow_local("Hickman_HALO - Hickman_PLAC (paired)", EST(e), mH - mP, ...
        median(data(iH).vals(:, e), "omitnan") - median(data(iP).vals(:, e), "omitnan"))]; %#ok<AGROW>
end
fprintf("\n(2)-(3) Differences in median beta_gen*: point estimate, percentile 95%% CI, two-sided bootstrap p\n");
disp(rowsD)

%% Predictions
K = rowsD(rowsD.estimator == "constellation", :);
q1 = ismember(K.contrast, ["Fraser - Cook_CTRL" "Fraser - Cook_ASD"]);
verdict_local("Q1 Cook CIs include zero against Fraser", all(K.includesZero(q1)), strjoin(K.contrast(q1) + " " + compose("%+.3f", K.diff(q1)) + " [" + compose("%+.3f", K.ciLo(q1)) + ", " + compose("%+.3f", K.ciHi(q1)) + "]", "; "));
q2 = ismember(K.contrast, ["Fraser - Hickman_PLAC" "Fraser - Hickman_HALO"]);
verdict_local("Q2 at least one Hickman arm is separated from Fraser", any(~K.includesZero(q2)), strjoin(K.contrast(q2) + " " + compose("%+.3f", K.diff(q2)) + " [" + compose("%+.3f", K.ciLo(q2)) + ", " + compose("%+.3f", K.ciHi(q2)) + "]", "; "));
q3 = K.halfWidth(K.contrast == "Cook_ASD - Cook_CTRL");
verdict_local("Q3 Cook ASD - CTRL CI half-width is 0.03-0.08", q3 >= 0.03 && q3 <= 0.08, sprintf("half-width %.3f", q3));

save(OUT_MAT, "rowsM", "rowsD", "cMed", "DATASETS", "EST", "B_BOOT", "RNG_SEED", "-v7.3");
fprintf("Saved: %s\n", OUT_MAT);

%% =========================================================================
function vs = bySubj_local(subj, v, ids)
% Finite beta_gen* values per subject id (empty where a subject has none).
vs = arrayfun(@(u) v(subj == u & isfinite(v)), ids, 'UniformOutput', false);
end

function m = bootMedian_local(vs, B, name)
ns = numel(vs);  m = zeros(B, 1);
for b = 1:B
    x = vertcat(vs{randi(ns, ns, 1)});
    if isempty(x), error("dsUnc:empty", "%s", "Bootstrap resample with no finite beta_gen* in " + name); end
    m(b) = median(x);
end
end

function [mA, mB] = bootMedianPair_local(vsA, vsB, B)
% Shared subject resampling for two arms measured on the same subjects.
ns = numel(vsA);  mA = zeros(B, 1);  mB = mA;
for b = 1:B
    ix = randi(ns, ns, 1);  xa = vertcat(vsA{ix});  xb = vertcat(vsB{ix});
    if isempty(xa) || isempty(xb), error("dsUnc:emptyPair", "%s", "Paired resample with no finite beta_gen* in one arm"); end
    mA(b) = median(xa);  mB(b) = median(xb);
end
end

function t = diffRow_local(name, est, d, pt)
q = prctile(d, [2.5 97.5]);  B = numel(d);
p = min(1, 2 * (1 + min(sum(d <= 0), sum(d >= 0))) / (B + 1));
t = table(name, est, pt, q(1), q(2), (q(2) - q(1)) / 2, p, q(1) <= 0 && q(2) >= 0, ...
    'VariableNames', ["contrast" "estimator" "diff" "ciLo" "ciHi" "halfWidth" "pBoot" "includesZero"]);
end

function verdict_local(tag, held, detail)
s = "NOT HELD";  if held, s = "HELD"; end
fprintf("%s %s: %s\n", tag, s, detail);
end
