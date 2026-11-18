%% checkDatasetMedianUncertainty_v004.m
% v004 (2026-11-16): VERSION "v028" (end-matched corpus via functions/loadV028_v001; workspace variables ARM (primary|budget1|shapeXu|harmonicsFc, default primary) and OBS_BASIS (window|whole, default window) select the arm and the beta_obs basis; output name carries _v028[_arm_basis]). The Table 6b anchor applies to v015 only, as before. Regression: VERSION v015 must reproduce checkDatasetMedianUncertainty_v003_cluster_v015.mat.
% (v003 header follows)
%% checkDatasetMedianUncertainty_v003.m
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

% v002 (2026-10-31, GD12): VERSION selects the corpus ("v015" reproduces v001 and is anchored; "v019" is the end-anchored rerun, read from src/ or
% the RDS src). Set VERSION in the workspace before running; default "v019". Outputs carry the version in their name.

% v003 (2026-10-31, GD5/GD13): CLUSTER basis. The statistic is the median of per-subject medians (Table 6b's primary estimate, "cluster_median.pointEst"),
% not the flat median over trials used in v001/v002; the bootstrap resamples subjects and recomputes that statistic. Subjects with no finite value
% drop out of the median. For v015 the five constellation point estimates must reproduce Table 6b (0.2749-0.3309; error if not). Everything else
% (seed, B, pipelines, contrasts) is as v002, so flat v002 and cluster v003 differ only in the statistic. SG-IRLS rows are the GD13 column.

%% CONFIG
if ~exist("VERSION", "var"), VERSION = "v019"; end
VERSION = string(VERSION);
if ~ismember(VERSION, ["v015" "v019" "v028"]), error("dsUnc:version", "VERSION must be v015, v019 or v028, got %s", VERSION); end
if ~exist("ARM", "var"), ARM = "primary"; end
if ~exist("OBS_BASIS", "var"), OBS_BASIS = "window"; end
ARM = string(ARM);  OBS_BASIS = string(OBS_BASIS);
VTAG = VERSION;  if VERSION == "v028" && (ARM ~= "primary" || OBS_BASIS ~= "window"), VTAG = "v028_" + ARM + "_" + OBS_BASIS; end
addpath(genpath(fullfile(fileparts(mfilename("fullpath")), "functions")));   % loadV028_v001
ROOT      = fileparts(fileparts(mfilename("fullpath")));
DATASETS  = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO"];
FILE_FMT  = "loopClosureResults_%s_all_shaped_xu_" + VERSION + ".mat";
OUT_MAT   = fullfile(ROOT, "results", "checkDatasetMedianUncertainty_v004_cluster_" + VTAG + ".mat");
EST       = ["constellation" "SG-IRLS"];
B_BOOT    = 5000;
RNG_SEED  = 1729;
ANCHOR    = [0.2749 0.3309];   % Table 6b cluster point estimates, v015 (min, max)
RDS_SRC   = "/Volumes/rdsprojects/f/fraserds-mpo-evaluation/2026_prereg/velocity-curvature-power-law-simulation-main/velocity-curvature-power-law-simulation-main/src";

%% Setup
if ~isfolder(fileparts(OUT_MAT)), error("dsUnc:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
rng(RNG_SEED, "twister");

%% Load
nD = numel(DATASETS);  data = repmat(struct(), 1, nD);
for di = 1:nD
    if VERSION == "v028"
        r = loadV028_v001(DATASETS(di), Arm=ARM, ObsBasis=OBS_BASIS);  r = r(:);  S = struct("pipelineLabels", {{"BWFD-OLS", "SG-OLS", "BWFD-LMLS", "SG-LMLS", "BWFD-IRLS", "SG-IRLS"}});
    else
    f = fullfile(ROOT, "src", sprintf(FILE_FMT, DATASETS(di)));
    if ~isfile(f)
        g = fullfile(RDS_SRC, sprintf(FILE_FMT, DATASETS(di)));
        if ~isfile(g), error("dsUnc:in", "%s", "Missing input in src/ and on the RDS: " + sprintf(FILE_FMT, DATASETS(di))); end
        f = fullfile(tempdir, sprintf(FILE_FMT, DATASETS(di)));  copyfile(g, f);
    end
    S = load(f, "results", "pipelineLabels");  r = S.results(:);
    end
    if string(S.pipelineLabels(6)) ~= "SG-IRLS", error("dsUnc:labels", "%s", "Pipeline 6 is not SG-IRLS in " + DATASETS(di)); end
    g = horzcat(r.betaGenStar)';
    data(di).subj = string({r.subjectID})';
    data(di).vals = [[r.betaGenStarMed]', g(:, 6)];
    data(di).ids  = unique(data(di).subj);
end

%% Anchor
cMed = arrayfun(@(d) clusterStat_local(bySubj_local(d.subj, d.vals(:, 1), d.ids)), data);
if VERSION == "v015" && (round(min(cMed), 4) ~= ANCHOR(1) || round(max(cMed), 4) ~= ANCHOR(2))
    error("dsUnc:anchor", "Dataset medians span %.4f-%.4f, not Table 6b cluster %.4f-%.4f.", min(cMed), max(cMed), ANCHOR(1), ANCHOR(2));
end
fprintf("%s: constellation CLUSTER medians span %.4f-%.4f (anchor 0.2749-0.3309 applies to v015 only)\n", VERSION, min(cMed), max(cMed));

%% (1) Own-variance bootstrap of each dataset's median
boot = cell(nD, numel(EST));  rowsM = table();
for di = 1:nD
    for e = 1:numel(EST)
        d = data(di);  vs = bySubj_local(d.subj, d.vals(:, e), d.ids);
        boot{di, e} = bootMedian_local(vs, B_BOOT, DATASETS(di));
        q = prctile(boot{di, e}, [2.5 97.5]);  m = clusterStat_local(vs);
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
            clusterStat_local(bySubj_local(data(1).subj, data(1).vals(:, e), data(1).ids)) - clusterStat_local(bySubj_local(data(di).subj, data(di).vals(:, e), data(di).ids)))]; %#ok<AGROW>
    end
end

%% (3) Clinical contrasts
iC = find(DATASETS == "Cook_CTRL");  iA = find(DATASETS == "Cook_ASD");
iP = find(DATASETS == "Hickman_PLAC");  iH = find(DATASETS == "Hickman_HALO");
ids = union(data(iP).ids, data(iH).ids);
fprintf("Hickman paired bootstrap: %d subjects in PLAC, %d in HALO, %d in the union, %d in both\n", numel(data(iP).ids), numel(data(iH).ids), numel(ids), numel(intersect(data(iP).ids, data(iH).ids)));
for e = 1:numel(EST)
    rowsD = [rowsD; diffRow_local("Cook_ASD - Cook_CTRL", EST(e), boot{iA, e} - boot{iC, e}, ...
        clusterStat_local(bySubj_local(data(iA).subj, data(iA).vals(:, e), data(iA).ids)) - clusterStat_local(bySubj_local(data(iC).subj, data(iC).vals(:, e), data(iC).ids)))]; %#ok<AGROW>
    [mH, mP] = bootMedianPair_local(bySubj_local(data(iH).subj, data(iH).vals(:, e), ids), bySubj_local(data(iP).subj, data(iP).vals(:, e), ids), B_BOOT);
    rowsD = [rowsD; diffRow_local("Hickman_HALO - Hickman_PLAC (paired)", EST(e), mH - mP, ...
        clusterStat_local(bySubj_local(data(iH).subj, data(iH).vals(:, e), ids)) - clusterStat_local(bySubj_local(data(iP).subj, data(iP).vals(:, e), ids)))]; %#ok<AGROW>
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
    x = vs(randi(ns, ns, 1));
    if isempty(vertcat(x{:})), error("dsUnc:empty", "%s", "Bootstrap resample with no finite beta_gen* in " + name); end
    m(b) = clusterStat_local(x);
end
end

function m = clusterStat_local(vs)
% Median of per-subject medians; subjects without a finite value are skipped. Error if none has one (fail loud).
sm = cellfun(@(x) median(x, "omitnan"), vs);  sm(cellfun(@isempty, vs)) = NaN;
if all(isnan(sm)), error("dsUnc:noSubject", "No subject has a finite beta_gen*"); end
m = median(sm, "omitnan");
end

function [mA, mB] = bootMedianPair_local(vsA, vsB, B)
% Shared subject resampling for two arms measured on the same subjects.
ns = numel(vsA);  mA = zeros(B, 1);  mB = mA;
for b = 1:B
    ix = randi(ns, ns, 1);  xa = vsA(ix);  xb = vsB(ix);
    if isempty(vertcat(xa{:})) || isempty(vertcat(xb{:})), error("dsUnc:emptyPair", "%s", "Paired resample with no finite beta_gen* in one arm"); end
    mA(b) = clusterStat_local(xa);  mB(b) = clusterStat_local(xb);
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
