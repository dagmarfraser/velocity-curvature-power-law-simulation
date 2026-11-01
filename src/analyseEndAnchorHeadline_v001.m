function out = analyseEndAnchorHeadline_v001(opts)
%ANALYSEENDANCHORHEADLINE_V001  v019 (end-anchored template, N_REPS 200) against v015: cccMed, valid trials, within-session tempo-slope attenuation (Sol's closure; Addendum 70).
%   out = analyseEndAnchorHeadline_v001()                 B = 2000 on 72 workers (run where the five v019 files are: the RDS src on BlueBEAR)
%   out = analyseEndAnchorHeadline_v001(SkipBoot=true)    metrics and point estimates only (minutes)
%   READING RULES, proposed before any v019 result exists (Dagmar to veto or change; thresholds are Claude's, not registered):
%     R1  cccMed: the end-anchored template is immaterial for the headline if every core dataset has |cccMed(v019) - cccMed(v015)| <= 0.03 and Cook CTRL is within 0.03 of the .776 baseline.
%         Otherwise the headline is re-quoted on v019 and the shift is reported as a sensitivity result.
%     R2  tempo attenuation (SG-IRLS, model V1): rule of spec 10.13 unchanged. A MAJORITY reduction stands if the lower 95% CI limit of redV1 > 0.50; the v019 point estimate is also compared
%         with the v015 interval (68% [51-90]; BWFD-OLS legacy 78% [49-113]).
%   Steps: (1) preflight (five v019 files, END_ANCHORED_TEMPLATE, N_REPS 200, trial IDs identical to v015); (2) constellationMetrics_v006 on v019 and on v015 (same datasets, same code, same seed), the v015 rerun
%   anchored to the saved constellationMetrics_v005 within 1e-9; (3) trial table from v019 as analyseTempoDecomposition_v002 (Dhieb absent) -> results/tempoDecomposition_v019.mat;
%   (4) bootstrapTempoReduction_v002 on that table (regression anchor skipped by design, disclosed) and comparison with bootstrapTempoReduction_v001.
%   Writes results/analyseEndAnchorHeadline_v001.mat (out) and .txt.   Fraser, D.S. (2026)  v001
arguments
    opts.Dir (1,1) string = ""
    opts.B (1,1) double {mustBePositive, mustBeInteger} = 2000
    opts.NWorkers (1,1) double {mustBeNonnegative, mustBeInteger} = 72
    opts.SkipBoot (1,1) logical = false
    opts.NPerm (1,1) double {mustBeInteger, mustBeNonnegative} = 99
end
srcDir = fileparts(mfilename('fullpath'));  root = fileparts(srcDir);  resDir = fullfile(root, 'results');  addpath(srcDir);
dataDir = opts.Dir;  if dataDir == "", dataDir = string(srcDir); end
ds = ["Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Fraser"];  PL = ["BWFD-OLS" "SG-OLS" "BWFD-LMLS" "SG-LMLS" "BWFD-IRLS" "SG-IRLS"];  PIPES = ["SG-IRLS" "BWFD-OLS"];  PAIRED = ["Hickman_PLAC" "Hickman_HALO"];
for d = {resDir, fullfile(resDir, 'tempoSim')}, if ~isfolder(d{1}), mkdir(d{1}); end, end
diary(fullfile(resDir, 'analyseEndAnchorHeadline_v001.txt'));  cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('analyseEndAnchorHeadline_v001 | %s | data %s\n', string(datetime('now')), dataDir);

%% 1 Preflight
for d = ds
    f19 = fullfile(dataDir, "loopClosureResults_" + d + "_all_shaped_xu_v019.mat");  f15 = fullfile(dataDir, "loopClosureResults_" + d + "_all_shaped_xu_v015.mat");
    for f = [f19 f15], if ~isfile(f), error('endAnchor:missing', 'Missing %s', f); end, end
    A = load(f19, 'config', 'results', 'pipelineLabels');  B = load(f15, 'results', 'pipelineLabels');
    if ~isfield(A.config, 'END_ANCHORED_TEMPLATE') || ~A.config.END_ANCHORED_TEMPLATE, error('endAnchor:flag', '%s is not an end-anchored run', f19); end
    if A.config.N_REPS ~= 200 || A.config.N_BETA ~= 25 || A.config.RNG_SEED ~= 1729, error('endAnchor:config', '%s: N_REPS/N_BETA/seed differ from 200/25/1729', f19); end
    if ~isequal(string({A.results.trialID}), string({B.results.trialID})), error('endAnchor:trials', '%s: trial IDs differ between v019 and v015', d); end
    if ~isequal(string(A.pipelineLabels(:))', PL), error('endAnchor:labels', '%s: pipeline order differs', d); end
    fprintf('  OK %-13s %5d trials, end-anchored, N_REPS 200\n', d, numel(A.results));
end

%% 2 Constellation metrics, v019 and v015 (same code)
[T19, PT19] = constellationMetrics_v006(Dir=dataDir, Datasets=ds, ShapedVer="v019", OutTag="_v019", NPerm=opts.NPerm);
[T15, PT15] = constellationMetrics_v006(Dir=dataDir, Datasets=ds, ShapedVer="v015", OutTag="_v015core", NPerm=opts.NPerm);
% regression anchor: saved constellationMetrics_v005 cccMed (shaped_xu, v015), embedded because results/ is git-excluded and absent on the RDS (the first HPC run failed here, 2026-10-29)
anchor = containers.Map(cellstr(["Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Fraser"]), num2cell([0.776131107206 0.68712559954 0.225184410568 0.284024016244 0.888405946777]));
for d = ds
    a = T15.cccMed(T15.dataset == d);  b = anchor(char(d));
    if abs(a - b) > 1e-9, error('endAnchor:anchor', 'v015 rerun cccMed for %s (%.9f) differs from the saved constellationMetrics_v005 value (%.9f)', d, a, b); end
end
fprintf('  v015 rerun reproduces constellationMetrics_v005 cccMed to 1e-9 for all five datasets\n');
C = table(ds', T15.Nvalid, T19.Nvalid, T15.cccMed, T19.cccMed, T19.cccMed - T15.cccMed, T15.nullCccMed, T19.nullCccMed, T15.medBetaGenStar, T19.medBetaGenStar, T15.loopCCC, T19.loopCCC, ...
    'VariableNames', {'dataset','Nvalid15','Nvalid19','ccc15','ccc19','dCcc','null15','null19','medGen15','medGen19','loopCCC15','loopCCC19'});
fprintf('\nTable E1: cccMed, v015 vs v019 (end-anchored)\n');  disp(C);
base = 0.776;  ccIdx = find(ds == "Cook_CTRL");
R1 = all(abs(C.dCcc) <= 0.03) && abs(C.ccc19(ccIdx) - base) <= 0.03;
fprintf('R1 (max |dCcc| %.3f; Cook CTRL v019 %.3f vs baseline %.3f): %s\n', max(abs(C.dCcc)), C.ccc19(ccIdx), base, iff_local(R1, 'immaterial for the headline', 'NOT immaterial: re-quote the headline on v019 and report the shift'));
fprintf('Pooled in-domain cccMed: v015 %.3f -> v019 %.3f\n', PT15.cccMed(PT15.scope == "inDomain"), PT19.cccMed(PT19.scope == "inDomain"));

%% 3 Tempo trial table from v019
rows = {};
for d = ds
    R = load(fullfile(dataDir, "loopClosureResults_" + d + "_all_shaped_xu_v019.mat"), 'results').results;  person = d;  if ismember(d, PAIRED), person = "Hickman"; end
    for i = 1:numel(R)
        r = R(i);
        for p = PIPES
            k = find(PL == p);
            rows{end+1, 1} = struct("dataset", d, "session", d + "_" + string(r.subjectID), "person", person + "_" + string(r.subjectID), "pipeline", p, "f0", r.f0, "sigma", r.sigmaMM, ...
                "a", r.a_mm, "alpha", r.alphaIRA, "betaObs", r.betaObs(k), "betaGen", r.betaGenStar(k), "invertible", r.invertStatus(k) == "rise"); %#ok<AGROW>
        end
    end
end
T = struct2table(vertcat(rows{:}));  T = T(isfinite(T.f0) & T.f0 > 0 & isfinite(T.sigma) & isfinite(T.a) & T.a > 0 & isfinite(T.alpha), :);
T.lf = log2(T.f0);  G = groupsummary(T(T.pipeline == PIPES(1), :), "session", "mean", "lf");  [~, loc] = ismember(T.session, G.session);  T.lfMean = G.mean_lf(loc);
T.lfw = T.lf - T.lfMean;  T.lfb = T.lfMean - mean(G.mean_lf);  zf = @(x) (x - mean(x, "omitnan")) ./ std(x, "omitnan");  T.lsa = zf(log2(T.sigma ./ T.a));  T.al = zf(T.alpha);
T.datasetC = categorical(T.dataset);  T.sessionC = categorical(T.session);  T.personC = categorical(T.person);
inMat = fullfile(resDir, 'tempoDecomposition_v019.mat');  save(inMat, 'T', 'PIPES', 'PL', '-v7.3');
fprintf('\nTempo table v019: %d trial x pipeline rows, %d sessions -> %s\n', height(T), numel(unique(T.session)), inMat);

%% 4 Tempo-slope attenuation
out = struct('C', C, 'PT15', PT15, 'PT19', PT19, 'R1', R1, 'runDate', datetime('now'), 'opts', opts);
if opts.SkipBoot
    fprintf('SkipBoot: tempo-slope bootstrap not run. Re-run without SkipBoot on the HPC.\n');
else
    b19 = bootstrapTempoReduction_v002(B=opts.B, NWorkers=opts.NWorkers, InMat=inMat, Tag="_v019");
    f1 = fullfile(resDir, 'tempoSim', 'bootstrapTempoReduction_v001.mat');
    if ~isfile(f1), error('endAnchor:noV001', 'Missing %s: cannot compare with the v015 interval', f1); end
    b15 = load(f1, 'out').out;
    E = table();
    for p = PIPES
        for n = ["bObsV1" "bGenV1" "redV1" "redV3"]
            g = @(o) o.R(o.R.pipeline == p & o.R.estimand == n, :);  a = g(b15);  b = g(b19);
            E = [E; table(p, n, a.estimate, a.ciLo, a.ciHi, b.estimate, b.ciLo, b.ciHi, b.estimate - a.estimate, 'VariableNames', {'pipeline','estimand','est15','lo15','hi15','est19','lo19','hi19','delta'})]; %#ok<AGROW>
    end, end
    fprintf('\nTable E2: within-session tempo slopes and attenuation, v015 (published bootstrap) vs v019\n');  disp(E);
    r = E(E.pipeline == "SG-IRLS" & E.estimand == "redV1", :);
    fprintf('R2 (SG-IRLS redV1 v019 %.3f [%.3f, %.3f]; v015 %.3f [%.3f, %.3f]): lower limit %s 0.50 -> %s; v019 point estimate %s the v015 interval\n', r.est19, r.lo19, r.hi19, r.est15, r.lo15, r.hi15, ...
        iff_local(r.lo19 > 0.5, '>', '<='), iff_local(r.lo19 > 0.5, 'a MAJORITY reduction is still supported', 'a majority is NOT supported'), iff_local(r.est19 >= r.lo15 && r.est19 <= r.hi15, 'lies inside', 'lies OUTSIDE'));
    out.E = E;  out.b19 = b19;  out.R2 = r.lo19 > 0.5;
end
save(fullfile(resDir, 'analyseEndAnchorHeadline_v001.mat'), 'out');
fprintf('Saved results/analyseEndAnchorHeadline_v001.mat\n=== COMPLETED OK ===\n');
end

function s = iff_local(c, a, b), if c, s = a; else, s = b; end, end
