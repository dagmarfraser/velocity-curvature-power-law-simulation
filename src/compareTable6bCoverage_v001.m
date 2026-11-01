%% compareTable6bCoverage_v001.m
% Table 6b's dataset intervals (cluster median of beta_gen*, t(M-1) CI from the stage-2 SE, inversion draw included) on the
% original template (v015) and the end-anchored template (v019), and whether each interval includes one-third.
% Reads src/loopClosureVarDecomp_v012_gated_all6_v015.mat and src/loopClosureVarDecomp_v013_gated_all6_v019.mat (same
% machinery, NBoot 2000, seed 20260619); the v015 rows must reproduce results/buildTable6bV015_v001.mat (N200) or it errors.
% Run from src/. Writes results/compareTable6bCoverage_v001.mat and .csv; log logs/compareTable6bCoverage_v001_<date>.log.
% Fraser, D.S. (2026)  v001
srcDir = fileparts(mfilename("fullpath"));  ROOT = fileparts(srcDir);  R = fullfile(ROOT, "results");  L = fullfile(ROOT, "logs");
if ~isfolder(L), mkdir(L); end
logFile = fullfile(L, "compareTable6bCoverage_v001_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log");
diary(logFile);  cleanupDiary = onCleanup(@() diary("off"));
ds  = ["Fraser" "Cook CTRL" "Cook ASD" "Hickman PLAC" "Hickman HALO"];
src = struct("v015", "loopClosureVarDecomp_v012_gated_all6_v015.mat", "v019", "loopClosureVarDecomp_v013_gated_all6_v019.mat");
rows = cell(0, 9);
for b = ["v015" "v019"]
    V = load(fullfile(srcDir, src.(b)));
    assert(V.meta.dataSource == b, "t6cov:src", "%s has dataSource %s, expected %s", src.(b), V.meta.dataSource, b);
    for k = 1:numel(ds)
        i = find(string({V.out.label}) == ds(k));  assert(isscalar(i), "t6cov:label", "dataset %s not in %s", ds(k), src.(b));
        cm = V.out(i).cluster_median;  ci = cm.ciT;
        rows(end+1, :) = {b, ds(k), cm.pointEst, ci(1), ci(2), cm.M, cm.N_trials, ci(1) <= 1/3 && ci(2) >= 1/3, 1/3 - ci(2)}; %#ok<AGROW>
    end
end
T = cell2table(rows, "VariableNames", ["template" "dataset" "clusterMedian" "ciLo" "ciHi" "M" "N_trials" "includesOneThird" "oneThirdMinusCiHi"]);
T.template = string(T.template);  T.dataset = string(T.dataset);
% regression anchor: the v015 rows are Table 6b's N200 column
B6 = load(fullfile(R, "buildTable6bV015_v001.mat")).T;  B6 = B6(B6.basis == "N200", :);
for k = 1:numel(ds)
    a = T(T.template == "v015" & T.dataset == ds(k), :);  p = B6(B6.dataset == ds(k), :);
    d = max(abs([a.clusterMedian - p.cluster, a.ciLo - p.ciLo, a.ciHi - p.ciHi]));
    if d > 1e-12, error("t6cov:anchor", "v015 row for %s does not reproduce Table 6b (max abs diff %.3g)", ds(k), d); end
end
fprintf("Anchor: v015 rows reproduce Table 6b (N_REPS = 200) exactly.\n\n");
disp(T);
for b = ["v015" "v019"]
    t = T(T.template == b, :);
    fprintf("%s: interval includes one-third in %d of %d (%s)\n", b, nnz(t.includesOneThird), height(t), strjoin(t.dataset(t.includesOneThird), ", "));
end
save(fullfile(R, "compareTable6bCoverage_v001.mat"), "T", "src");  writetable(T, fullfile(R, "compareTable6bCoverage_v001.csv"));
fprintf("Saved results/compareTable6bCoverage_v001.mat\n=== compareTable6bCoverage_v001 DONE ===\n");
