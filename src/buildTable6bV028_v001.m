function T = buildTable6bV028_v001(opts)
%BUILDTABLE6BV028_V001 Table 6b (validity gradient) on the v028 end-matched corpus, with the frozen one-third band rule (2026-11-16).
%   Fork of buildTable6bV015_v001 (its N_REPS = 20 block and printed-table anchor are dropped; the N200 row is kept as the regression anchor).
%   Basis "v015": loopClosureVarDecomp_v014_gated_all6_v015.mat + loopClosureResults_<ds>_all_shaped_xu_v015.mat; must reproduce the N200 rows of
%                 results/buildTable6bV015_v001.mat to 1e-12 or the script errors.
%   Basis "v028": loopClosureVarDecomp_v014_gated_all6_v028[_<arm>_<basis>].mat + loadV028_v001 (Arm, ObsBasis as given).
%   flat median, cluster median and the t(M-1) interval come from the varDecomp file (out(k).flat_median / cluster_median; the interval includes the
%   inversion draw); mean and SD are over trials with a finite betaGenStarMed (row-wise NaN-omitting median of the six pipelines).
%   Band rule (docs/V028_PreCommitment_FROZEN_20261116_v001.md, section 6.2): band = 1/3 +/- MDC (0.03); status "indistinguishable at the
%   registered resolution" if the interval lies within the band, "different" if wholly outside it, otherwise "inconclusive". Whether the interval
%   includes 1/3 is also a column (reported, not a headline count). Rows: the five core datasets (Zarandi and Dhieb are reported separately).
%   Saves results/buildTable6bV028_v001_<basis>.mat and .csv.
arguments
    opts.Basis (1,1) string {mustBeMember(opts.Basis, ["v015", "v028"])} = "v028"
    opts.Arm (1,1) string {mustBeMember(opts.Arm, ["primary", "budget1", "shapeXu", "harmonicsFc"])} = "primary"
    opts.ObsBasis (1,1) string {mustBeMember(opts.ObsBasis, ["window", "whole"])} = "window"
end
srcDir = fileparts(mfilename('fullpath'));  resDir = fullfile(fileparts(srcDir), 'results');  addpath(genpath(fullfile(srcDir, 'functions')));
ds  = ["Fraser" "Cook ASD" "Cook CTRL" "Hickman HALO" "Hickman PLAC"];
tag = ["Fraser" "Cook_ASD" "Cook_CTRL" "Hickman_HALO" "Hickman_PLAC"];
MDC = 0.03;  third = 1/3;  lo0 = third - MDC;  hi0 = third + MDC;
suffix = opts.Basis;
if opts.Basis == "v028" && (opts.Arm ~= "primary" || opts.ObsBasis ~= "window"), suffix = "v028_" + opts.Arm + "_" + opts.ObsBasis; end
vdFile = fullfile(srcDir, "loopClosureVarDecomp_v014_gated_all6_" + suffix + ".mat");
assert(isfile(vdFile), 'tbl6bV028:src', 'varDecomp file missing: %s (run loopClosureVarDecomp_v014_gated with DataSource ''%s'')', vdFile, opts.Basis);
V = load(vdFile);  assert(V.meta.dataSource == opts.Basis, 'tbl6bV028:src', 'varDecomp file is not the %s source', opts.Basis);
rows = cell(0, 13);
for k = 1:numel(ds)
    i = find(string({V.out.label}) == ds(k));  assert(isscalar(i), 'tbl6bV028:label', 'dataset %s not in %s', ds(k), vdFile);
    cm = V.out(i).cluster_median;  fl = V.out(i).flat_median;
    if opts.Basis == "v015"
        R = load(fullfile(srcDir, sprintf('loopClosureResults_%s_all_shaped_xu_v015.mat', tag(k))), 'results').results;
    else
        R = loadV028_v001(ds(k), Arm=opts.Arm, ObsBasis=opts.ObsBasis);
    end
    x = arrayfun(@(r) double(r.betaGenStarMed), R(:));  x = x(isfinite(x));
    ci = cm.ciT;  inBand = ci(1) >= lo0 && ci(2) <= hi0;  outBand = ci(2) < lo0 || ci(1) > hi0;
    status = "inconclusive";  if inBand, status = "indistinguishable"; elseif outBand, status = "different"; end
    rows(end+1, :) = {ds(k), fl.median, cm.pointEst, ci(1), ci(2), mean(x), std(x), cm.N_trials, numel(x), ci(1) <= third && third <= ci(2), status, ...
        cm.pointEst - third, opts.Basis}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'dataset', 'flat', 'cluster', 'ciLo', 'ciHi', 'mean', 'sd', 'Nfin', 'nFiniteMed', 'ciIncludesOneThird', 'bandStatus', 'clusterMinusOneThird', 'basis'});
T.dataset = string(T.dataset);  T.bandStatus = string(T.bandStatus);  T.basis = string(T.basis);
if opts.Basis == "v015"   % regression anchor: the N200 rows of buildTable6bV015_v001
    P = load(fullfile(resDir, 'buildTable6bV015_v001.mat'), 'T').T;  P = P(P.basis == "N200", :);
    [~, j] = ismember(T.dataset, P.dataset);  assert(all(j > 0), 'tbl6bV028:anchor', 'anchor rows missing');
    err = max(abs([T.flat - P.flat(j), T.cluster - P.cluster(j), T.ciLo - P.ciLo(j), T.ciHi - P.ciHi(j), T.mean - P.mean(j), T.sd - P.sd(j), T.Nfin - P.Nfin(j)]), [], 'all');
    fprintf('regression anchor: max abs difference to buildTable6bV015_v001 N200 rows %.3g\n', err);
    assert(err < 1e-12, 'tbl6bV028:anchor', 'v015 rebuild does not reproduce buildTable6bV015_v001 (max diff %.3g)', err);
end
disp(T);
save(fullfile(resDir, "buildTable6bV028_v001_" + suffix + ".mat"), 'T');  writetable(T, fullfile(resDir, "buildTable6bV028_v001_" + suffix + ".csv"));
end
