function T = buildTable6bV015_v001()
% Table 6b (validity gradient) on both bases, with a regression anchor.
%   N_REPS = 20  : loopClosureVarDecomp_v012_gated_all6_legacy.mat + loopClosureResults_<ds>_all_shaped_xu_{v007,v008}.mat
%                  -> must reproduce the printed Table 6b (flat, cluster [t CI], mean+-SD, N_fin) or the script errors.
%   N_REPS = 200 : loopClosureVarDecomp_v012_gated_all6_v015.mat + loopClosureResults_<ds>_all_shaped_xu_v015.mat (new Table 6b).
% flat median, cluster median and t(M-1) CI come from the varDecomp files (out(k).flat_median / cluster_median);
% mean +- SD are over the trials with finite betaGenStarMed (row-wise NaN-omitting median of the six pipelines).
% Saves results/buildTable6bV015_v001.mat (T) and CSV.
srcDir = fileparts(mfilename('fullpath'));  resDir = fullfile(fileparts(srcDir), 'results');
ds   = ["Fraser" "Cook ASD" "Cook CTRL" "Hickman HALO" "Hickman PLAC"];
tag  = ["Fraser" "Cook_ASD" "Cook_CTRL" "Hickman_HALO" "Hickman_PLAC"];
v0   = ["v008" "v007" "v007" "v007" "v007"];
pub  = struct('flat', [0.3292 0.294 0.284 0.299 0.298], 'cl', [0.3296 0.2904 0.2815 0.2823 0.2985], ...
              'lo', [0.317 0.233 0.209 0.257 0.283], 'hi', [0.342 0.348 0.354 0.307 0.314], ...
              'mean', [0.3424 0.276 0.283 0.293 0.303], 'sd', [0.1074 0.130 0.130 0.114 0.086], 'N', [2744 81 72 301 343]);
vd   = struct('N20', "loopClosureVarDecomp_v012_gated_all6_legacy.mat", 'N200', "loopClosureVarDecomp_v012_gated_all6_v015.mat");
rows = cell(0, 10);
for b = ["N20" "N200"]
    V = load(fullfile(srcDir, vd.(b)));
    if b == "N200", assert(V.meta.dataSource == "v015", 'tbl6b:src', 'varDecomp file is not the v015 source'); end
    for k = 1:numel(ds)
        i = find(string({V.out.label}) == ds(k));  assert(isscalar(i), 'tbl6b:label', 'dataset %s not in %s', ds(k), vd.(b));
        o = V.out(i);  cm = o.cluster_median;  fl = o.flat_median;
        ver = v0(k);  if b == "N200", ver = "v015"; end
        R = load(fullfile(srcDir, sprintf('loopClosureResults_%s_all_shaped_xu_%s.mat', tag(k), ver)), 'results', 'config').results;
        x = arrayfun(@(r) double(r.betaGenStarMed), R(:));  x = x(isfinite(x));
        rows(end+1, :) = {b, ds(k), fl.median, cm.pointEst, cm.ciT(1), cm.ciT(2), mean(x), std(x), cm.N_trials, numel(x)}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'basis','dataset','flat','cluster','ciLo','ciHi','mean','sd','Nfin','nFiniteMed'});
T.basis = string(T.basis);  T.dataset = string(T.dataset);
% regression anchor: N20 rows against the printed Table 6b
A = T(T.basis == "N20", :);  tol = @(x, p, d) max(abs(round(x(:), d) - p(:)));
err = [tol(A.flat, pub.flat', 3) tol(A.cluster, pub.cl', 4) tol(A.ciLo, pub.lo', 3) tol(A.ciHi, pub.hi', 3) tol(A.sd, pub.sd', 3) max(abs(A.Nfin - pub.N'))];
fprintf('anchor max abs diff (flat, cluster, ciLo, ciHi, sd, N): %s\n', mat2str(err, 3));
tolMean = max(abs(round(A.mean, 3) - round(pub.mean', 3)));  fprintf('anchor mean max diff %.4f\n', tolMean);
sdBad = A.dataset(abs(round(A.sd, 3) - pub.sd') > 0.0011);
if ~isempty(sdBad), warning('tbl6b:sdAnchor', 'printed Table 6b SD does not reproduce for %s (printed %s, rebuilt %s): the printed SD predates this basis and has no matching source', strjoin(sdBad, ', '), mat2str(pub.sd(ismember(ds, sdBad))), mat2str(round(A.sd(ismember(A.dataset, sdBad)), 3)')); end
assert(all(err([1 2 3 4 6]) <= [0.0011 0.00011 0.0011 0.0011 0]) && tolMean <= 0.0011, 'tbl6b:anchor', 'N20 rebuild does not reproduce the printed Table 6b');
disp(T);
save(fullfile(resDir, 'buildTable6bV015_v001.mat'), 'T');  writetable(T, fullfile(resDir, 'buildTable6bV015_v001.csv'));
end
