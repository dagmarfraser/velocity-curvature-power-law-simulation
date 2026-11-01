function T = analyseA60pPattern_v001(opts)
%ANALYSEA60PPATTERN_V001  Where on the beta grid is the A60p forward map NaN, and does it depend on the trial? (SESSION_LOG Addendum 40, follows analyseA60pSupport_v001)
%   For the A60p (_sub4) run, profile SG-IRLS and BWFD-OLS: share of trials NaN at each of the 25 grid points (beta 0 to 0.75), the number of trials with any NaN, the number whose own g falls on or beside a NaN point, and the
%   NaN share by tercile of the trial's cycles count (trial length x f0 is not saved, so the tercile is of the reference run's observed beta_obs as a coarse proxy is NOT used): instead by the trial's betaObs tercile and gStar tercile.
%   Descriptive; no data modified. Saves results/analyseA60pPattern_v001.mat.   Fraser, D.S. (2026)  v001
arguments
    opts.Dir (1,1) string = "/Volumes/rdsprojects/f/fraserds-mpo-evaluation/2026_prereg/velocity-curvature-power-law-simulation-main/velocity-curvature-power-law-simulation-main/src"
end
srcDir = fileparts(mfilename('fullpath'));  base = "loopClosureResults_Fraser_all_shaped_xu";  PROF = ["SG-IRLS" "BWFD-OLS"];  P = [6 1];
A = load_local(fullfile(opts.Dir, base + "_v015.mat"));  B = load_local(fullfile(opts.Dir, base + "_sub4_v017.mat"));  [~, ia, ib] = intersect(A.ids, B.ids, 'stable');  n = numel(ia);  vec = linspace(0, 0.75, size(A.curve, 2));
T = struct('profile', {}, 'nanByGrid', {}, 'anyNaN', {}, 'gOnNaN', {}, 'nanByGTercile', {}, 'nanByObsTercile', {});
for q = 1:2
    p = P(q);  cb = squeeze(B.curve(p, :, ib));  nanMask = ~isfinite(cb);  g = A.gStar(ia, p);  ob = B.obs(ib, p);  anyN = any(nanMask, 1)';
    onG = false(n, 1);  for i = find(isfinite(g))', [~, k] = min(abs(vec - g(i)));  lo = max(1, k - 1);  hi = min(numel(vec), k + 1);  onG(i) = any(nanMask(lo:hi, i)); end
    gt = tercile_local(g);  ot = tercile_local(ob);
    r = struct('profile', PROF(q), 'nanByGrid', 100 * mean(nanMask, 2)', 'anyNaN', nnz(anyN), 'gOnNaN', nnz(onG), 'nanByGTercile', arrayfun(@(b) 100 * mean(anyN(gt == b)), 1:3), 'nanByObsTercile', arrayfun(@(b) 100 * mean(anyN(ot == b)), 1:3));
    T(end + 1) = r; %#ok<AGROW>
    fprintf('\n%s (n %d): trials with any NaN on the grid %d; g on or beside a NaN point %d\n  %% NaN by grid point (beta 0 to 0.75 in 25 steps):', PROF(q), n, r.anyNaN, r.gOnNaN);  fprintf(' %.0f', r.nanByGrid);
    fprintf('\n  any-NaN %% by gStar tercile (low|mid|high): %.0f | %.0f | %.0f ; by betaObs(A60p) tercile: %.0f | %.0f | %.0f\n', r.nanByGTercile, r.nanByObsTercile);
end
save(fullfile(srcDir, '..', 'results', 'analyseA60pPattern_v001.mat'), 'T', 'opts');  fprintf('\nSaved results/analyseA60pPattern_v001.mat\n=== A60P PATTERN DONE ===\n');
end

function b = tercile_local(v)
    e = quantile(v(isfinite(v)), [0 1/3 2/3 1]);  b = discretize(v, e);  b(~isfinite(v)) = NaN;
end

function M = load_local(f)
    if ~isfile(f), error('a60p:missing', '%s not found (RDS mounted?)', f), end
    R = load(f, 'results').results;  n = numel(R);
    M.ids = string({R.trialID})';  M.gStar = reshape([R.betaGenStar], 6, n)';  M.obs = reshape([R.betaObs], 6, n)';  M.curve = reshape([R.betaRecCurveMed], 6, [], n);
end
