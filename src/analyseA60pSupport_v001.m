function T = analyseA60pSupport_v001(opts)
%ANALYSEA60PSUPPORT_V001  Why does R2's predicted change lose about half the A60p trials? (SESSION_LOG Addendum 40)
%   Replays analyseD16S2_v001's support rule for each rate arm (A60p vs A240, A120p vs A240) and profile (SG-IRLS, BWFD-OLS): dPred is NaN when the reference run's own-pipeline beta_gen* (g) is NaN, outside the grid
%   0 to 0.75, or when either forward-map curve is NaN at g. Counts each cause, the share of curve grid points that are NaN, and the observed beta_obs NaN counts. Descriptive; no data modified.
%   T = analyseA60pSupport_v001(Dir=<folder with the v015/v017 Fraser files>). Saves results/analyseA60pSupport_v001.mat.   Fraser, D.S. (2026)  v001
arguments
    opts.Dir (1,1) string = "/Volumes/rdsprojects/f/fraserds-mpo-evaluation/2026_prereg/velocity-curvature-power-law-simulation-main/velocity-curvature-power-law-simulation-main/src"
end
srcDir = fileparts(mfilename('fullpath'));  base = "loopClosureResults_Fraser_all_shaped_xu";  PROF = ["SG-IRLS" "BWFD-OLS"];  P = [6 1];
ST = struct('name', {"A60p" "A120p"}, 'tag', {"_sub4" "_sub2"});
A = load_local(fullfile(opts.Dir, base + "_v015.mat"));  T = struct('setting', {}, 'profile', {}, 'n', {}, 'gNaN', {}, 'gLow', {}, 'gHigh', {}, 'refCurveNaN', {}, 'setCurveNaN', {}, 'usable', {}, 'medGridNaN', {}, 'medG', {});
fprintf('%-6s %-9s %5s | %5s %5s %5s | %6s %6s | %6s | %8s %7s\n', 'sett', 'profile', 'n', 'g NaN', 'g<0', 'g>.75', 'refNaN', 'setNaN', 'usable', 'gridNaN', 'med g');
for s = ST
    B = load_local(fullfile(opts.Dir, base + s.tag + "_v017.mat"));  [~, ia, ib] = intersect(A.ids, B.ids, 'stable');  n = numel(ia);  vec = linspace(0, 0.75, size(A.curve, 2));
    for q = 1:2
        p = P(q);  g = A.gStar(ia, p);  gNaN = ~isfinite(g);  gLow = isfinite(g) & g < vec(1);  gHigh = isfinite(g) & g > vec(end);  gin = isfinite(g) & ~gLow & ~gHigh;
        fa = nan(n, 1);  fb = fa;
        for i = find(gin)', fa(i) = interp1(vec, squeeze(A.curve(p, :, ia(i))), g(i));  fb(i) = interp1(vec, squeeze(B.curve(p, :, ib(i))), g(i)); end
        cb = squeeze(B.curve(p, :, ib));
        r = struct('setting', s.name, 'profile', PROF(q), 'n', n, 'gNaN', nnz(gNaN), 'gLow', nnz(gLow), 'gHigh', nnz(gHigh), 'refCurveNaN', nnz(gin & ~isfinite(fa)), 'setCurveNaN', nnz(gin & ~isfinite(fb)), ...
            'usable', nnz(isfinite(fa) & isfinite(fb)), 'medGridNaN', median(mean(~isfinite(cb), 1)), 'medG', median(g(isfinite(g))));
        T(end + 1) = r; %#ok<AGROW>
        fprintf('%-6s %-9s %5d | %5d %5d %5d | %6d %6d | %6d | %8.2f %7.3f\n', r.setting, r.profile, n, r.gNaN, r.gLow, r.gHigh, r.refCurveNaN, r.setCurveNaN, r.usable, r.medGridNaN, r.medG);
    end
end
save(fullfile(srcDir, '..', 'results', 'analyseA60pSupport_v001.mat'), 'T', 'opts');  fprintf('\nSaved results/analyseA60pSupport_v001.mat\n=== A60P SUPPORT DONE ===\n');
end

function M = load_local(f)
    if ~isfile(f), error('a60p:missing', '%s not found (RDS mounted?)', f), end
    R = load(f, 'results').results;  n = numel(R);
    M.ids = string({R.trialID})';  M.gStar = reshape([R.betaGenStar], 6, n)';  M.curve = reshape([R.betaRecCurveMed], 6, [], n);
end
