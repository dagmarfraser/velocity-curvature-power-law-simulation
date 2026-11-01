function T = analyseD16Checks_v001(opts)
%ANALYSED16CHECKS_V001  Saves two D16 numbers that previously existed only in logs or log text (SESSION_LOG Addendum 45), so the manuscript ledger can check them.
%   (1) PAIRING CHECK of the rep-paired run: median over trials and six profiles of corr(betaRecAll(:,p,gi)) between S240fc10 (reference) and S240fc05, gi = grid node nearest the trial's betaGenStarMed
%       (identical to pairingCheck_local in runD16PairedHPC_v001; the HPC log printed 0.658).
%   (2) RATE-ARM SPREADS: for A60p and A120p against A240 and each profile (SG-IRLS, BWFD-OLS), on trials with finite dObs and dPred (analyseD16S2_v001's rule, own-pipeline g, subset all):
%       SD(dPred), SD(dObs), n, and the OLS slope of dObs on dPred. Registered R2 text quoted these; the S2 .mat held only the CCC rows.
%   T = analyseD16Checks_v001(Dir=<folder with the Fraser v015/v017/v018 files>, default = this script's folder). Saves results/analyseD16Checks_v001.mat.   Fraser, D.S. (2026)  v001
arguments
    opts.Dir (1,1) string = ""
end
srcDir = fileparts(mfilename('fullpath'));  if opts.Dir == "", opts.Dir = string(srcDir); end
base = "loopClosureResults_Fraser_all_shaped_xu";  vec = linspace(0, 0.75, 25);  PROF = ["SG-IRLS" "BWFD-OLS"];  P = [6 1];   % Dir default = this script's folder (the RDS src on BlueBEAR); on the Mac pass Dir explicitly
% (1) pairing check
A = load_local(fullfile(opts.Dir, base + "_clip41_v018.mat"), true);  B = load_local(fullfile(opts.Dir, base + "_bw5_sg34_clip41_v018.mat"), true);
[~, ia, ib] = intersect(A.ids, B.ids, 'stable');  c = nan(numel(ia), 6);
for i = 1:numel(ia)
    gi = 13;  if isfinite(A.gMed(ia(i))), [~, gi] = min(abs(vec - A.gMed(ia(i)))); end
    for p = 1:6, c(i, p) = corr(double(A.all(:, p, gi, ia(i))), double(B.all(:, p, gi, ib(i))), 'Rows', 'complete'); end
end
T.pairing = struct('median', median(c(:), 'omitnan'), 'nTrials', numel(ia));
fprintf('PAIRING CHECK: median rep-wise correlation fc10 vs fc05 = %.3f over %d trials\n', T.pairing.median, numel(ia));
% (2) rate-arm spreads
A = load_local(fullfile(opts.Dir, base + "_v015.mat"), false);  T.rate = struct('setting', {}, 'profile', {}, 'n', {}, 'sdPred', {}, 'sdObs', {}, 'slope', {});
fprintf('%-6s %-9s %5s | %8s %8s | %7s\n', 'sett', 'profile', 'n', 'SD dPred', 'SD dObs', 'slope');
for s = struct('name', {"A60p" "A120p"}, 'tag', {"_sub4" "_sub2"})
    B = load_local(fullfile(opts.Dir, base + s.tag + "_v017.mat"), false);  [~, ia, ib] = intersect(A.ids, B.ids, 'stable');  n = numel(ia);
    for q = 1:2
        p = P(q);  g = A.gStar(ia, p);  dObs = B.obs(ib, p) - A.obs(ia, p);  dPred = nan(n, 1);
        for i = find(isfinite(g) & g >= vec(1) & g <= vec(end))'
            dPred(i) = interp1(vec, squeeze(B.curve(p, :, ib(i))), g(i)) - interp1(vec, squeeze(A.curve(p, :, ia(i))), g(i));
        end
        use = isfinite(dObs) & isfinite(dPred);  b = [ones(nnz(use), 1) dPred(use)] \ dObs(use);
        r = struct('setting', s.name, 'profile', PROF(q), 'n', nnz(use), 'sdPred', std(dPred(use)), 'sdObs', std(dObs(use)), 'slope', b(2));
        T.rate(end + 1) = r; %#ok<AGROW>
        fprintf('%-6s %-9s %5d | %8.4f %8.4f | %7.3f\n', r.setting, r.profile, r.n, r.sdPred, r.sdObs, r.slope);
    end
end
save(fullfile(srcDir, '..', 'results', 'analyseD16Checks_v001.mat'), 'T', 'opts');  fprintf('\nSaved results/analyseD16Checks_v001.mat\n=== D16 CHECKS DONE ===\n');
end

function M = load_local(f, withAll)
    if ~isfile(f), error('d16chk:missing', '%s not found (RDS mounted?)', f), end
    R = load(f, 'results').results;  n = numel(R);
    M.ids = string({R.trialID})';  M.gStar = reshape([R.betaGenStar], 6, n)';  M.gMed = [R.betaGenStarMed]';  M.obs = reshape([R.betaObs], 6, n)';  M.curve = reshape([R.betaRecCurveMed], 6, [], n);
    if withAll, M.all = cat(4, R.betaRecAll); end
end
