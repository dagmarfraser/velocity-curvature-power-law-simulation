%% checkNodeInterpolationSpec_v001  Deep Review GD11: does snapping each trial's predicted constellation to the nearest beta_gen node (spacing 0.03125) matter?
% For every trial that enters cccMed (v005's mask: finite betaObs and slice, finite betaGenStarMed, max|dObs| >= DEG, not all |dObs| < MDC) the predicted
% pipeline outputs are (a) the median of betaRecSlice, the registered nearest-node prediction (checked here to equal betaRecCurveMed at the nearest node), and
% (b) betaRecCurveMed linearly interpolated to betaGenStarMed itself. Reported per dataset: the median / 95th percentile / maximum absolute change in the 15
% predicted contrasts, the share of trial-contrast cells moving by more than 0.0108 and by more than the MDC, and the median per-trial CCC with each
% prediction, its median paired change and the number of trials that change by more than 0.01. Run from src/. Writes results/checkNodeInterpolationSpec_v001.mat.
% Fraser, D.S. (2026)  v001
ROOT = fileparts(fileparts(mfilename("fullpath")));  addpath(genpath(fullfile(ROOT, "src", "functions")));
cfg = load(fullfile(ROOT, "results", "constellationMetrics_v005.mat"), "results").results.config;  MDC = cfg.MDC;  DEG = cfg.DEG;  TOL = 0.0108;
DS = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Dhieb"];  nP = 6;  pP = [];  pQ = [];
for p = 2:nP, for q = 1:p-1, pP(end+1) = p;  pQ(end+1) = q; end, end %#ok<AGROW>
out = struct("dataset", DS, "N", NaN(1, numel(DS)), "medAbs", NaN(1, numel(DS)), "p95Abs", NaN(1, numel(DS)), "maxAbs", NaN(1, numel(DS)), "fracOverTol", NaN(1, numel(DS)), ...
    "fracOverMDC", NaN(1, numel(DS)), "cccMedNearest", NaN(1, numel(DS)), "cccMedInterp", NaN(1, numel(DS)), "medDeltaCcc", NaN(1, numel(DS)), "nMoved001", NaN(1, numel(DS)), ...
    "maxSliceVsCurve", 0);
for k = 1:numel(DS)
    f = fullfile(ROOT, "src", "loopClosureResults_" + DS(k) + "_all_shaped_xu_v015.mat");  if ~isfile(f), error("nis:file", "%s", "FAILED PATH: " + f); end
    S = load(f, "results", "betaGenVec");  r = S.results;  g = S.betaGenVec(:)';  N = numel(r);
    dO = NaN(N, numel(pP));  dN = dO;  dI = dO;  ok = false(N, 1);
    for t = 1:N
        bo = r(t).betaObs(:)';  sl = r(t).betaRecSlice;  gs = r(t).betaGenStarMed;  cm = r(t).betaRecCurveMed;
        if numel(bo) ~= nP || any(~isfinite(bo)) || isempty(sl) || size(sl, 2) ~= nP || ~isfinite(gs), continue, end
        bn = median(sl, 1, "omitnan");  if any(~isfinite(bn)), continue, end
        [~, j] = min(abs(g - gs));  out.maxSliceVsCurve = max(out.maxSliceVsCurve, max(abs(bn - cm(:, j)')));
        bi = NaN(1, nP);  for p = 1:nP, bi(p) = interp1(g, cm(p, :), gs, "linear"); end
        if any(~isfinite(bi)), continue, end
        dO(t, :) = bo(pP) - bo(pQ);  dN(t, :) = bn(pP) - bn(pQ);  dI(t, :) = bi(pP) - bi(pQ);  ok(t) = true;
    end
    vi = find(ok & ~(max(abs(dO), [], 2) < DEG) & ~all(abs(dO) < MDC, 2));  n = numel(vi);
    dc = abs(dI(vi, :) - dN(vi, :));  cN = NaN(n, 1);  cI = cN;
    for i = 1:n, cN(i) = linCCC_v001(dN(vi(i), :)', dO(vi(i), :)').ccc;  cI(i) = linCCC_v001(dI(vi(i), :)', dO(vi(i), :)').ccc; end
    out.N(k) = n;  out.medAbs(k) = median(dc(:));  out.p95Abs(k) = prctile(dc(:), 95);  out.maxAbs(k) = max(dc(:));  out.fracOverTol(k) = mean(dc(:) > TOL);  out.fracOverMDC(k) = mean(dc(:) > MDC);
    out.cccMedNearest(k) = median(cN, "omitnan");  out.cccMedInterp(k) = median(cI, "omitnan");  out.medDeltaCcc(k) = median(cI - cN, "omitnan");  out.nMoved001(k) = sum(abs(cI - cN) > 0.01);
end
if out.maxSliceVsCurve > 1e-12, error("nis:slice", "median(betaRecSlice) differs from betaRecCurveMed at the nearest node by %.3g", out.maxSliceVsCurve); end
disp(struct2table(rmfield(out, "maxSliceVsCurve")));
save(fullfile(ROOT, "results", "checkNodeInterpolationSpec_v001.mat"), "-struct", "out");
