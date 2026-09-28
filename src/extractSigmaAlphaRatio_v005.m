function results = extractSigmaAlphaRatio_v005(opts)
%EXTRACTSIGMAALPHARATIO_V005  Sigma vs alpha channel decomposition of beta-bias (#44), corrected.
%   Partial derivatives of predicted Delta-beta (beta_gen - beta_rec) with respect to physical
%   sigma (per mm) and alpha (per unit), from the Stage 1 L9 LMM at each dataset's operating
%   point; ratio = |dDb/dSigma| / |dDb/dAlpha| (< 1: alpha dominates).
%
%   v005 changes from v004. Operating points are v004's, unchanged: per-trial loop-closure
%   mats (loopClosureResults_<name>_all_shaped_xu_v007.mat), valid = finite betaGenStarMed,
%   alpha = mean(alphaMaj, alphaMin), betaGen = median betaGenStarMed, VGF = vgfObsMM at the
%   pipeline's column, native fs. Three corrections, each attributed by a variant:
%     A  v004 as run: L9_coefficients_v004.csv, hard-coded scaling, VGF in mm. Reproduces v004.
%     B  A with the fitted model's 192 fixed effects (results/L9ModelInspect_v002.mat).
%        The CSV's interactions do not match the model (Finding #237). A -> B isolates it.
%     C  B with (i) scaling from the data (nReps-weighted row moments of
%        perCoordinateSEM_v2_001, as checkLMMTempoConfound_v002, whose gate reproduced all
%        192 fixed effects) and (ii) VGF in the grid's units. The grid's VGF is canvas-px
%        (kappa in 1/px; Finding #236), so VGF_px = VGF_mm * pixelScale^(1 - beta_obs) per
%        trial, beta_obs being the same pipeline's fitted exponent that produced vgfObs.
%        B -> C isolates the scaling and units fix. C is the corrected result.
%   Also: paths resolve from this file's location (run from src/); several pipelines per
%   call; fail loud on missing fields (the beta_obs field is looked up, and the error lists
%   what exists); results saved to results/sigmaAlphaRatio_v005.mat.
%   v004's caveats carry over: Zarandi is Scenario A (saturation-consensus coordinate); fs
%   is asserted per dataset, not measured.
%   Dagmar Fraser, 2026.  v005

arguments
    opts.DatasetNames (1,:) string = ["Zarandi","Hickman_PLAC","Pilot","Cook_ASD","Cook_CTRL"]
    opts.Pipelines (1,:) string {mustBeMember(opts.Pipelines, ...
        ["BWFD-OLS","BWFD-LMLS","BWFD-IRLS","SG-OLS","SG-LMLS","SG-IRLS"])} = ["SG-IRLS","BWFD-OLS"]
    opts.PrintTable (1,1) logical = true
end
BETA_OBS_FIELDS = ["betaObs" "betaObsAll" "betaRec" "beta"];   % candidates, first found is used

%% Paths
SRC  = string(fileparts(mfilename("fullpath")));  ROOT = fileparts(SRC);
MODEL_MAT = fullfile(ROOT, "results", "L9ModelInspect_v002.mat");
CSV       = fullfile(SRC, "L9_coefficients_v004.csv");
SEM_MAT   = fullfile(SRC, "perCoordinateSEM_v2_001.mat");
LC_DIR    = SRC;
OUT_MAT   = fullfile(ROOT, "results", "sigmaAlphaRatio_v005.mat");
for f = [MODEL_MAT CSV SEM_MAT]
    if ~isfile(f), error("sar5:input", "%s", "Missing: " + f); end
end
if ~isfolder(fileparts(OUT_MAT)), error("sar5:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end

%% Coefficients
I = load(MODEL_MAT, "out");
if ~strcmp(I.out.verdict, "fitted object is clean; the CSV anomaly arose downstream")
    error("sar5:model", "%s", "Model inspection verdict not clean: " + I.out.verdict);
end
coef.model = table(I.out.fixed.name, I.out.fixed.est, 'VariableNames', ["Name" "Estimate"]);
if height(coef.model) ~= 192, error("sar5:nCoef", "%s", sprintf("%d coefficients, expected 192", height(coef.model))); end
coef.csv = readtable(CSV, "TextType", "string");

%% Scaling: v004 hard-coded, and from the data
scale.v004 = struct("betaGenerated", struct("mean", 0.3402, "std", 0.2072), ...
                    "VGF",            struct("mean", 184.3,  "std", 73.37), ...
                    "noiseMagnitude", struct("mean", 4.122,  "std", 5.549), ...
                    "noiseColor",     struct("mean", 3.000,  "std", 1.789));
S = load(SEM_MAT, "coordTable");  T = S.coordTable;  w = double(T.nReps);
wm = @(x) sum(w .* x) / sum(w);  ws = @(x) sqrt(sum(w .* (x - wm(x)).^2) / (sum(w) - 1));
src = struct("betaGenerated", double(T.betaGen), "VGF", double(T.VGF), ...
             "noiseMagnitude", double(T.sigma), "noiseColor", double(T.alpha));
fprintf("Scaling, data-derived (v004 hard-coded):\n");
for fn = string(fieldnames(src))'
    scale.data.(fn) = struct("mean", wm(src.(fn)), "std", ws(src.(fn)));
    fprintf("  %-15s mean %9.4f (%9.4f)  std %8.4f (%8.4f)\n", fn, scale.data.(fn).mean, ...
        scale.v004.(fn).mean, scale.data.(fn).std, scale.v004.(fn).std);
end
clear S T src
G = gridShapeGeometry_v002(ROOT);  pxs = G.pixelScale;
fprintf("pixelScale %.2f px/mm (Toolchain_caller_v058 L213, via gridShapeGeometry_v002)\n", pxs);

%% Variants x pipelines x datasets
VAR = struct("A", struct("coef", "csv",   "scale", "v004", "vgfPx", false), ...
             "B", struct("coef", "model", "scale", "v004", "vgfPx", false), ...
             "C", struct("coef", "model", "scale", "data", "vgfPx", true));
R = table();
for pl = opts.Pipelines
    [f6, r4, r5] = pipelineToDummies_local(pl);
    for d = opts.DatasetNames
        P = aggregateDataset_local(d, LC_DIR, pl, pxs, BETA_OBS_FIELDS);
        row = table(pl, d, P.betaGen, P.VGFmm, P.VGFpx, P.fs, P.sigma, P.alpha, P.nValid, ...
            'VariableNames', ["pipeline" "dataset" "betaGen" "VGFmm" "VGFpx" "fs" "sigma" "alpha" "nValid"]);
        for v = string(fieldnames(VAR))'
            c = VAR.(v);  sc = scale.(c.scale);
            vgf = ifelse_local(c.vgfPx, P.VGFpx, P.VGFmm);
            pv = dictionary("betaGenerated", (P.betaGen - sc.betaGenerated.mean) / sc.betaGenerated.std, ...
                            "VGF",           (vgf - sc.VGF.mean) / sc.VGF.std, ...
                            "samplingRate",  P.fs, ...
                            "filterType_6", f6, "regressionType_4", r4, "regressionType_5", r5, ...
                            "noiseMagnitude", (P.sigma - sc.noiseMagnitude.mean) / sc.noiseMagnitude.std, ...
                            "noiseColor",     (P.alpha - sc.noiseColor.mean) / sc.noiseColor.std);
            dS = partialDerivWrt_local(coef.(c.coef), "noiseMagnitude", pv) / sc.noiseMagnitude.std;
            dA = partialDerivWrt_local(coef.(c.coef), "noiseColor", pv) / sc.noiseColor.std;
            row.("dSig_" + v) = dS;  row.("dAlp_" + v) = dA;  row.("ratio_" + v) = abs(dS) / abs(dA);
        end
        R = [R; row]; %#ok<AGROW>
    end
end

%% Report
if opts.PrintTable
    fprintf("\n=== extractSigmaAlphaRatio_v005 ===\n");
    fprintf("ratio = |dDb/dSigma per mm| / |dDb/dAlpha per unit|; < 1: alpha dominates.\n");
    fprintf("A = v004 as run | B = A + fitted-model coefficients (#237) | C = B + data scaling + VGF in px (#236)\n\n");
    disp(R(:, ["pipeline" "dataset" "betaGen" "VGFmm" "VGFpx" "sigma" "alpha" "nValid" ...
               "ratio_A" "ratio_B" "ratio_C"]));
    for pl = opts.Pipelines
        m = R.pipeline == pl;
        fprintf("%-9s mean ratio  A %.3f | B %.3f | C %.3f;  alpha dominates at %d/%d points (C)\n", pl, ...
            mean(R.ratio_A(m)), mean(R.ratio_B(m)), mean(R.ratio_C(m)), sum(R.ratio_C(m) < 1), sum(m));
    end
end
results = struct("table", R, "scaling", scale, "pixelScale", pxs, "variants", VAR, ...
    "coefSource", MODEL_MAT, "timestamp", datetime("now"));
save(OUT_MAT, "results");
fprintf("Saved %s\n", OUT_MAT);
end

%% ------------------------------------------------------------------------
function P = aggregateDataset_local(dsName, lcDir, pipeline, pxs, betaFields)
lcFile = fullfile(lcDir, sprintf("loopClosureResults_%s_all_shaped_xu_v007.mat", dsName));
if ~isfile(lcFile), error("sar5:lc", "%s", "Missing: " + lcFile); end
S = load(lcFile, "results", "pipelineLabels");  r = S.results;
pIdx = find(strcmpi(strrep(string(S.pipelineLabels), "_", "-"), pipeline), 1);
if isempty(pIdx)
    error("sar5:pipeline", "%s", pipeline + " not in pipelineLabels of " + dsName + ": " + strjoin(string(S.pipelineLabels), ", "));
end
fns = string(fieldnames(r));
need = ["betaGenStarMed" "alphaMaj" "alphaMin" "sigmaMM" "vgfObsMM"];
if ~all(ismember(need, fns))
    error("sar5:fields", "%s", dsName + " lacks: " + strjoin(need(~ismember(need, fns)), ", ") + ...
        " (vgfObsMM: run correctVgfObsUnits_v001)");
end
bf = betaFields(ismember(betaFields, fns));
if isempty(bf)
    error("sar5:betaObs", "%s", dsName + ": no beta_obs field among " + strjoin(betaFields, ", ") + ...
        ". Fields present: " + strjoin(fns', ", ") + ". Add the right name to BETA_OBS_FIELDS.");
end
bf = bf(1);
vg = vertcat(r.vgfObsMM);  bo = vertcat(r.(bf));
if size(vg, 2) ~= numel(S.pipelineLabels) || ~isequal(size(bo), size(vg))
    error("sar5:shape", "%s", sprintf("%s: vgfObsMM %s and %s %s must both be nTrials x %d", ...
        dsName, mat2str(size(vg)), bf, mat2str(size(bo)), numel(S.pipelineLabels)));
end
bgs = [r.betaGenStarMed];  ok = isfinite(bgs);
if ~any(ok), error("sar5:noValid", "%s", dsName + ": no valid (finite betaGenStarMed) trials"); end
alpha = mean([[r.alphaMaj]; [r.alphaMin]], 1, "omitnan");
vgfMM = vg(:, pIdx)';  bObs = bo(:, pIdx)';
vgfPx = vgfMM .* pxs .^ (1 - bObs);               % #236: VGF_mm = VGF_px * pxs^(beta - 1)
P = struct("betaGen", median(bgs(ok), "omitnan"), "VGFmm", median(vgfMM(ok), "omitnan"), ...
    "VGFpx", median(vgfPx(ok), "omitnan"), "fs", nativeSamplingRate_local(dsName), ...
    "sigma", median([r(ok).sigmaMM], "omitnan"), "alpha", median(alpha(ok), "omitnan"), "nValid", sum(ok));
fprintf("  %-13s %-9s %4d valid trials; beta_obs field '%s'\n", dsName, pipeline, sum(ok), bf);
end

function fs = nativeSamplingRate_local(dsName)
% Asserted per dataset (v004; docs/EMPIRICAL_DATASETS.md), not measured here.
switch dsName
    case "Zarandi",                                            fs = 100;
    case {"Cook_CTRL", "Cook_ASD", "Hickman_PLAC", "Hickman_HALO"}, fs = 133;
    case "Pilot",                                              fs = 240;
    case "Dhieb",                                              fs = 100;
    otherwise, error("sar5:fs", "%s", "No native fs on file for " + dsName);
end
end

function deriv = partialDerivWrt_local(coeffTbl, target, pv)
% Sum of every term containing target once, times the other predictors' values.
deriv = 0;
for r = 1:height(coeffTbl)
    name = coeffTbl.Name(r);
    if name == "(Intercept)", continue, end
    parts = split(name, ":");
    n = sum(parts == target);
    if n == 0, continue, end
    if n > 1, error("sar5:term", "%s", "Term repeats " + target + ": " + name); end
    f = 1;
    for o = parts(parts ~= target)'
        if ~isKey(pv, o), error("sar5:predictor", "%s", "Unknown predictor " + o + " in " + name); end
        f = f * pv(o);
    end
    deriv = deriv + coeffTbl.Estimate(r) * f;
end
end

function [f6, r4, r5] = pipelineToDummies_local(pl)
parts = split(pl, "-");
f6 = double(parts(1) == "SG");  r4 = double(parts(2) == "LMLS");  r5 = double(parts(2) == "IRLS");
end

function out = ifelse_local(c, a, b)
if c, out = a; else, out = b; end
end
