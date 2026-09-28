%% checkLMMTempoConfound_v002.m
% P-1 (docs/TODO_v003_Rewrite_v001.md): is the Stage 1 LMM's reading of beta_gen as the
% dominant bias predictor a tempo effect? On the grid, realised tempo is
% f0 = VGF / I(beta), I(beta) = closed integral of kappa^beta ds over the grid path
% (#224, #236), so tempo is a function of (beta_gen, VGF) and cannot enter the full
% factorial as an extra covariate. We compare parameterisations of the same cells:
%   M0  VGF (z)      reproduction gate against the fitted L9 model's 192 fixed effects
%   M1  log VGF (z)  same information on a log scale (isolates the transform)
%   M2  log f0 (z)   beta_gen's effect is now at fixed tempo, not at fixed VGF
% Route: OLS on per-cell mean bias (perCoordinateSEM_v2_001.mat, 3.49M cells) instead of
% fitlme on 17.5M rows (ModelAdequacy_Stage1_KitchenSink_v2_001 L451, L958). The pooled
% factorial with filterType x regressionType equals six per-pipeline 5-way factorials
% (same fitted values); pooled reference-coded coefficients are rebuilt by differences
% for the gate. Scaling as the LMM (L925-939): beta_gen, VGF-term, sigma, alpha z-scored
% with row-level (nReps-weighted) moments; samplingRate raw. If the gate fails, escalate
% to a BlueBEAR fitlme refit; nothing downstream runs.
% v002: gate target is results/L9ModelInspect_v002.mat (fixed effects read from the fitted
% lme; inspectL9Model_v002 verdict: object clean). v001 gated on L9_coefficients_v004.csv,
% whose interaction terms do not match the fitted model (Finding pending), and failed.
% Run from src/:  checkLMMTempoConfound_v002
% Writes results/lmmTempoConfound_v002.mat.
% Fraser, D.S. (2026)  v002

%% CONFIG
ROOT    = fileparts(fileparts(mfilename("fullpath")));
SEM_MAT = fullfile(ROOT, "src", "perCoordinateSEM_v2_001.mat");
L9_MAT  = fullfile(ROOT, "results", "L9ModelInspect_v002.mat");
OUT_MAT = fullfile(ROOT, "results", "lmmTempoConfound_v002.mat");
N_LMM   = 17458535;                       % Finding #189 (post Step 5b)
GATE_Z  = 2;  GATE_REL = 0.05;  GATE_R = 0.999;   % per term: |d| <= 2 SE or <= 5% |est|; r(est) >= .999
REF     = [2 3];                          % nominal reference levels (lowest codes), as fitlme
PARAMS  = ["VGF" "logVGF" "logF0"];

%% Load
addpath(genpath(fullfile(ROOT, "src", "functions")));
for f = [SEM_MAT L9_MAT]
    if ~isfile(f), error("p1:input", "%s", "Missing: " + f); end
end
if ~isfolder(fileparts(OUT_MAT)), error("p1:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
S = load(SEM_MAT, "coordTable");  T = S.coordTable;  nAll = height(T);
T = T(isfinite(T.meanBias), :);
D = table(double(T.betaGen), double(T.VGF), double(T.fs), double(T.sigma), double(T.alpha), ...
    double(T.meanBias), double(T.nReps), double(T.filterType), double(T.regressType), string(T.pipeline), ...
    'VariableNames', ["betaGen" "VGF" "fs" "sigma" "alpha" "bias" "n" "filt" "reg" "pipeline"]);
fprintf("Cells: %d (%d dropped for non-finite meanBias); rows represented: %d vs LMM N %d\n", ...
    height(D), nAll - height(D), sum(D.n), N_LMM);
if sum(D.n) ~= N_LMM
    fprintf("*** ROW-COUNT MISMATCH: cell table represents %d rows, LMM used %d. The gate decides. ***\n", sum(D.n), N_LMM);
end

%% Tempo per cell from the true grid path (checkGridShapeUnits_v001: median error 0.14%)
G  = gridShapeGeometry_v002(ROOT);
ds = vecnorm(diff([G.pathPx; G.pathPx(1, :)]), 2, 2);  dsV = (ds + circshift(ds, 1)) / 2;
bU = unique(D.betaGen);  I = arrayfun(@(b) sum(G.kappaPx .^ b .* dsV), bU);
[~, ib] = ismember(D.betaGen, bU);  D.f0 = D.VGF ./ I(ib);
for bt = [1/3 2/3]
    [~, k] = min(abs(bU - bt));  f = D.f0(D.betaGen == bU(k));
    fprintf("f0 at beta_gen %.4f: %.2f-%.2f Hz (#224 realised: %s)\n", bU(k), min(f), max(f), ...
        ifelse_local(bt < 0.5, "0.49-1.79", "2.76-10.71"));
end

%% Scaling (row-level moments, as the LMM)
w  = D.n;
zw = @(x) (x - sum(w .* x) / sum(w)) ./ sqrt(sum(w .* (x - sum(w .* x) / sum(w)).^2) / (sum(w) - 1));
D.betaGenerated  = zw(D.betaGen);
D.noiseMagnitude = zw(D.sigma);
D.noiseColor     = zw(D.alpha);
D.samplingRate   = D.fs;
XV = struct("VGF", zw(D.VGF), "logVGF", zw(log(D.VGF)), "logF0", zw(log(D.f0)));
for p = PARAMS, fprintf("r(beta_gen, %s) = %.3f\n", p, corr(D.betaGenerated, XV.(p))); end

%% Pipelines
[pk, ~, pIdx] = unique(D(:, ["filt" "reg"]), "rows");
pk.label = strings(height(pk), 1);
for i = 1:height(pk)
    lab = unique(D.pipeline(pIdx == i));
    if numel(lab) ~= 1, error("p1:pipelineMap", "%s", "Codes map to several labels"); end
    pk.label(i) = lab;
end
if height(pk) ~= 6, error("p1:pipelines", "%s", sprintf("Expected 6 pipelines, found %d", height(pk))); end
disp(pk);

%% Fits: per parameterisation, per pipeline; full, drop-beta, drop-X
FORMS = struct("full", "bias ~ betaGenerated*X*samplingRate*noiseMagnitude*noiseColor", ...
               "dropBeta", "bias ~ X*samplingRate*noiseMagnitude*noiseColor", ...
               "dropX", "bias ~ betaGenerated*samplingRate*noiseMagnitude*noiseColor");
tss = sum((D.bias - mean(D.bias)).^2);  N = height(D);
R = struct();
for p = PARAMS
    rss = struct("full", 0, "dropBeta", 0, "dropX", 0);  coef = cell(6, 1);
    ame = table(pk.label, zeros(6, 1), zeros(6, 1), 'VariableNames', ["pipeline" "ameBeta" "ameX"]);
    for i = 1:6
        Tp = D(pIdx == i, ["betaGenerated" "samplingRate" "noiseMagnitude" "noiseColor" "bias"]);
        Tp.X = XV.(p)(pIdx == i);
        for m = string(fieldnames(FORMS))'
            mdl = fitlm(Tp, FORMS.(m));
            if any(~isfinite(mdl.Coefficients.Estimate))
                error("p1:rank", "%s", sprintf("%s %s %s: non-finite coefficient (rank deficient)", p, pk.label(i), m));
            end
            rss.(m) = rss.(m) + mdl.SSE;
            if m == "full"
                coef{i} = table(canon_local(mdl.CoefficientNames', p), mdl.Coefficients.Estimate, ...
                    'VariableNames', ["key" "est"]);
                y0 = mdl.Fitted;  Tb = Tp;  Tb.betaGenerated = Tb.betaGenerated + 1;  Tx = Tp;  Tx.X = Tx.X + 1;
                ame.ameBeta(i) = mean(predict(mdl, Tb) - y0);  ame.ameX(i) = mean(predict(mdl, Tx) - y0);
            end
        end
    end
    k = 6 * 32;
    R.(p) = struct("R2", 1 - rss.full / tss, "AIC", N * log(rss.full / N) + 2 * k, ...
        "dR2beta", (rss.dropBeta - rss.full) / tss, "dR2X", (rss.dropX - rss.full) / tss, ...
        "ame", ame, "coef", {coef}, "pooled", pooled_local(coef, pk, REF));
end

%% Gate: M0 pooled coefficients vs L9
I9 = load(L9_MAT, "out");
if ~strcmp(I9.out.verdict, "fitted object is clean; the CSV anomaly arose downstream")
    error("p1:inspect", "%s", "L9 inspection verdict is not clean: " + I9.out.verdict);
end
L9 = table(I9.out.fixed.name, I9.out.fixed.key, I9.out.fixed.est, I9.out.fixed.SE, ...
    'VariableNames', ["Name" "key" "Estimate" "SE"]);
P0 = R.VGF.pooled;
[ok, loc] = ismember(L9.key, P0.key);
if ~all(ok)
    error("p1:gateNames", "%s", "L9 terms not reconstructed: " + strjoin(L9.Name(~ok), ", "));
end
L9.cellMean = P0.est(loc);  L9.d = L9.cellMean - L9.Estimate;  L9.zd = abs(L9.d) ./ L9.SE;
L9.pass = L9.zd <= GATE_Z | abs(L9.d) <= GATE_REL * abs(L9.Estimate);
rEst = corr(L9.cellMean, L9.Estimate);
fprintf("\nGATE: %d/%d terms pass (|d| <= %g SE or <= %g%% |est|); r(est) = %.5f; median |d|/SE = %.2f\n", ...
    sum(L9.pass), height(L9), GATE_Z, 100 * GATE_REL, rEst, median(L9.zd));
if height(P0) ~= height(L9), error("p1:gateCount", "%s", sprintf("Reconstructed %d terms, model has %d", height(P0), height(L9))); end
worst = sortrows(L9(:, ["Name" "Estimate" "cellMean" "SE" "zd" "pass"]), "zd", "descend");  disp(worst(1:8, :));
if ~all(L9.pass) || rEst < GATE_R
    save(OUT_MAT, "R", "L9", "pk");
    error("p1:gateFail", "%s", "Cell-mean route does not reproduce L9; escalate to a BlueBEAR fitlme refit. " + ...
        "Diagnostics saved to " + OUT_MAT);
end

%% Report
fprintf("\n%-7s %8s %12s %10s %10s\n", "param", "R2", "AIC", "dR2 beta", "dR2 X");
for p = PARAMS
    fprintf("%-7s %8.4f %12.0f %10.4f %10.4f\n", p, R.(p).R2, R.(p).AIC, R.(p).dR2beta, R.(p).dR2X);
end
fprintf("\nAverage marginal effect per SD (cells weighted equally), bias = beta_gen - beta_rec:\n");
A = table(pk.label, 'VariableNames', "pipeline");
for p = PARAMS
    A.("beta_" + p) = R.(p).ame.ameBeta;  A.("X_" + p) = R.(p).ame.ameX;
end
disp(A);

out = struct("R", R, "gate", L9, "gateR", rEst, "pipelines", pk, "ame", A, "f0Integral", table(bU, I), ...
    "nCells", N, "nRows", sum(D.n), "runDate", string(datetime("now")));
save(OUT_MAT, "out");
fprintf("Saved %s\n", OUT_MAT);

%% Local functions
function k = canon_local(names, xname)
% Order-free term key: tokens sorted; the parameterisation's column "X" renamed to xname.
    names = string(names);  k = strings(size(names));
    for i = 1:numel(names)
        if names(i) == "(Intercept)", k(i) = names(i); continue, end
        t = split(names(i), ":");  t(t == "X") = xname;
        k(i) = strjoin(sort(t), ":");
    end
end

function P = pooled_local(coef, pk, ref)
% Reference-coded pooled coefficients from per-pipeline fits (exact for OLS).
    g = @(f, r) coef{pk.filt == f & pk.reg == r};
    b0 = g(ref(1), ref(2));  keys = b0.key;
    val = @(f, r) g(f, r).est(match_local(keys, g(f, r).key));
    fOther = setdiff(unique(pk.filt), ref(1));  rOther = setdiff(unique(pk.reg), ref(2));
    K = keys;  E = b0.est;
    for f = fOther'
        K = [K; addTok_local(keys, "filterType_" + f)];  E = [E; val(f, ref(2)) - b0.est]; %#ok<AGROW>
    end
    for r = rOther'
        K = [K; addTok_local(keys, "regressionType_" + r)];  E = [E; val(ref(1), r) - b0.est]; %#ok<AGROW>
        for f = fOther'
            K = [K; addTok_local(keys, ["filterType_" + f "regressionType_" + r])]; %#ok<AGROW>
            E = [E; val(f, r) - val(f, ref(2)) - val(ref(1), r) + b0.est]; %#ok<AGROW>
        end
    end
    P = table(K, E, 'VariableNames', ["key" "est"]);
end

function k = addTok_local(keys, tok)
    k = strings(size(keys));
    for i = 1:numel(keys)
        t = tok(:);
        if keys(i) ~= "(Intercept)", t = [split(keys(i), ":"); t]; end %#ok<AGROW>
        k(i) = strjoin(sort(t), ":");
    end
end

function idx = match_local(want, have)
    [ok, idx] = ismember(want, have);
    if ~all(ok), error("p1:terms", "%s", "Term sets differ across pipelines"); end
end

function out = ifelse_local(c, a, b)
    if c, out = a; else, out = b; end
end
