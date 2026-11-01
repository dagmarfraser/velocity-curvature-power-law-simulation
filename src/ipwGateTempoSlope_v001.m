%% ipwGateTempoSlope_v001.m
% D-T6 (TODO_PreSubmission_v001; PreKnife T6): how much of the residual beta_gen* tempo slope could the
% invertibility gate's selection on tempo produce? (v004 §2.5a: the movement is a collider of tempo and
% beta_gen, so selecting trials by the gate can open a tempo-beta_gen association.)
% Inverse-probability weighting: model P(passing the gate) on the five-dataset domain, weight each analysed
% trial by its stabilised inverse probability, and refit
%   (1) the pooled within/between tempo model of checkTempoDecompositionDomain_v001 (V1, Dhieb out):
%       betaGen ~ 1 + lfw + lfb + lsa + al + datasetC + (1|sessionC), per profile;
%   (2) the Cook ASD-CTRL beta_gen* contrast: all trials, and matched on log2 f0 at 0.5 doublings
%       (powerGroupContrast_v001 / checkCookMatchBins_v001), match weights times inverse-probability weights.
% Propensity models (both reported): P1 the published selection GLMM (invertible ~ lfw + lfb + lsa + al +
% datasetC + (1|sessionC), binomial, Laplace); P2 = P1 + lfw^2 + lfb^2 (Part 2 found a quadratic tempo term).
% Probabilities include the session effects (conditional). Weights: stabilised (dataset pass rate / p), on
% the analysed sample (finite beta_gen*, as published); untrimmed, and trimmed at the 1st/99th percentile.
% LIMITS (stated with the result): passing also depends on beta_obs, which carries the exponent itself, so
% weighting on tempo and noise corrects only the selection they explain; fitlme weights rescale the residual
% variance, so the standard errors are model-based, not design-based.
% Anchors (fail loud): the unweighted refits reproduce the saved V1 lfw slopes (checkTempoDecompositionDomain_v001,
% 4 d.p.) and the unweighted Cook contrasts reproduce powerGroupContrast_v001 (4 d.p.).
% Reads results/tempoDecomposition_v002.mat, checkTempoDecompositionDomain_v001.mat, powerGroupContrast_v001.mat.
% Writes results/ipwGateTempoSlope_v001.mat and a dated log in logs/. Run from src/.

%% CONFIG
ROOT    = fileparts(fileparts(mfilename("fullpath")));
IN_T    = fullfile(ROOT, "results", "tempoDecomposition_v002.mat");
IN_D    = fullfile(ROOT, "results", "checkTempoDecompositionDomain_v001.mat");
IN_P    = fullfile(ROOT, "results", "powerGroupContrast_v001.mat");
OUT_MAT = fullfile(ROOT, "results", "ipwGateTempoSlope_v001.mat");
OUT_DOM = "Dhieb";  PIPES = ["SG-IRLS" "BWFD-OLS"];  TOL = 1e-4;  TRIM = [1 99];
PROP = struct("P1", "invertible ~ 1 + lfw + lfb + lsa + al + datasetC + (1|sessionC)", ...
              "P2", "invertible ~ 1 + lfw + lfb + lfw2 + lfb2 + lsa + al + datasetC + (1|sessionC)");
FTEMPO = "betaGen ~ 1 + lfw + lfb + lsa + al + datasetC + (1|sessionC)";
FCOOK  = "betaGen ~ 1 + group + (1|personC)";

%% Setup
logFile = fullfile(ROOT, "logs", "ipwGateTempoSlope_v001_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log");
LG = fopen(logFile, "w");  if LG < 0, error("ipw:log", "%s", "Cannot open log: " + logFile); end
say_local(LG, "ipwGateTempoSlope_v001 run %s\n", string(datetime("now")));
for f = [IN_T IN_D IN_P], if ~isfile(f), error("ipw:input", "%s", "Missing: " + f); end, end
T0 = load(IN_T, "T").T;  D1 = load(IN_D, "D1").D1{2};  Rp = load(IN_P, "R").R;

%% Domain table, built exactly as checkTempoDecompositionDomain_v001 (V1)
T = T0(T0.dataset ~= OUT_DOM, :);
G = groupsummary(T(T.pipeline == PIPES(1), :), "session", "mean", "lf");
T.lfb = T.lfMean - mean(G.mean_lf);
z = @(x) (x - mean(x, "omitnan")) ./ std(x, "omitnan");
T.lsa = z(log2(T.sigma ./ T.a));  T.al = z(T.alpha);
T.datasetC = removecats(categorical(string(T.datasetC)));  T.sessionC = removecats(T.sessionC);
T.lfw2 = T.lfw .^ 2;  T.lfb2 = T.lfb .^ 2;  T.group = double(T.dataset == "Cook_ASD");

%% Propensity, weights and refits
S = table();  C = table();  W = table();
for p = PIPES
    V = T(T.pipeline == p, :);
    for pk = ["none" "P1" "P2"]
        if pk == "none", w = ones(height(V), 1); wt = w;
        else
            ms = fitglme(V, PROP.(pk), "Distribution", "Binomial", "FitMethod", "Laplace");
            ph = predict(ms);                                         % conditional on the session effects
            [~, ~, di] = unique(V.dataset);  rate = accumarray(di, double(V.invertible), [], @mean);
            w = rate(di) ./ ph;  a = isfinite(V.betaGen);
            lim = prctile(w(a), TRIM);  wt = min(max(w, lim(1)), lim(2));
            W = [W; table(p, pk, min(ph(a)), max(w(a)), prctile(w(a), 99), sum(w(a)) ^ 2 / sum(w(a) .^ 2), sum(a), ...
                'VariableNames', ["pipeline" "propensity" "minP" "maxW" "p99W" "ESS" "nAnalysed"])]; %#ok<AGROW>
        end
        for tr = ["untrimmed" "trimmed"]
            if pk == "none" && tr == "trimmed", continue, end
            ww = w; if tr == "trimmed", ww = wt; end
            A = V(isfinite(V.betaGen), :);  wa = ww(isfinite(V.betaGen));
            m = fitlme(A, FTEMPO, "Weights", wa);
            c = m.Coefficients;  k = strcmp(c.Name, "lfw");  kb = strcmp(c.Name, "lfb");
            S = [S; table(p, pk, tr, height(A), c.Estimate(k), c.SE(k), c.Lower(k), c.Upper(k), c.Estimate(kb), c.SE(kb), ...
                'VariableNames', ["pipeline" "propensity" "trim" "nTrials" "lfw" "lfwSE" "lfwLo" "lfwHi" "lfb" "lfbSE"])]; %#ok<AGROW>
            % Cook contrast: all trials and matched at 0.5 doublings
            K = A(ismember(A.dataset, ["Cook_CTRL" "Cook_ASD"]), :);  wk = wa(ismember(A.dataset, ["Cook_CTRL" "Cook_ASD"]));
            for an = ["all" "matched0.5"]
                [Km, wm, keep] = match_local(K, an);
                mc = fitlme(Km, FCOOK, "Weights", wm .* wk(keep));
                cc = mc.Coefficients;  kg = strcmp(cc.Name, "group");
                C = [C; table(p, pk, tr, an, height(Km), cc.Estimate(kg), cc.Lower(kg), cc.Upper(kg), ...
                    'VariableNames', ["pipeline" "propensity" "trim" "analysis" "nTrials" "b" "CIlo" "CIhi"])]; %#ok<AGROW>
            end
        end
    end
end

%% Anchors: unweighted refits reproduce the published values
for p = PIPES
    got = round(S.lfw(S.pipeline == p & S.propensity == "none"), 4);
    want = D1.b(D1.pipeline == p & D1.outcome == "beta_gen*" & D1.term == "lfw");
    if abs(got - want) > TOL, error("ipw:anchorSlope", "%s", sprintf("%s unweighted lfw %.4f, saved %.4f", p, got, want)); end
    for an = ["all" "matched0.5"]
        got = round(C.b(C.pipeline == p & C.propensity == "none" & C.analysis == an), 4);
        want = round(Rp.b(Rp.pipeline == p & Rp.outcome == "betaGen" & Rp.analysis == an), 4);
        if abs(got - want) > TOL, error("ipw:anchorCook", "%s", sprintf("%s Cook %s unweighted %.4f, published %.4f", p, an, got, want)); end
    end
end
say_local(LG, "Anchors passed: unweighted refits reproduce the saved V1 slopes and the published Cook contrasts.\n");

%% Report
rnd = @(X, v) varfun(@(x) round(x, 4), X, "InputVariables", v);
say_local(LG, "\nWeights on the analysed trials (stabilised by dataset pass rate):\n%s", evalc("disp(W)"));
S.shift = S.lfw - S.lfw(1);  for p = PIPES, k = S.pipeline == p; S.shift(k) = S.lfw(k) - S.lfw(k & S.propensity == "none"); end
say_local(LG, "\nWithin-session beta_gen* tempo slope (lfw, per doubling), unweighted and inverse-probability weighted:\n%s", evalc("disp(S)"));
say_local(LG, "\nCook ASD minus CTRL, beta_gen*:\n%s", evalc("disp(C)"));
save(OUT_MAT, "S", "C", "W", "PROP", "TRIM");
say_local(LG, "Saved: %s\nLog: %s\n", OUT_MAT, logFile);  fclose(LG);

%% =========================================================================
function [K, w, keep] = match_local(K, an)
% As powerGroupContrast_v001 / checkCookMatchBins_v001 (L32-52); "all": unit weights. keep indexes the input rows.
    keep = (1:height(K))';
    if an == "all", w = ones(height(K), 1); return, end
    bw = str2double(extractAfter(an, "matched"));
    K.bin = floor(K.lf / bw);
    cnt = groupsummary(K, ["bin" "group"]);
    ok = intersect(cnt.bin(cnt.group == 0), cnt.bin(cnt.group == 1));
    k = ismember(K.bin, ok);  K = K(k, :);  keep = keep(k);
    w = zeros(height(K), 1);
    for gq = 0:1
        Kg = K.group == gq;  nG = sum(Kg);
        for b = ok(:)'
            kk = Kg & K.bin == b;  w(kk) = (sum(K.bin == b) / height(K)) / (sum(kk) / nG);
        end
    end
end

function say_local(fid, fmt, varargin)
    fprintf(fmt, varargin{:});  fprintf(fid, fmt, varargin{:});
end
