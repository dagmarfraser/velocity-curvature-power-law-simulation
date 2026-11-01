%% powerGroupContrast_v001.m
% D-T7 (TODO_PreSubmission_v001; PreKnife T7): how large a Cook ASD-CTRL difference could this design detect?
% Part 4 reports no group difference in beta_gen* under any analysis, with intervals of about +-0.05 that
% cannot establish equivalence at the MDC (0.03). This gives the design's sensitivity, so the null can be read
% as "differences smaller than X would usually be missed".
% For each profile (SG-IRLS vetted, BWFD-OLS legacy), outcome (beta_gen*, beta_obs) and analysis (all trials
% unadjusted; matched on log2 f0 at 0.5 and 0.1 doublings, coarsened exact matching as checkCookMatchBins_v001,
% L32-52): fit outcome ~ 1 + group + (1|person) (weights for the matched fits), then
%   MDE80 = (t(0.975, df) + t(0.80, df)) * SE, df by Satterthwaite (about the number of subjects, not trials);
%   90% interval (Satterthwaite) for a two one-sided tests reading at +-0.03.
% Anchors (fail loud): (1) the matched fits reproduce checkCookMatchBins_v001.mat (b and the published 95% CI,
% 4 d.p.); (2) simulating from the fitted SG-IRLS beta_gen* unadjusted model with the group effect set to its
% MDE80 gives power within 0.05 of 0.80 (N_SIM draws, Satterthwaite p).
% Reads results/tempoDecomposition_v002.mat (T) and results/checkCookMatchBins_v001.mat (M).
% Writes results/powerGroupContrast_v001.mat and a dated log in logs/. Run from src/.

%% CONFIG
ROOT    = fileparts(fileparts(mfilename("fullpath")));
IN_T    = fullfile(ROOT, "results", "tempoDecomposition_v002.mat");
IN_M    = fullfile(ROOT, "results", "checkCookMatchBins_v001.mat");
OUT_MAT = fullfile(ROOT, "results", "powerGroupContrast_v001.mat");
PIPES   = ["SG-IRLS" "BWFD-OLS"];  OUTS = ["betaGen" "betaObs"];  ANALYSES = ["all" "matched0.5" "matched0.1"];
ALPHA   = 0.05;  POWER = 0.80;  EQ = 0.03;  TOL = 1e-4;
N_SIM   = 500;  TOL_SIM = 0.05;  SEED = 20261024;

%% Setup
logFile = fullfile(ROOT, "logs", "powerGroupContrast_v001_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log");
LG = fopen(logFile, "w");  if LG < 0, error("pgc:log", "%s", "Cannot open log: " + logFile); end
say_local(LG, "powerGroupContrast_v001 run %s\n", string(datetime("now")));
for f = [IN_T IN_M], if ~isfile(f), error("pgc:input", "%s", "Missing: " + f); end, end
T = load(IN_T, "T").T;  Mpub = load(IN_M, "M").M;

%% Fits
R = table();  fits = struct();
for p = PIPES
    V = T(T.pipeline == p & ismember(T.dataset, ["Cook_CTRL" "Cook_ASD"]), :);  V.group = double(V.dataset == "Cook_ASD");
    for out = OUTS
        W0 = V(isfinite(V.(out)), :);
        for an = ANALYSES
            [W, w] = match_local(W0, an);
            if an == "all", lme = fitlme(W, out + " ~ 1 + group + (1|personC)");
            else, lme = fitlme(W, out + " ~ 1 + group + (1|personC)", "Weights", w); end
            [~, ~, st] = fixedEffects(lme, "DFMethod", "satterthwaite", "Alpha", ALPHA);  k = st.Name == "group";
            [~, ~, st90] = fixedEffects(lme, "DFMethod", "satterthwaite", "Alpha", 2 * ALPHA);
            ci = coefCI(lme);  kd = strcmp(lme.CoefficientNames, "group");          % default (residual) df, as published
            b = st.Estimate(k);  se = st.SE(k);  df = st.DF(k);
            mde = (tinv(1 - ALPHA / 2, df) + tinv(POWER, df)) * se;
            lo90 = st90.Lower(k);  hi90 = st90.Upper(k);
            R = [R; table(p, out, an, height(W), numel(unique(W.personC(W.group == 0))), numel(unique(W.personC(W.group == 1))), ...
                b, ci(kd, 1), ci(kd, 2), se, df, st.Lower(k), st.Upper(k), mde, lo90, hi90, lo90 > -EQ && hi90 < EQ, ...
                'VariableNames', ["pipeline" "outcome" "analysis" "nTrials" "nSubjCTRL" "nSubjASD" "b" "CIlo_pub" "CIhi_pub" "SE" "dfSatt" ...
                "CIlo_satt" "CIhi_satt" "MDE80" "CI90lo" "CI90hi" "equivalentAtMDC"])]; %#ok<AGROW>
            fits.(matlab.lang.makeValidName(p + "_" + out + "_" + an)) = lme;
        end
    end
end

%% Anchor 1: matched fits reproduce the published table
for r = find(R.analysis ~= "all")'
    bw = str2double(extractAfter(R.analysis(r), "matched"));
    j = find(Mpub.binWidth == bw & Mpub.pipeline == R.pipeline(r) & Mpub.outcome == R.outcome(r));
    got = round([R.b(r) R.CIlo_pub(r) R.CIhi_pub(r)], 4);  want = round([Mpub.b_matched(j) Mpub.CIlo_m(j) Mpub.CIhi_m(j)], 4);
    if numel(j) ~= 1 || any(abs(got - want) > TOL)
        error("pgc:anchor1", "%s", sprintf("%s %s %s: %s against published %s", R.pipeline(r), R.outcome(r), R.analysis(r), mat2str(got), mat2str(want)));
    end
end
say_local(LG, "Anchor 1 passed: the matched fits reproduce checkCookMatchBins_v001 (%d rows).\n", sum(R.analysis ~= "all"));

%% Anchor 2: simulated power at the MDE (SG-IRLS beta_gen*, all trials)
r0 = find(R.pipeline == "SG-IRLS" & R.outcome == "betaGen" & R.analysis == "all");
lme = fits.(matlab.lang.makeValidName("SG-IRLS_betaGen_all"));  W = lme.Variables(~lme.ObservationInfo.Missing & ~lme.ObservationInfo.Excluded, :);
[psi, mse] = covarianceParameters(lme);  sdU = sqrt(psi{1});  sdE = sqrt(mse);
b0 = fixedEffects(lme);  [g, ~, gi] = unique(W.personC);  isA = accumarray(gi, W.group, [], @max);
rng(SEED);  hit = false(N_SIM, 1);
for s = 1:N_SIM
    u = sdU * randn(numel(g), 1);
    y = b0(1) + R.MDE80(r0) * W.group + u(gi) + sdE * randn(height(W), 1);
    D = table(y, W.group, W.personC, 'VariableNames', ["y" "group" "personC"]);
    m = fitlme(D, "y ~ 1 + group + (1|personC)");
    [~, ~, ss] = fixedEffects(m, "DFMethod", "satterthwaite");  hit(s) = ss.pValue(ss.Name == "group") < ALPHA;
end
pw = mean(hit);
say_local(LG, "Anchor 2: simulated power at MDE80 = %.4f (SG-IRLS beta_gen*, all trials, %d subjects: %d CTRL, %d ASD): %.3f over %d draws (target %.2f +- %.2f)\n", ...
    R.MDE80(r0), numel(g), sum(isA == 0), sum(isA == 1), pw, N_SIM, POWER, TOL_SIM);
if abs(pw - POWER) > TOL_SIM, error("pgc:anchor2", "%s", sprintf("Simulated power %.3f is not within %.2f of %.2f", pw, TOL_SIM, POWER)); end

%% Report
show = R;  vn = ["b" "CIlo_pub" "CIhi_pub" "SE" "dfSatt" "CIlo_satt" "CIhi_satt" "MDE80" "CI90lo" "CI90hi"];
for v = vn, show.(v) = round(show.(v), 4); end
say_local(LG, "\nCook ASD minus CTRL: estimate, 95%% CI (published, residual df; and Satterthwaite), MDE at 80%% power, 90%% CI for equivalence at +-%.2f:\n%s", EQ, evalc("disp(show)"));
save(OUT_MAT, "R", "pw", "N_SIM", "SEED", "ALPHA", "POWER", "EQ");
say_local(LG, "Saved: %s\nLog: %s\n", OUT_MAT, logFile);  fclose(LG);

%% =========================================================================
function [W, w] = match_local(W, an)
% Coarsened exact matching on log2 f0, exactly as checkCookMatchBins_v001 (L32-52): keep bins with both groups,
% weight each trial so each group's bin shares equal the pooled bin shares. "all": no matching, unit weights.
    if an == "all", w = ones(height(W), 1); return, end
    bw = str2double(extractAfter(an, "matched"));
    W.bin = floor(W.lf / bw);
    cnt = groupsummary(W, ["bin" "group"]);
    ok = intersect(cnt.bin(cnt.group == 0), cnt.bin(cnt.group == 1));
    W = W(ismember(W.bin, ok), :);
    w = zeros(height(W), 1);
    for gq = 0:1
        Wg = W.group == gq;  nG = sum(Wg);
        for b = ok(:)'
            k = Wg & W.bin == b;  w(k) = (sum(W.bin == b) / height(W)) / (sum(k) / nG);
        end
    end
end

function say_local(fid, fmt, varargin)
    fprintf(fmt, varargin{:});  fprintf(fid, fmt, varargin{:});
end
