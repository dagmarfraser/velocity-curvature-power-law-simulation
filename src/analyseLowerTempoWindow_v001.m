%% analyseLowerTempoWindow_v001.m
% The LOWER edge of the usable tempo window (exploratory, post hoc; not pre-registered).
% v008 states that invertibility rises with tempo (Finding #227/#245, GLMM) and gives per-dataset
% pass rates in four coarse bins, but has no lower-edge statement analogous to the upper rule
% (cutoff >= 5.5 x f0). This script derives one descriptively from the SAVED per-trial table of
% analyseTrialTempoWindow_v002.m. It reads that file only; nothing is re-simulated, differentiated
% or regressed on kinematics.
%
% Three edges, each defined on the in-domain (Finding #242) trials:
%   E1 can invert      lowest bin lower edge from which the pass rate (status == "rise") is
%                      >= 90% in that bin and every faster bin
%   E2 informative     lowest bin lower edge from which the median local map slope of invertible
%                      trials is >= 0.7 in that bin and every faster bin
%   E3 precise         lowest bin lower edge from which the median generator-scale CI width of
%                      invertible trials is <= 0.05 in that bin and every faster bin
% Bins with fewer than MIN_N trials are skipped (shown as "-") and do not break a run.
%
% Inputs:  results/trialTempoWindow_v002.mat (read only; T = 23,166 trial x pipeline rows)
% Outputs: results/lowerTempoWindow_v001.mat, logs/analyseLowerTempoWindow_v001_<stamp>.log
% Anchor:  reproduces t3 of trialTempoWindow_v002.mat (Finding #245 per-dataset SG-IRLS pass
%          rates in the bins <0.35, 0.35-0.7, 0.7-1.4, 1.4-2.8, >2.8 Hz) to 0.05 points, or errors.
% Claims:  the final block prints each number quoted to Dagmar on 2026-11-04 against the value
%          computed here (MATCH / DIFF); a DIFF is reported, not hidden.
% USAGE:   from the project root: analyseLowerTempoWindow_v001   (seconds; no corpora loaded)
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT    = fileparts(fileparts(mfilename("fullpath")));
IN_MAT  = fullfile(ROOT, "results", "trialTempoWindow_v002.mat");
OUT_MAT = fullfile(ROOT, "results", "lowerTempoWindow_v001.mat");
LOGDIR  = fullfile(ROOT, "logs");
OUT_DOMAIN = ["Dhieb" "Zarandi"];                       % Finding #242
DATASETS   = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO"];
PIPES      = ["SG-IRLS" "BWFD-OLS"];                     % vetted, legacy profile
FINE       = [0.1 0.2 0.3 0.4 0.5 0.7 1.0 1.4 Inf];      % pass-rate bins (Hz)
FINE_S     = [0.1 0.2 0.3 0.4 0.45 0.5 0.55 0.6 0.7 0.85 1.0 1.4 2.0 Inf];  % slope / CI bins
NOISE_F0   = [0 0.3 0.5 0.7 Inf];                        % tempo bins for the noise x tempo table
EXPO       = [0.3 0.5 0.7];                              % exposure thresholds (Hz)
MIN_N      = 5;
PASS_MIN   = 0.90;  SLOPE_MIN = 0.70;  CIW_MAX = 0.05;
REQ_T      = ["dataset" "subj" "f0" "sigma" "alpha" "a" "M" "pipeline" "status" "slope" "ciW"];

if ~isfile(IN_MAT), error("lowerTempo:input", "%s", "Missing: " + IN_MAT); end
if ~isfolder(LOGDIR), error("lowerTempo:logdir", "%s", "Missing folder: " + LOGDIR); end
stamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
logf  = fullfile(LOGDIR, "analyseLowerTempoWindow_v001_" + stamp + ".log");
diary(logf);  cleanupObj = onCleanup(@() diary("off"));
fprintf("analyseLowerTempoWindow_v001  run %s\ninput %s\n", stamp, IN_MAT);

%% LOAD AND CHECK
S = load(IN_MAT);
for f = ["T" "t3" "EDGES"]
    if ~isfield(S, f), error("lowerTempo:field", "%s", "Input lacks variable " + f); end
end
T = S.T;
miss = setdiff(REQ_T, string(T.Properties.VariableNames));
if ~isempty(miss), error("lowerTempo:cols", "%s", "T lacks columns: " + strjoin(miss, ", ")); end
T.dataset = string(T.dataset);  T.pipeline = string(T.pipeline);  T.status = string(T.status);
T = T(~ismember(T.dataset, OUT_DOMAIN), :);
if ~isequal(sort(unique(T.dataset))', sort(DATASETS)), error("lowerTempo:ds", "%s", "Unexpected in-domain datasets"); end
for v = ["f0" "sigma" "a" "alpha"]
    nb = sum(~isfinite(T.(v)));
    if nb > 0, error("lowerTempo:nan", "%s", sprintf("%d non-finite %s values in-domain", nb, v)); end
end
T.relNoise = T.sigma ./ T.a;
T.subjArm  = T.dataset + "_" + string(T.subj);
% Denominator = every trial x pipeline with a forward map; a trial with no beta_obs counts as a
% failure, exactly as Finding #245 Table 3 (t3) does (checked by the anchor below).
EV = true(height(T), 1);
fprintf("in-domain rows %d; per pipeline %d; no_beta_obs rows %d (counted as failures)\n", ...
    height(T), height(T)/6, sum(T.status == "no_beta_obs"));

%% ANCHOR: Finding #245 Table 3 (SG-IRLS, coarse bins)
U = T(T.pipeline == "SG-IRLS" & EV, :);
b = discretize(U.f0, S.EDGES);  t3 = S.t3;  maxd = 0;
pc = ["pct_x_0_35" "pct_x0_35_0_7" "pct_x0_7_1_4" "pct_x1_4_2_8" "pct_x_2_8"];
for d = DATASETS
    r = t3(string(t3.dataset) == d, :);
    for k = 1:numel(S.EDGES)-1
        m = U.dataset == d & b == k;
        if ~any(m), continue; end
        maxd = max(maxd, abs(100*mean(U.status(m) == "rise") - r.(pc(k))));
    end
end
if maxd > 0.05, error("lowerTempo:anchor", "%s", sprintf("Table 3 not reproduced (max diff %.3f points)", maxd)); end
fprintf("ANCHOR Finding #245 Table 3 reproduced (max diff %.3g points)\n", maxd);

R = struct();  R.runStamp = stamp;  R.input = IN_MAT;

%% 1. Pass rate by fine tempo bin, per dataset and pooled
for p = PIPES
    U = T(T.pipeline == p & EV, :);  b = discretize(U.f0, FINE);
    P = nan(numel(DATASETS)+1, numel(FINE)-1);  N = zeros(size(P));
    for i = 1:numel(DATASETS)+1
        if i <= numel(DATASETS), mi = U.dataset == DATASETS(i); else, mi = true(height(U), 1); end
        for k = 1:numel(FINE)-1
            m = mi & b == k;  N(i,k) = sum(m);
            if N(i,k) >= MIN_N, P(i,k) = 100*mean(U.status(m) == "rise"); end
        end
    end
    R.(fld(p)).passFine = array2table(P, "RowNames", [DATASETS "ALL"], "VariableNames", binNames(FINE));
    R.(fld(p)).passFineN = array2table(N, "RowNames", [DATASETS "ALL"], "VariableNames", binNames(FINE));
    printGrid("1. Pass rate (%) by tempo, " + p, P, N, [DATASETS "ALL"], FINE);
end

%% 2. Pass rate by relative-noise third x tempo (pooled in-domain)
for p = PIPES
    U = T(T.pipeline == p & EV, :);
    q = quantile(U.relNoise, [1/3 2/3]);  tb = discretize(U.relNoise, [-Inf q Inf]);  b = discretize(U.f0, NOISE_F0);
    P = nan(3, numel(NOISE_F0)-1);  N = zeros(size(P));  lab = strings(1, 3);
    for t = 1:3
        lab(t) = sprintf("sigma/a %.3f-%.3f", min(U.relNoise(tb == t)), max(U.relNoise(tb == t)));
        for k = 1:numel(NOISE_F0)-1
            m = tb == t & b == k;  N(t,k) = sum(m);
            if N(t,k) >= MIN_N, P(t,k) = 100*mean(U.status(m) == "rise"); end
        end
    end
    R.(fld(p)).noiseCuts = q;
    R.(fld(p)).passNoise = array2table(P, "RowNames", lab, "VariableNames", binNames(NOISE_F0));
    R.(fld(p)).passNoiseN = N;
    printGrid("2. Pass rate (%) by relative-noise third x tempo, " + p, P, N, lab, NOISE_F0);
end

%% 3-4. Median local map slope and generator-scale CI width (invertible trials)
for p = PIPES
    U = T(T.pipeline == p & T.status == "rise", :);
    q = R.(fld(p)).noiseCuts;  lowNoise = U.relNoise <= q(1);
    fprintf("\n%s invertible trials %d; slope NaN %d; ciW NaN %d (%.1f%%)\n", p, height(U), ...
        sum(isnan(U.slope)), sum(isnan(U.ciW)), 100*mean(isnan(U.ciW)));
    groups = {"ALL", true(height(U),1); "Fraser", U.dataset=="Fraser"; "Hickman_PLAC", U.dataset=="Hickman_PLAC"; ...
              "Cook_ASD", U.dataset=="Cook_ASD"; "lowest-noise third", lowNoise};
    b = discretize(U.f0, FINE_S);
    SL = nan(size(groups,1), numel(FINE_S)-1);  CW = SL;  NN = zeros(size(SL));
    for g = 1:size(groups,1)
        for k = 1:numel(FINE_S)-1
            m = groups{g,2} & b == k;  NN(g,k) = sum(m);
            if NN(g,k) >= MIN_N
                SL(g,k) = median(U.slope(m), "omitnan");  CW(g,k) = median(U.ciW(m), "omitnan");
            end
        end
    end
    gl = string(groups(:,1))';
    R.(fld(p)).slope = array2table(SL, "RowNames", gl, "VariableNames", binNames(FINE_S));
    R.(fld(p)).ciW   = array2table(CW, "RowNames", gl, "VariableNames", binNames(FINE_S));
    R.(fld(p)).slopeN = NN;
    printGrid("3. Median local map slope, " + p, SL, NN, gl, FINE_S, "%6.2f");
    printGrid("4. Median generator-scale CI width, " + p, CW, NN, gl, FINE_S, "%6.3f");
end

%% 5. Exposure: share of in-domain trials below each threshold (SG-IRLS rows = one per trial)
U = T(T.pipeline == "SG-IRLS", :);  E = nan(numel(DATASETS)+1, numel(EXPO)+2);
fprintf("\n5. Exposure (SG-IRLS rows, one per trial with a map)\n%-14s %6s %8s %8s %8s %8s %10s\n", "dataset", "n", "medf0", "<0.3", "<0.5", "<0.7", "med s/a");
for i = 1:numel(DATASETS)+1
    if i <= numel(DATASETS), m = U.dataset == DATASETS(i); nm = DATASETS(i); else, m = true(height(U),1); nm = "ALL"; end
    x = U.f0(m);  E(i,:) = [median(x), 100*mean(x < EXPO(1)), 100*mean(x < EXPO(2)), 100*mean(x < EXPO(3)), median(U.relNoise(m))];
    fprintf("%-14s %6d %8.2f %7.1f%% %7.1f%% %7.1f%% %10.3f\n", nm, sum(m), E(i,:));
end
R.exposure = array2table(E, "RowNames", [DATASETS "ALL"], "VariableNames", ["medF0" "pctBelow0p3" "pctBelow0p5" "pctBelow0p7" "medRelNoise"]);

%% 6. The three edges
fprintf("\n6. Edges (lowest bin lower edge from which the criterion holds in every faster bin; Hz)\n");
for p = PIPES
    e1 = runEdge(R.(fld(p)).passFine{"ALL",:} >= 100*PASS_MIN, R.(fld(p)).passFine{"ALL",:}, FINE);
    pn = R.(fld(p)).passNoise{:,:};  e1n = nan(1,3);
    for t = 1:3, e1n(t) = runEdge(pn(t,:) >= 100*PASS_MIN, pn(t,:), NOISE_F0); end
    sl = R.(fld(p)).slope{"ALL",:};  cw = R.(fld(p)).ciW{"ALL",:};
    e2 = runEdge(sl >= SLOPE_MIN, sl, FINE_S);  e3 = runEdge(cw <= CIW_MAX, cw, FINE_S);
    R.(fld(p)).edges = struct("canInvertPooled", e1, "canInvertByNoiseThird", e1n, "gain0p7", e2, "ciW0p05", e3);
    fprintf("%-9s E1 can invert (pass >= %.0f%%): pooled %.2f; by noise third low/mid/high %.2f / %.2f / %.2f\n", ...
        p, 100*PASS_MIN, e1, e1n);
    fprintf("%-9s E2 informative (median slope >= %.2f): %.2f   E3 precise (median CI width <= %.2f): %.2f\n", ...
        p, SLOPE_MIN, e2, CIW_MAX, e3);
end

%% 7. Exploratory: does a crude noise-at-3f0 index summarise pass/fail? (fixed-effects logistic)
% index = -(2 log10(sigma/a) - alpha log10(3 f0)): the relative noise power implied at the third
% harmonic for a power-law spectrum, sign flipped so larger = cleaner. Descriptive only (no random
% effects); compared by AIC with separate terms and with tempo alone.
for p = PIPES
    U = T(T.pipeline == p & EV, :);
    D = table(double(U.status == "rise"), log10(U.relNoise), log10(U.f0), U.alpha, log10(U.M), ...
        'VariableNames', ["y" "rn" "lf" "al" "lM"]);
    D.snr = -(2*D.rn - D.al .* log10(3*U.f0));
    a1 = fitglm(D, "y ~ rn + lf + al + lM", "Distribution", "binomial").ModelCriterion.AIC;
    a2 = fitglm(D, "y ~ snr + lM", "Distribution", "binomial").ModelCriterion.AIC;
    a3 = fitglm(D, "y ~ lf + lM", "Distribution", "binomial").ModelCriterion.AIC;
    R.(fld(p)).aic = struct("separate", a1, "index3f0", a2, "tempoOnly", a3);
    fprintf("7. %s AIC: separate terms %.1f | 3f0 index %.1f | tempo only %.1f\n", p, a1, a2, a3);
end

%% CLAIMS CHECK: numbers quoted to Dagmar on 2026-11-04 (tolerance = rounding of the quote)
C = { ...
 "SG-IRLS pooled pass 0.1-0.2 Hz", 45,  R.SG_IRLS.passFine{"ALL",1}, 0.5;
 "SG-IRLS pooled pass 0.2-0.3 Hz", 89,  R.SG_IRLS.passFine{"ALL",2}, 0.5;
 "SG-IRLS pooled pass 0.3-0.4 Hz", 96,  R.SG_IRLS.passFine{"ALL",3}, 0.5;
 "SG-IRLS pooled pass 0.4-0.5 Hz", 97,  R.SG_IRLS.passFine{"ALL",4}, 0.5;
 "SG-IRLS pooled pass 0.5-0.7 Hz", 95,  R.SG_IRLS.passFine{"ALL",5}, 0.5;
 "SG-IRLS pooled pass 0.7-1.0 Hz", 95,  R.SG_IRLS.passFine{"ALL",6}, 0.5;
 "SG-IRLS pooled pass 1.0-1.4 Hz", 97,  R.SG_IRLS.passFine{"ALL",7}, 0.5;
 "SG-IRLS pooled pass >1.4 Hz",    99,  R.SG_IRLS.passFine{"ALL",8}, 0.5;
 "BWFD-OLS pooled pass 0.1-0.2 Hz", 52, R.BWFD_OLS.passFine{"ALL",1}, 0.5;
 "BWFD-OLS pooled pass >1.4 Hz",    89, R.BWFD_OLS.passFine{"ALL",8}, 0.5;
 "SG-IRLS high-noise third pass <0.3 Hz", 41, R.SG_IRLS.passNoise{3,1}, 0.5;
 "SG-IRLS high-noise third pass >0.7 Hz", 93, R.SG_IRLS.passNoise{3,4}, 0.5;
 "SG-IRLS mid-noise third pass <0.3 Hz",  86, R.SG_IRLS.passNoise{2,1}, 0.5;
 "SG-IRLS slope 0.2-0.3 Hz", 0.28, R.SG_IRLS.slope{"ALL",2}, 0.005;
 "SG-IRLS slope 0.3-0.4 Hz", 0.32, R.SG_IRLS.slope{"ALL",3}, 0.005;
 "SG-IRLS slope 0.6-0.7 Hz", 0.74, R.SG_IRLS.slope{"ALL",8}, 0.005;
 "SG-IRLS slope 1.0-1.4 Hz (peak)", 0.87, R.SG_IRLS.slope{"ALL",11}, 0.005;
 "SG-IRLS slope >2 Hz", 0.75, R.SG_IRLS.slope{"ALL",13}, 0.005;
 "Fraser slope 1.0-1.4 Hz", 0.85, R.SG_IRLS.slope{"Fraser",11}, 0.005;
 "SG-IRLS CI width 0.85-1.0 Hz", 0.10, R.SG_IRLS.ciW{"ALL",10}, 0.005;
 "SG-IRLS CI width 1.0-1.4 Hz",  0.05, R.SG_IRLS.ciW{"ALL",11}, 0.005;
 "Exposure: all trials below 0.5 Hz (%)", 43.3, R.exposure{"ALL","pctBelow0p5"}, 0.05;
 "Exposure: all trials below 0.7 Hz (%)", 59.6, R.exposure{"ALL","pctBelow0p7"}, 0.05;
 "Cook CTRL median f0 (Hz)", 0.26, R.exposure{"Cook_CTRL","medF0"}, 0.005;
 "AIC separate terms (SG-IRLS)", 599.4, R.SG_IRLS.aic.separate, 0.05;
 "AIC 3f0 index (SG-IRLS)",      941.9, R.SG_IRLS.aic.index3f0, 0.05};
fprintf("\nCLAIMS CHECK\n");  nd = 0;
for i = 1:size(C,1)
    ok = abs(C{i,3} - C{i,2}) <= C{i,4};  nd = nd + ~ok;
    fprintf("%-44s quoted %8.3f  computed %8.3f  %s\n", C{i,1}, C{i,2}, C{i,3}, ifelseStr(ok, "MATCH", "DIFF"));
end
fprintf("%d of %d quoted numbers match; %d differ\n", size(C,1)-nd, size(C,1), nd);
R.claims = cell2table(C, "VariableNames", ["claim" "quoted" "computed" "tol"]);

save(OUT_MAT, "R", "-v7");
fprintf("saved %s\nlog %s\n", OUT_MAT, logf);

%% ---- local functions ----
function s = fld(p), s = strrep(p, "-", "_"); end
function s = ifelseStr(c, a, b), if c, s = a; else, s = b; end, end
function n = binNames(E)
    n = strings(1, numel(E)-1);
    for k = 1:numel(E)-1, n(k) = matlab.lang.makeValidName(sprintf("f%g_%g", E(k), E(k+1))); end
end
function e = runEdge(ok, vals, E)
    % lowest lower edge from which ok holds in every evaluable faster bin; NaN if never
    e = NaN;  ev = ~isnan(vals);
    for k = 1:numel(ok)
        if ev(k) && all(ok(k:end) | ~ev(k:end)), e = E(k); return; end
    end
end
function printGrid(ttl, P, N, rows, E, fmt)
    if nargin < 6, fmt = "%6.0f"; end
    fprintf("\n%s\n%-20s", ttl, "");
    for k = 1:numel(E)-1, fprintf("%12s", sprintf("%g-%g", E(k), E(k+1))); end
    fprintf("\n");
    for i = 1:numel(rows)
        fprintf("%-20s", rows(i));
        for k = 1:numel(E)-1
            if isnan(P(i,k)), fprintf("%12s", "-"); else, fprintf("%12s", sprintf(fmt + "/%d", P(i,k), N(i,k))); end
        end
        fprintf("\n");
    end
end
