%% analyseTrialTempoWindow_v001.m
% Does each trial's own tempo predict whether its forward map is invertible, and how
% precise inversion is? Read-only analysis of the v015 per-trial forward maps (each
% built at the trial's own f0), 2026-09-26. No simulation.
%   Table 1  per dataset: curves folding (median drop > FOLD_DROP), status mix, where
%            failed trials sit (below the map / above its top / within)
%   Table 2  in-domain, by f0 bin: invertible, on a fold branch, failure position, map top
%   Table 3  per dataset x f0 bin: % invertible for PIPE_TAB (and n)
%   Table 4  LMLS family: % on a fold branch by f0 bin
%   Models   per MODEL_PIPES, unit = trial, in-domain datasets:
%            GLMM  rise ~ lf + lf2 + lsa + al + lM + dataset + (1|subject)      (binomial, Laplace)
%            LMM   log(CI width) ~ same, invertible trials only
%            lf = z(log f0), lsa = z(log sigma/a), al = z(alpha), lM = z(log M)
%   Bins     median CI width and local map slope at beta_gen* by f0 bin
% Inputs: src/loopClosureResults_<dataset>_all_shaped_xu_v015.mat
% Outputs: results/trialTempoWindow_v001.mat, figures/trialTempoWindow_v001.png

%% CONFIG
ROOT      = fileparts(fileparts(mfilename("fullpath")));
DATASETS  = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Dhieb" "Zarandi"];
OUT_DOMAIN = "Zarandi";                          % D20: outside the validity domain
PL        = ["BWFD-OLS" "SG-OLS" "BWFD-LMLS" "SG-LMLS" "BWFD-IRLS" "SG-IRLS"];  % runner order (v012 L118)
EDGES     = [0 0.35 0.7 1.4 2.8 Inf];
BIN_NAMES = ["<0.35" "0.35-0.7" "0.7-1.4" "1.4-2.8" ">2.8"];
FOLD_DROP = 0.02;                                % descriptive; production check is invertStatus
PIPE_TAB  = "SG-IRLS";
MODEL_PIPES = ["SG-IRLS" "BWFD-OLS"];
OUT_MAT   = fullfile(ROOT, "results", "trialTempoWindow_v001.mat");
OUT_PNG   = fullfile(ROOT, "figures", "trialTempoWindow_v001.png");

%% Build the trial x pipeline table
rows = {};
for d = DATASETS
    f = fullfile(ROOT, "src", "loopClosureResults_" + d + "_all_shaped_xu_v015.mat");
    if ~isfile(f), error("trialTempo:input", "%s", "Missing: " + f); end
    L = load(f, "results", "betaGenVec");  Rr = L.results;  bgv = L.betaGenVec;
    for i = 1:numel(Rr)
        Mc = Rr(i).betaRecCurveMed;
        if isempty(Mc), continue; end
        for p = 1:numel(PL)
            c = Mc(p, :);
            if all(isnan(c)), continue; end
            b = Rr(i).betaObs(p);
            if ~isfinite(b), pos = "nan"; elseif b < min(c), pos = "below"; elseif b > max(c), pos = "above"; else, pos = "within"; end
            gr = gradient(c, bgv);  bs = Rr(i).betaGenStar(p);  sl = NaN;
            if isfinite(bs), sl = interp1(bgv, gr, bs); end
            rows(end+1, :) = {d, d + "_" + string(Rr(i).subjectID), Rr(i).f0, Rr(i).sigmaMM, Rr(i).alphaIRA, ...
                Rr(i).a_mm, Rr(i).M, PL(p), Rr(i).invertStatus(p), pos, c(end), c(1), max(c) - c(end), ...
                Rr(i).ciHi(p) - Rr(i).ciLo(p), sl}; %#ok<SAGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', ["dataset" "subj" "f0" "sigma" "alpha" "a" "M" "pipeline" ...
    "status" "pos" "top" "floor" "drop" "ciW" "slope"]);
T.rise = T.status == "rise";  T.fail = T.status == "neither";  T.foldBranch = ismember(T.status, ["ambiguous" "desc"]);
T.bin = discretize(T.f0, EDGES, 'categorical', BIN_NAMES);
fprintf("Trial x pipeline cells: %d (%d datasets)\n", height(T), numel(DATASETS));

%% Table 1: per dataset
t1 = table();
for d = DATASETS
    x = T(T.dataset == d, :);  fx = x(x.fail, :);
    t1 = [t1; table(d, height(x), 100*mean(x.drop > FOLD_DROP), 100*mean(x.fail), 100*mean(fx.pos == "below"), ...
        100*mean(fx.pos == "above"), 100*mean(fx.pos == "within"), 100*mean(x.status == "ambiguous"), ...
        'VariableNames', ["dataset" "n" "pctFolding" "pctNeither" "neitherBelow" "neitherAbove" "neitherWithin" "pctAmbiguous"])]; %#ok<AGROW>
end
t1{:, 3:end} = round(t1{:, 3:end}, 1);
fprintf("\nTable 1: per dataset (all pipelines)\n");  disp(t1)

%% Table 2: in-domain by f0 bin
D = T(T.dataset ~= OUT_DOMAIN, :);
t2 = table();
for b = BIN_NAMES
    x = D(D.bin == b, :);  fx = x(x.fail, :);
    t2 = [t2; table(b, height(x), 100*mean(x.rise), 100*mean(x.foldBranch), 100*mean(fx.pos == "below"), ...
        100*mean(fx.pos == "above"), median(x.top, "omitnan"), 'VariableNames', ...
        ["f0bin" "n" "pctRise" "pctFoldBranch" "failBelowPct" "failAbovePct" "medianTop"])]; %#ok<AGROW>
end
t2{:, 3:end} = round(t2{:, 3:end}, 3);
fprintf("\nTable 2: in-domain (excl. %s), all pipelines, by trial f0 bin\n", OUT_DOMAIN);  disp(t2)

%% Table 3: per dataset x f0 bin, PIPE_TAB
s = D(D.pipeline == PIPE_TAB, :);
[t3pct, t3n] = deal(array2table(NaN(numel(unique(s.dataset)), numel(BIN_NAMES)), 'VariableNames', BIN_NAMES, ...
    'RowNames', cellstr(unique(s.dataset))));
for d = string(t3pct.Properties.RowNames)'
    for b = BIN_NAMES
        x = s(s.dataset == d & s.bin == b, :);
        t3n{char(d), b} = height(x);
        if height(x) > 0, t3pct{char(d), b} = round(100*mean(x.rise), 1); end
    end
end
fprintf("\nTable 3: %s %% invertible by dataset x f0 bin\n", PIPE_TAB);  disp(t3pct);  fprintf("n trials:\n");  disp(t3n)

%% Table 4: LMLS family on a fold branch by f0 bin
t4 = groupsummary(D(contains(D.pipeline, "LMLS"), :), "bin", "mean", "foldBranch");
t4.pct = round(100*t4.mean_foldBranch, 2);
fprintf("\nTable 4: LMLS pipelines, %% on a fold branch (ambiguous/desc) by f0 bin\n");  disp(t4(:, ["bin" "GroupCount" "pct"]))

%% Models, unit = trial, in-domain
D.subj = categorical(D.subj);  D.datasetC = categorical(D.dataset);
z = @(x) (x - mean(x, "omitnan")) ./ std(x, "omitnan");
D.lf = z(log(D.f0));  D.lf2 = D.lf.^2;  D.lsa = z(log(D.sigma ./ D.a));  D.al = z(D.alpha);  D.lM = z(log(D.M));
mu = mean(log(D.f0));  sd = std(log(D.f0));
terms = ["lf" "lf2" "lsa" "al" "lM"];  Models = struct();
for p = MODEL_PIPES
    V = D(D.pipeline == p, :);
    g = fitglme(V, "rise ~ 1 + lf + lf2 + lsa + al + lM + datasetC + (1|subj)", "Distribution", "Binomial", "FitMethod", "Laplace");
    W = V(V.rise & isfinite(V.ciW) & V.ciW > 0, :);  W.lciW = log(W.ciW);
    m = fitlme(W, "lciW ~ 1 + lf + lf2 + lsa + al + lM + datasetC + (1|subj)");
    cg = coefTab_local(g, terms);  cm = coefTab_local(m, terms);
    b1 = cg.b(cg.term == "lf");  b2 = cg.b(cg.term == "lf2");  vtx = NaN;
    if b2 < 0, vtx = exp(mu + sd * (-b1 / (2*b2))); end
    W.bin = discretize(W.f0, EDGES, 'categorical', BIN_NAMES);
    pb = groupsummary(W, "bin", "median", ["ciW" "slope"]);
    fprintf("\n===== %s  (unit = trial; %d trials, %d subjects)\n", p, height(V), numel(unique(V.subj)));
    fprintf("GLMM invertible: binomial, Laplace; covariates z-scored; dataset fixed, subject random intercept\n");  disp(cg)
    fprintf("  fitted peak of invertibility: %s\n", ifelse_local(isnan(vtx), "none (quadratic >= 0)", sprintf("%.2f Hz (edge of data if > 2.8 Hz)", vtx)));
    fprintf("LMM log(CI width, generator scale), invertible trials only (n = %d)\n", height(W));  disp(cm)
    fprintf("Median CI width and local slope at beta_gen* by f0 bin:\n");  disp(pb)
    Models.(erase(p, "-")) = struct("glme", g, "lme", m, "coefGLMM", cg, "coefLMM", cm, "peakHz", vtx, "bins", pb);
end

%% Figure: invertibility and slope by f0 bin
if ~isfolder(fileparts(OUT_PNG)), error("trialTempo:outDir", "%s", "Missing folder: " + fileparts(OUT_PNG)); end
if ~isfolder(fileparts(OUT_MAT)), error("trialTempo:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
fg = figure("Color", "w", "Position", [80 80 1000 400]);  tl = tiledlayout(1, 2, "TileSpacing", "compact");
nexttile; hold on;
for d = string(t3pct.Properties.RowNames)'
    plot(1:numel(BIN_NAMES), t3pct{char(d), :}, "-o", "LineWidth", 1.3, "DisplayName", strrep(d, "_", " "));
end
xticks(1:numel(BIN_NAMES)); xticklabels(BIN_NAMES); xlabel("trial f_0 (Hz)"); ylabel("% invertible (" + PIPE_TAB + ")");
ylim([0 100]); box on; legend("Location", "southeast");
nexttile; pb = Models.(erase(PIPE_TAB, "-")).bins;
yyaxis left;  plot(1:height(pb), pb.median_slope, "-o", "LineWidth", 1.5); ylabel("median local slope at \beta_{gen}*");
yyaxis right; plot(1:height(pb), pb.median_ciW, "-s", "LineWidth", 1.5); ylabel("median CI width (generator scale)");
xticks(1:height(pb)); xticklabels(string(pb.bin)); xlabel("trial f_0 (Hz)"); box on;
title(tl, "Per-trial forward maps (v015): tempo window in real data, in-domain datasets");
set(findall(fg, "Type", "axes"), "Toolbar", []);  exportgraphics(fg, OUT_PNG, "Resolution", 200);

save(OUT_MAT, "T", "t1", "t2", "t3pct", "t3n", "t4", "Models", "EDGES", "FOLD_DROP", "OUT_DOMAIN", "-v7.3");
fprintf("\nSaved: %s\nFigure: %s\n", OUT_MAT, OUT_PNG);

%% =========================================================================
function t = coefTab_local(mdl, terms)
    c = mdl.Coefficients;  k = ismember(string(c.Name), terms);
    t = table(string(c.Name(k)), round(c.Estimate(k), 3), round(c.SE(k), 3), round(c.tStat(k), 2), c.DF(k), ...
        c.pValue(k), round(c.Lower(k), 3), round(c.Upper(k), 3), 'VariableNames', ["term" "b" "SE" "t" "df" "p" "CIlo" "CIhi"]);
end

function v = ifelse_local(c, a, b)
    if c, v = a; else, v = b; end
end
