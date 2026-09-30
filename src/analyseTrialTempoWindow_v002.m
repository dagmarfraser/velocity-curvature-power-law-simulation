%% analyseTrialTempoWindow_v002.m
% Finding #227 re-run on the validity domain of Finding #242. v001 (2026-09-26) predates #242
% and excluded only Zarandi (OUT_DOMAIN = "Zarandi"), so Dhieb's trials entered every pooled
% table and both models: the GLMM's 3,728 trials are the five in-domain datasets plus Dhieb, not
% the 3,722 in-domain trials (coherence item A4.3). v002 builds the same trial x pipeline table
% once and runs two domain definitions on it, so the change is attributable:
%   V0 v001 DOMAIN   Zarandi out                regression anchor: must reproduce #227's GLMM
%                                               (SG-IRLS n = 3,728, b(lf) = 2.01, SE 0.201,
%                                               df 3717) and LMM (n = 3,166, b(lf) = -0.775).
%   V1 #242 DOMAIN   Zarandi and Dhieb out      the result to cite.
% Everything else is v001 unchanged: v015 per-trial maps, bins, covariates, model formulas,
% subject random intercept, dataset fixed. Per-dataset invertibility by f0 bin (Table 3) is
% computed for all seven datasets and flagged by domain, so an out-of-domain row is reported
% beside the in-domain ones, not pooled. Trials that carry no forward map (empty curve) are
% counted per dataset rather than silently dropped.
% Inputs:  src/loopClosureResults_<dataset>_all_shaped_xu_v015.mat (anchor = #227's published numbers, in ANCHOR)
% Outputs: results/trialTempoWindow_v002.mat, figures/trialTempoWindow_v002.png
% USAGE:   from the project root: analyseTrialTempoWindow_v002   (loads seven corpora: run as a batch job)
% Fraser, D.S. (2026)  v002

%% CONFIG
ROOT      = fileparts(fileparts(mfilename("fullpath")));
DATASETS  = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Dhieb" "Zarandi"];
DOMAINS   = {["Zarandi"], ["Zarandi" "Dhieb"]};                     % V0 (v001), V1 (#242)
VNAME     = ["V0 v001 domain (Zarandi out)", "V1 #242 domain (Zarandi, Dhieb out)"];
N_6A      = [2829 94 102 359 338 NaN NaN];                         % Table 6a all-trial Ns
PL        = ["BWFD-OLS" "SG-OLS" "BWFD-LMLS" "SG-LMLS" "BWFD-IRLS" "SG-IRLS"];  % runner order (v012 L118)
EDGES     = [0 0.35 0.7 1.4 2.8 Inf];
BIN_NAMES = ["<0.35" "0.35-0.7" "0.7-1.4" "1.4-2.8" ">2.8"];
PIPE_TAB  = "SG-IRLS";
MODEL_PIPES = ["SG-IRLS" "BWFD-OLS"];
ANCHOR    = struct("nGLMM", 3728, "b", 2.01, "SE", 0.201, "df", 3717, "nLMM", 3166, "bLMM", -0.775);
OUT_MAT   = fullfile(ROOT, "results", "trialTempoWindow_v002.mat");
OUT_PNG   = fullfile(ROOT, "figures", "trialTempoWindow_v002.png");
for f = [fileparts(OUT_MAT) fileparts(OUT_PNG)]
    if ~isfolder(f), error("trialTempo2:outDir", "%s", "Missing folder: " + f); end
end

%% Build the trial x pipeline table (v001 L31-57)
rows = {};  cover = table();
for k = 1:numel(DATASETS)
    d = DATASETS(k);
    f = fullfile(ROOT, "src", "loopClosureResults_" + d + "_all_shaped_xu_v015.mat");
    if ~isfile(f), error("trialTempo2:input", "%s", "Missing: " + f); end
    L = load(f, "results", "betaGenVec");  Rr = L.results;  bgv = L.betaGenVec;
    if ~isnan(N_6A(k)) && numel(Rr) ~= N_6A(k)
        error("trialTempo2:Count", "%s", sprintf("%s: %d trials, Table 6a says %d", d, numel(Rr), N_6A(k)));
    end
    nNoCurve = 0;  nTab = 0;
    for i = 1:numel(Rr)
        Mc = Rr(i).betaRecCurveMed;
        if isempty(Mc), nNoCurve = nNoCurve + 1; continue; end
        c6 = Mc(strcmp(PL, PIPE_TAB), :);
        if all(isnan(c6)), nTab = nTab + 1; end
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
    cover = [cover; table(d, numel(Rr), nNoCurve, nTab, 'VariableNames', ["dataset" "nTrials" "nNoMap" "nNoMapSGIRLS"])]; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', ["dataset" "subj" "f0" "sigma" "alpha" "a" "M" "pipeline" ...
    "status" "pos" "top" "floor" "drop" "ciW" "slope"]);
T.rise = T.status == "rise";  T.fail = T.status == "neither";  T.foldBranch = ismember(T.status, ["ambiguous" "desc"]);
T.bin = discretize(T.f0, EDGES, 'categorical', BIN_NAMES);
fprintf("Trial x pipeline cells: %d (%d datasets)\nTrials without a forward map (not in any model):\n", height(T), numel(DATASETS));
disp(cover)

%% Table 3 for all seven datasets, flagged (PIPE_TAB)
s = T(T.pipeline == PIPE_TAB, :);
t3 = table();
for d = DATASETS
    r = table(d, ~ismember(d, DOMAINS{2}), 'VariableNames', ["dataset" "inDomain242"]);
    for b = BIN_NAMES
        x = s(s.dataset == d & s.bin == b, :);
        r.("pct_" + matlab.lang.makeValidName(b)) = ifelse_local(height(x) > 0, round(100 * mean(x.rise), 1), NaN);
        r.("n_" + matlab.lang.makeValidName(b)) = height(x);
    end
    t3 = [t3; r]; %#ok<AGROW>
end
fprintf("\nTable 3: %s %% invertible by dataset x f0 bin (all datasets; inDomain242 flags #242)\n", PIPE_TAB);  disp(t3)

%% Each domain definition
Res = cell(1, numel(DOMAINS));
for v = 1:numel(DOMAINS)
    D = T(~ismember(T.dataset, DOMAINS{v}), :);
    fprintf("\n################ %s ################\n", VNAME(v));
    Res{v} = fitDomain_local(D, MODEL_PIPES, EDGES, BIN_NAMES);
end

%% Regression anchor: V0 must reproduce Finding #227
g0 = Res{1}.(erase(PIPE_TAB, "-"));
c0 = g0.coefGLMM(g0.coefGLMM.term == "lf", :);  l0 = g0.coefLMM(g0.coefLMM.term == "lf", :);
got = [g0.nGLMM, c0.b, c0.SE, c0.df, g0.nLMM, l0.b];
want = [ANCHOR.nGLMM, ANCHOR.b, ANCHOR.SE, ANCHOR.df, ANCHOR.nLMM, ANCHOR.bLMM];
if any(abs(got - want) > [0 0.005 0.0005 0 0 0.0005])
    error("trialTempo2:Anchor", "%s", "V0 does not reproduce #227: got " + mat2str(got, 5) + ", want " + mat2str(want, 5));
end
fprintf("\nREGRESSION ANCHOR passed: V0 reproduces #227 (GLMM n = %d, b = %.2f, SE %.3f, df %d; LMM n = %d, b = %.3f).\n", got);

%% Side by side
fprintf("\n%-38s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n", "", "nGLMM", "nSubj", "b(lf)", "SE", "t", "df", "b(lf2)", "peakHz", "nLMM");
for v = 1:numel(DOMAINS)
    for p = MODEL_PIPES
        g = Res{v}.(erase(p, "-"));  c = g.coefGLMM;  a = c(c.term == "lf", :);  q = c(c.term == "lf2", :);
        fprintf("%-38s %8d %8d %8.3f %8.3f %8.2f %8d %8.3f %8.2f %8d\n", VNAME(v) + " " + p, g.nGLMM, g.nSubj, ...
            a.b, a.SE, a.t, a.df, q.b, g.peakHz, g.nLMM);
    end
end

%% Figure: in-domain (#242) invertibility and slope by f0 bin
fg = figure("Color", "w", "Position", [80 80 1000 400]);  tl = tiledlayout(1, 2, "TileSpacing", "compact");
nexttile; hold on;
for r = 1:height(t3)
    y = t3{r, startsWith(t3.Properties.VariableNames, "pct_")};
    ls = ifelse_local(t3.inDomain242(r), "-o", ":x");
    if t3.dataset(r) == "Zarandi", continue; end
    plot(1:numel(BIN_NAMES), y, ls, "LineWidth", 1.3, "DisplayName", strrep(t3.dataset(r), "_", " ") + ifelse_local(t3.inDomain242(r), "", " (outside domain)"));
end
xticks(1:numel(BIN_NAMES)); xticklabels(BIN_NAMES); xlabel("trial f_0 (Hz)"); ylabel("% invertible (" + PIPE_TAB + ")");
ylim([0 100]); box on; legend("Location", "southeast");
nexttile; pb = Res{2}.(erase(PIPE_TAB, "-")).bins;
yyaxis left;  plot(1:height(pb), pb.median_slope, "-o", "LineWidth", 1.5); ylabel("median local slope at \beta_{gen}*");
yyaxis right; plot(1:height(pb), pb.median_ciW, "-s", "LineWidth", 1.5); ylabel("median CI width (generator scale)");
xticks(1:height(pb)); xticklabels(string(pb.bin)); xlabel("trial f_0 (Hz)"); box on;
title(tl, "Per-trial forward maps (v015): tempo window, validity domain of #242");
set(findall(fg, "Type", "axes"), "Toolbar", []);  exportgraphics(fg, OUT_PNG, "Resolution", 200);

runDate = string(datetime("now", "Format", "yyyy-MM-dd"));
save(OUT_MAT, "T", "t3", "cover", "Res", "DOMAINS", "VNAME", "EDGES", "ANCHOR", "runDate", "-v7.3");
fprintf("\nSaved: %s\nFigure: %s\n", OUT_MAT, OUT_PNG);

%% =========================================================================
function R = fitDomain_local(D, pipes, edges, binNames)
% v001 L70-122 on one domain definition.
    t2 = table();
    for b = binNames
        x = D(D.bin == b, :);  fx = x(x.fail, :);
        t2 = [t2; table(b, height(x), 100*mean(x.rise), 100*mean(x.foldBranch), 100*mean(fx.pos == "below"), ...
            100*mean(fx.pos == "above"), median(x.top, "omitnan"), 'VariableNames', ...
            ["f0bin" "n" "pctRise" "pctFoldBranch" "failBelowPct" "failAbovePct" "medianTop"])]; %#ok<AGROW>
    end
    t2{:, 3:end} = round(t2{:, 3:end}, 3);
    fprintf("Table 2: by trial f0 bin, all pipelines\n");  disp(t2)
    D.subj = categorical(D.subj);  D.datasetC = categorical(D.dataset);
    z = @(x) (x - mean(x, "omitnan")) ./ std(x, "omitnan");
    D.lf = z(log(D.f0));  D.lf2 = D.lf.^2;  D.lsa = z(log(D.sigma ./ D.a));  D.al = z(D.alpha);  D.lM = z(log(D.M));
    mu = mean(log(D.f0));  sd = std(log(D.f0));
    terms = ["lf" "lf2" "lsa" "al" "lM"];  R = struct("t2", t2);
    for p = pipes
        V = D(D.pipeline == p, :);
        g = fitglme(V, "rise ~ 1 + lf + lf2 + lsa + al + lM + datasetC + (1|subj)", "Distribution", "Binomial", "FitMethod", "Laplace");
        W = V(V.rise & isfinite(V.ciW) & V.ciW > 0, :);  W.lciW = log(W.ciW);
        m = fitlme(W, "lciW ~ 1 + lf + lf2 + lsa + al + lM + datasetC + (1|subj)");
        cg = coefTab_local(g, terms);  cm = coefTab_local(m, terms);
        b1 = cg.b(cg.term == "lf");  b2 = cg.b(cg.term == "lf2");  vtx = NaN;
        if b2 < 0, vtx = exp(mu + sd * (-b1 / (2*b2))); end
        W.bin = discretize(W.f0, edges, 'categorical', binNames);
        pb = groupsummary(W, "bin", "median", ["ciW" "slope"]);
        byDs = groupcounts(V, "dataset");
        fprintf("\n===== %s  (unit = trial; %d trials, %d subjects; above 2.8 Hz: %d)\n", p, height(V), numel(unique(V.subj)), nnz(V.f0 > 2.8));
        disp(byDs(:, ["dataset" "GroupCount"]))
        fprintf("GLMM invertible (binomial, Laplace)\n");  disp(cg)
        fprintf("  fitted peak: %s\n", ifelse_local(isnan(vtx), "none (quadratic >= 0)", sprintf("%.2f Hz", vtx)));
        fprintf("LMM log(CI width), invertible trials only (n = %d)\n", height(W));  disp(cm)
        fprintf("Median CI width and local slope at beta_gen* by f0 bin:\n");  disp(pb)
        R.(erase(p, "-")) = struct("glme", g, "lme", m, "coefGLMM", cg, "coefLMM", cm, "peakHz", vtx, "bins", pb, ...
            "nGLMM", height(V), "nSubj", numel(unique(V.subj)), "nLMM", height(W), "nAbove28", nnz(V.f0 > 2.8), "byDataset", byDs);
    end
end

function t = coefTab_local(mdl, terms)
    c = mdl.Coefficients;  k = ismember(string(c.Name), terms);
    t = table(string(c.Name(k)), round(c.Estimate(k), 3), round(c.SE(k), 3), round(c.tStat(k), 2), c.DF(k), ...
        c.pValue(k), round(c.Lower(k), 3), round(c.Upper(k), 3), 'VariableNames', ["term" "b" "SE" "t" "df" "p" "CIlo" "CIhi"]);
end

function v = ifelse_local(c, a, b)
    if c, v = a; else, v = b; end
end
