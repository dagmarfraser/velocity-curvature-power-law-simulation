%% plotTempoSlope_v001.m
% Figure for the within-session tempo-slope result (v008 Part 4, 'How much of the beta_obs tempo slope
% can the pipeline alone produce?' and 'Tempo through the pipeline: what inversion attenuates').
% Reads saved results only; nothing is re-simulated or re-fitted.
%
%   Panel A  the analysis alone: for each of the five real-noise synthetic arms (AXr, AXrs2, Sm1Xr,
%            Sm2Xr, Sp1Xr), the within-session slope of the measured exponent in the synthetic copies
%            (true exponent independent of speed) against the real slope on the same matched trials;
%            vetted (SG-IRLS) and legacy (BWFD-OLS). Source: results/tempoSim/rerunSharesSameTrials_v001.mat
%            (out.Rows; shareOld = synthetic / empirical-matched, the share v008 quotes).
%   Panel B  pooled real slope in beta_obs against beta_gen*, both profiles, with model 95% CIs.
%            Source: results/checkTempoDecompositionDomain_v001.mat (D1{2} = #242 domain, term lfw).
%            The person-level bootstrap CIs v008 quotes (68%, 51-90; 78%, 49-113) come from
%            bootstrapTempoReduction_v001.mat, which is not on this machine, so they are NOT plotted.
%   Panel C  per dataset (vetted): beta_obs slope against beta_gen* slope with model CIs, annotated with
%            invertible / all trials. Source: results/tempoDecomposition_v002.mat (P1: per-dataset model,
%            session random intercept, no noise covariates; Sel: invertible shares).
%
% Basis: v015 corpus (Fraser N = 45-era files); Fraser values may move with the N = 48 rerun.
% Outputs: figures/tempoSlope_v001.png (300 dpi) and .pdf; results/plotTempoSlope_v001.mat;
%          logs/plotTempoSlope_v001_<stamp>.log. Ends with a claims check against v008 Part 4.
% USAGE: from src/ or anywhere: plotTempoSlope_v001
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT   = fileparts(fileparts(mfilename("fullpath")));
F_SH   = fullfile(ROOT, "results", "tempoSim", "rerunSharesSameTrials_v001.mat");
F_D1   = fullfile(ROOT, "results", "checkTempoDecompositionDomain_v001.mat");
F_TD   = fullfile(ROOT, "results", "tempoDecomposition_v002.mat");
OUT_PNG = fullfile(ROOT, "figures", "tempoSlope_v001.png");
OUT_PDF = fullfile(ROOT, "figures", "tempoSlope_v001.pdf");
OUT_MAT = fullfile(ROOT, "results", "plotTempoSlope_v001.mat");
LOGDIR  = fullfile(ROOT, "logs");
ARMS  = ["AXr" "AXrs2" "Sm1Xr" "Sm2Xr" "Sp1Xr"];
PIPES = ["SG-IRLS" "BWFD-OLS"];  PLAB = ["Vetted (SG-IRLS)" "Legacy (BWFD-OLS)"];
DS    = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO"];
DSLAB = ["Fraser" "Cook CTRL" "Cook ASD" "Hickman PLAC" "Hickman HALO"];
C_OBS = [0 114 178]/255;  C_SYN = [230 159 0]/255;  C_GEN = [0 158 115]/255;   % Okabe-Ito blue, orange, bluish green

for f = [F_SH F_D1 F_TD]
    if ~isfile(f), error("tempoSlopeFig:input", "%s", "Missing: " + f); end
end
for f = [fileparts(OUT_PNG) fileparts(OUT_MAT) LOGDIR]
    if ~isfolder(f), error("tempoSlopeFig:dir", "%s", "Missing folder: " + f); end
end
stamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
logf = fullfile(LOGDIR, "plotTempoSlope_v001_" + stamp + ".log");  diary(logf);  cu = onCleanup(@() diary("off"));
fprintf("plotTempoSlope_v001 run %s\n", stamp);

%% LOAD AND CHECK
S = load(F_SH);  if ~isfield(S, "out") || ~isfield(S.out, "Rows"), error("tempoSlopeFig:sh", "%s", "rerunShares lacks out.Rows"); end
A = S.out.Rows;  A.arm = string(A.arm);  A.pipeline = string(A.pipeline);
needA = ["arm" "pipeline" "nSyn" "nMatched" "synObs" "empMatched" "shareOld"];
if ~all(ismember(needA, string(A.Properties.VariableNames))), error("tempoSlopeFig:shcols", "%s", "Rows lacks a needed column"); end
if ~all(ismember(ARMS, A.arm)) || height(A) ~= numel(ARMS)*numel(PIPES), error("tempoSlopeFig:arms", "%s", "Unexpected arm rows"); end

C = load(F_D1);  if ~isfield(C, "D1") || numel(C.D1) < 2, error("tempoSlopeFig:d1", "%s", "D1{2} missing"); end
D = C.D1{2};  D.pipeline = string(D.pipeline);  D.outcome = string(D.outcome);  D.term = string(D.term);
if ~contains(string(C.VNAME(2)), "#242"), error("tempoSlopeFig:dom", "%s", "D1{2} is not the #242 domain"); end
D = D(D.term == "lfw", :);
pool = @(p, o) D(D.pipeline == p & D.outcome == o, :);

T = load(F_TD, "P1", "Sel");  P = T.P1;  P.pipeline = string(P.pipeline);  P.dataset = string(P.dataset);  P.outcome = string(P.outcome);
Sel = T.Sel;  Sel.pipeline = string(Sel.pipeline);  Sel.dataset = string(Sel.dataset);
pd = @(d, o) P(P.pipeline == "SG-IRLS" & P.dataset == d & P.outcome == o, :);
vals = [];
for d = DS
    for o = ["betaObs" "betaGen"]
        r = pd(d, o);  if height(r) ~= 1, error("tempoSlopeFig:p1", "%s", "P1 row missing: " + d + " " + o); end
        vals = [vals; r.b r.CIlo r.CIhi]; %#ok<AGROW>
    end
end
nNaN = sum(isnan(vals(:)));  fprintf("Panel C values: %d NaN\n", nNaN);
if nNaN > 0, error("tempoSlopeFig:nan", "%s", "NaN in panel C values"); end

%% FIGURE
fig = figure("Units", "centimeters", "Position", [2 2 30 11], "Color", "w");
tl = tiledlayout(fig, 1, 3, "TileSpacing", "compact", "Padding", "compact");
title(tl, "Within-session tempo slope of the measured exponent (change in \beta per doubling of drawing speed)", "FontSize", 11);
subtitle(tl, "v015 corpus (Fraser N = 45-era files; Fraser values may move with the N = 48 rerun)", "FontSize", 8);

% A: the analysis alone
ax = nexttile(tl);  hold(ax, "on");
xg = [1 2];  off = linspace(-0.28, 0.28, numel(ARMS));
for k = 1:numel(PIPES)
    for j = 1:numel(ARMS)
        r = A(A.arm == ARMS(j) & A.pipeline == PIPES(k), :);
        x = xg(k) + off(j);
        plot(ax, [x x], [r.empMatched r.synObs], "-", "Color", [0.6 0.6 0.6]);
        h1 = plot(ax, x, r.empMatched, "o", "MarkerFaceColor", C_OBS, "MarkerEdgeColor", C_OBS, "MarkerSize", 6);
        h2 = plot(ax, x, r.synObs, "s", "MarkerFaceColor", "w", "MarkerEdgeColor", C_SYN, "LineWidth", 1.5, "MarkerSize", 6);
    end
    sh = A.shareOld(A.pipeline == PIPES(k));
    text(ax, xg(k), 0.098, sprintf("synthetic / real\n%.0f-%.0f%% (mean %.0f%%)", 100*min(sh), 100*max(sh), 100*mean(sh)), ...
        "HorizontalAlignment", "center", "FontSize", 8);
end
set(ax, "XTick", xg, "XTickLabel", PLAB, "XLim", [0.5 2.5], "YLim", [0 0.11], "TickLabelInterpreter", "none");
ylabel(ax, "Slope (\Delta\beta per doubling of speed)");
legend(ax, [h1 h2], ["Real trials (same matched trials)", "Synthetic copies, true \beta unrelated to speed"], "Location", "southoutside", "FontSize", 8, "Box", "off");
title(ax, "A. The analysis alone produces most of it", "FontSize", 10);
text(ax, 0.55, 0.004, "each pair = one simulation arm (two seeds, three level shifts)", "FontSize", 7, "Color", [0.35 0.35 0.35]);

% B: pooled real vs beta_gen*
ax = nexttile(tl);  hold(ax, "on");
red = nan(1, 2);
for k = 1:numel(PIPES)
    o = pool(PIPES(k), "beta_obs");  g = pool(PIPES(k), "beta_gen*");
    if height(o) ~= 1 || height(g) ~= 1, error("tempoSlopeFig:pool", "%s", "Pooled row missing: " + PIPES(k)); end
    errorbar(ax, k - 0.12, o.b, o.b - o.CIlo, o.CIhi - o.b, "o", "Color", C_OBS, "MarkerFaceColor", C_OBS, "LineWidth", 1.2);
    errorbar(ax, k + 0.12, g.b, g.b - g.CIlo, g.CIhi - g.b, "d", "Color", C_GEN, "MarkerFaceColor", C_GEN, "LineWidth", 1.2);
    red(k) = 100 * (1 - g.b / o.b);
    text(ax, k, max(o.CIhi) + 0.008, sprintf("%.0f%% smaller", red(k)), "HorizontalAlignment", "center", "FontSize", 8);
end
yline(ax, 0, ":", "Color", [0.5 0.5 0.5]);
set(ax, "XTick", [1 2], "XTickLabel", PLAB, "XLim", [0.5 2.5], "YLim", [-0.01 0.11], "TickLabelInterpreter", "none");
legend(ax, ["Measured \beta (\beta_{obs})", "Corrected \beta (\beta_{gen}*, pointwise inversion)"], "Location", "southoutside", "FontSize", 8, "Box", "off");
title(ax, "B. Correction removes most of it (pooled)", "FontSize", 10);
text(ax, 0.55, -0.0075, "bars: model 95% CI; five in-domain datasets", "FontSize", 7, "Color", [0.35 0.35 0.35]);

% C: per dataset, vetted
ax = nexttile(tl);  hold(ax, "on");
for i = 1:numel(DS)
    o = pd(DS(i), "betaObs");  g = pd(DS(i), "betaGen");
    errorbar(ax, i - 0.12, o.b, o.b - o.CIlo, o.CIhi - o.b, "o", "Color", C_OBS, "MarkerFaceColor", C_OBS, "LineWidth", 1.1);
    errorbar(ax, i + 0.12, g.b, g.b - g.CIlo, g.CIhi - g.b, "d", "Color", C_GEN, "MarkerFaceColor", C_GEN, "LineWidth", 1.1);
    s = Sel(Sel.pipeline == "SG-IRLS" & Sel.dataset == DS(i), :);
    if height(s) ~= 1, error("tempoSlopeFig:sel", "%s", "Sel row missing: " + DS(i)); end
    text(ax, i, -0.035, sprintf("%d/%d\ninvertible", round(s.shareInvertible * s.nTrials), s.nTrials), "HorizontalAlignment", "center", "FontSize", 7);
end
yline(ax, 0, ":", "Color", [0.5 0.5 0.5]);
set(ax, "XTick", 1:numel(DS), "XTickLabel", DSLAB, "XLim", [0.4 numel(DS) + 0.6], "YLim", [-0.05 0.27], "TickLabelInterpreter", "none");
legend(ax, ["Measured \beta", "Corrected \beta (\beta_{gen}*)"], "Location", "southoutside", "FontSize", 8, "Box", "off");
title(ax, "C. By dataset (vetted)", "FontSize", 10);
text(ax, 0.5, 0.258, "per-dataset model, no noise covariates; Cook: few, slow trials", "FontSize", 7, "Color", [0.35 0.35 0.35]);

exportgraphics(fig, OUT_PNG, "Resolution", 300);  exportgraphics(fig, OUT_PDF, "ContentType", "vector");
fprintf("saved %s and %s\n", OUT_PNG, OUT_PDF);

%% CLAIMS CHECK against v008 Part 4 (quoted value, computed, tolerance = rounding of the prose)
vs = A(A.pipeline == "SG-IRLS", :);  vl = A(A.pipeline == "BWFD-OLS", :);
Q = { ...
 "pooled vetted beta_obs slope +0.084",  0.084,  pool("SG-IRLS", "beta_obs").b,  0.0005;
 "pooled vetted beta_gen* slope +0.0265", 0.0265, pool("SG-IRLS", "beta_gen*").b, 0.00005;
 "pooled legacy beta_obs slope +0.065",  0.065,  pool("BWFD-OLS", "beta_obs").b, 0.0005;
 "pooled legacy beta_gen* slope +0.0145", 0.0145, pool("BWFD-OLS", "beta_gen*").b, 0.00005;
 "vetted reduction 68%", 68, red(1), 0.5;
 "legacy reduction 78%", 78, red(2), 0.5;
 "vetted share min 57%", 57, 100*min(vs.shareOld), 0.5;
 "vetted share max 82%", 82, 100*max(vs.shareOld), 0.5;
 "vetted share mean 69%", 69, 100*mean(vs.shareOld), 0.5;
 "legacy share min 72%", 72, 100*min(vl.shareOld), 0.5;
 "legacy share max 120%", 120, 100*max(vl.shareOld), 0.5;
 "legacy share mean 96%", 96, 100*mean(vl.shareOld), 0.5;
 "vetted matched real slope min +0.057", 0.057, min(vs.empMatched), 0.0005;
 "vetted matched real slope max +0.083", 0.083, max(vs.empMatched), 0.0005;
 "Fraser vetted beta_gen* slope +0.012", 0.012, pd("Fraser", "betaGen").b, 0.0005};
fprintf("\nCLAIMS CHECK\n");  nd = 0;
for i = 1:size(Q, 1)
    ok = abs(Q{i,3} - Q{i,2}) <= Q{i,4};  nd = nd + ~ok;
    fprintf("%-72s quoted %8.4f computed %8.4f  %s\n", Q{i,1}, Q{i,2}, Q{i,3}, tern(ok, "MATCH", "DIFF"));
end
fprintf("%d of %d match\n", size(Q,1) - nd, size(Q,1));
R = struct("arms", A, "pooled", D, "perDataset", P(P.pipeline == "SG-IRLS", :), "sel", Sel, "reductionPct", red, ...
    "claims", cell2table(Q, "VariableNames", ["claim" "quoted" "computed" "tol"]), "runStamp", stamp);
save(OUT_MAT, "R", "-v7");  fprintf("saved %s\nlog %s\n", OUT_MAT, logf);
diary("off");   % a script's onCleanup lives in the base workspace, so close the log explicitly

function s = tern(c, a, b), if c, s = a; else, s = b; end, end
