%% extractGainTarget_v002.m
% Gain and target of every forward map (Session 123 vocabulary, claude.md v003 row).
%   gain   = how strongly beta_rec is pulled: slope g of the line beta_rec = c + g * beta_gen
%            fitted over FIT_WIN (plus a local gain near 1/3, LOC_WIN)
%   target = the generator value the map returns unchanged; departures from it are scaled by
%            the gain: beta_rec - T = g * (beta_gen - T). Estimated two ways:
%              line-fit target  T = c / (1 - g)  (fixed point of the fitted line)
%              crossing target  where the curve itself crosses the identity line (model-free;
%                               the crossing nearest the line-fit value)
% v002 (vs v001): a target is RELIABLE only if the pull is >= 10% (g <= GAIN_MAX), the map is
%   straight enough for a line (residual RMS <= RMSE_MAX; R2 is unusable for fully pulled, flat
%   maps), the crossing target exists and agrees with the line-fit target within AGREE_TOL,
%   the target lies inside the swept range +/- EXT (no extrapolation), and the cell is not
%   survivor-conditioned (failed fits <= FAIL_MAX %). Adds: sensitivity of the reliable share
%   to RMSE_MAX; within-dataset target vs each trial's own f0 (bins, Spearman per dataset, LMM
%   with relative noise, colour and subject); terminology "line-fit" and "crossing" targets.
% Predictions: (1) sigma = 0 at fixed tempo -> target 1/3 (the pull to one-third);
%   (2) noise at fixed tempo -> target beta_noise, varying with colour and pipeline (the pull to
%   the noise exponent); (3) real per-trial maps (v015): where their targets sit, and whether
%   they move towards 1/3 with the trial's tempo.
% Reads (local results/, then the RDS mount): tempoFoldSweep_v001, tempoFoldSweepC_v001,
%   tempoFoldSweepE_v001; src/loopClosureResults_<dataset>_all_shaped_xu_v015.mat
% Writes: results/gainTarget_v002.mat, figures/gainTarget_v002.png
% Edited in place 2026-09-26 (figure only): dataset tick labels had underscores rendered as
% TeX subscripts; now cleaned, with TickLabelInterpreter 'none' on that axis. Panel 3 labels
% carry n per dataset (Zarandi has no reliable SG-IRLS target). Results unchanged.

%% CONFIG
ROOT      = fileparts(fileparts(mfilename("fullpath")));
MOUNT     = "/Volumes/rdsprojects/f/fraserds-mpo-evaluation/2026_prereg/velocity-curvature-power-law-simulation-main/velocity-curvature-power-law-simulation-main";
OUT_MAT   = fullfile(ROOT, "results", "gainTarget_v002.mat");
OUT_PNG   = fullfile(ROOT, "figures", "gainTarget_v002.png");
FIT_WIN   = [0 0.75];
LOC_WIN   = [0.2 0.5];
GAIN_MAX  = 0.9;
RMSE_MAX  = 0.02;   RMSE_SENS = [0.01 0.02 0.04];
AGREE_TOL = 0.02;
EXT       = 0.05;
FAIL_MAX  = 10;
EDGES     = [0 0.35 0.7 1.4 2.8 Inf];  BIN_NAMES = ["<0.35" "0.35-0.7" "0.7-1.4" "1.4-2.8" ">2.8"];
OUT_DOMAIN = "Zarandi";                                                         % D20
LAB_ENG   = ["BWFD-OLS" "BWFD-LMLS" "BWFD-IRLS" "SG-OLS" "SG-LMLS" "SG-IRLS"];   % engine order
LAB_RUN   = ["BWFD-OLS" "SG-OLS" "BWFD-LMLS" "SG-LMLS" "BWFD-IRLS" "SG-IRLS"];   % v015 runner order (v012 L118)
V015      = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Dhieb" "Zarandi"];
for f = [OUT_MAT OUT_PNG], if ~isfolder(fileparts(f)), error("gainTarget:outDir", "%s", "Missing folder: " + fileparts(f)); end, end

%% Collect every curve (structs, one table at the end)
P = {};
R = loadFirst_local("tempoFoldSweep_v001.mat", ROOT, MOUNT).R;                   % Block A, sigma = 0
R = R(startsWith(R.block, "A") & ~R.lengthExcluded, :);
K = unique(R(:, ["block" "ab" "f0" "FS" "edgeClip" "nCycles"]), "rows");
for k = 1:height(K)
    c = sortrows(R(R.block == K.block(k) & R.ab == K.ab(k) & R.f0 == K.f0(k) & R.FS == K.FS(k) & ...
        R.edgeClip == K.edgeClip(k) & R.nCycles == K.nCycles(k), :), "betaGen");
    cond = sprintf("ab%.2f_fs%d_clip%d_ncyc%d", K.ab(k), K.FS(k), K.edgeClip(k), K.nCycles(k));
    for p = 1:6
        P{end+1, 1} = row_local("A", K.block(k), cond, "", 0, NaN, NaN, K.f0(k), LAB_ENG(p), c.betaGen, c.betaMean(:, p), 0, FIT_WIN, LOC_WIN); %#ok<SAGROW>
    end
end
R = loadFirst_local("tempoFoldSweepC_v001.mat", ROOT, MOUNT).R;                  % Blocks C/D
K = unique(R(:, ["block" "dataset" "sigma" "alpha" "f0"]), "rows");
for k = 1:height(K)
    c = sortrows(R(R.block == K.block(k) & R.dataset == K.dataset(k) & R.sigma == K.sigma(k) & ...
        R.alpha == K.alpha(k) & R.f0 == K.f0(k), :), "betaGen");
    br = c.betaRec;  if ~iscell(br), br = num2cell(br, 2); end
    for p = 1:6
        med = cellfun(@(x) median(x(:, p), "omitnan"), br);
        fp = 100 * sum(c.nFail(:, p)) / sum(c.nReps);
        P{end+1, 1} = row_local(K.block(k), K.dataset(k), K.dataset(k), "", K.sigma(k), K.alpha(k), c.a(1), K.f0(k), LAB_ENG(p), c.betaGen, med, fp, FIT_WIN, LOC_WIN); %#ok<SAGROW>
    end
end
R = loadFirst_local("tempoFoldSweepE_v001.mat", ROOT, MOUNT).R;                  % Block E
R.f0key = R.f0;  R.f0key(isnan(R.f0key)) = -1;                                  % replay: f0 realised, not set
K = unique(R(:, ["block" "sigmaMM" "alpha" "f0key"]), "rows");
K.f0 = K.f0key;  K.f0(K.f0 < 0) = NaN;
for k = 1:height(K)
    m = R.block == K.block(k) & R.sigmaMM == K.sigmaMM(k) & R.alpha == K.alpha(k) & R.f0key == K.f0key(k);
    c = sortrows(R(m & ~R.lengthExcluded, :), "betaGen");
    for p = 1:6
        fp = 100 * sum(c.nFail(:, p)) / sum(c.nReps);
        P{end+1, 1} = row_local("E_" + K.block(k), "grid ellipse", "E", "", K.sigmaMM(k), K.alpha(k), NaN, K.f0(k), LAB_ENG(p), c.betaGen, c.betaMean(:, p), fp, FIT_WIN, LOC_WIN); %#ok<SAGROW>
    end
end
nSim = numel(P);
for d = V015                                                                    % v015 per-trial maps
    f = fullfile(ROOT, "src", "loopClosureResults_" + d + "_all_shaped_xu_v015.mat");
    if ~isfile(f), error("gainTarget:input", "%s", "Missing: " + f); end
    L = load(f, "results", "betaGenVec");  bg = L.betaGenVec(:);
    for i = 1:numel(L.results)
        r = L.results(i);  if isempty(r.betaRecCurveMed), continue; end
        for p = 1:6
            y = r.betaRecCurveMed(p, :)';  if all(isnan(y)), continue; end
            P{end+1, 1} = row_local("v015", d, d, d + "_" + string(r.subjectID), r.sigmaMM, r.alphaIRA, r.a_mm, r.f0, LAB_RUN(p), bg, y, 0, FIT_WIN, LOC_WIN); %#ok<SAGROW>
        end
    end
end
P = struct2table(vertcat(P{:}));
P.masked = P.failPct > FAIL_MAX;
rel = @(t, rm) t.gain <= GAIN_MAX & t.rmse <= rm & t.nCross >= 1 & abs(t.targetLine - t.targetCross) <= AGREE_TOL & ...
               t.targetLine >= FIT_WIN(1) - EXT & t.targetLine <= FIT_WIN(2) + EXT & ~t.masked;
P.reliable = rel(P, RMSE_MAX);
P.target = NaN(height(P), 1);  P.target(P.reliable) = P.targetLine(P.reliable);
fprintf("Curves: %d simulated, %d v015 per-trial; masked (fit failures > %d%%): %d; reliable targets: %d\n", ...
    nSim, height(P) - nSim, FAIL_MAX, sum(P.masked), sum(P.reliable));
fprintf("Reliable share of pulled curves (gain <= %.2f) by RMSE_MAX:", GAIN_MAX);
pulled = P.gain <= GAIN_MAX & ~P.masked;
for rm = RMSE_SENS, fprintf("  %.2f -> %.1f%%", rm, 100 * sum(rel(P, rm)) / sum(pulled)); end
fprintf("\n");

%% (1) sigma = 0: the pull to one-third
A = P(P.source == "A" & P.reliable, :);
fprintf("\n(1) sigma = 0, fixed tempo: reliable targets %d of %d curves (%d pulled)\n", height(A), sum(P.source == "A"), sum(P.source == "A" & pulled));
fprintf("    target median %.4f, IQR %.4f-%.4f, range %.4f-%.4f; residual RMS median %.4f\n", median(A.target), ...
    prctile(A.target, 25), prctile(A.target, 75), min(A.target), max(A.target), median(A.rmse));
disp(groupsummary(A, "pipeline", ["median" "min" "max"], "target"))

%% (2) noise at fixed tempo: the pull to the noise exponent
E = P(P.source == "E_FIXED", :);
fprintf("\n(2) Block E, fixed tempo, sigma 10 mm: reliable target (NaN = unreliable or masked)\n");
for f0 = unique(E.f0)'
    t = E(E.sigma == 10 & E.f0 == f0, ["alpha" "pipeline" "target"]);
    w = unstack(t, "target", "pipeline", "VariableNamingRule", "preserve");  w{:, 2:end} = round(w{:, 2:end}, 3);
    fprintf("  f0 = %g Hz\n", f0);  disp(w)
end
fprintf("The Maoz coincidence: legacy OLS (BWFD-OLS) under white noise, fixed tempo <= 1 Hz:\n");
disp(P(P.source == "E_FIXED" & P.alpha == 0 & P.pipeline == "BWFD-OLS" & P.f0 <= 1 & P.sigma > 0, ...
    ["sigma" "f0" "gain" "targetLine" "targetCross" "rmse" "reliable"]))
fprintf("Blocks C and D (noisy): median gain; reliable targets (n, median)\n");
CD = P(ismember(P.source, ["C" "D"]) & P.sigma > 0, :);
disp(groupsummary(CD, ["source" "dataset" "alpha" "f0"], {"median", @(x) sum(isfinite(x))}, ["gain" "target"]))

%% (3) real per-trial maps
V = P(P.source == "v015", :);
S3 = groupsummary(V(ismember(V.pipeline, ["SG-IRLS" "BWFD-OLS"]), :), ["dataset" "pipeline"], ...
    {"median", @(x) mean(x)}, ["gain" "rmse" "reliable"]);
fprintf("\n(3) v015 per-trial maps (fun1 = mean): gain, residual RMS, share with a reliable target\n");  disp(S3)
T3 = groupsummary(V(V.reliable & ismember(V.pipeline, ["SG-IRLS" "BWFD-OLS"]), :), ["dataset" "pipeline"], ...
    {"median", @(x) prctile(x, 25), @(x) prctile(x, 75)}, "target");
fprintf("Reliable targets by dataset (fun1/fun2 = IQR):\n");  disp(T3)

%% (3b) does a trial's target move with its own tempo?
V.bin = discretize(V.f0, EDGES, 'categorical', BIN_NAMES);
Vs = V(V.pipeline == "SG-IRLS", :);
B3 = groupsummary(Vs, ["dataset" "bin"], {"median", @(x) sum(isfinite(x))}, ["gain" "target"]);
fprintf("\n(3b) SG-IRLS by dataset x trial-f0 bin: median gain; reliable targets (n, median)\n");  disp(B3)
fprintf("Spearman correlation of reliable target with trial f0, per dataset (unit = trial; trials are\n");
fprintf("nested in subjects, so p is approximate; the LMM below models the nesting):\n");
C3 = table();
for d = V015
    for p = ["SG-IRLS" "BWFD-OLS"]
        x = V(V.dataset == d & V.pipeline == p & V.reliable, :);
        if height(x) >= 10, [rho, pv] = corr(x.f0, x.target, "Type", "Spearman"); else, [rho, pv] = deal(NaN); end
        C3 = [C3; table(d, p, height(x), rho, pv, 'VariableNames', ["dataset" "pipeline" "n" "rho" "p"])]; %#ok<AGROW>
    end
end
disp(C3)
D = V(V.reliable & V.pipeline == "SG-IRLS" & V.dataset ~= OUT_DOMAIN & isfinite(V.a) & V.a > 0, :);
z = @(x) (x - mean(x, "omitnan")) ./ std(x, "omitnan");
D.lf = z(log(D.f0));  D.lsa = z(log(D.sigma ./ D.a));  D.al = z(D.alpha);
D.datasetC = categorical(D.dataset);  D.subjC = categorical(D.subject);
m3 = fitlme(D, "target ~ 1 + lf + lsa + al + datasetC + (1|subjC)");
cm = m3.Coefficients;  k = ismember(string(cm.Name), ["lf" "lsa" "al"]);
fprintf("\nLMM: reliable target (SG-IRLS, in-domain) ~ z(log f0) + z(log sigma/a) + z(alpha) + dataset + (1|subject)\n");
fprintf("unit = trial, n = %d trials, %d subjects\n", height(D), numel(unique(D.subjC)));
disp(table(string(cm.Name(k)), round(cm.Estimate(k), 4), round(cm.SE(k), 4), round(cm.tStat(k), 2), cm.DF(k), cm.pValue(k), ...
    round(cm.Lower(k), 4), round(cm.Upper(k), 4), 'VariableNames', ["term" "b" "SE" "t" "df" "p" "CIlo" "CIhi"]))

%% Figure
fg = figure("Color", "w", "Position", [60 60 1300 800]);  tl = tiledlayout(2, 2, "TileSpacing", "compact");
nexttile; hold on;
a = P(P.source == "A" & P.cond == "ab2.16_fs120_clip20_ncyc10", :);
for p = ["BWFD-OLS" "SG-IRLS"]
    x = sortrows(a(a.pipeline == p, :), "f0");
    plot(x.f0, x.target, "-o", "DisplayName", p + " target");  plot(x.f0, x.gain, ":", "DisplayName", p + " gain");
end
yline(1/3, "k--", "HandleVisibility", "off");  set(gca, "XScale", "log");  ylim([0 1.05]);  box on;
xlabel("f_0 (Hz)");  title("\sigma = 0: the pull to one-third");  legend("Location", "southwest");
nexttile; hold on;
for p = ["BWFD-OLS" "SG-OLS" "BWFD-IRLS" "SG-IRLS"]
    for f0 = [0.25 1]
        x = sortrows(E(E.sigma == 10 & E.f0 == f0 & E.pipeline == p, :), "alpha");
        plot(x.alpha, x.target, ifelse_local(f0 == 1, "-o", "--s"), "DisplayName", sprintf("%s, %g Hz", p, f0));
    end
end
yline(1/3, "k--", "HandleVisibility", "off");  ylim([-0.05 0.4]);  box on;
xlabel("noise colour \alpha");  ylabel("target \beta_{noise}");  title("Noise, \sigma = 10 mm: the noise exponent");  legend("Location", "eastoutside");
nexttile; hold on;
v = V(V.pipeline == "SG-IRLS" & V.reliable, :);
nBox = arrayfun(@(d) sum(v.dataset == d), V015);                  % n per dataset: an empty box must read as "no data"
boxchart(categorical(v.dataset, V015, strrep(V015, "_", " ") + " (n = " + nBox + ")"), v.target);  yline(1/3, "k--");  box on;
set(gca, "TickLabelInterpreter", "none");                    % data-derived tick labels: no TeX
ylabel("reliable target");  title("Real per-trial maps (v015), SG-IRLS");
nexttile; hold on;
for d = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO"]
    x = B3(B3.dataset == d, :);  plot(double(x.bin), x.median_target, "-o", "DisplayName", strrep(d, "_", " "));
end
yline(1/3, "k--", "HandleVisibility", "off");  xticks(1:numel(BIN_NAMES));  xticklabels(BIN_NAMES);  box on;
xlabel("trial f_0 (Hz)");  ylabel("median reliable target");  title("Target vs the trial's tempo (SG-IRLS; Dhieb and Zarandi omitted, see panel 3)");  legend("Location", "best");
title(tl, "Gain and target of every forward map (v002: reliable targets only)");
set(findall(fg, "Type", "axes"), "Toolbar", []);  exportgraphics(fg, OUT_PNG, "Resolution", 200);

save(OUT_MAT, "P", "S3", "T3", "B3", "C3", "m3", "FIT_WIN", "LOC_WIN", "GAIN_MAX", "RMSE_MAX", "AGREE_TOL", "EXT", "FAIL_MAX", "-v7.3");
fprintf("\nSaved: %s\nFigure: %s\n", OUT_MAT, OUT_PNG);

%% =========================================================================
function r = row_local(src, dataset, cond, subj, sigma, alpha, a, f0, pipe, b, y, failPct, fitWin, locWin)
    b = b(:);  y = y(:);
    ok = isfinite(y) & b >= fitWin(1) & b <= fitWin(2);  b = b(ok);  y = y(ok);
    [g, c, R2, rm, gl, Tl, Tc] = deal(NaN);  nC = 0;
    if numel(b) >= 4
        p = polyfit(b, y, 1);  g = p(1);  c = p(2);  res = y - polyval(p, b);
        rm = sqrt(mean(res.^2));  R2 = 1 - sum(res.^2) / sum((y - mean(y)).^2);
        w = b >= locWin(1) & b <= locWin(2);
        if sum(w) >= 3, q = polyfit(b(w), y(w), 1); gl = q(1); end
        if g < 1, Tl = c / (1 - g); end
        d = y - b;  s = find(d(1:end-1) .* d(2:end) < 0 | d(1:end-1) == 0);  nC = numel(s);
        if nC >= 1
            xs = zeros(nC, 1);
            for k = 1:nC
                j = s(k);
                if d(j) == 0, xs(k) = b(j); else, xs(k) = b(j) - d(j) * (b(j+1) - b(j)) / (d(j+1) - d(j)); end
            end
            if isfinite(Tl), [~, kk] = min(abs(xs - Tl)); else, kk = 1; end
            Tc = xs(kk);
        end
    end
    r = struct("source", string(src), "dataset", string(dataset), "cond", string(cond), "subject", string(subj), ...
        "sigma", sigma, "alpha", alpha, "a", a, "f0", f0, "pipeline", string(pipe), "nNodes", numel(b), "gain", g, ...
        "intercept", c, "targetLine", Tl, "targetCross", Tc, "nCross", nC, "R2", R2, "rmse", rm, "gainLocal", gl, "failPct", failPct);
end

function s = loadFirst_local(name, root, mount)
    cands = [fullfile(root, "results", name); fullfile(mount, "results", name)];
    i = find(isfile(cands), 1);
    if isempty(i), error("gainTarget:input", "%s", "Missing: " + strjoin(cands, " | ")); end
    fprintf("Input: %s\n", cands(i));  s = load(cands(i), "R");
end

function v = ifelse_local(c, a, b)
    if c, v = a; else, v = b; end
end
