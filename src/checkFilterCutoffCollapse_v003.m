%% checkFilterCutoffCollapse_v003.m
% The pull to one-third at zero noise, zoomed to the tempi people actually draw at.
% v001/v002 swept f0 from 0.25 to 16 Hz to show that gain collapses on f0/f_c (Finding #233,
% #249); manual drawing never reaches most of that range. Here f0 runs 0.50-2.00 Hz in 0.05 Hz
% steps (the registration intended ~0.5-2 Hz at beta = 1/3, prereg 3.1), with the same six
% filter settings, geometry (a/b 2.16), fs 240 Hz, sigma = 0, and v002's edge-safe Savitzky-Golay
% clip (half-window + 1; Butterworth clip 20). A third panel shows where the five in-domain
% datasets' trials sit on the same axis (per-trial f0 from the v015 corpora).
% Anchors (fail loud): (1) jobs at 0.5, 1 and 2 Hz, which are also v002 nodes, reproduce v002's
% betas in all six filters (OLS to 1e-12, IRLS to 1e-9; see TOL_A1); (2) the production Butterworth-OLS onset (f_c 10 Hz) and the
% edge-safe production Savitzky-Golay-OLS onset (0.171 s) lie within 3% of #249's fine-grid
% values, 1.83 and 1.841 Hz.
% Run from src/. Reads results/filterCutoffCollapse_v002.mat and the v015 corpora in src/.
% Writes results/filterCutoffZoom_v003.mat, figures/filterCutoffZoom_v003.png and a dated log
% in logs/ (fprintf, not diary: diary does not capture output run through the MCP).

%% CONFIG
ROOT     = fileparts(fileparts(mfilename("fullpath")));
IN_MAT   = fullfile(ROOT, "results", "filterCutoffCollapse_v002.mat");
OUT_MAT  = fullfile(ROOT, "results", "filterCutoffZoom_v003.mat");
OUT_PNG  = fullfile(ROOT, "figures", "filterCutoffZoom_v003.png");
FS       = 240;  A_MM = 50;  AB = 2.16;  NCYC = 10;  CLIP_BW = 20;
F0       = 0.50:0.05:2.00;
BETA     = 0:1/30:0.75;
GAIN_WIN = [0.2 0.5];  GAIN_REF = 0.9;
REG      = [3 5];  REG_NAMES = ["OLS" "IRLS"];
ONSET_249 = struct("BW10_OLS", 1.83, "SG41_OLS", 1.841);  TOL_249 = 0.03;   % relative
DS       = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO"];  VER = "v015";
SHOW_B   = [0 1/3 2/3];   % beta_gen values tabulated at production settings in the log
TOL_A1   = [1e-12 1e-9];  % anchor 1, [OLS IRLS]. OLS is closed-form; IRLS is an iterative fit, whose
                          % re-runs differed from v002 by up to 4.5e-11 in the first run of this script
                          % (2026-10-23, which used 1e-12 for both and stopped); 1e-9 is still far below
                          % any precision reported.

%% Setup
SRC = fullfile(ROOT, "src");
addpath(genpath(fullfile(SRC, "functions")));  addpath(genpath(fullfile(SRC, "req")));
logFile = fullfile(ROOT, "logs", "checkFilterCutoffCollapse_v003_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log");
LG = fopen(logFile, "w");  if LG < 0, error("zoom3:log", "%s", "Cannot open log: " + logFile); end
say_local(LG, "checkFilterCutoffCollapse_v003 run %s\n", string(datetime("now")));
if ~isfile(IN_MAT), error("zoom3:in", "%s", "Missing input: " + IN_MAT); end
V2 = load(IN_MAT, "J", "flt", "F0", "BETA", "FS");
if V2.FS ~= FS || ~isequal(V2.BETA, BETA), error("zoom3:cfg", "%s", "v002 configuration differs"); end
flt = V2.flt;
[fi, ki, bi] = ndgrid(1:height(flt), 1:numel(F0), 1:numel(BETA));
J = table(fi(:), F0(ki(:))', BETA(bi(:))', 'VariableNames', ["filt" "f0" "betaGen"]);
say_local(LG, "Jobs: %d (6 filters x %d tempos x %d beta_gen), sigma = 0\n", height(J), numel(F0), numel(BETA));

%% Anchor 1: v002 nodes at 0.5, 1 and 2 Hz, all six filters
kNode = find(abs(V2.F0 - 0.5) < 1e-9 | abs(V2.F0 - 1) < 1e-9 | abs(V2.F0 - 2) < 1e-9);
if numel(kNode) ~= 3, error("zoom3:anchor1", "%s", "0.5, 1 and 2 Hz are not all v002 nodes"); end
JA = V2.J(ismember(V2.J.f0, V2.F0(kNode)), :);
bA = NaN(height(JA), numel(REG));
for f = 1:height(flt)
    k = find(JA.filt == f);
    bA(k, :) = runJobs_local(JA(k, :), flt, clip_local(flt, f, FS, CLIP_BW), FS, A_MM, AB, NCYC, REG);
end
dOLS  = max(abs(bA(:, 1) - JA.betaOLS), [], "omitnan");
dIRLS = max(abs(bA(:, 2) - JA.betaIRLS), [], "omitnan");
nanMismatch = nnz(isnan(bA(:, 1)) ~= isnan(JA.betaOLS)) + nnz(isnan(bA(:, 2)) ~= isnan(JA.betaIRLS));
say_local(LG, "Anchor 1: max |diff| from v002, OLS %.3g (tolerance %.0g), IRLS %.3g (tolerance %.0g); NaN pattern mismatches %d\n", ...
    dOLS, TOL_A1(1), dIRLS, TOL_A1(2), nanMismatch);
if ~(dOLS <= TOL_A1(1)) || ~(dIRLS <= TOL_A1(2)) || nanMismatch > 0
    error("zoom3:anchor1", "%s", "v002 nodes not reproduced (see the line above)");
end
say_local(LG, "Anchor 1 passed: %d jobs at 0.5, 1 and 2 Hz reproduce v002.\n", height(JA));

%% Run the fine grid
bRec = NaN(height(J), numel(REG));
for f = 1:height(flt)
    k = find(J.filt == f);
    bRec(k, :) = runJobs_local(J(k, :), flt, clip_local(flt, f, FS, CLIP_BW), FS, A_MM, AB, NCYC, REG);
    say_local(LG, "%s %s: %d jobs, %d non-finite fits\n", flt.family(f), setting_local(flt, f), numel(k), nnz(~isfinite(bRec(k, :))));
end
J.betaOLS = bRec(:, 1);  J.betaIRLS = bRec(:, 2);

%% Gain against tempo and onsets
G = table();
for f = 1:height(flt)
    for k = 1:numel(F0)
        c = sortrows(J(J.filt == f & abs(J.f0 - F0(k)) < 1e-12, :), "betaGen");
        w = c.betaGen >= GAIN_WIN(1) & c.betaGen <= GAIN_WIN(2);  g = NaN(1, numel(REG));
        for q = 1:numel(REG)
            y = c.("beta" + REG_NAMES(q))(w);  bw = c.betaGen(w);  ok = isfinite(y);
            if sum(ok) >= 3, p = polyfit(bw(ok), y(ok), 1); g(q) = p(1); end
        end
        G = [G; table(f, flt.family(f), flt.fc(f), flt.Tw(f), F0(k), g(1), g(2), ...
            'VariableNames', ["filt" "family" "fc" "Tw" "f0" "gainOLS" "gainIRLS"])]; %#ok<AGROW>
    end
end
O = table();
for f = 1:height(flt)
    g = G(G.filt == f, :);
    for q = 1:numel(REG)
        O = [O; table(flt.family(f), setting_local(flt, f), REG_NAMES(q), crossing_local(g.f0, g.("gain" + REG_NAMES(q)), GAIN_REF), ...
            min(g.("gain" + REG_NAMES(q))), 'VariableNames', ["family" "setting" "regression" "onsetHz" "minGainInRange"])]; %#ok<AGROW>
    end
end
say_local(LG, "\nOnset (first tempo with gain below %.2f) within 0.5-2 Hz; NaN = gain stays above it:\n%s", GAIN_REF, evalc("disp(O)"));

%% Anchor 2: production onsets against #249's fine-grid values
oBW = O.onsetHz(O.family == "BW" & O.setting == "f_c 10 Hz" & O.regression == "OLS");
oSG = O.onsetHz(O.family == "SG" & O.setting == "T_w 0.171 s" & O.regression == "OLS");
say_local(LG, "Anchor 2: BW-OLS f_c 10 Hz onset %.3f Hz (#249: %.2f); SG-OLS 0.171 s onset %.3f Hz (#249: %.3f)\n", oBW, ONSET_249.BW10_OLS, oSG, ONSET_249.SG41_OLS);
if ~(abs(oBW / ONSET_249.BW10_OLS - 1) <= TOL_249) || ~(abs(oSG / ONSET_249.SG41_OLS - 1) <= TOL_249)
    error("zoom3:anchor2", "%s", "Production onsets do not reproduce Finding #249 within 3%");
end
say_local(LG, "Anchor 2 passed.\n");

%% What the production settings return, in exponent units (OLS)
fProd = [find(flt.family == "BW" & flt.fc == 10), find(flt.family == "SG" & abs(flt.Tw - 41/FS) < 1e-12)];
P = table();
for f = fProd
    for b = SHOW_B
        r = J(J.filt == f & abs(J.betaGen - b) < 1e-9 & ismember(round(J.f0, 4), [0.5 1 1.5 2]), :);
        r = sortrows(r, "f0");
        P = [P; table(setting_local(flt, f), b, r.betaOLS', 'VariableNames', ["setting" "betaGen" "betaRec_at_0p5_1_1p5_2Hz"])]; %#ok<AGROW>
    end
end
say_local(LG, "\nbeta_rec (OLS, zero noise) at 0.5, 1, 1.5 and 2 Hz, production settings:\n%s", evalc("disp(P)"));

%% Where real trials sit (per-trial f0, v015)
Tq = table();  f0all = cell(numel(DS), 1);
for d = 1:numel(DS)
    fm = fullfile(SRC, sprintf("loopClosureResults_%s_all_shaped_xu_%s.mat", DS(d), VER));
    if ~isfile(fm), error("zoom3:corpus", "%s", "Missing corpus: " + fm); end
    S = load(fm, "results");  v = [S.results.f0]';  v = v(isfinite(v));
    if isempty(v), error("zoom3:f0", "%s", "No finite f0 in " + DS(d)); end
    f0all{d} = v;
    Tq = [Tq; table(DS(d), numel(v), prctile(v, 10), median(v), prctile(v, 90), 100 * mean(v < 0.5), 100 * mean(v > 2), ...
        'VariableNames', ["dataset" "nTrials" "p10" "median" "p90" "pctBelow0p5" "pctAbove2"])]; %#ok<AGROW>
end
say_local(LG, "\nPer-trial tempo (Hz) in the five in-domain datasets, v015:\n%s", evalc("disp(Tq)"));

%% Figure
fg = figure("Color", "w", "Position", [80 80 1150 860]);  tl = tiledlayout(3, 2, "TileSpacing", "compact");
cols = [0.52 0.72 0.92; 0.13 0.47 0.80; 0.05 0.22 0.42];
oProd = [oBW oSG];
for fmy = ["BW" "SG"]
    nexttile([2 1]); hold on;
    sets = find(flt.family == fmy);
    for s = 1:numel(sets)
        g = G(G.filt == sets(s), :);  isProd = ismember(sets(s), fProd);
        plot(g.f0, g.gainOLS, "-", "Color", cols(s, :), "LineWidth", 1.2 + 1.2 * isProd, "DisplayName", setting_local(flt, sets(s)) + ", OLS" + ifelse_local(isProd, " (production)", ""));
        plot(g.f0, g.gainIRLS, "--", "Color", cols(s, :), "LineWidth", 1.0 + 0.8 * isProd, "DisplayName", setting_local(flt, sets(s)) + ", IRLS");
    end
    yline(GAIN_REF, ":", "gain 0.9", "LabelHorizontalAlignment", "left", "HandleVisibility", "off");
    xline(oProd(1 + (fmy == "SG")), "--k", "HandleVisibility", "off");
    xlim([0.5 2]);  ylim([0.6 1.02]);  box on;  grid on;
    xlabel("f_0 (Hz)");  ylabel("gain d\beta_{rec}/d\beta_{gen} near 1/3");
    title(ifelse_local(fmy == "BW", "Butterworth (zero-lag, 2nd order)", "Savitzky-Golay (order 4, edge-safe)"));
    legend("Location", "southwest", "FontSize", 8);
end
nexttile([1 2]); hold on;
for d = 1:numel(DS)
    q = Tq(d, :);  y = numel(DS) + 1 - d;
    plot([q.p10 q.p90], [y y], "-", "Color", [0.4 0.4 0.4], "LineWidth", 1.5, "HandleVisibility", "off");
    plot(prctile(f0all{d}, [25 75]), [y y], "-", "Color", [0.15 0.15 0.15], "LineWidth", 5, "HandleVisibility", "off");
    plot(q.median, y, "o", "MarkerFaceColor", "w", "MarkerEdgeColor", "k", "HandleVisibility", "off");
end
labs = strings(numel(DS), 1);
for d = 1:numel(DS)
    labs(d) = sprintf("%s (%.0f%% below 0.5 Hz, %.0f%% above 2 Hz)", strrep(DS(d), "_", " "), Tq.pctBelow0p5(d), Tq.pctAbove2(d));
end
xline(oProd(1), "--k", "BW onset", "LabelVerticalAlignment", "bottom", "HandleVisibility", "off");
xline(oProd(2), ":k", "SG onset", "LabelVerticalAlignment", "top", "HandleVisibility", "off");
set(gca, "YTick", 1:numel(DS), "YTickLabel", flip(labs), "TickLabelInterpreter", "none");
xlim([0.5 2]);  ylim([0.4 numel(DS) + 0.6]);  box on;  grid on;
xlabel("per-trial f_0 (Hz): 10th-90th percentile, interquartile range, median");
title("Where real trials sit (five in-domain datasets, v015)");
title(tl, "The pull to one-third over the tempi people draw at (zero noise, a/b 2.16, 240 Hz)");
set(findall(fg, "Type", "axes"), "Toolbar", []);  exportgraphics(fg, OUT_PNG, "Resolution", 200);
save(OUT_MAT, "J", "G", "O", "P", "Tq", "flt", "F0", "BETA", "FS", "GAIN_WIN", "GAIN_REF", "CLIP_BW", "-v7.3");
say_local(LG, "Saved: %s\nFigure: %s\nLog: %s\n", OUT_MAT, OUT_PNG, logFile);  fclose(LG);

%% =========================================================================
function b = runJobs_local(Jk, flt, clip, FS, A_MM, AB, NCYC, REG)
    % Zero-noise ellipse jobs through the production pipeline (as v001/v002).
    b = NaN(height(Jk), numel(REG));  ef = false(height(Jk), 1);
    fam = flt.family;  fc = flt.fc;  tw = flt.Tw;
    nW = 0;  if ~isempty(gcp("nocreate")), nW = gcp("nocreate").NumWorkers; end
    parfor (k = 1:height(Jk), nW)
        r = Jk(k, :);  bb = NaN(1, numel(REG));
        [x, y] = generatePowerLawEllipse_v001(A_MM, A_MM / AB, r.f0, FS, r.betaGen, 'nCycles', NCYC);
        if fam(r.filt) == "BW", ft = 2; fp = [2 fc(r.filt) 1];
        else, fr = round(tw(r.filt) * FS); fr = fr + mod(fr + 1, 2); ft = 4; fp = [4 fr]; end
        try
            [dx, dy] = differentiateKinematicsEBR(x(:), y(:), ft, fp, FS);
            vx = dx(clip:end-clip, 2); vy = dy(clip:end-clip, 2); ax = dx(clip:end-clip, 3); ay = dy(clip:end-clip, 3);
            sp = hypot(vx, vy);  kp = curvatureKinematicEBR(vx, vy, ax, ay);  ok = isfinite(sp) & isfinite(kp) & kp > 0;
            lm = [1, -1/3];
            for q = 1:numel(REG)
                [bq, vq] = regressDataEBR(sp(ok), kp(ok), REG(q), lm, 0, 0);
                bb(q) = bq;
                if q == 1 && isfinite(bq) && isfinite(vq), lm = [vq, bq]; end
            end
        catch
            ef(k) = true;
        end
        b(k, :) = bb;
    end
    if any(ef), fprintf("  %d trajectories threw an error (left NaN)\n", sum(ef)); end
end

function c = clip_local(flt, f, FS, clipBW)
    % Butterworth: production clip. Savitzky-Golay: half-window + 1 (edge-safe, #249).
    c = clipBW;
    if flt.family(f) == "SG", fr = round(flt.Tw(f) * FS); fr = fr + mod(fr + 1, 2); c = (fr - 1) / 2 + 1; end
end

function s = setting_local(flt, f)
    if flt.family(f) == "BW", s = sprintf("f_c %g Hz", flt.fc(f)); else, s = sprintf("T_w %.3f s", flt.Tw(f)); end
    s = string(s);
end

function f = crossing_local(f0, g, ref)
    % First tempo at which gain falls below ref (linear interpolation on the fine grid); NaN if never.
    i = find(g < ref, 1);
    if isempty(i), f = NaN; return, end
    if i == 1, error("zoom3:crossing", "%s", "Gain already below the reference at the lowest tempo"); end
    t = (ref - g(i-1)) / (g(i) - g(i-1));  f = f0(i-1) + t * (f0(i) - f0(i-1));
end

function v = ifelse_local(c, a, b)
    if c, v = a; else, v = b; end
end

function say_local(fid, fmt, varargin)
    fprintf(fmt, varargin{:});  fprintf(fid, fmt, varargin{:});
end
