%% checkFilterCutoffCollapse_v004.m
% v003 (the pull to one-third at zero noise, 0.50-2.00 Hz) widened to the empirical tempo range,
% with 0.5-2 Hz as the central, shaded band (Dagmar, 2026-10-23).
% Range: from the lowest 1st percentile to the highest 99th percentile of per-trial f0 in the
% five in-domain datasets (v015), rounded outwards. Tempi inside 0.5-2 Hz are v003's (0.05 Hz
% steps, reused); tempi outside it are new, log-spaced at 2^(1/8). Same six filter settings,
% geometry (a/b 2.16), fs 240 Hz, sigma = 0, edge-safe Savitzky-Golay clip (half-window + 1).
% Gain is now left undefined at any tempo where a fit in the beta_gen window failed (v003's IRLS
% kinks came partly from slopes fitted to fewer points); the number masked is logged.
% Anchors (fail loud): (1) a re-run at 0.5, 1.25 and 2 Hz reproduces v003 (OLS 1e-12, IRLS 1e-9);
% (2) production onsets within 3% of #249 (BW-OLS f_c 10 Hz 1.83 Hz; SG-OLS 0.171 s 1.841 Hz).
% Run from src/. Reads results/filterCutoffZoom_v003.mat and the v015 corpora in src/. Writes
% results/filterCutoffZoom_v004.mat, figures/filterCutoffZoom_v004.png and a dated log in logs/.

%% CONFIG
ROOT     = fileparts(fileparts(mfilename("fullpath")));
IN_MAT   = fullfile(ROOT, "results", "filterCutoffZoom_v003.mat");
OUT_MAT  = fullfile(ROOT, "results", "filterCutoffZoom_v004.mat");
OUT_PNG  = fullfile(ROOT, "figures", "filterCutoffZoom_v004.png");
FS       = 240;  A_MM = 50;  AB = 2.16;  NCYC = 10;  CLIP_BW = 20;
BAND     = [0.5 2];                 % central band (v003's fine grid)
STEP_OUT = 2^(1/8);                 % log step outside the band
PCT      = [1 99];                  % empirical range: lowest 1st and highest 99th percentile
BETA     = 0:1/30:0.75;
GAIN_WIN = [0.2 0.5];  GAIN_REF = 0.9;
REG      = [3 5];  REG_NAMES = ["OLS" "IRLS"];
TOL_A1   = [1e-12 1e-9];            % [OLS IRLS], as v003
ONSET_249 = struct("BW10_OLS", 1.83, "SG41_OLS", 1.841);  TOL_249 = 0.03;
DS       = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO"];  VER = "v015";

%% Setup
SRC = fullfile(ROOT, "src");
addpath(genpath(fullfile(SRC, "functions")));  addpath(genpath(fullfile(SRC, "req")));
logFile = fullfile(ROOT, "logs", "checkFilterCutoffCollapse_v004_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log");
LG = fopen(logFile, "w");  if LG < 0, error("zoom4:log", "%s", "Cannot open log: " + logFile); end
say_local(LG, "checkFilterCutoffCollapse_v004 run %s\n", string(datetime("now")));
if ~isfile(IN_MAT), error("zoom4:in", "%s", "Missing input: " + IN_MAT); end
V3 = load(IN_MAT, "J", "flt", "F0", "BETA", "FS");
if V3.FS ~= FS || ~isequal(V3.BETA, BETA), error("zoom4:cfg", "%s", "v003 configuration differs"); end
flt = V3.flt;

%% Empirical tempo range (per-trial f0, v015)
Tq = table();  f0all = cell(numel(DS), 1);
for d = 1:numel(DS)
    fm = fullfile(SRC, sprintf("loopClosureResults_%s_all_shaped_xu_%s.mat", DS(d), VER));
    if ~isfile(fm), error("zoom4:corpus", "%s", "Missing corpus: " + fm); end
    S = load(fm, "results");  v = [S.results.f0]';  v = v(isfinite(v));
    if isempty(v), error("zoom4:f0", "%s", "No finite f0 in " + DS(d)); end
    f0all{d} = v;
    Tq = [Tq; table(DS(d), numel(v), min(v), prctile(v, PCT(1)), prctile(v, 10), median(v), prctile(v, 90), prctile(v, PCT(2)), max(v), ...
        100 * mean(v < BAND(1)), 100 * mean(v > BAND(2)), ...
        'VariableNames', ["dataset" "nTrials" "min" "p1" "p10" "median" "p90" "p99" "max" "pctBelowBand" "pctAboveBand"])]; %#ok<AGROW>
end
lo = floor(min(Tq.p1) * 20) / 20;  hi = ceil(max(Tq.p99) * 4) / 4;
if ~(lo > 0 && lo < BAND(1) && hi > BAND(2)), error("zoom4:range", "%s", sprintf("Unexpected empirical range %.3f-%.3f Hz", lo, hi)); end
say_local(LG, "Per-trial tempo (Hz), five in-domain datasets, v015:\n%s", evalc("disp(Tq)"));
say_local(LG, "Swept range: %.2f-%.2f Hz (lowest %dst and highest %dth percentile, rounded outwards); band %.1f-%.1f Hz\n", lo, hi, PCT(1), PCT(2), BAND(1), BAND(2));

%% Tempo grid: v003 inside the band, new log-spaced tempi outside it
fLow  = BAND(1) ./ STEP_OUT .^ (1:ceil(log(BAND(1) / lo) / log(STEP_OUT)));  fLow = sort(fLow(fLow >= lo * 0.999));
fHigh = BAND(2) .* STEP_OUT .^ (1:ceil(log(hi / BAND(2)) / log(STEP_OUT)));  fHigh = fHigh(fHigh <= hi * 1.001);
if isempty(fLow) || isempty(fHigh), error("zoom4:grid", "%s", "Empty outer grid"); end
F0new = [fLow fHigh];  F0 = sort([F0new V3.F0]);
[fi, ki, bi] = ndgrid(1:height(flt), 1:numel(F0new), 1:numel(BETA));
Jn = table(fi(:), F0new(ki(:))', BETA(bi(:))', 'VariableNames', ["filt" "f0" "betaGen"]);
say_local(LG, "New tempi: %d below the band (%.3f-%.3f Hz), %d above (%.3f-%.3f Hz); %d new jobs; %d reused from v003\n", ...
    numel(fLow), fLow(1), fLow(end), numel(fHigh), fHigh(1), fHigh(end), height(Jn), height(V3.J));

%% Anchor 1: re-run three v003 tempi in all six filters
fA = V3.F0(ismember(round(V3.F0, 4), [0.5 1.25 2]));
if numel(fA) ~= 3, error("zoom4:anchor1", "%s", "0.5, 1.25 and 2 Hz are not all v003 tempi"); end
JA = V3.J(ismember(V3.J.f0, fA), :);  bA = NaN(height(JA), numel(REG));
for f = 1:height(flt)
    k = find(JA.filt == f);  bA(k, :) = runJobs_local(JA(k, :), flt, clip_local(flt, f, FS, CLIP_BW), FS, A_MM, AB, NCYC, REG);
end
dOLS = max(abs(bA(:, 1) - JA.betaOLS), [], "omitnan");  dIRLS = max(abs(bA(:, 2) - JA.betaIRLS), [], "omitnan");
nanMismatch = nnz(isnan(bA(:, 1)) ~= isnan(JA.betaOLS)) + nnz(isnan(bA(:, 2)) ~= isnan(JA.betaIRLS));
say_local(LG, "Anchor 1: max |diff| from v003, OLS %.3g (tolerance %.0g), IRLS %.3g (tolerance %.0g); NaN pattern mismatches %d\n", ...
    dOLS, TOL_A1(1), dIRLS, TOL_A1(2), nanMismatch);
if ~(dOLS <= TOL_A1(1)) || ~(dIRLS <= TOL_A1(2)) || nanMismatch > 0, error("zoom4:anchor1", "%s", "v003 tempi not reproduced"); end
say_local(LG, "Anchor 1 passed (%d jobs).\n", height(JA));

%% Run the new tempi
bN = NaN(height(Jn), numel(REG));
for f = 1:height(flt)
    k = find(Jn.filt == f);  bN(k, :) = runJobs_local(Jn(k, :), flt, clip_local(flt, f, FS, CLIP_BW), FS, A_MM, AB, NCYC, REG);
    say_local(LG, "%s %s: %d new jobs, %d non-finite fits\n", flt.family(f), setting_local(flt, f), numel(k), nnz(~isfinite(bN(k, :))));
end
Jn.betaOLS = bN(:, 1);  Jn.betaIRLS = bN(:, 2);
J = [V3.J(:, ["filt" "f0" "betaGen" "betaOLS" "betaIRLS"]); Jn];

%% Gain, masked where any fit in the window failed
G = table();  nMask = zeros(height(flt), numel(REG));
for f = 1:height(flt)
    for k = 1:numel(F0)
        c = sortrows(J(J.filt == f & abs(J.f0 - F0(k)) < 1e-12, :), "betaGen");
        w = c.betaGen >= GAIN_WIN(1) - 1e-12 & c.betaGen <= GAIN_WIN(2) + 1e-12;  g = NaN(1, numel(REG));
        for q = 1:numel(REG)
            y = c.("beta" + REG_NAMES(q))(w);
            if all(isfinite(y)), p = polyfit(c.betaGen(w), y, 1); g(q) = p(1); else, nMask(f, q) = nMask(f, q) + 1; end
        end
        G = [G; table(f, flt.family(f), F0(k), g(1), g(2), 'VariableNames', ["filt" "family" "f0" "gainOLS" "gainIRLS"])]; %#ok<AGROW>
    end
end
M = table(string(arrayfun(@(f) setting_local(flt, f), (1:height(flt))', "UniformOutput", false)), nMask(:, 1), nMask(:, 2), 'VariableNames', ["setting" "maskedOLS" "maskedIRLS"]);
say_local(LG, "\nTempi with gain left undefined (a fit in the window failed), of %d:\n%s", numel(F0), evalc("disp(M)"));
nAbove1 = nnz(G.gainIRLS > 1.005);
say_local(LG, "IRLS gains above 1.005 after masking: %d of %d (zero-noise bisquare weighting; not a filter effect)\n", nAbove1, nnz(isfinite(G.gainIRLS)));

%% Onsets and anchor 2
O = table();
for f = 1:height(flt)
    g = G(G.filt == f, :);
    for q = 1:numel(REG)
        O = [O; table(flt.family(f), setting_local(flt, f), REG_NAMES(q), crossing_local(g.f0, g.("gain" + REG_NAMES(q)), GAIN_REF), ...
            'VariableNames', ["family" "setting" "regression" "onsetHz"])]; %#ok<AGROW>
    end
end
say_local(LG, "\nOnset (first tempo with gain below %.2f; NaN = never within %.2f-%.2f Hz):\n%s", GAIN_REF, lo, hi, evalc("disp(O)"));
oBW = O.onsetHz(O.family == "BW" & O.setting == "f_c 10 Hz" & O.regression == "OLS");
oSG = O.onsetHz(O.family == "SG" & O.setting == "T_w 0.171 s" & O.regression == "OLS");
if ~(abs(oBW / ONSET_249.BW10_OLS - 1) <= TOL_249) || ~(abs(oSG / ONSET_249.SG41_OLS - 1) <= TOL_249)
    error("zoom4:anchor2", "%s", sprintf("Production onsets %.3f and %.3f Hz do not reproduce Finding #249 within 3%%", oBW, oSG));
end
say_local(LG, "Anchor 2 passed: BW-OLS f_c 10 Hz %.3f Hz, SG-OLS 0.171 s %.3f Hz.\n", oBW, oSG);

%% Figure
fg = figure("Color", "w", "Position", [80 80 1150 860]);  tl = tiledlayout(3, 2, "TileSpacing", "compact");
cols = [0.52 0.72 0.92; 0.13 0.47 0.80; 0.05 0.22 0.42];
fProd = [find(flt.family == "BW" & flt.fc == 10), find(flt.family == "SG" & abs(flt.Tw - 41/FS) < 1e-12)];
oProd = [oBW oSG];  xl = [lo hi];  xt = [0.1 0.2 0.3 0.5 1 1.5 2 3 4];  xt = xt(xt >= lo & xt <= hi);
for fmy = ["BW" "SG"]
    nexttile([2 1]); hold on;
    band_local(BAND, [0 1.05]);
    sets = find(flt.family == fmy);
    for s = 1:numel(sets)
        g = G(G.filt == sets(s), :);  isProd = ismember(sets(s), fProd);
        plot(g.f0, g.gainOLS, "-", "Color", cols(s, :), "LineWidth", 1.2 + 1.2 * isProd, "DisplayName", setting_local(flt, sets(s)) + ", OLS" + ifelse_local(isProd, " (production)", ""));
        plot(g.f0, g.gainIRLS, "--", "Color", cols(s, :), "LineWidth", 1.0 + 0.8 * isProd, "DisplayName", setting_local(flt, sets(s)) + ", IRLS");
    end
    yline(GAIN_REF, ":", "gain 0.9", "LabelHorizontalAlignment", "left", "HandleVisibility", "off");
    xline(oProd(1 + (fmy == "SG")), "--k", "HandleVisibility", "off");
    set(gca, "XScale", "log", "XTick", xt);  xlim(xl);  ylim([0 1.05]);  box on;  grid on;
    xlabel("f_0 (Hz)");  ylabel("gain d\beta_{rec}/d\beta_{gen} near 1/3");
    title(ifelse_local(fmy == "BW", "Butterworth (zero-lag, 2nd order)", "Savitzky-Golay (order 4, edge-safe)"));
    legend("Location", "southwest", "FontSize", 8);
end
nexttile([1 2]); hold on;
band_local(BAND, [0.4 numel(DS) + 0.6]);
labs = strings(numel(DS), 1);
for d = 1:numel(DS)
    q = Tq(d, :);  y = numel(DS) + 1 - d;
    plot([q.p10 q.p90], [y y], "-", "Color", [0.4 0.4 0.4], "LineWidth", 1.5, "HandleVisibility", "off");
    plot(prctile(f0all{d}, [25 75]), [y y], "-", "Color", [0.15 0.15 0.15], "LineWidth", 5, "HandleVisibility", "off");
    plot(q.median, y, "o", "MarkerFaceColor", "w", "MarkerEdgeColor", "k", "HandleVisibility", "off");
    labs(d) = sprintf("%s (%.0f%% below %.1f Hz, %.0f%% above %.0f Hz)", strrep(DS(d), "_", " "), q.pctBelowBand, BAND(1), q.pctAboveBand, BAND(2));
end
xline(oProd(1), "--k", "production onsets", "LabelVerticalAlignment", "bottom", "HandleVisibility", "off");
set(gca, "XScale", "log", "XTick", xt, "YTick", 1:numel(DS), "YTickLabel", flip(labs), "TickLabelInterpreter", "none");
xlim(xl);  ylim([0.4 numel(DS) + 0.6]);  box on;  grid on;
xlabel("per-trial f_0 (Hz): 10th-90th percentile, interquartile range, median");
title(sprintf("Where real trials sit (five in-domain datasets, v015); shaded: %.1f-%.0f Hz", BAND(1), BAND(2)));
title(tl, "The pull to one-third across the empirical tempo range (zero noise, a/b 2.16, 240 Hz)");
set(findall(fg, "Type", "axes"), "Toolbar", []);  exportgraphics(fg, OUT_PNG, "Resolution", 200);
save(OUT_MAT, "J", "G", "O", "M", "Tq", "flt", "F0", "BETA", "FS", "GAIN_WIN", "GAIN_REF", "CLIP_BW", "BAND", "lo", "hi", "-v7.3");
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
    if i == 1, error("zoom4:crossing", "%s", "Gain already below the reference at the lowest tempo"); end
    t = (ref - g(i-1)) / (g(i) - g(i-1));  f = f0(i-1) + t * (f0(i) - f0(i-1));
end

function v = ifelse_local(c, a, b)
    if c, v = a; else, v = b; end
end

function say_local(fid, fmt, varargin)
    fprintf(fmt, varargin{:});  fprintf(fid, fmt, varargin{:});
end

function band_local(band, yl)
    % Shade the central tempo band.
    patch([band(1) band(2) band(2) band(1)], [yl(1) yl(1) yl(2) yl(2)], [0.93 0.93 0.93], "EdgeColor", "none", "HandleVisibility", "off");
end
