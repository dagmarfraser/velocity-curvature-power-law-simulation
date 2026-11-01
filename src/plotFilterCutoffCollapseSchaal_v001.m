%% plotFilterCutoffCollapseSchaal_v001.m
% Fig 1B-2 with Schaal & Sternad's (2001) demonstration placed on the collapse curve.
% Their case (Exp Brain Res 136:60-72, p. 71): a constant-speed ellipse (true beta = 0), 0.2 x 0.1 m
% axes (a/b 2), period 0.6 s (f0 = 1.667 Hz), passed through a fourth-order Butterworth at 2.5 Hz,
% read beta = 0.328. On the scaled axis that is f0/f_c = 0.667, with an effective gain
% (beta_rec - 1/3)/(beta_gen - 1/3) = 0.016. This script plots it beside our zero-noise curves
% (zero-lag second-order Butterworth, the production filter) and prints, from the saved v002 run:
% our gain near 1/3 and our effective gain at beta_gen = 0 at their f0/f_c, and the f0/f_c at which
% Butterworth-OLS gain falls to 0.9, 0.8 and 0.5. No simulation is run.
% Their filter's pass (single or zero-lag) is not stated in the paper; magnitudes at their third
% harmonic are close either way (|H| 0.062 single-pass fourth order, 0.059 zero-lag second order).
% Run from src/. Reads results/filterCutoffCollapse_v002.mat. Writes
% figures/filterCutoffCollapseSchaal_v001.png, results/plotFilterCutoffCollapseSchaal_v001.mat and a
% dated log in logs/ (fprintf, not diary).

%% CONFIG
ROOT    = fileparts(fileparts(mfilename("fullpath")));
IN_MAT  = fullfile(ROOT, "results", "filterCutoffCollapse_v002.mat");
OUT_PNG = fullfile(ROOT, "figures", "filterCutoffCollapseSchaal_v001.png");
OUT_MAT = fullfile(ROOT, "results", "plotFilterCutoffCollapseSchaal_v001.mat");
SS      = struct("period", 0.6, "fc", 2.5, "order", 4, "betaGen", 0, "betaRec", 0.328, "ab", 2);   % p. 71
LEVELS  = [0.9 0.8 0.5];
FC_PROD = 10;

%% Load and derive
logFile = fullfile(ROOT, "logs", "plotFilterCutoffCollapseSchaal_v001_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log");
LG = fopen(logFile, "w");  if LG < 0, error("ssplot:log", "%s", "Cannot open log: " + logFile); end
say_local(LG, "plotFilterCutoffCollapseSchaal_v001 run %s\n", string(datetime("now")));
if ~isfile(IN_MAT), fclose(LG); error("ssplot:in", "%s", "Missing input: " + IN_MAT); end
S = load(IN_MAT, "J", "G", "flt", "F0", "GAIN_REF");
G = S.G;  J = S.J;  flt = S.flt;
xSS = (1 / SS.period) / SS.fc;
gSS = (SS.betaRec - 1/3) / (SS.betaGen - 1/3);
say_local(LG, "Schaal & Sternad (2001): f0 %.3f Hz, f_c %.1f Hz, f0/f_c %.3f; beta %g -> %.3f, effective gain %.4f\n", ...
    1 / SS.period, SS.fc, xSS, SS.betaGen, SS.betaRec, gSS);
% magnitude response at their third harmonic, for the two filter forms (stated in the header)
r = 3 * xSS;
say_local(LG, "|H| at 3 f0: single-pass 4th order %.4f; zero-lag 2nd order %.4f\n", 1 / sqrt(1 + r^8), 1 / (1 + r^4));

fBW = find(flt.family == "BW")';
T = table();
for f = fBW
    g = G(G.family == "BW" & G.fc == flt.fc(f), :);  x = g.f0 / flt.fc(f);
    if ~(min(x) < xSS && max(x) > xSS), fclose(LG); error("ssplot:range", "%s", "Their f0/f_c lies outside the swept range"); end
    gNear = interp1(log(x), g.gainOLS, log(xSS));
    c = J(J.filt == f & abs(J.betaGen) < 1e-12, :);  c = sortrows(c, "f0");
    b0 = interp1(log(c.f0 / flt.fc(f)), c.betaOLS, log(xSS));
    cr = arrayfun(@(L) crossing_local(x, g.gainOLS, L), LEVELS);
    T = [T; table(flt.fc(f), gNear, b0, 1 - b0 / (1/3), cr(1), cr(2), cr(3), ...
        'VariableNames', ["fc" "gainNear13_atSS" "betaRec_at_beta0_atSS" "effGain_beta0_atSS" "x_gain09" "x_gain08" "x_gain05"])]; %#ok<AGROW>
end
say_local(LG, "\nOur zero-noise Butterworth-OLS curves (zero-lag 2nd order, a/b 2.16, 240 Hz), from v002:\n%s", evalc("disp(T)"));
p = T(T.fc == FC_PROD, :);
say_local(LG, "Production (f_c %d Hz): gain 0.9, 0.8, 0.5 at f0 = %.2f, %.2f, %.2f Hz\n", FC_PROD, p.x_gain09 * FC_PROD, p.x_gain08 * FC_PROD, p.x_gain05 * FC_PROD);

%% Figure (v002 layout; their point on the scaled Butterworth panel)
fg = figure("Color", "w", "Position", [80 80 1100 750]);  tl = tiledlayout(2, 2, "TileSpacing", "compact");
cols = [0.52 0.72 0.92; 0.22 0.54 0.87; 0.05 0.27 0.49];
for row = 1:2
    for fmy = ["BW" "SG"]
        ax = nexttile; hold on;
        sets = flt(flt.family == fmy, :);
        for s = 1:height(sets)
            g = G(G.family == fmy & eqn_local(G.fc, sets.fc(s)) & eqn_local(G.Tw, sets.Tw(s)), :);
            xv = g.f0;  if row == 2, xv = g.x; end
            lab = sprintf("f_c %g Hz", sets.fc(s));  if fmy == "SG", lab = sprintf("T_w %.3f s", sets.Tw(s)); end
            plot(xv, g.gainOLS, "-o", "Color", cols(s, :), "LineWidth", 1.4, "MarkerSize", 4, "DisplayName", lab);
        end
        yline(S.GAIN_REF, ":", "HandleVisibility", "off");  set(ax, "XScale", "log");  ylim([0 1.1]);  box on;
        if row == 1, xlabel("f_0 (Hz)"); elseif fmy == "BW", xlabel("f_0 / f_c"); else, xlabel("f_0 \times T_w"); end
        ylabel("gain d\beta_{rec}/d\beta_{gen} (OLS, near 1/3)");
        ttl = "Butterworth";  if fmy == "SG", ttl = "Savitzky-Golay (edge-safe clip)"; end
        if row == 1, ttl = ttl + ", raw tempo"; else, ttl = ttl + ", tempo scaled by the filter"; end
        if row == 2 && fmy == "BW"
            xline(xSS, ":", "HandleVisibility", "off");
            plot(xSS, gSS, "kd", "MarkerFaceColor", [0.85 0.33 0.10], "MarkerSize", 9, ...
                "DisplayName", "Schaal & Sternad (2001)");
            text(xSS * 1.15, 0.30, {"Schaal & Sternad (2001)", "1.67 Hz ellipse, 2.5 Hz", "4th-order Butterworth", "\beta 0 \rightarrow 0.328"}, "FontSize", 8);
        end
        loc = "southwest";
        title(ttl);  legend("Location", loc);
    end
end
title(tl, "The pull to one-third at zero noise: gain collapses when tempo is scaled by the filter");
set(findall(fg, "Type", "axes"), "Toolbar", []);  exportgraphics(fg, OUT_PNG, "Resolution", 200);
save(OUT_MAT, "T", "SS", "xSS", "gSS", "LEVELS", "FC_PROD");
say_local(LG, "Figure: %s\nSaved: %s\nLog: %s\n", OUT_PNG, OUT_MAT, logFile);  fclose(LG);

%% =========================================================================
function x0 = crossing_local(x, g, ref)
    % First f0/f_c at which gain falls below ref (log-linear interpolation); NaN if never.
    i = find(g < ref, 1);
    if isempty(i) || i == 1, x0 = NaN; return, end
    t = (ref - g(i-1)) / (g(i) - g(i-1));  x0 = exp(log(x(i-1)) + t * (log(x(i)) - log(x(i-1))));
end

function m = eqn_local(v, s)
    m = (isnan(v) & isnan(s)) | v == s;
end

function say_local(fid, fmt, varargin)
    fprintf(fmt, varargin{:});  fprintf(fid, fmt, varargin{:});
end
