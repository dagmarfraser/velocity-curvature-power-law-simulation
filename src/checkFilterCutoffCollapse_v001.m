%% checkFilterCutoffCollapse_v001.m
% Is the zero-noise pull towards 1/3 (Finding #225) smoothing attenuation, i.e. the
% derivation stage's filter removing the movement's own harmonics? If so, the gain of the
% forward map against tempo should depend on tempo only relative to the filter:
%   Butterworth (BWFD, filtfilt 2nd order): collapse against f0 / fc
%   Savitzky-Golay (order 4): collapse against f0 * Tw (window length in seconds)
% Test: f0_90 = tempo at which the local gain first falls below GAIN_REF; the scaled value
% (f0_90/fc, or f0_90*Tw) should be constant across the three settings of each family.
% sigma = 0 (no noise anywhere), exact ellipse, fixed tempo, as Block A.
% Gain = OLS/IRLS slope of beta_rec on beta_gen over GAIN_WIN (local, around 1/3).
% Pipeline: differentiateKinematicsEBR -> curvatureKinematicEBR -> regressDataEBR
% (limitBreak = 0; OLS-first seeding as the runner). Runs locally (parfor if a pool exists).
% Writes results/filterCutoffCollapse_v001.mat, figures/filterCutoffCollapse_v001.png

%% CONFIG
ROOT     = fileparts(fileparts(mfilename("fullpath")));
OUT_MAT  = fullfile(ROOT, "results", "filterCutoffCollapse_v001.mat");
OUT_PNG  = fullfile(ROOT, "figures", "filterCutoffCollapse_v001.png");
FS       = 240;
A_MM     = 50;  AB = 2.16;                 % Block A geometry (sigma = 0 is scale-invariant)
NCYC     = 10;  CLIP = 20;
F0       = logspace(log10(0.25), log10(16), 19);
BETA     = 0:1/30:0.75;
BW_FC    = [5 10 20];                      % Hz; production is 10
SG_TW    = [21 41 81] / FS;                % s;  production is 0.17 s (41 samples at 240 Hz)
GAIN_WIN = [0.2 0.5];                      % beta_gen range for the local gain, around 1/3
GAIN_REF = 0.9;
REG      = [3 5];  REG_NAMES = ["OLS" "IRLS"];

%% Setup
addpath(genpath(fullfile(ROOT, "src", "functions")));  addpath(genpath(fullfile(ROOT, "src", "req")));
for f = [OUT_MAT OUT_PNG], if ~isfolder(fileparts(f)), error("filterCollapse:outDir", "%s", "Missing folder: " + fileparts(f)); end, end
flt = [table(repmat("BW", 3, 1), BW_FC(:), NaN(3, 1), 'VariableNames', ["family" "fc" "Tw"]); ...
       table(repmat("SG", 3, 1), NaN(3, 1), SG_TW(:), 'VariableNames', ["family" "fc" "Tw"])];
[fi, ki, bi] = ndgrid(1:height(flt), 1:numel(F0), 1:numel(BETA));
J = table(fi(:), F0(ki(:))', BETA(bi(:))', 'VariableNames', ["filt" "f0" "betaGen"]);
fprintf("Jobs: %d (6 filters x %d tempos x %d beta_gen), sigma = 0\n", height(J), numel(F0), numel(BETA));

%% Run
bRec = NaN(height(J), numel(REG));  nFail = zeros(height(J), 1);  errFlag = false(height(J), 1);
fam = flt.family;  fc = flt.fc;  tw = flt.Tw;
nW = 0;  if ~isempty(gcp("nocreate")), nW = gcp("nocreate").NumWorkers; end
parfor (k = 1:height(J), nW)
    r = J(k, :);  b = NaN(1, numel(REG));  ef = false;
    [x, y] = generatePowerLawEllipse_v001(A_MM, A_MM / AB, r.f0, FS, r.betaGen, 'nCycles', NCYC);
    if fam(r.filt) == "BW", ft = 2; fp = [2 fc(r.filt) 1];
    else, fr = round(tw(r.filt) * FS); fr = fr + mod(fr + 1, 2); ft = 4; fp = [4 fr]; end
    try
        [dx, dy] = differentiateKinematicsEBR(x(:), y(:), ft, fp, FS);
        vx = dx(CLIP:end-CLIP, 2); vy = dy(CLIP:end-CLIP, 2); ax = dx(CLIP:end-CLIP, 3); ay = dy(CLIP:end-CLIP, 3);
        sp = hypot(vx, vy);  kp = curvatureKinematicEBR(vx, vy, ax, ay);  ok = isfinite(sp) & isfinite(kp) & kp > 0;
        lm = [1, -1/3];
        for q = 1:numel(REG)
            [bb, vv] = regressDataEBR(sp(ok), kp(ok), REG(q), lm, 0, 0);
            b(q) = bb;
            if q == 1 && isfinite(bb) && isfinite(vv), lm = [vv, bb]; end
        end
    catch
        ef = true;       % counted and reported below; never silently absorbed
    end
    bRec(k, :) = b;  nFail(k) = sum(~isfinite(b));  errFlag(k) = ef;
end
J.betaOLS = bRec(:, 1);  J.betaIRLS = bRec(:, 2);  J.nFail = nFail;  J.err = errFlag;
fprintf("Non-finite estimates: %d of %d regressions; trajectories that threw an error: %d\n", ...
    sum(nFail), height(J) * numel(REG), sum(errFlag));

%% Gain against tempo, and the collapse test
G = table();
for f = 1:height(flt)
    for k = 1:numel(F0)
        c = sortrows(J(J.filt == f & J.f0 == F0(k), :), "betaGen");
        w = c.betaGen >= GAIN_WIN(1) & c.betaGen <= GAIN_WIN(2);
        g = NaN(1, numel(REG));
        for q = 1:numel(REG)
            bw = c.betaGen(w);  y = c.("beta" + REG_NAMES(q))(w);  ok = isfinite(y);
            if sum(ok) >= 3, p = polyfit(bw(ok), y(ok), 1); g(q) = p(1); end
        end
        [~, i23] = min(abs(c.betaGen - 2/3));
        G = [G; table(flt.family(f), flt.fc(f), flt.Tw(f), F0(k), g(1), g(2), c.betaOLS(i23), ...
            'VariableNames', ["family" "fc" "Tw" "f0" "gainOLS" "gainIRLS" "bOLSat23"])]; %#ok<AGROW>
    end
end
G.x = G.f0 ./ G.fc;  G.x(G.family == "SG") = G.f0(G.family == "SG") .* G.Tw(G.family == "SG");

C = table();
for f = 1:height(flt)
    g = G(G.family == flt.family(f) & isequaln_local(G.fc, flt.fc(f)) & isequaln_local(G.Tw, flt.Tw(f)), :);
    for q = 1:numel(REG)
        f90 = crossing_local(g.f0, g.("gain" + REG_NAMES(q)), GAIN_REF);
        sc = f90 / flt.fc(f);  if flt.family(f) == "SG", sc = f90 * flt.Tw(f); end
        C = [C; table(flt.family(f), flt.fc(f), flt.Tw(f), REG_NAMES(q), f90, sc, ...
            'VariableNames', ["family" "fc" "Tw" "regression" "f0_at_gain90" "scaled"])]; %#ok<AGROW>
    end
end
fprintf("\nTempo at which the local gain falls to %.2f, and its scaled value (f0/fc for BW, f0*Tw for SG):\n", GAIN_REF);
disp(C)
S = groupsummary(C, ["family" "regression"], ["mean" "std"], "scaled");  S.cv = S.std_scaled ./ S.mean_scaled;
fprintf("Collapse test (coefficient of variation of the scaled value across the three settings; small = collapse):\n");
disp(S(:, ["family" "regression" "mean_scaled" "std_scaled" "cv"]))
fprintf("Unscaled, for contrast: CV of f0_at_gain90 itself:\n");
U = groupsummary(C, ["family" "regression"], ["mean" "std"], "f0_at_gain90");  U.cv = U.std_f0_at_gain90 ./ U.mean_f0_at_gain90;
disp(U(:, ["family" "regression" "cv"]))

%% Figure: raw (top) and scaled (bottom), OLS gain
fg = figure("Color", "w", "Position", [80 80 1100 750]);  tl = tiledlayout(2, 2, "TileSpacing", "compact");
cols = [0.52 0.72 0.92; 0.22 0.54 0.87; 0.05 0.27 0.49];
for row = 1:2
    for fmy = ["BW" "SG"]
        nexttile; hold on;
        sets = flt(flt.family == fmy, :);
        for s = 1:height(sets)
            g = G(G.family == fmy & isequaln_local(G.fc, sets.fc(s)) & isequaln_local(G.Tw, sets.Tw(s)), :);
            xv = g.f0;  if row == 2, xv = g.x; end
            lab = sprintf("f_c %g Hz", sets.fc(s));  if fmy == "SG", lab = sprintf("T_w %.3f s", sets.Tw(s)); end
            plot(xv, g.gainOLS, "-o", "Color", cols(s, :), "LineWidth", 1.4, "MarkerSize", 4, "DisplayName", lab);
        end
        yline(GAIN_REF, ":", "HandleVisibility", "off");  set(gca, "XScale", "log");  ylim([0 1.1]);  box on;
        if row == 1, xlabel("f_0 (Hz)");
        elseif fmy == "BW", xlabel("f_0 / f_c");
        else, xlabel("f_0 \times T_w"); end
        ylabel("gain d\beta_{rec}/d\beta_{gen} (OLS, near 1/3)");
        title(sprintf("%s, %s", ifelse_local(fmy == "BW", "Butterworth", "Savitzky-Golay"), ifelse_local(row == 1, "raw tempo", "tempo scaled by the filter")));
        legend("Location", "southwest");
    end
end
title(tl, "Smoothing attenuation at zero noise: does gain collapse when tempo is scaled by the filter?");
set(findall(fg, "Type", "axes"), "Toolbar", []);  exportgraphics(fg, OUT_PNG, "Resolution", 200);

save(OUT_MAT, "J", "G", "C", "S", "U", "flt", "F0", "BETA", "FS", "GAIN_WIN", "GAIN_REF", "-v7.3");
fprintf("Saved: %s\nFigure: %s\n", OUT_MAT, OUT_PNG);

%% =========================================================================
function f = crossing_local(f0, g, ref)
    % First tempo at which gain falls below ref (log-linear interpolation); NaN if never.
    i = find(g < ref, 1);
    if isempty(i) || i == 1, f = NaN; return, end
    t = (ref - g(i-1)) / (g(i) - g(i-1));
    f = exp(log(f0(i-1)) + t * (log(f0(i)) - log(f0(i-1))));
end

function m = isequaln_local(v, s)
    m = (isnan(v) & isnan(s)) | v == s;
end

function v = ifelse_local(c, a, b)
    if c, v = a; else, v = b; end
end
