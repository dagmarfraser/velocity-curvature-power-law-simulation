%% checkFilterCutoffCollapse_v002.m
% v001 (Finding #233) with the Savitzky-Golay rows re-run EDGE-SAFE, for Fig 1B-2.
% Finding #249: v001's fixed 20-sample clip keeps samples whose SG window reaches the zero
% padding of conv(...,'same'); at the long window (81 samples) that produced a spurious
% SG-OLS "exception". Here the SG clip is half-window + 1 (11, 21, 41 samples); Butterworth
% rows are unchanged and reused from v001 after an anchor re-run reproduces them exactly.
% Anchors (fail loud): (1) v001's config equals this one's; (2) a re-run of the f_c = 10 Hz
% Butterworth rows at three tempos reproduces v001 to 1e-12 (sigma = 0, deterministic);
% (3) edge-safe SG-OLS onsets reproduce #249's node-grid values 0.311, 0.309, 0.308 (+-0.003).
% Run from src/. Reads results/filterCutoffCollapse_v001.mat. Writes
% results/filterCutoffCollapse_v002.mat, figures/filterCutoffCollapse_v002.png, and a dated
% log in logs/ (written with fprintf, not diary, so it also works when run through the MCP).

%% CONFIG
ROOT     = fileparts(fileparts(mfilename("fullpath")));
IN_MAT   = fullfile(ROOT, "results", "filterCutoffCollapse_v001.mat");
OUT_MAT  = fullfile(ROOT, "results", "filterCutoffCollapse_v002.mat");
OUT_PNG  = fullfile(ROOT, "figures", "filterCutoffCollapse_v002.png");
LOG_DIR  = fullfile(ROOT, "logs");
FS       = 240;  A_MM = 50;  AB = 2.16;  NCYC = 10;  CLIP_BW = 20;
F0       = logspace(log10(0.25), log10(16), 19);
BETA     = 0:1/30:0.75;
BW_FC    = [5 10 20];
SG_TW    = [21 41 81] / FS;
GAIN_WIN = [0.2 0.5];  GAIN_REF = 0.9;
REG      = [3 5];  REG_NAMES = ["OLS" "IRLS"];
SG_OLS_249 = [0.311 0.309 0.308];  TOL_249 = 0.003;   % #249, edge-safe, #233's node grid

%% Setup and anchor 1
addpath(genpath(fullfile(ROOT, "src", "functions")));  addpath(genpath(fullfile(ROOT, "src", "req")));
for f = [OUT_MAT OUT_PNG], if ~isfolder(fileparts(f)), error("filterCollapse2:outDir", "%s", "Missing folder: " + fileparts(f)); end, end
if ~isfolder(LOG_DIR), mkdir(LOG_DIR); end
logFile = fullfile(LOG_DIR, "checkFilterCutoffCollapse_v002_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log");
LG = fopen(logFile, "w");  if LG < 0, error("filterCollapse2:log", "%s", "Cannot open log: " + logFile); end
say_local(LG, "checkFilterCutoffCollapse_v002 run %s\n", string(datetime("now")));
if ~isfile(IN_MAT), error("filterCollapse2:in", "%s", "Missing input: " + IN_MAT); end
V1 = load(IN_MAT, "J", "flt", "F0", "BETA", "FS", "GAIN_WIN", "GAIN_REF");
if V1.FS ~= FS || ~isequal(V1.F0, F0) || ~isequal(V1.BETA, BETA) || ~isequal(V1.GAIN_WIN, GAIN_WIN) || V1.GAIN_REF ~= GAIN_REF
    error("filterCollapse2:anchor1", "%s", "v001 configuration differs from this script's");
end
flt = V1.flt;
if ~isequaln(flt.fc(flt.family == "BW")', BW_FC) || ~isequal(flt.Tw(flt.family == "SG")', SG_TW)
    error("filterCollapse2:anchor1", "%s", "v001 filter table differs");
end
say_local(LG, "Anchor 1 passed: v001 configuration matches.\n");

%% Anchor 2: re-run three Butterworth tempos at f_c = 10 Hz, compare with v001
fBW10 = find(flt.family == "BW" & flt.fc == 10);
kA = find(V1.J.filt == fBW10 & ismember(V1.J.f0, F0([5 10 15])));
bA = runJobs_local(V1.J(kA, :), flt, CLIP_BW, FS, A_MM, AB, NCYC, REG);
d2 = max(abs([bA(:, 1) - V1.J.betaOLS(kA); bA(:, 2) - V1.J.betaIRLS(kA)]), [], "omitnan");
if ~(d2 <= 1e-12), error("filterCollapse2:anchor2", "%s", sprintf("BW re-run differs from v001 by %.3g", d2)); end
say_local(LG, "Anchor 2 passed: %d Butterworth jobs reproduce v001 (max |diff| %.2g).\n", numel(kA), d2);

%% Re-run the Savitzky-Golay rows edge-safe
J = V1.J;
for f = find(flt.family == "SG")'
    fr = round(flt.Tw(f) * FS);  fr = fr + mod(fr + 1, 2);  clipSG = (fr - 1) / 2 + 1;
    k = find(J.filt == f);
    b = runJobs_local(J(k, :), flt, clipSG, FS, A_MM, AB, NCYC, REG);
    J.betaOLS(k) = b(:, 1);  J.betaIRLS(k) = b(:, 2);  J.nFail(k) = sum(~isfinite(b), 2);
    say_local(LG, "SG window %d samples: clip %d, %d jobs, %d non-finite fits\n", fr, clipSG, numel(k), sum(~isfinite(b(:))));
end

%% Gain against tempo and the collapse test (as v001)
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
        G = [G; table(flt.family(f), flt.fc(f), flt.Tw(f), F0(k), g(1), g(2), ...
            'VariableNames', ["family" "fc" "Tw" "f0" "gainOLS" "gainIRLS"])]; %#ok<AGROW>
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
say_local(LG, "\nTempo at which the local gain falls to %.2f, and its scaled value (f0/fc for BW, f0*Tw for SG):\n", GAIN_REF);
say_local(LG, "%s", evalc("disp(C)"));
S = groupsummary(C, ["family" "regression"], ["mean" "std"], "scaled");  S.cv = S.std_scaled ./ S.mean_scaled;
say_local(LG, "Collapse test (CV of the scaled value across the three settings):\n");
say_local(LG, "%s", evalc("disp(S(:, [""family"" ""regression"" ""mean_scaled"" ""std_scaled"" ""cv""]))"));

%% Anchor 3: edge-safe SG-OLS onsets against #249
sgOLS = C.scaled(C.family == "SG" & C.regression == "OLS")';
say_local(LG, "Anchor 3: SG-OLS onsets %s against #249's %s\n", mat2str(sgOLS, 4), mat2str(SG_OLS_249, 3));
if any(~isfinite(sgOLS)) || any(abs(sgOLS - SG_OLS_249) > TOL_249)
    error("filterCollapse2:anchor3", "%s", "Edge-safe SG-OLS onsets do not reproduce Finding #249");
end
say_local(LG, "Anchor 3 passed.\n");

%% Figure (as v001; bottom row is the collapse)
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
        if row == 1, xlabel("f_0 (Hz)"); elseif fmy == "BW", xlabel("f_0 / f_c"); else, xlabel("f_0 \times T_w"); end
        ylabel("gain d\beta_{rec}/d\beta_{gen} (OLS, near 1/3)");
        ttl = "Butterworth";  if fmy == "SG", ttl = "Savitzky-Golay (edge-safe clip)"; end
        if row == 1, ttl = ttl + ", raw tempo"; else, ttl = ttl + ", tempo scaled by the filter"; end
        title(ttl);  legend("Location", "southwest");
    end
end
title(tl, "The pull to one-third at zero noise: gain collapses when tempo is scaled by the filter");
set(findall(fg, "Type", "axes"), "Toolbar", []);  exportgraphics(fg, OUT_PNG, "Resolution", 200);
save(OUT_MAT, "J", "G", "C", "S", "flt", "F0", "BETA", "FS", "GAIN_WIN", "GAIN_REF", "CLIP_BW", "-v7.3");
say_local(LG, "Saved: %s\nFigure: %s\nLog: %s\n", OUT_MAT, OUT_PNG, logFile);
fclose(LG);

%% =========================================================================
function b = runJobs_local(Jk, flt, clip, FS, A_MM, AB, NCYC, REG)
    % Zero-noise ellipse jobs through the production pipeline (as v001); errors are counted, not absorbed.
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

function say_local(fid, fmt, varargin)
    % Print to the console and to the dated log.
    fprintf(fmt, varargin{:});  fprintf(fid, fmt, varargin{:});
end

function f = crossing_local(f0, g, ref)
    i = find(g < ref, 1);
    if isempty(i) || i == 1, f = NaN; return, end
    t = (ref - g(i-1)) / (g(i) - g(i-1));
    f = exp(log(f0(i-1)) + t * (log(f0(i)) - log(f0(i-1))));
end

function m = isequaln_local(v, s)
    m = (isnan(v) & isnan(s)) | v == s;
end
