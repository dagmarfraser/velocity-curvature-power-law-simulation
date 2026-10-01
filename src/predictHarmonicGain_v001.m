%% predictHarmonicGain_v001.m
% First-order harmonic theory of the zero-noise pull towards 1/3 (Finding #233), tested
% against #233's measured gain curves.
%
% Theory. On an ellipse, beta = 1/3 is harmonic motion: x and y are pure sinusoids at f0.
% At beta = 1/3 + eps the traversal phase is modulated at even multiples of f0, so to first
% order the trajectory gains eps-proportional energy at odd harmonics k*f0 (k = 1, 3, 5, ...).
% A linear analysis chain scales each harmonic by its response relative to an ideal
% differentiator, R_c(f) = |H_c(f) / D_c(f)|, for channel c (velocity, acceleration). A
% uniform amplitude change per channel leaves beta unchanged (scale invariance), so the
% forward map's gain at 1/3 is
%     g(f0) = sum_c sum_k w_ck * R_c(k f0) / R_c(f0)
% where the weights w_ck depend on geometry and estimator only, not on tempo or filter.
% The theory is magnitude-only: the half-sample delay of BWFD's forward-difference velocity
% is ignored by it and captured by the emulation below.
%
% Three levels, zero noise, exact ellipse, #233's settings (fs 240 Hz, 10 cycles, clip 20):
%   measured  : #233's pipeline gains (results/filterCutoffCollapse_v001.mat), the anchor
%   emulation : the pipeline's exact complex frequency response (filter and differences)
%               applied to a periodic record; tests "the effect is the linear filter"
%   theory    : the first-order sum above, and a 3f0-only variant (all non-fundamental
%               weight placed at 3f0), the back-of-envelope form
% Eccentricity test: weights shift with a/b, so the theory predicts geometry-dependent onsets;
% the real pipeline (#233's code path) is run at a/b 1.25, 2.16 and 4.1 for the production
% filters. At 2.16 it must reproduce #233's gains exactly (anchor).
% Gain as #233: slope of beta_rec on beta_gen over GAIN_WIN (local, around 1/3); the theory
% is the local derivative at 1/3.
% Linearity gate (Block A): OLS must pass homogeneity, channel and harmonic additivity (error
% otherwise). IRLS's bisquare scale comes from O(eps) residuals at zero noise, so it may be
% homogeneous but non-additive; if so its theory values are flagged APPROX in every table and
% a DEGRADED banner is printed (first run, 2026-10-01: a/b 4.1 IRLS sum off by 0.027).
% Run from src/. Writes results/predictHarmonicGain_v001.mat, figures/predictHarmonicGain_v001.png

%% CONFIG
ROOT    = fileparts(fileparts(mfilename("fullpath")));
IN_233  = fullfile(ROOT, "results", "filterCutoffCollapse_v001.mat");
OUT_MAT = fullfile(ROOT, "results", "predictHarmonicGain_v001.mat");
OUT_PNG = fullfile(ROOT, "figures", "predictHarmonicGain_v001.png");
A_MM    = 50;  AB_REF = 2.16;  AB_SET = [1.25 2.16 4.1];   % #233 geometry; eccentricity test
NCYC    = 10;  CLIP = 20;                                   % as #233
EPS_B   = 0.01;                                             % beta step for the weights
KMAX    = 21;                                               % highest odd harmonic kept
NPHI    = 1e5;                                              % generator grid for the weights
W_TOL   = 0.02;                                             % |sum(w) - total sensitivity|
W_MIN   = 1e-3;                                             % weights above this must not alias
F_FINE  = logspace(log10(0.1), log10(20), 3000);            % theory grid (Hz)
PROD    = [2 5];                                            % #233 flt rows: BW 10 Hz, SG 41 samples
REG     = [3 5];  REG_NAMES = ["OLS" "IRLS"];
K       = 1:2:KMAX;

%% Setup
addpath(genpath(fullfile(ROOT, "src", "functions")));  addpath(genpath(fullfile(ROOT, "src", "req")));
if ~isfile(IN_233), error("harmGain:in233", "%s", "Missing #233 results: " + IN_233); end
for f = [OUT_MAT OUT_PNG], if ~isfolder(fileparts(f)), error("harmGain:outDir", "%s", "Missing folder: " + fileparts(f)); end, end
R233 = load(IN_233, "G", "C", "flt", "F0", "BETA", "FS", "GAIN_WIN", "GAIN_REF");
FS = R233.FS;  F0 = R233.F0;  BETA = R233.BETA;  GAIN_WIN = R233.GAIN_WIN;  GAIN_REF = R233.GAIN_REF;
flt = R233.flt;  G233 = R233.G;  C233 = R233.C;
iRef = find(AB_SET == AB_REF);
if FS ~= 240 || height(flt) ~= 6 || flt.fc(PROD(1)) ~= 10 || abs(flt.Tw(PROD(2)) * FS - 41) > 1e-9 || isempty(iRef)
    error("harmGain:config", "%s", "#233 settings differ from what this script assumes");
end

%% Anchor 1: #233's onsets recomputed from its own gain table with its own crossing rule
for r = 1:height(C233)
    gm = gainRows_local(G233, C233.family(r), C233.fc(r), C233.Tw(r));
    f90 = crossing_local(gm.f0, gm.("gain" + C233.regression(r)), GAIN_REF);
    if ~(isequaln(f90, C233.f0_at_gain90(r)) || abs(f90 - C233.f0_at_gain90(r)) <= 1e-12)
        error("harmGain:anchor233", "%s", sprintf("Row %d: %.6f vs #233 %.6f", r, f90, C233.f0_at_gain90(r)));
    end
end
bwS = C233.scaled(C233.family == "BW");  sgS = C233.scaled(C233.family == "SG" & C233.regression == "IRLS");
fprintf("Anchor 1 passed: #233 onsets reproduced (BW f0/fc %.3f-%.3f; SG-IRLS f0*Tw %.3f-%.3f)\n", ...
    min(bwS), max(bwS), min(sgS), max(sgS));

%% Anchor 2: modelled chain responses equal an ideal differentiator at low frequency
for r = 1:height(flt)
    [Hv, Ha] = chainResponse_local(flt.family(r), flt.fc(r), flt.Tw(r), FS, 0.2);
    D = 1i * 2 * pi * 0.2;
    if abs(Hv / D - 1) > 1e-2 || abs(Ha / D^2 - 1) > 1e-2
        error("harmGain:response", "%s", sprintf("Filter %d: H/D = %.4f%+.4fi (v), %.4f%+.4fi (a) at 0.2 Hz", ...
            r, real(Hv / D), imag(Hv / D), real(Ha / D^2), imag(Ha / D^2)));
    end
end
fprintf("Anchor 2 passed: all six chain responses within 1%% of ideal derivatives at 0.2 Hz\n");

%% Block A: first-order harmonic weights per geometry (ideal spectral derivatives)
% Gate. OLS is linear in log space, so its first-order weights must add up to its total
% sensitivity: a hard requirement. IRLS (bisquare, fitnlm) takes its robust scale from the
% residuals, which at zero noise are themselves O(eps), so its response may be homogeneous
% but not additive across harmonics. Diagnosed, not assumed: homogeneity (eps vs 2 eps),
% channel additivity (velocity-only + acceleration-only vs full) and harmonic additivity are
% printed for every geometry. IRLS failing additivity only -> visible DEGRADED flag on its
% theory values; anything else -> error.
W = cell(numel(AB_SET), 1);  Wtot = NaN(numel(AB_SET), numel(REG));  info = struct([]);  DG = table();
for i = 1:numel(AB_SET)
    [W{i}, Wtot(i, :), inf_i] = harmonicWeights_local(AB_SET(i), A_MM, FS, NCYC, KMAX, EPS_B, NPHI, REG);
    info = [info; inf_i]; %#ok<AGROW>
    sumW = reshape(sum(W{i}, [1 2]), 1, []);
    for q = 1:numel(REG)
        DG = [DG; table(AB_SET(i), REG_NAMES(q), Wtot(i, q), inf_i.tot2(q), abs(Wtot(i, q) - inf_i.tot2(q)), ...
            abs(inf_i.totV(q) + inf_i.totA(q) - Wtot(i, q)), abs(sumW(q) - Wtot(i, q)), inf_i.energyResidual, ...
            'VariableNames', ["ab" "regression" "totalSens" "totalSens2eps" "homogeneityDev" "channelAddDev" ...
            "harmonicAddDev" "energyNotOddK"])]; %#ok<AGROW>
    end
end
fprintf("\nBlock A diagnostics (tolerance %.3f):\n", W_TOL);  disp(DG)
isOLS = DG.regression == "OLS";  nonAdd = DG.harmonicAddDev > W_TOL;
if any(isOLS & (nonAdd | DG.channelAddDev > W_TOL | DG.homogeneityDev > W_TOL))
    error("harmGain:olsNonlinear", "%s", "OLS weights fail a linearity check: the first-order theory does not hold for OLS");
end
if any(~isOLS & DG.homogeneityDev > W_TOL)
    error("harmGain:irlsNonHomogeneous", "%s", "IRLS total sensitivity depends on eps: not a first-order response");
end
irlsApprox = false(numel(AB_SET), 1);
for i = 1:numel(AB_SET), irlsApprox(i) = any(~isOLS & nonAdd & DG.ab == AB_SET(i)); end
if any(irlsApprox)
    msg = sprintf("DEGRADED: IRLS weights do not add up at a/b %s (harmonic additivity dev %s). IRLS theory " + ...
        "values there are approximate and flagged in every table; OLS is the test of the theory.", ...
        mat2str(AB_SET(irlsApprox)), mat2str(DG.harmonicAddDev(~isOLS & nonAdd)', 3));
    warning("harmGain:irlsNonAdditive", "%s", msg);
    fprintf("\n*** %s ***\n", msg);
end
W2 = harmonicWeights_local(AB_REF, A_MM, FS, NCYC, KMAX, 2 * EPS_B, NPHI, REG);
fprintf("\nBlock A: harmonic weights (OLS | IRLS), velocity + acceleration channels summed\n");
TW = table();
for i = 1:numel(AB_SET)
    for q = 1:numel(REG)
        wk = sum(W{i}(:, :, q), 1);
        TW = [TW; table(AB_SET(i), REG_NAMES(q), Wtot(i, q), sum(wk), wk(1), wk(2), wk(3), sum(wk(4:end)), ...
            sum(W{i}(1, :, q)), sum(W{i}(2, :, q)), info(i).energyResidual, ...
            'VariableNames', ["ab" "regression" "totalSens" "sumW" "w1" "w3" "w5" "w7up" "wVel" "wAcc" "energyNotOddK"])]; %#ok<AGROW>
    end
end
disp(TW)
fprintf("Linearity (a/b %.2f): max |w(eps) - w(2 eps)| = %.4f\n", AB_REF, max(abs(W{iRef} - W2), [], "all"));

%% Jobs: emulation (six filters, a/b 2.16) and pipeline (production filters, three a/b)
[fi, ki, bi] = ndgrid(1:height(flt), 1:numel(F0), 1:numel(BETA));
JE = table(repmat("emul", numel(fi), 1), fi(:), repmat(AB_REF, numel(fi), 1), F0(ki(:))', BETA(bi(:))', ...
    'VariableNames', ["block" "filt" "ab" "f0" "betaGen"]);
[pf, ai, ki, bi] = ndgrid(PROD, AB_SET, 1:numel(F0), 1:numel(BETA));
JP = table(repmat("pipe", numel(pf), 1), pf(:), ai(:), F0(ki(:))', BETA(bi(:))', ...
    'VariableNames', ["block" "filt" "ab" "f0" "betaGen"]);
J = [JE; JP];
fprintf("\nJobs: %d emulation + %d pipeline trajectories, sigma = 0\n", height(JE), height(JP));
bRec = NaN(height(J), numel(REG));  errFlag = false(height(J), 1);
fam = flt.family;  fcv = flt.fc;  twv = flt.Tw;
blk = J.block;  fl = J.filt;  abv = J.ab;  f0v = J.f0;  bgv = J.betaGen;
nW = 0;  if ~isempty(gcp("nocreate")), nW = gcp("nocreate").NumWorkers; end
parfor (k = 1:height(J), nW)
    r = fl(k);  b = NaN(1, numel(REG));  ef = false;
    if blk(k) == "pipe"
        [x, y] = generatePowerLawEllipse_v001(A_MM, A_MM / abv(k), f0v(k), FS, bgv(k), 'nCycles', NCYC);
        try
            b = pipelineBeta_local(x, y, fam(r), fcv(r), twv(r), FS, CLIP, REG);
        catch
            ef = true;   % counted and reported below, as #233; never silently absorbed
        end
    else
        b = emulateBeta_local(A_MM, abv(k), f0v(k), bgv(k), fam(r), fcv(r), twv(r), FS, NCYC, REG);
    end
    bRec(k, :) = b;  errFlag(k) = ef;
end
J.bOLS = bRec(:, 1);  J.bIRLS = bRec(:, 2);  J.err = errFlag;
fprintf("Non-finite estimates: %d emulation, %d pipeline; pipeline trajectories that threw: %d\n", ...
    sum(~isfinite(bRec(J.block == "emul", :)), "all"), sum(~isfinite(bRec(J.block == "pipe", :)), "all"), sum(errFlag));

%% Gains
[grp, GG] = findgroups(J(:, ["block" "filt" "ab" "f0"]));
GG.gainOLS  = splitapply(@(b, y) gainAt_local(b, y, GAIN_WIN), J.betaGen, J.bOLS, grp);
GG.gainIRLS = splitapply(@(b, y) gainAt_local(b, y, GAIN_WIN), J.betaGen, J.bIRLS, grp);

%% Anchor 3: the pipeline at a/b 2.16 reproduces #233's gains
for r = PROD
    gp = GG(GG.block == "pipe" & GG.filt == r & GG.ab == AB_REF, :);
    gm = gainRows_local(G233, flt.family(r), flt.fc(r), flt.Tw(r));
    if height(gp) ~= height(gm), error("harmGain:anchorPipe", "%s", sprintf("Filter %d: %d vs %d tempi", r, height(gp), height(gm))); end
    d = [gp.gainOLS - gm.gainOLS; gp.gainIRLS - gm.gainIRLS];
    if ~isequal(isnan([gp.gainOLS; gp.gainIRLS]), isnan([gm.gainOLS; gm.gainIRLS])) || max(abs(d), [], "omitnan") > 1e-10
        error("harmGain:anchorPipe", "%s", sprintf("Filter %d: pipeline gains differ from #233 (max %.3g)", r, max(abs(d), [], "omitnan")));
    end
end
fprintf("Anchor 3 passed: pipeline gains at a/b %.2f reproduce #233 for both production filters\n", AB_REF);

%% Theory curves
TH = cell(height(flt), 1);  TH3 = TH;  THF0 = TH;
for r = 1:height(flt)
    [TH{r}, TH3{r}] = theoryGain_local(W{iRef}, K, flt.family(r), flt.fc(r), flt.Tw(r), FS, F_FINE, W_MIN);
    THF0{r} = theoryGain_local(W{iRef}, K, flt.family(r), flt.fc(r), flt.Tw(r), FS, F0, W_MIN);
end
THecc = cell(numel(PROD), numel(AB_SET));
for p = 1:numel(PROD)
    for i = 1:numel(AB_SET)
        r = PROD(p);  THecc{p, i} = theoryGain_local(W{i}, K, flt.family(r), flt.fc(r), flt.Tw(r), FS, F_FINE, W_MIN);
    end
end

%% Onsets: measured, emulated, theory, theory (3f0 only)
T_on = table();  T_rms = table();
for r = 1:height(flt)
    gm = gainRows_local(G233, flt.family(r), flt.fc(r), flt.Tw(r));
    ge = GG(GG.block == "emul" & GG.filt == r, :);
    sc = scale_local(flt(r, :));
    for q = 1:numel(REG)
        nm = "gain" + REG_NAMES(q);
        fM = C233.f0_at_gain90(C233.family == flt.family(r) & isequaln_local(C233.fc, flt.fc(r)) & ...
            isequaln_local(C233.Tw, flt.Tw(r)) & C233.regression == REG_NAMES(q));
        fE = crossing_local(ge.f0, ge.(nm), GAIN_REF);
        fT = crossing_local(F_FINE(:), TH{r}(:, q), GAIN_REF);
        f3 = crossing_local(F_FINE(:), TH3{r}(:, q), GAIN_REF);
        T_on = [T_on; table(flt.family(r), flt.fc(r), flt.Tw(r), REG_NAMES(q), fM * sc, fE * sc, fT * sc, f3 * sc, ...
            theoryFlag_local(REG_NAMES(q), irlsApprox(iRef)), ...
            'VariableNames', ["family" "fc" "Tw" "regression" "measured" "emulated" "theory" "theory3f0" "theoryStatus"])]; %#ok<AGROW>
        ok = isfinite(gm.(nm)) & gm.(nm) >= 0.5;
        T_rms = [T_rms; table(flt.family(r), flt.fc(r), flt.Tw(r), REG_NAMES(q), sum(ok), ...
            rms_local(ge.(nm)(ok) - gm.(nm)(ok)), rms_local(THF0{r}(ok, q) - gm.(nm)(ok)), ...
            'VariableNames', ["family" "fc" "Tw" "regression" "nTempi" "rmsEmulVsMeasured" "rmsTheoryVsMeasured"])]; %#ok<AGROW>
    end
end
fprintf("\nOnset of gain %.2f, scaled (f0/fc for BW, f0*Tw for SG). Measured = #233.\n", GAIN_REF);
disp(T_on)
bw = T_on(T_on.family == "BW", :);
fprintf("BW cutoff multiple implied (1 / scaled onset): measured %.2f-%.2f, theory %.2f-%.2f, 3f0-only %.2f-%.2f\n", ...
    1 / max(bw.measured), 1 / min(bw.measured), 1 / max(bw.theory), 1 / min(bw.theory), 1 / max(bw.theory3f0), 1 / min(bw.theory3f0));
fprintf("\nRMS gain error over tempi where measured gain >= 0.5:\n");
disp(T_rms)

%% Eccentricity: pipeline onsets against theory at each a/b
T_ecc = table();
for p = 1:numel(PROD)
    r = PROD(p);
    for i = 1:numel(AB_SET)
        gp = GG(GG.block == "pipe" & GG.filt == r & GG.ab == AB_SET(i), :);
        for q = 1:numel(REG)
            T_ecc = [T_ecc; table(flt.family(r), AB_SET(i), REG_NAMES(q), crossing_local(gp.f0, gp.("gain" + REG_NAMES(q)), GAIN_REF), ...
                crossing_local(F_FINE(:), THecc{p, i}(:, q), GAIN_REF), theoryFlag_local(REG_NAMES(q), irlsApprox(i)), ...
                'VariableNames', ["family" "ab" "regression" "pipelineOnsetHz" "theoryOnsetHz" "theoryStatus"])]; %#ok<AGROW>
        end
    end
end
fprintf("Eccentricity test (production filters; onset of gain %.2f, Hz):\n", GAIN_REF);
disp(T_ecc)

%% Figure
fg = figure("Color", "w", "Position", [80 80 1150 800]);  tl = tiledlayout(2, 2, "TileSpacing", "compact");
cols = [0.52 0.72 0.92; 0.22 0.54 0.87; 0.05 0.27 0.49];
for fmy = ["BW" "SG"]
    nexttile; hold on;  rows = find(flt.family == fmy);
    for s = 1:numel(rows)
        r = rows(s);  sc = scale_local(flt(r, :));
        gm = gainRows_local(G233, fmy, flt.fc(r), flt.Tw(r));  ge = GG(GG.block == "emul" & GG.filt == r, :);
        lab = sprintf("f_c %g Hz", flt.fc(r));  if fmy == "SG", lab = sprintf("T_w %.3f s", flt.Tw(r)); end
        plot(gm.f0 * sc, gm.gainOLS, "o", "Color", cols(s, :), "MarkerSize", 5, "DisplayName", lab);
        plot(ge.f0 * sc, ge.gainOLS, "-", "Color", cols(s, :), "HandleVisibility", "off");
        plot(F_FINE * sc, TH{r}(:, 1), "--", "Color", cols(s, :), "LineWidth", 1.2, "HandleVisibility", "off");
    end
    yline(GAIN_REF, ":", "HandleVisibility", "off");  set(gca, "XScale", "log");  ylim([0 1.1]);  box on;
    if fmy == "BW", xlabel("f_0 / f_c"); title("Butterworth (filtfilt) + finite differences");
    else, xlabel("f_0 \times T_w"); title("Savitzky-Golay, order 4"); end
    ylabel("gain d\beta_{rec}/d\beta_{gen} (OLS)");  legend("Location", "southwest");
    text(0.98, 0.95, "o measured (#233)   - emulation   -- first-order theory", "Units", "normalized", ...
        "HorizontalAlignment", "right", "FontSize", 8);
end
nexttile;  wk = zeros(numel(K), numel(AB_SET));
for i = 1:numel(AB_SET), wk(:, i) = sum(W{i}(:, :, 1), 1)'; end
bar(K, wk);  xlabel("harmonic k (of f_0)");  ylabel("first-order weight w_k (OLS)");  box on;
legend(compose("a/b %.2f", AB_SET), "Location", "northeast");  title("Where the \beta information sits");
nexttile; hold on;  mk = ["o" "s"];
for p = 1:numel(PROD)
    r = PROD(p);
    for i = 1:numel(AB_SET)
        gp = GG(GG.block == "pipe" & GG.filt == r & GG.ab == AB_SET(i), :);
        plot(gp.f0, gp.gainOLS, mk(p), "Color", cols(i, :), "MarkerSize", 5, ...
            "DisplayName", sprintf("%s, a/b %.2f", flt.family(r), AB_SET(i)));
        plot(F_FINE, THecc{p, i}(:, 1), "--", "Color", cols(i, :), "HandleVisibility", "off");
    end
end
yline(GAIN_REF, ":", "HandleVisibility", "off");  set(gca, "XScale", "log");  ylim([0 1.1]);  box on;
xlabel("f_0 (Hz)");  ylabel("gain (OLS)");  title("Eccentricity: pipeline (markers) against theory (dashed)");
legend("Location", "southwest", "FontSize", 7);
title(tl, "Zero-noise pull towards 1/3: first-order harmonic theory against measured gain");
set(findall(fg, "Type", "axes"), "Toolbar", []);  exportgraphics(fg, OUT_PNG, "Resolution", 200);

save(OUT_MAT, "J", "GG", "W", "W2", "Wtot", "info", "TW", "DG", "irlsApprox", "T_on", "T_rms", "T_ecc", "TH", "TH3", "THecc", ...
    "F_FINE", "K", "AB_SET", "AB_REF", "EPS_B", "KMAX", "W_MIN", "-v7.3");
fprintf("Saved: %s\nFigure: %s\n", OUT_MAT, OUT_PNG);

%% =========================================================================
function [W, tot, info] = harmonicWeights_local(ab, A_MM, FS, NCYC, KMAX, epsB, NPHI, REG)
% First-order weights of beta_rec on each odd harmonic of d(x,y)/d(beta) at beta = 1/3, per
% channel (1 velocity, 2 acceleration), with ideal spectral derivatives on a periodic record.
M = NCYC * FS;                                   % f0 = 1 Hz: harmonic k sits at exactly k Hz
gen = @(bt) genPeriodic_local(A_MM, A_MM / ab, 1, FS, bt, M, NPHI);
[x0, y0] = gen(1/3);  [xp, yp] = gen(1/3 + epsB);  [xm, ym] = gen(1/3 - epsB);
f = signedBins_local(M, FS);  Dv = 1i * 2 * pi * f;  Da = Dv.^2;
nyq = abs(abs(f) - FS/2) < 1e-9;  Dv(nyq) = 0;  Da(nyq) = 0;
d = @(s, D) real(ifft(D .* fft(s)));
v0 = d([x0 y0], Dv);  a0 = d([x0 y0], Da);
[xp2, yp2] = gen(1/3 + 2 * epsB);  [xm2, ym2] = gen(1/3 - 2 * epsB);
tot  = (betaKin_local(d([xp yp], Dv), d([xp yp], Da), REG) - betaKin_local(d([xm ym], Dv), d([xm ym], Da), REG)) / (2 * epsB);
tot2 = (betaKin_local(d([xp2 yp2], Dv), d([xp2 yp2], Da), REG) - betaKin_local(d([xm2 ym2], Dv), d([xm2 ym2], Da), REG)) / (4 * epsB);
dP = ([xp yp] - [xm ym]) / (2 * epsB);  dvF = d(dP, Dv);  daF = d(dP, Da);   % d(x,y)/d(beta) and its derivatives
totV = (betaKin_local(v0 + epsB * dvF, a0, REG) - betaKin_local(v0 - epsB * dvF, a0, REG)) / (2 * epsB);
totA = (betaKin_local(v0, a0 + epsB * daF, REG) - betaKin_local(v0, a0 - epsB * daF, REG)) / (2 * epsB);
Fp = fft(dP);                                    % spectrum of d(x,y)/d(beta)
E = sum(abs(Fp).^2, "all");  binw = FS / M;
K = 1:2:KMAX;  W = NaN(2, numel(K), numel(REG));  eK = NaN(1, numel(K));
for j = 1:numel(K)
    m = abs(abs(f) - K(j)) < binw / 2;
    eK(j) = sum(abs(Fp(m, :)).^2, "all") / E;
    dk = real(ifft(m .* Fp));  dv = d(dk, Dv);  da = d(dk, Da);
    W(1, j, :) = (betaKin_local(v0 + epsB * dv, a0, REG) - betaKin_local(v0 - epsB * dv, a0, REG)) / (2 * epsB);
    W(2, j, :) = (betaKin_local(v0, a0 + epsB * da, REG) - betaKin_local(v0, a0 - epsB * da, REG)) / (2 * epsB);
end
if any(~isfinite(W), "all") || any(~isfinite(tot))
    error("harmGain:weightNaN", "%s", sprintf("a/b %.2f: a weight regression returned NaN", ab));
end
info = struct("ab", ab, "energyK", eK, "energyResidual", 1 - sum(eK), "tot2", tot2, "totV", totV, "totA", totA);
end

function [x, y] = genPeriodic_local(a, b, f0, FS, bt, M, NPHI)
[x, y] = generatePowerLawEllipse_v001(a, b, f0, FS, bt, 'M', M, 'nPhi', NPHI);  x = x(:);  y = y(:);
end

function f = signedBins_local(M, FS)
f = (0:M-1)' * FS / M;  f(f > FS/2) = f(f > FS/2) - FS;
end

function [Hv, Ha] = chainResponse_local(fam, fc, Tw, FS, f)
% Complex responses of velocity and acceleration as differentiateKinematicsEBR computes them:
% case 2 (2nd-order Butterworth via filtfilt, then forward differences) or case 4 (SG order 4).
w = 2 * pi * f(:) / FS;  z = @(c) exp(-1i * w * (0:numel(c)-1)) * c(:);
if fam == "BW"
    [b, a] = butter(2, fc / (FS/2));
    Hb = abs(z(b) ./ z(a)).^2;                   % filtfilt: zero phase, magnitude squared
    Hv = Hb .* FS .* (1 - exp(-1i * w));         % v(i) = fs (x(i) - x(i-1)): half-sample delay
    Ha = Hb .* FS^2 .* (2 * cos(w) - 2);         % a(i) = fs (v(i+1) - v(i)): centred
else
    L = sgFrame_local(Tw, FS);  [~, g] = sgolay(4, L);  E = exp(-1i * w * ((1:L) - (L+1)/2));
    Hv = E * (factorial(1) / (-1/FS)^1 * g(:, 2));   % conv(x, c, 'same') with the code's scaling
    Ha = E * (factorial(2) / (-1/FS)^2 * g(:, 3));
end
end

function L = sgFrame_local(Tw, FS)
L = round(Tw * FS);  L = L + mod(L + 1, 2);      % as #233
end

function [g, g3] = theoryGain_local(Wab, K, fam, fc, Tw, FS, f0, wMin)
% First-order gain from the weights and the chain's magnitude response relative to ideal.
f0 = f0(:);  nq = size(Wab, 3);  g = NaN(numel(f0), nq);  g3 = g;
kUse = K(max(abs(Wab), [], [1 3]) > wMin);       % harmonics carrying weight
for i = 1:numel(f0)
    if max(kUse) * f0(i) >= FS/2, continue, end  % would alias: theory not defined, left NaN
    fq = [f0(i); K(:) * f0(i)];  [Hv, Ha] = chainResponse_local(fam, fc, Tw, FS, fq);
    D = 1i * 2 * pi * fq;  Rv = abs(Hv ./ D);  Ra = abs(Ha ./ D.^2);
    rv = Rv(2:end) / Rv(1);  ra = Ra(2:end) / Ra(1);
    for q = 1:nq
        wv = Wab(1, :, q)';  wa = Wab(2, :, q)';
        g(i, q)  = sum(wv .* rv) + sum(wa .* ra);
        g3(i, q) = wv(1) + wa(1) + sum(wv(2:end)) * rv(K == 3) + sum(wa(2:end)) * ra(K == 3);
    end
end
end

function b = emulateBeta_local(A_MM, ab, f0, bt, fam, fc, Tw, FS, NCYC, REG)
% The pipeline's exact linear response applied to a whole-cycle (exactly periodic) record.
M = round(NCYC * FS / f0);  f0e = NCYC * FS / M;
[x, y] = generatePowerLawEllipse_v001(A_MM, A_MM / ab, f0e, FS, bt, 'M', M);
f = signedBins_local(M, FS);  [Hv, Ha] = chainResponse_local(fam, fc, Tw, FS, f);
nyq = abs(abs(f) - FS/2) < 1e-9;  Hv(nyq) = 0;  Ha(nyq) = 0;
X = fft([x(:) y(:)]);
b = betaKin_local(real(ifft(Hv .* X)), real(ifft(Ha .* X)), REG);
end

function b = pipelineBeta_local(x, y, fam, fc, tw, FS, CLIP, REG)
% #233's code path, verbatim in effect.
if fam == "BW", ft = 2; fp = [2 fc 1];
else, fr = sgFrame_local(tw, FS); ft = 4; fp = [4 fr]; end
[dx, dy] = differentiateKinematicsEBR(x(:), y(:), ft, fp, FS);
b = betaKin_local([dx(CLIP:end-CLIP, 2) dy(CLIP:end-CLIP, 2)], [dx(CLIP:end-CLIP, 3) dy(CLIP:end-CLIP, 3)], REG);
end

function b = betaKin_local(v, a, REG)
% Speed and curvature, then #233's seeded regression sequence (limitBreak = 0).
sp = hypot(v(:, 1), v(:, 2));  kp = curvatureKinematicEBR(v(:, 1), v(:, 2), a(:, 1), a(:, 2));
ok = isfinite(sp) & isfinite(kp) & kp > 0;
b = NaN(1, numel(REG));  lm = [1, -1/3];
for q = 1:numel(REG)
    [bb, vv] = regressDataEBR(sp(ok), kp(ok), REG(q), lm, 0, 0);
    b(q) = bb;
    if q == 1 && isfinite(bb) && isfinite(vv), lm = [vv, bb]; end
end
end

function g = gainAt_local(betaGen, bRec, win)
w = betaGen >= win(1) & betaGen <= win(2) & isfinite(bRec);
g = NaN;  if sum(w) >= 3, p = polyfit(betaGen(w), bRec(w), 1); g = p(1); end
end

function gm = gainRows_local(G, fam, fc, Tw)
gm = sortrows(G(G.family == fam & isequaln_local(G.fc, fc) & isequaln_local(G.Tw, Tw), :), "f0");
end

function s = scale_local(row)
if row.family == "BW", s = 1 / row.fc; else, s = row.Tw; end
end

function s = theoryFlag_local(reg, approx)
s = "exact first order";  if reg == "IRLS" && approx, s = "APPROX (IRLS non-additive)"; end
end

function r = rms_local(e)
e = e(isfinite(e));  r = NaN;  if ~isempty(e), r = sqrt(mean(e.^2)); end
end

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
