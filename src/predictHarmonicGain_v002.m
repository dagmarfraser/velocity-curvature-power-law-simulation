%% predictHarmonicGain_v002.m
% Follow-up to predictHarmonicGain_v001 (first-order harmonic theory of the zero-noise pull
% towards 1/3, Finding #233). v001 confirmed the theory for OLS but left three residuals; v002
% tests one explanation for each. Predictions and thresholds fixed here before the run
% (2026-10-01); reported as HELD / NOT HELD, never gated:
%  P1 phase   : Butterworth onsets drift with f_c (emulated f0/f_c 0.1805 -> 0.1756, OLS) while
%               the magnitude-only theory is flat (0.1844-0.1846). Proposed cause: the
%               half-sample delay of the forward-difference velocity. Test: emulation with that
%               delay removed, magnitude unchanged ("emulZP"). Prediction: BW-OLS emulZP scaled
%               onsets spread <= P1_SPREAD across the three cutoffs (emulated spread 0.0049).
%  P2 window  : theory sits above emulation even where phase cannot act (SG derivatives are
%               zero phase: SG-OLS 0.3103 emulated against 0.3159). Proposed cause: theory is
%               the local derivative at 1/3, the gains are #233's secant over beta_gen 0.2-0.5.
%               Test: gains also over LOCAL_WIN (nodes 0.300, 1/3, 0.367, a central secant).
%               Prediction: OLS local-window onsets (emulZP for BW, emul for SG) within P2_ABS
%               (scaled) of theory for all six filters.
%  P3 edges   : #233's long SG window (OLS onset 0.219 against 0.308 emulated) and the production
%               SG at a/b 1.25 (OLS 1.57 Hz against theory 1.92) are record-end artefacts: conv
%               'same' zero-pads, and CLIP = 20 keeps samples whose window reaches the padding
%               (half-windows 20 and 40). Test: the pipeline with CLIP = half-window + 1
%               ("pipeSafe"). Prediction: pipeSafe matches emulation (RMS <= P3_RMS, OLS and IRLS,
%               tempi with emulated gain >= 0.5) in all five SG conditions, and the clip-20
%               pipeline does not, for the long window and for a/b 1.25 OLS.
%  P4 geometry: production filters, OLS, a/b 1.25, 2.16, 4.1: local-window onsets of the
%               corrected emulation (emulZP for BW, emul for SG) within P4_REL of theory computed
%               from each geometry's own weights.
% Anchors (error if broken): 1-3 as v001; 4, Block A weights and every gain v001 computed are
% reproduced to ANCHOR_TOL (results/predictHarmonicGain_v001.mat).
% Zero noise, exact ellipse, #233 settings (fs 240 Hz, 10 cycles). About 13,000 trajectories
% (OLS and IRLS fits): roughly twice v001. Run from src/.
% Writes results/predictHarmonicGain_v002.mat, figures/predictHarmonicGain_v002.png

%% CONFIG
ROOT       = fileparts(fileparts(mfilename("fullpath")));
IN_233     = fullfile(ROOT, "results", "filterCutoffCollapse_v001.mat");
IN_V001    = fullfile(ROOT, "results", "predictHarmonicGain_v001.mat");
OUT_MAT    = fullfile(ROOT, "results", "predictHarmonicGain_v002.mat");
OUT_PNG    = fullfile(ROOT, "figures", "predictHarmonicGain_v002.png");
A_MM       = 50;  AB_REF = 2.16;  AB_SET = [1.25 2.16 4.1];
NCYC       = 10;  CLIP = 20;                                % as #233
EPS_B      = 0.01;  KMAX = 21;  NPHI = 1e5;                 % as v001
W_TOL      = 0.02;  W_MIN = 1e-3;                           % as v001
F_FINE     = logspace(log10(0.1), log10(20), 3000);
PROD       = [2 5];                                         % #233 flt rows: BW 10 Hz, SG 41 samples
LOCAL_WIN  = [0.29 0.37];                                   % nodes 0.300, 0.333, 0.367
REG        = [3 5];  REG_NAMES = ["OLS" "IRLS"];
K          = 1:2:KMAX;
ANCHOR_TOL = 1e-10;                                         % as v001's anchor 3 (passed across machines)
P1_SPREAD  = 0.001;  P2_ABS = 0.002;  P3_RMS = 0.005;  P4_REL = 0.01;

%% Setup
addpath(genpath(fullfile(ROOT, "src", "functions")));  addpath(genpath(fullfile(ROOT, "src", "req")));
for f = [IN_233 IN_V001], if ~isfile(f), error("harmGain2:in", "%s", "Missing input: " + f); end, end
for f = [OUT_MAT OUT_PNG], if ~isfolder(fileparts(f)), error("harmGain2:outDir", "%s", "Missing folder: " + fileparts(f)); end, end
R233 = load(IN_233, "G", "C", "flt", "F0", "BETA", "FS", "GAIN_WIN", "GAIN_REF");
FS = R233.FS;  F0 = R233.F0;  BETA = R233.BETA;  GAIN_WIN = R233.GAIN_WIN;  GAIN_REF = R233.GAIN_REF;
flt = R233.flt;  G233 = R233.G;  C233 = R233.C;
V1 = load(IN_V001, "W", "GG");
iRef = find(AB_SET == AB_REF);  abOther = AB_SET(AB_SET ~= AB_REF);
iBW = find(flt.family == "BW")';  iSG = find(flt.family == "SG")';
[~, j] = max(flt.Tw(iSG));  rLong = iSG(j);
if FS ~= 240 || height(flt) ~= 6 || flt.fc(PROD(1)) ~= 10 || abs(flt.Tw(PROD(2)) * FS - 41) > 1e-9 || isempty(iRef)
    error("harmGain2:config", "%s", "#233 settings differ from what this script assumes");
end
if nnz(abs(BETA - 1/3) < 1e-12) ~= 1 || nnz(BETA >= LOCAL_WIN(1) & BETA <= LOCAL_WIN(2)) ~= 3
    error("harmGain2:localWin", "%s", "LOCAL_WIN must hold exactly three BETA nodes centred on 1/3");
end

%% Anchor 1: #233's onsets recomputed from its own gain table with its own crossing rule
for r = 1:height(C233)
    gm = gainRows_local(G233, C233.family(r), C233.fc(r), C233.Tw(r));
    f90 = crossing_local(gm.f0, gm.("gain" + C233.regression(r)), GAIN_REF);
    if ~(isequaln(f90, C233.f0_at_gain90(r)) || abs(f90 - C233.f0_at_gain90(r)) <= 1e-12)
        error("harmGain2:anchor233", "%s", sprintf("Row %d: %.6f vs #233 %.6f", r, f90, C233.f0_at_gain90(r)));
    end
end
fprintf("Anchor 1 passed: #233 onsets reproduced\n");

%% Anchor 2: modelled chain responses equal an ideal differentiator at low frequency
for r = 1:height(flt)
    [Hv, Ha] = chainResponse_local(flt.family(r), flt.fc(r), flt.Tw(r), FS, 0.2);
    D = 1i * 2 * pi * 0.2;
    if abs(Hv / D - 1) > 1e-2 || abs(Ha / D^2 - 1) > 1e-2
        error("harmGain2:response", "%s", sprintf("Filter %d: chain response off ideal at 0.2 Hz", r));
    end
end
fprintf("Anchor 2 passed: chain responses within 1%% of ideal derivatives at 0.2 Hz\n");

%% Block A: harmonic weights and linearity gate (as v001)
W = cell(numel(AB_SET), 1);  Wtot = NaN(numel(AB_SET), numel(REG));  info = struct([]);  DG = table();
for i = 1:numel(AB_SET)
    [W{i}, Wtot(i, :), inf_i] = harmonicWeights_local(AB_SET(i), A_MM, FS, NCYC, KMAX, EPS_B, NPHI, REG);
    info = [info; inf_i]; %#ok<AGROW>
    sumW = reshape(sum(W{i}, [1 2]), 1, []);
    for q = 1:numel(REG)
        DG = [DG; table(AB_SET(i), REG_NAMES(q), Wtot(i, q), abs(Wtot(i, q) - inf_i.tot2(q)), ...
            abs(inf_i.totV(q) + inf_i.totA(q) - Wtot(i, q)), abs(sumW(q) - Wtot(i, q)), ...
            'VariableNames', ["ab" "regression" "totalSens" "homogeneityDev" "channelAddDev" "harmonicAddDev"])]; %#ok<AGROW>
    end
end
disp(DG)
isOLS = DG.regression == "OLS";  nonAdd = DG.harmonicAddDev > W_TOL;
if any(isOLS & (nonAdd | DG.channelAddDev > W_TOL | DG.homogeneityDev > W_TOL))
    error("harmGain2:olsNonlinear", "%s", "OLS weights fail a linearity check");
end
if any(~isOLS & DG.homogeneityDev > W_TOL)
    error("harmGain2:irlsNonHomogeneous", "%s", "IRLS total sensitivity depends on eps");
end
irlsApprox = false(numel(AB_SET), 1);
for i = 1:numel(AB_SET), irlsApprox(i) = any(~isOLS & nonAdd & DG.ab == AB_SET(i)); end
if any(irlsApprox)
    msg = sprintf("DEGRADED: IRLS weights do not add up at a/b %s; IRLS theory values there are flagged APPROX.", ...
        mat2str(AB_SET(irlsApprox)));
    warning("harmGain2:irlsNonAdditive", "%s", msg);  fprintf("\n*** %s ***\n", msg);
end

%% Anchor 4a: Block A weights reproduce v001
for i = 1:numel(AB_SET)
    dW = max(abs(W{i} - V1.W{i}), [], "all");
    if ~(dW <= ANCHOR_TOL), error("harmGain2:anchorW", "%s", sprintf("a/b %.2f: weights differ from v001 by %.3g", AB_SET(i), dW)); end
end
fprintf("Anchor 4a passed: Block A weights reproduce v001\n");

%% Jobs: v001's blocks plus the three tests
J = [jobs_local("emul", 1:height(flt), AB_REF, F0, BETA); jobs_local("emul", PROD, abOther, F0, BETA); ...
     jobs_local("emulZP", iBW, AB_REF, F0, BETA); jobs_local("emulZP", PROD(1), abOther, F0, BETA); ...
     jobs_local("pipe", PROD, AB_SET, F0, BETA); ...
     jobs_local("pipeSafe", iSG, AB_REF, F0, BETA); jobs_local("pipeSafe", PROD(2), abOther, F0, BETA)];
if any((J.block == "pipeSafe" & flt.family(J.filt) ~= "SG") | (J.block == "emulZP" & flt.family(J.filt) ~= "BW"))
    error("harmGain2:jobs", "%s", "pipeSafe is defined for SG only and emulZP for BW only");
end
[gb, gBlk] = findgroups(J.block);
fprintf("\nJobs: %d trajectories (%s), sigma = 0\n", height(J), countList_local(gBlk, splitapply(@numel, J.filt, gb)));
bRec = NaN(height(J), numel(REG));  errFlag = false(height(J), 1);
fam = flt.family;  fcv = flt.fc;  twv = flt.Tw;
blk = J.block;  fl = J.filt;  abv = J.ab;  f0v = J.f0;  bgv = J.betaGen;
nW = 0;  if ~isempty(gcp("nocreate")), nW = gcp("nocreate").NumWorkers; end
parfor (k = 1:height(J), nW)
    r = fl(k);  b = NaN(1, numel(REG));  ef = false;
    if blk(k) == "pipe" || blk(k) == "pipeSafe"
        [x, y] = generatePowerLawEllipse_v001(A_MM, A_MM / abv(k), f0v(k), FS, bgv(k), 'nCycles', NCYC);
        cl = CLIP;  if blk(k) == "pipeSafe", cl = (sgFrame_local(twv(r), FS) - 1) / 2 + 1; end
        try
            b = pipelineBeta_local(x, y, fam(r), fcv(r), twv(r), FS, cl, REG);
        catch
            ef = true;   % counted and reported below, as #233; never silently absorbed
        end
    else
        b = emulateBeta_local(A_MM, abv(k), f0v(k), bgv(k), fam(r), fcv(r), twv(r), FS, NCYC, REG, blk(k) == "emulZP");
    end
    bRec(k, :) = b;  errFlag(k) = ef;
end
J.bOLS = bRec(:, 1);  J.bIRLS = bRec(:, 2);  J.err = errFlag;
fprintf("Non-finite estimates by block: %s; pipeline trajectories that threw: %d\n", ...
    countList_local(gBlk, splitapply(@(b) nnz(~isfinite(b)), [J.bOLS J.bIRLS], gb)), sum(errFlag));

%% Gains: #233's window and the local window
[grp, GG] = findgroups(J(:, ["block" "filt" "ab" "f0"]));
GG.gainOLS     = splitapply(@(b, y) gainAt_local(b, y, GAIN_WIN),  J.betaGen, J.bOLS,  grp);
GG.gainIRLS    = splitapply(@(b, y) gainAt_local(b, y, GAIN_WIN),  J.betaGen, J.bIRLS, grp);
GG.gainOLSloc  = splitapply(@(b, y) gainAt_local(b, y, LOCAL_WIN), J.betaGen, J.bOLS,  grp);
GG.gainIRLSloc = splitapply(@(b, y) gainAt_local(b, y, LOCAL_WIN), J.betaGen, J.bIRLS, grp);

%% Anchor 3: pipeline at a/b 2.16 reproduces #233; anchor 4b: every v001 gain reproduced
for r = PROD
    gp = sel_local(GG, "pipe", r, AB_REF);  gm = gainRows_local(G233, flt.family(r), flt.fc(r), flt.Tw(r));
    checkSame_local([gp.gainOLS; gp.gainIRLS], [gm.gainOLS; gm.gainIRLS], ANCHOR_TOL, sprintf("anchor 3, filter %d", r));
end
fprintf("Anchor 3 passed: pipeline gains at a/b %.2f reproduce #233\n", AB_REF);
A1 = renamevars(V1.GG(:, ["block" "filt" "ab" "f0" "gainOLS" "gainIRLS"]), ["gainOLS" "gainIRLS"], ["v1OLS" "v1IRLS"]);
M1 = innerjoin(A1, GG(:, ["block" "filt" "ab" "f0" "gainOLS" "gainIRLS"]), "Keys", ["block" "filt" "ab" "f0"]);
if height(M1) ~= height(A1), error("harmGain2:anchorV1", "%s", sprintf("%d of %d v001 gain rows found", height(M1), height(A1))); end
checkSame_local([M1.gainOLS; M1.gainIRLS], [M1.v1OLS; M1.v1IRLS], ANCHOR_TOL, "anchor 4b");
fprintf("Anchor 4b passed: all %d v001 gain rows reproduced\n", height(A1));

%% Theory curves for every filter and geometry
TH = cell(height(flt), numel(AB_SET));
for r = 1:height(flt)
    for i = 1:numel(AB_SET), TH{r, i} = theoryGain_local(W{i}, K, flt.family(r), flt.fc(r), flt.Tw(r), FS, F_FINE, W_MIN); end
end

%% Butterworth: phase and window (P1, P2), scaled f0/f_c
T_bw = table();
for r = iBW
    sc = 1 / flt.fc(r);
    for q = 1:numel(REG)
        nm = "gain" + REG_NAMES(q);
        fM = C233.f0_at_gain90(C233.family == "BW" & C233.fc == flt.fc(r) & C233.regression == REG_NAMES(q));
        ge = sel_local(GG, "emul", r, AB_REF);  gz = sel_local(GG, "emulZP", r, AB_REF);
        T_bw = [T_bw; table(flt.fc(r), REG_NAMES(q), fM * sc, onset_local(ge, nm, GAIN_REF) * sc, ...
            onset_local(gz, nm, GAIN_REF) * sc, onset_local(ge, nm + "loc", GAIN_REF) * sc, ...
            onset_local(gz, nm + "loc", GAIN_REF) * sc, thOn_local(TH, F_FINE, r, iRef, q, GAIN_REF) * sc, theoryFlag_local(REG_NAMES(q), irlsApprox(iRef)), ...
            'VariableNames', ["fc" "regression" "measured" "emul" "emulZP" "emulLoc" "emulZPLoc" "theory" "theoryStatus"])]; %#ok<AGROW>
    end
end
fprintf("\nButterworth onsets of gain %.2f, f0/f_c (a/b %.2f). Loc = local window %s.\n", GAIN_REF, AB_REF, mat2str(LOCAL_WIN));
disp(T_bw)

%% Savitzky-Golay: edges and window (P2, P3), scaled f0*T_w
condSG = [iSG(:) repmat(AB_REF, numel(iSG), 1); repmat(PROD(2), numel(abOther), 1) abOther(:)];
T_sg = table();
for c = 1:height(condSG)
    r = condSG(c, 1);  ab = condSG(c, 2);  i = find(AB_SET == ab);  sc = flt.Tw(r);
    gp = pipeRef_local(GG, G233, flt, r, ab, AB_REF);  gs = sel_local(GG, "pipeSafe", r, ab);  ge = sel_local(GG, "emul", r, ab);
    if height(gp) ~= height(ge) || height(gs) ~= height(ge), error("harmGain2:sgRows", "%s", sprintf("Condition %d: tempo grids differ", c)); end
    for q = 1:numel(REG)
        nm = "gain" + REG_NAMES(q);  ok = ge.(nm) >= 0.5;
        T_sg = [T_sg; table(flt.Tw(r), ab, REG_NAMES(q), onset_local(gp, nm, GAIN_REF) * sc, onset_local(gs, nm, GAIN_REF) * sc, ...
            onset_local(ge, nm, GAIN_REF) * sc, onset_local(ge, nm + "loc", GAIN_REF) * sc, thOn_local(TH, F_FINE, r, i, q, GAIN_REF) * sc, ...
            rms_local(gp.(nm)(ok) - ge.(nm)(ok)), rms_local(gs.(nm)(ok) - ge.(nm)(ok)), theoryFlag_local(REG_NAMES(q), irlsApprox(i)), ...
            'VariableNames', ["Tw" "ab" "regression" "pipe" "pipeSafe" "emul" "emulLoc" "theory" "rmsPipeVsEmul" "rmsSafeVsEmul" "theoryStatus"])]; %#ok<AGROW>
    end
end
fprintf("\nSavitzky-Golay onsets of gain %.2f, f0*T_w. pipe = clip %d (a/b %.2f: #233); pipeSafe = clip half-window + 1.\n", GAIN_REF, CLIP, AB_REF);
disp(T_sg)

%% Geometry (P4), production filters, Hz
T_geo = table();
for r = PROD
    for i = 1:numel(AB_SET)
        gp = sel_local(GG, "pipe", r, AB_SET(i));  ge = sel_local(GG, "emul", r, AB_SET(i));
        if flt.family(r) == "BW", gc = sel_local(GG, "emulZP", r, AB_SET(i)); gs = [];
        else, gc = ge; gs = sel_local(GG, "pipeSafe", r, AB_SET(i)); end
        for q = 1:numel(REG)
            nm = "gain" + REG_NAMES(q);  fC = onset_local(gc, nm + "loc", GAIN_REF);  fT = thOn_local(TH, F_FINE, r, i, q, GAIN_REF);
            fS = NaN;  if ~isempty(gs), fS = onset_local(gs, nm, GAIN_REF); end
            T_geo = [T_geo; table(flt.family(r), AB_SET(i), REG_NAMES(q), onset_local(gp, nm, GAIN_REF), fS, ...
                onset_local(ge, nm, GAIN_REF), fC, fT, (fC - fT) / fT, theoryFlag_local(REG_NAMES(q), irlsApprox(i)), ...
                'VariableNames', ["family" "ab" "regression" "pipe" "pipeSafe" "emul" "corrected" "theory" "relDiff" "theoryStatus"])]; %#ok<AGROW>
        end
    end
end
fprintf("\nProduction filters, onsets in Hz. corrected = local-window onset of emulZP (BW) or emul (SG).\n");
disp(T_geo)

%% Predictions
P = struct();
bwO = T_bw(T_bw.regression == "OLS", :);
P.P1 = max(bwO.emulZP) - min(bwO.emulZP) <= P1_SPREAD;
verdict_local("P1 phase", P.P1, sprintf("BW-OLS f0/f_c spread across f_c: emulZP %.4f (threshold %.3f); emul %.4f; theory %.4f", ...
    max(bwO.emulZP) - min(bwO.emulZP), P1_SPREAD, max(bwO.emul) - min(bwO.emul), max(bwO.theory) - min(bwO.theory)));
sgO = T_sg(T_sg.regression == "OLS" & T_sg.ab == AB_REF, :);
d2 = [bwO.emulZPLoc - bwO.theory; sgO.emulLoc - sgO.theory];
P.P2 = all(abs(d2) <= P2_ABS);
verdict_local("P2 window", P.P2, sprintf("OLS local-window onset minus theory (scaled), six filters: %s (threshold %.3f); #233 window: %s", ...
    mat2str(d2', 3), P2_ABS, mat2str([bwO.emulZP - bwO.theory; sgO.emul - sgO.theory]', 3)));
pLong = T_sg.rmsPipeVsEmul(T_sg.Tw == flt.Tw(rLong) & T_sg.ab == AB_REF & T_sg.regression == "OLS");
p125  = T_sg.rmsPipeVsEmul(T_sg.Tw == flt.Tw(PROD(2)) & T_sg.ab == AB_SET(1) & T_sg.regression == "OLS");
P.P3 = all(T_sg.rmsSafeVsEmul <= P3_RMS) && pLong > P3_RMS && p125 > P3_RMS;
verdict_local("P3 edges", P.P3, sprintf("pipeSafe RMS vs emulation max %.4f (threshold %.3f); clip-%d pipe RMS: long window %.4f, a/b 1.25 %.4f", ...
    max(T_sg.rmsSafeVsEmul), P3_RMS, CLIP, pLong, p125));
geoO = T_geo(T_geo.regression == "OLS", :);
P.P4 = all(abs(geoO.relDiff) <= P4_REL);
verdict_local("P4 geometry", P.P4, sprintf("OLS corrected-vs-theory relative difference: %s (threshold %.2f)", mat2str(geoO.relDiff', 3), P4_REL));

%% Figure
fg = figure("Color", "w", "Position", [80 80 1150 800]);  tl = tiledlayout(2, 2, "TileSpacing", "compact");
cB = [0.05 0.27 0.49];  cS = [0.80 0.40 0.00];  cT = [0.45 0.45 0.45];
nexttile; hold on;
plot(bwO.fc, bwO.measured, "o", "Color", cB, "DisplayName", "measured (#233)");
plot(bwO.fc, bwO.emul, "-x", "Color", cB, "DisplayName", "emulation");
plot(bwO.fc, bwO.emulZP, "-s", "Color", cS, "DisplayName", "emulation, delay removed");
plot(bwO.fc, bwO.emulZPLoc, "-d", "Color", cS, "MarkerFaceColor", cS, "DisplayName", "delay removed, local window");
plot(bwO.fc, bwO.theory, "--", "Color", cT, "LineWidth", 1.2, "DisplayName", "first-order theory");
set(gca, "XScale", "log", "XTick", flt.fc(iBW));  box on;  xlabel("f_c (Hz)");  ylabel("onset f_0 / f_c (OLS)");
title("P1/P2: Butterworth onset against cutoff");  legend("Location", "best", "FontSize", 7);
sgPanel_local(sel_local(GG, "emul", rLong, AB_REF), pipeRef_local(GG, G233, flt, rLong, AB_REF, AB_REF), ...
    sel_local(GG, "pipeSafe", rLong, AB_REF), F_FINE, TH{rLong, iRef}(:, 1), flt.Tw(rLong), GAIN_REF, ...
    sprintf("P3: SG T_w %.3f s, a/b %.2f (OLS)", flt.Tw(rLong), AB_REF), "f_0 \times T_w");
sgPanel_local(sel_local(GG, "emul", PROD(2), AB_SET(1)), sel_local(GG, "pipe", PROD(2), AB_SET(1)), ...
    sel_local(GG, "pipeSafe", PROD(2), AB_SET(1)), F_FINE, TH{PROD(2), 1}(:, 1), 1, GAIN_REF, ...
    sprintf("P3: SG production, a/b %.2f (OLS)", AB_SET(1)), "f_0 (Hz)");
nexttile; hold on;
for fmy = ["BW" "SG"]
    g = geoO(geoO.family == fmy, :);  cc = cB;  if fmy == "SG", cc = cS; end
    plot(g.ab, g.corrected, "-o", "Color", cc, "MarkerFaceColor", cc, "DisplayName", fmy + " corrected emulation");
    plot(g.ab, g.theory, "--", "Color", cc, "LineWidth", 1.2, "DisplayName", fmy + " theory");
    plot(g.ab, g.pipe, "x", "Color", cc, "MarkerSize", 8, "DisplayName", fmy + " pipeline (clip 20)");
end
box on;  xlabel("a/b");  ylabel("onset f_0 (Hz, OLS)");  title("P4: geometry, production filters");  legend("Location", "best", "FontSize", 7);
title(tl, "Zero-noise pull towards 1/3: the three residuals of v001");
set(findall(fg, "Type", "axes"), "Toolbar", []);  exportgraphics(fg, OUT_PNG, "Resolution", 200);

save(OUT_MAT, "J", "GG", "W", "DG", "irlsApprox", "TH", "T_bw", "T_sg", "T_geo", "P", "F_FINE", "K", "AB_SET", ...
    "AB_REF", "LOCAL_WIN", "CLIP", "P1_SPREAD", "P2_ABS", "P3_RMS", "P4_REL", "-v7.3");
fprintf("Saved: %s\nFigure: %s\n", OUT_MAT, OUT_PNG);

%% =========================================================================
function T = jobs_local(blk, filts, abv, F0, BETA)
[fi, ai, ki, bi] = ndgrid(filts, abv, 1:numel(F0), 1:numel(BETA));
T = table(repmat(blk, numel(fi), 1), fi(:), ai(:), F0(ki(:))', BETA(bi(:))', ...
    'VariableNames', ["block" "filt" "ab" "f0" "betaGen"]);
end

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
dP = ([xp yp] - [xm ym]) / (2 * epsB);  dvF = d(dP, Dv);  daF = d(dP, Da);
totV = (betaKin_local(v0 + epsB * dvF, a0, REG) - betaKin_local(v0 - epsB * dvF, a0, REG)) / (2 * epsB);
totA = (betaKin_local(v0, a0 + epsB * daF, REG) - betaKin_local(v0, a0 - epsB * daF, REG)) / (2 * epsB);
Fp = fft(dP);  E = sum(abs(Fp).^2, "all");  binw = FS / M;
K = 1:2:KMAX;  W = NaN(2, numel(K), numel(REG));  eK = NaN(1, numel(K));
for j = 1:numel(K)
    m = abs(abs(f) - K(j)) < binw / 2;
    eK(j) = sum(abs(Fp(m, :)).^2, "all") / E;
    dk = real(ifft(m .* Fp));  dv = d(dk, Dv);  da = d(dk, Da);
    W(1, j, :) = (betaKin_local(v0 + epsB * dv, a0, REG) - betaKin_local(v0 - epsB * dv, a0, REG)) / (2 * epsB);
    W(2, j, :) = (betaKin_local(v0, a0 + epsB * da, REG) - betaKin_local(v0, a0 - epsB * da, REG)) / (2 * epsB);
end
if any(~isfinite(W), "all") || any(~isfinite(tot))
    error("harmGain2:weightNaN", "%s", sprintf("a/b %.2f: a weight regression returned NaN", ab));
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
    Hv = E * (factorial(1) / (-1/FS)^1 * g(:, 2));
    Ha = E * (factorial(2) / (-1/FS)^2 * g(:, 3));
end
end

function L = sgFrame_local(Tw, FS)
L = round(Tw * FS);  L = L + mod(L + 1, 2);      % as #233
end

function g = theoryGain_local(Wab, K, fam, fc, Tw, FS, f0, wMin)
% First-order gain from the weights and the chain's magnitude response relative to ideal.
f0 = f0(:);  nq = size(Wab, 3);  g = NaN(numel(f0), nq);
kUse = K(max(abs(Wab), [], [1 3]) > wMin);
for i = 1:numel(f0)
    if max(kUse) * f0(i) >= FS/2, continue, end  % would alias: theory not defined, left NaN
    fq = [f0(i); K(:) * f0(i)];  [Hv, Ha] = chainResponse_local(fam, fc, Tw, FS, fq);
    D = 1i * 2 * pi * fq;  Rv = abs(Hv ./ D);  Ra = abs(Ha ./ D.^2);
    rv = Rv(2:end) / Rv(1);  ra = Ra(2:end) / Ra(1);
    for q = 1:nq, g(i, q) = sum(Wab(1, :, q)' .* rv) + sum(Wab(2, :, q)' .* ra); end
end
end

function b = emulateBeta_local(A_MM, ab, f0, bt, fam, fc, Tw, FS, NCYC, REG, zeroPhase)
% The pipeline's exact linear response applied to a whole-cycle (exactly periodic) record.
% zeroPhase: remove the forward difference's half-sample delay (BW only); magnitude unchanged.
M = round(NCYC * FS / f0);  f0e = NCYC * FS / M;
[x, y] = generatePowerLawEllipse_v001(A_MM, A_MM / ab, f0e, FS, bt, 'M', M);
f = signedBins_local(M, FS);  [Hv, Ha] = chainResponse_local(fam, fc, Tw, FS, f);
if zeroPhase
    if fam ~= "BW", error("harmGain2:zpFamily", "%s", "Zero-phase variant is defined for BW only"); end
    Hv = Hv .* exp(1i * pi * f / FS);
end
nyq = abs(abs(f) - FS/2) < 1e-9;  Hv(nyq) = 0;  Ha(nyq) = 0;
X = fft([x(:) y(:)]);
b = betaKin_local(real(ifft(Hv .* X)), real(ifft(Ha .* X)), REG);
end

function b = pipelineBeta_local(x, y, fam, fc, tw, FS, clip, REG)
% #233's code path, with the clip as a parameter.
if fam == "BW", ft = 2; fp = [2 fc 1];
else, fr = sgFrame_local(tw, FS); ft = 4; fp = [4 fr]; end
[dx, dy] = differentiateKinematicsEBR(x(:), y(:), ft, fp, FS);
b = betaKin_local([dx(clip:end-clip, 2) dy(clip:end-clip, 2)], [dx(clip:end-clip, 3) dy(clip:end-clip, 3)], REG);
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

function g = sel_local(GG, blk, r, ab)
g = sortrows(GG(GG.block == blk & GG.filt == r & GG.ab == ab, :), "f0");
if isempty(g), error("harmGain2:sel", "%s", sprintf("No gains for %s, filter %d, a/b %.2f", blk, r, ab)); end
end

function g = pipeRef_local(GG, G233, flt, r, ab, abRef)
% Clip-20 pipeline gains: #233's own table at the reference geometry, this run's elsewhere.
if ab == abRef, g = gainRows_local(G233, flt.family(r), flt.fc(r), flt.Tw(r));
else, g = sel_local(GG, "pipe", r, ab); end
end

function gm = gainRows_local(G, fam, fc, Tw)
gm = sortrows(G(G.family == fam & isequaln_local(G.fc, fc) & isequaln_local(G.Tw, Tw), :), "f0");
end

function f = thOn_local(TH, F_FINE, r, i, q, ref)
f = crossing_local(F_FINE(:), TH{r, i}(:, q), ref);
end

function f = onset_local(g, col, ref)
f = crossing_local(g.f0, g.(col), ref);
end

function checkSame_local(a, b, tol, tag)
if numel(a) ~= numel(b) || ~isequal(isnan(a), isnan(b)) || max(abs(a - b), [], "omitnan") > tol
    error("harmGain2:anchor", "%s", sprintf("%s: gains differ (max %.3g)", tag, max(abs(a - b), [], "omitnan")));
end
end

function sgPanel_local(ge, gp, gs, F_FINE, th, sc, ref, ttl, xl)
nexttile; hold on;
plot(gp.f0 * sc, gp.gainOLS, "o", "Color", [0.05 0.27 0.49], "DisplayName", "pipeline, clip 20");
plot(gs.f0 * sc, gs.gainOLS, "s", "Color", [0.80 0.40 0.00], "DisplayName", "pipeline, clip half-window + 1");
plot(ge.f0 * sc, ge.gainOLS, "-", "Color", [0.2 0.2 0.2], "DisplayName", "emulation (no record ends)");
plot(F_FINE * sc, th, "--", "Color", [0.45 0.45 0.45], "LineWidth", 1.2, "DisplayName", "first-order theory");
yline(ref, ":", "HandleVisibility", "off");  set(gca, "XScale", "log");  ylim([0 1.1]);  box on;
xlabel(xl);  ylabel("gain (OLS, #233 window)");  title(ttl);  legend("Location", "southwest", "FontSize", 7);
end

function s = countList_local(names, counts)
% "name n, name n, ..." from a string array and matching counts.
s = strjoin(arrayfun(@(i) sprintf("%s %d", names(i), counts(i)), 1:numel(names)), ", ");
end

function verdict_local(tag, held, detail)
s = "NOT HELD";  if held, s = "HELD"; end
fprintf("%s %s: %s\n", tag, s, detail);
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
