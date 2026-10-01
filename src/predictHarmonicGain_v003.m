%% predictHarmonicGain_v003.m
% Settles the residual left by predictHarmonicGain_v002 (P2 and P4 NOT HELD: theory onsets 1-2%
% above emulation). An exploratory read of v002's saved results (2026-10-01) traced it to the
% comparison, not the theory: emulated and measured onsets were interpolated between #233's
% tempo nodes (a factor of 1.26 apart), where the chord of the concave gain curve crosses 0.9
% early, while the theory was crossed on a fine grid. Sampled at the same nodes, the OLS residual
% fell to 0.08-0.30%. This script makes that a scripted result and removes the node grid:
%  Block 1  the node-sampled comparison, from v002's saved gains and theory.
%  Block 2  fine tempo grids (N_FINE log-spaced tempi over SPAN x the theory's OLS onset, plus the
%           #233 nodes inside that range) for #233's six filters at a/b 2.16 and the production
%           filters at a/b 1.25 and 4.1: the registered pipeline (clip 20, as #233), the edge-safe
%           pipeline (SG only) and the corrected emulation (BW: half-sample delay removed; SG as is).
% Predictions, fixed before the run (2026-10-01), reported HELD / NOT HELD:
%  Q1 theory: on the fine grid, OLS local-window onsets of the corrected emulation lie within
%             Q1_REL of the first-order theory in all ten conditions.
%  Q2 #233  : #233's node-interpolated onsets are biased low: for all six filters and both
%             regressions the fine-grid pipeline onset (#233's definition: clip 20, gain over
%             beta_gen 0.2-0.5) exceeds #233's, and the OLS biases lie within Q2_BAND percent.
% Anchors (error if broken): (1) v002's theory curves reproduced; (2) at the #233 nodes inside each
% fine grid, pipeline gains equal #233's (a/b 2.16) or v002's (other a/b), and edge-safe and
% corrected-emulation gains equal v002's, within TOL_GAIN; (0, on resume only) two jobs recomputed
% from scratch match the checkpoint within TOL_GAIN. TOL_GAIN is tiered: OLS is closed form, IRLS
% stops at fitnlm's iteration tolerance, and #233 ran on the iMac while this runs on the lab PC.
% First run (2026-10-01): a single 1e-10 bound stopped it at IRLS 1.56e-10 (condition 4, SG
% T_w 0.0875, against #233), after the compute. For the same machine pair, zero-noise IRLS betas
% differed by at most 3.9e-9 (checkBWFDHalfSampleNoise_v001, anchor 2).
% Checkpoint: the raw run is saved to CKPT straight after the parfor, before any anchor; if CKPT
% exists and its job table matches exactly, the run is not repeated. The first run's checkpoint
% was saved by hand from the workspace (J only); anchor 0 is what verifies it.
% Zero noise, exact ellipse, #233 settings (fs 240 Hz, 10 cycles). Roughly 9,000 trajectories.
% Run from src/. Writes results/predictHarmonicGain_v003.mat, figures/predictHarmonicGain_v003.png

%% CONFIG
ROOT       = fileparts(fileparts(mfilename("fullpath")));
IN_233     = fullfile(ROOT, "results", "filterCutoffCollapse_v001.mat");
IN_V002    = fullfile(ROOT, "results", "predictHarmonicGain_v002.mat");
OUT_MAT    = fullfile(ROOT, "results", "predictHarmonicGain_v003.mat");
OUT_PNG    = fullfile(ROOT, "figures", "predictHarmonicGain_v003.png");
A_MM       = 50;  AB_REF = 2.16;  NCYC = 10;  CLIP = 20;    % as #233
PROD       = [2 5];                                         % #233 flt rows: BW 10 Hz, SG 41 samples
LOCAL_WIN  = [0.29 0.37];                                   % as v002
W_MIN      = 1e-3;                                          % as v001/v002
N_FINE     = 33;  SPAN = [0.85 1.6];                        % about 2% tempo spacing
REG        = [3 5];  REG_NAMES = ["OLS" "IRLS"];
TOL_GAIN   = [1e-12 1e-7];                              % OLS, IRLS (see header)
CKPT       = fullfile(ROOT, "results", "predictHarmonicGain_v003_runCheckpoint.mat");
Q1_REL     = 0.005;  Q2_BAND = [0.5 3];

%% Setup
addpath(genpath(fullfile(ROOT, "src", "functions")));  addpath(genpath(fullfile(ROOT, "src", "req")));
for f = [IN_233 IN_V002], if ~isfile(f), error("harmGain3:in", "%s", "Missing input: " + f); end, end
for f = [OUT_MAT OUT_PNG], if ~isfolder(fileparts(f)), error("harmGain3:outDir", "%s", "Missing folder: " + fileparts(f)); end, end
R233 = load(IN_233, "G", "C", "flt", "F0", "BETA", "FS", "GAIN_WIN", "GAIN_REF");
FS = R233.FS;  F0 = R233.F0(:);  BETA = R233.BETA;  GAIN_WIN = R233.GAIN_WIN;  GAIN_REF = R233.GAIN_REF;
flt = R233.flt;  G233 = R233.G;  C233 = R233.C;
V2 = load(IN_V002, "W", "TH", "GG", "AB_SET", "F_FINE", "K", "irlsApprox");
AB_SET = V2.AB_SET;  F_FINE = V2.F_FINE(:);  K = V2.K;
iRef = find(AB_SET == AB_REF);  abOther = AB_SET(AB_SET ~= AB_REF);
if FS ~= 240 || height(flt) ~= 6 || isempty(iRef) || flt.fc(PROD(1)) ~= 10 || abs(flt.Tw(PROD(2)) * FS - 41) > 1e-9
    error("harmGain3:config", "%s", "#233 or v002 settings differ from what this script assumes");
end
BSUB = BETA(BETA >= GAIN_WIN(1) - 1e-12 & BETA <= GAIN_WIN(2) + 1e-12);
if nnz(BSUB >= LOCAL_WIN(1) & BSUB <= LOCAL_WIN(2)) ~= 3, error("harmGain3:localWin", "%s", "LOCAL_WIN must hold three nodes"); end
COND = [(1:height(flt))' repmat(AB_REF, height(flt), 1); repmat(PROD(:), numel(abOther), 1) repelem(abOther(:), numel(PROD))];

%% Anchor 1: v002's theory curves reproduced
for r = 1:height(flt)
    for i = 1:numel(AB_SET)
        th = theoryGain_local(V2.W{i}, K, flt.family(r), flt.fc(r), flt.Tw(r), FS, F_FINE, W_MIN);
        checkSame_local(th(:), V2.TH{r, i}(:), 1e-12, sprintf("anchor 1, filter %d, a/b %.2f", r, AB_SET(i)));
    end
end
fprintf("Anchor 1 passed: v002 theory curves reproduced\n");

%% Block 1: theory sampled at #233's nodes, crossed by the same rule (from v002's saved gains)
T_nodes = table();
for c = 1:height(COND)
    r = COND(c, 1);  i = find(AB_SET == COND(c, 2));
    g = sel_local(V2.GG, corrBlk_local(flt, r), r, COND(c, 2));
    for q = 1:numel(REG)
        th = V2.TH{r, i}(:, q);  ok = isfinite(th);
        thN = interp1(log(F_FINE(ok)), th(ok), log(F0), "pchip", NaN);
        fF = crossing_local(F_FINE, th, GAIN_REF);  fN = crossing_local(F0, thN, GAIN_REF);
        fE = crossing_local(g.f0, g.("gain" + REG_NAMES(q) + "loc"), GAIN_REF);
        T_nodes = [T_nodes; table(flt.family(r), r, COND(c, 2), REG_NAMES(q), fF, fN, fE, 100 * (fE - fF) / fF, ...
            100 * (fE - fN) / fN, theoryFlag_local(REG_NAMES(q), V2.irlsApprox(i)), 'VariableNames', ["family" "filt" "ab" ...
            "regression" "theoryFine" "theoryNodes" "emulNodes" "residFinePct" "residNodesPct" "theoryStatus"])]; %#ok<AGROW>
    end
end
fprintf("\nBlock 1: onsets (Hz) of the corrected emulation's local-window gain on #233's nodes, against theory\n");
fprintf("crossed on the fine grid and on the same nodes (residual, %%).\n");  disp(T_nodes)

%% Block 2 jobs: fine tempo grids
J = table();
for c = 1:height(COND)
    r = COND(c, 1);  ab = COND(c, 2);  i = find(AB_SET == ab);
    fT = crossing_local(F_FINE, V2.TH{r, i}(:, 1), GAIN_REF);
    if ~isfinite(fT), error("harmGain3:fT", "%s", sprintf("No theory onset for condition %d", c)); end
    tf = fT * logspace(log10(SPAN(1)), log10(SPAN(2)), N_FINE)';
    tf = unique([tf; F0(F0 >= min(tf) & F0 <= max(tf))]);
    blks = ["pipe" "emulC"];  if flt.family(r) == "SG", blks = ["pipe" "pipeSafe" "emulC"]; end
    for b = blks
        [ti, bi] = ndgrid(1:numel(tf), 1:numel(BSUB));
        J = [J; table(repmat(b, numel(ti), 1), repmat(c, numel(ti), 1), repmat(r, numel(ti), 1), repmat(ab, numel(ti), 1), ...
            tf(ti(:)), BSUB(bi(:))', 'VariableNames', ["block" "cond" "filt" "ab" "f0" "betaGen"])]; %#ok<AGROW>
    end
end
[gb, gBlk] = findgroups(J.block);
fprintf("\nBlock 2 jobs: %d trajectories (%s), sigma = 0\n", height(J), countList_local(gBlk, splitapply(@numel, J.cond, gb)));
KEYS = ["block" "cond" "filt" "ab" "f0" "betaGen"];
fam = flt.family;  fcv = flt.fc;  twv = flt.Tw;
resumed = isfile(CKPT);
if resumed
    S = load(CKPT, "J");
    if ~all(ismember([KEYS "bOLS" "bIRLS" "err"], string(S.J.Properties.VariableNames))) || ~isequaln(S.J(:, KEYS), J(:, KEYS))
        error("harmGain3:ckpt", "%s", "Checkpoint jobs differ from this configuration; move it aside to rerun: " + CKPT);
    end
    J = S.J;
    fprintf("Resumed from checkpoint (run not repeated): %s\n", CKPT);
    for k = round(linspace(1, height(J), 2))                % anchor 0: first and last job, from scratch
        [b, ef] = jobBeta_local(J.block(k), J.filt(k), J.ab(k), J.f0(k), J.betaGen(k), fam, fcv, twv, A_MM, FS, NCYC, CLIP, REG);
        if ef ~= J.err(k), error("harmGain3:anchor0", "%s", sprintf("Job %d: error flag differs from the checkpoint", k)); end
        checkSame_local(b(:), [J.bOLS(k); J.bIRLS(k)], TOL_GAIN(:), sprintf("anchor 0, job %d", k));
    end
    fprintf("Anchor 0 passed: two jobs recomputed and match the checkpoint\n");
else
    bRec = NaN(height(J), numel(REG));  errFlag = false(height(J), 1);
    blk = J.block;  fl = J.filt;  abv = J.ab;  f0v = J.f0;  bgv = J.betaGen;
    nW = 0;  if ~isempty(gcp("nocreate")), nW = gcp("nocreate").NumWorkers; end
    parfor (k = 1:height(J), nW)
        [b, ef] = jobBeta_local(blk(k), fl(k), abv(k), f0v(k), bgv(k), fam, fcv, twv, A_MM, FS, NCYC, CLIP, REG);
        bRec(k, :) = b;  errFlag(k) = ef;
    end
    J.bOLS = bRec(:, 1);  J.bIRLS = bRec(:, 2);  J.err = errFlag;
    save(CKPT, "J", "-v7.3");
    fprintf("Checkpoint written: %s\n", CKPT);
end
fprintf("Non-finite estimates by block: %s; pipeline trajectories that threw: %d\n", ...
    countList_local(gBlk, splitapply(@(b) nnz(~isfinite(b)), [J.bOLS J.bIRLS], gb)), sum(errFlag));

%% Gains on the fine grids
[grp, GF] = findgroups(J(:, ["block" "cond" "filt" "ab" "f0"]));
for q = 1:numel(REG)
    y = J.("b" + REG_NAMES(q));
    GF.("gain" + REG_NAMES(q))         = splitapply(@(b, v) gainAt_local(b, v, GAIN_WIN),  J.betaGen, y, grp);
    GF.("gain" + REG_NAMES(q) + "loc") = splitapply(@(b, v) gainAt_local(b, v, LOCAL_WIN), J.betaGen, y, grp);
end

%% Anchor 2: gains at #233's nodes reproduce #233 and v002
nAnch = 0;  dev2 = zeros(1, numel(REG));
for c = 1:height(COND)
    r = COND(c, 1);  ab = COND(c, 2);
    for b = unique(GF.block(GF.cond == c))'
        g = GF(GF.block == b & GF.cond == c & ismember(GF.f0, F0), :);
        if isempty(g), continue, end
        switch b
            case "pipe",     ref = pipeRef_local(V2.GG, G233, flt, r, ab, AB_REF);  sfx = "";
            case "pipeSafe", ref = sel_local(V2.GG, "pipeSafe", r, ab);              sfx = "";
            otherwise,       ref = sel_local(V2.GG, corrBlk_local(flt, r), r, ab);  sfx = "loc";
        end
        ref = ref(ismember(ref.f0, g.f0), :);  g = sortrows(g, "f0");
        if height(ref) ~= height(g), error("harmGain3:anchor2", "%s", sprintf("Condition %d %s: node rows differ", c, b)); end
        for q = 1:numel(REG)
            nm = "gain" + REG_NAMES(q) + sfx;
            dev2(q) = max(dev2(q), max(abs(g.(nm) - ref.(nm)), [], "omitnan"));
            checkSame_local(g.(nm), ref.(nm), TOL_GAIN(q), sprintf("anchor 2, condition %d, %s, %s", c, b, nm));
        end
        nAnch = nAnch + height(g);
    end
end
fprintf("Anchor 2 passed: %d node gains reproduce #233 and v002 (max |difference| %s)\n", nAnch, ...
    strjoin(arrayfun(@(q) sprintf("%s %.3g, tol %.0e", REG_NAMES(q), dev2(q), TOL_GAIN(q)), 1:numel(REG)), "; "));

%% Fine-grid onsets
T_fine = table();
for c = 1:height(COND)
    r = COND(c, 1);  ab = COND(c, 2);  i = find(AB_SET == ab);
    for q = 1:numel(REG)
        nm = "gain" + REG_NAMES(q);
        fT = crossing_local(F_FINE, V2.TH{r, i}(:, q), GAIN_REF);
        fP = onset_local(GF, "pipe", c, nm, GAIN_REF);
        fS = NaN;  if flt.family(r) == "SG", fS = onset_local(GF, "pipeSafe", c, nm, GAIN_REF); end
        fC = onset_local(GF, "emulC", c, nm + "loc", GAIN_REF);
        T_fine = [T_fine; table(flt.family(r), r, flt.fc(r), flt.Tw(r), ab, REG_NAMES(q), fP, fS, fC, fT, (fC - fT) / fT, ...
            theoryFlag_local(REG_NAMES(q), V2.irlsApprox(i)), 'VariableNames', ["family" "filt" "fc" "Tw" "ab" "regression" ...
            "pipe" "pipeSafe" "emulCLoc" "theory" "relDiff" "theoryStatus"])]; %#ok<AGROW>
    end
end
fprintf("\nFine-grid onsets (Hz). pipe: clip %d, #233 window; pipeSafe: clip half-window + 1; emulCLoc: corrected emulation, local window.\n", CLIP);
disp(T_fine)

%% #233 corrected: node-interpolated against fine-grid onsets, scaled as #233
T_233 = T_fine(T_fine.ab == AB_REF, ["family" "filt" "fc" "Tw" "regression" "pipe" "pipeSafe"]);
T_233.node233 = NaN(height(T_233), 1);
for k = 1:height(T_233)
    T_233.node233(k) = C233.f0_at_gain90(C233.family == T_233.family(k) & isequaln_local(C233.fc, T_233.fc(k)) & ...
        isequaln_local(C233.Tw, T_233.Tw(k)) & C233.regression == T_233.regression(k));
end
T_233.biasPct = 100 * (T_233.pipe - T_233.node233) ./ T_233.node233;
sc = 1 ./ T_233.fc;  sc(T_233.family == "SG") = T_233.Tw(T_233.family == "SG");
T_233.scaledNode = T_233.node233 .* sc;  T_233.scaledFine = T_233.pipe .* sc;  T_233.scaledFineSafe = T_233.pipeSafe .* sc;
fprintf("\n#233 onsets: node-interpolated (#233) against fine grid; scaled f0/f_c (BW) or f0*T_w (SG).\n");  disp(T_233)
S233 = groupsummary(T_233, ["family" "regression"], ["mean" "std"], ["scaledNode" "scaledFine" "scaledFineSafe"]);
S233.cvNode = S233.std_scaledNode ./ S233.mean_scaledNode;  S233.cvFine = S233.std_scaledFine ./ S233.mean_scaledFine;
S233.cvFineSafe = S233.std_scaledFineSafe ./ S233.mean_scaledFineSafe;
disp(S233(:, ["family" "regression" "mean_scaledNode" "cvNode" "mean_scaledFine" "cvFine" "mean_scaledFineSafe" "cvFineSafe"]))
bwF = T_233(T_233.family == "BW" & T_233.regression == "OLS", :);
fprintf("BW-OLS cutoff multiple (1 / scaled onset): node-interpolated %.2f-%.2f, fine grid %.2f-%.2f\n", ...
    1 / max(bwF.scaledNode), 1 / min(bwF.scaledNode), 1 / max(bwF.scaledFine), 1 / min(bwF.scaledFine));

%% Predictions
P = struct();
o = T_fine(T_fine.regression == "OLS", :);
P.Q1 = all(abs(o.relDiff) <= Q1_REL);
verdict_local("Q1 theory", P.Q1, sprintf("OLS fine-grid corrected emulation vs theory, relative: %s (threshold %.3f)", mat2str(o.relDiff', 3), Q1_REL));
bO = T_233.biasPct(T_233.regression == "OLS");
P.Q2 = all(T_233.biasPct > 0) && all(bO >= Q2_BAND(1) & bO <= Q2_BAND(2));
verdict_local("Q2 #233", P.Q2, sprintf("fine minus node onset, %%: OLS %s, IRLS %s (OLS band %s)", mat2str(bO', 3), ...
    mat2str(T_233.biasPct(T_233.regression == "IRLS")', 3), mat2str(Q2_BAND)));

%% Figure
fg = figure("Color", "w", "Position", [80 80 1150 800]);  tl = tiledlayout(2, 2, "TileSpacing", "compact");
cols = [0.52 0.72 0.92; 0.22 0.54 0.87; 0.05 0.27 0.49];
for fmy = ["BW" "SG"]
    nexttile; hold on;  rows = find(flt.family == fmy)';
    for s = 1:numel(rows)
        r = rows(s);  c = find(COND(:, 1) == r & COND(:, 2) == AB_REF);  i = iRef;  scl = scale_local(flt(r, :));
        gp = sortrows(GF(GF.block == "pipe" & GF.cond == c, :), "f0");  ge = sortrows(GF(GF.block == "emulC" & GF.cond == c, :), "f0");
        gm = gainRows_local(G233, fmy, flt.fc(r), flt.Tw(r));
        lab = sprintf("f_c %g Hz", flt.fc(r));  if fmy == "SG", lab = sprintf("T_w %.3f s", flt.Tw(r)); end
        plot(gm.f0 * scl, gm.gainOLS, "o", "Color", cols(s, :), "MarkerSize", 6, "DisplayName", lab + " #233 nodes");
        plot(gp.f0 * scl, gp.gainOLS, "-", "Color", cols(s, :), "HandleVisibility", "off");
        plot(ge.f0 * scl, ge.gainOLSloc, ":", "Color", cols(s, :), "LineWidth", 1.2, "HandleVisibility", "off");
        plot(F_FINE * scl, V2.TH{r, i}(:, 1), "--", "Color", [0.5 0.5 0.5], "HandleVisibility", "off");
    end
    yline(GAIN_REF, ":", "HandleVisibility", "off");  set(gca, "XScale", "log");  box on;
    xr = [0.2 0.5];  if fmy == "BW", xr = [0.1 0.3]; end
    xlim(xr);  ylim([0.6 1.05]);
    if fmy == "BW", xlabel("f_0 / f_c"); title("Butterworth: fine grid around the onset");
    else, xlabel("f_0 \times T_w"); title("Savitzky-Golay: fine grid around the onset"); end
    ylabel("gain (OLS)");  legend("Location", "southwest", "FontSize", 7);
    text(0.98, 0.97, "o #233   - pipeline (fine)   : corrected emulation   -- theory", "Units", "normalized", ...
        "HorizontalAlignment", "right", "FontSize", 7);
end
nexttile;  t2 = T_233(T_233.regression == "OLS", :);
bar(categorical(t2.family + " " + string(t2.filt)), t2.biasPct);  box on;
ylabel("fine minus node onset (%)");  title("#233's node-interpolation bias (OLS)");
nexttile;  hold on;  tn = T_nodes(T_nodes.regression == "OLS", :);  xi = 1:height(tn);
plot(xi, tn.residFinePct, "o-", "DisplayName", "v002: emulation on nodes vs fine theory");
plot(xi, tn.residNodesPct, "s-", "DisplayName", "v002: emulation on nodes vs theory on nodes");
plot(xi, 100 * o.relDiff, "d-", "DisplayName", "v003: both on the fine grid");
yline(0, ":", "HandleVisibility", "off");  box on;  xticks(xi);
xticklabels(tn.family + string(tn.filt) + " a/b " + compose("%.2f", tn.ab));  xtickangle(45);
ylabel("emulation minus theory (%)");  title("The residual was the comparison");  legend("Location", "southwest", "FontSize", 7);
title(tl, "First-order harmonic theory on fine tempo grids");
set(findall(fg, "Type", "axes"), "Toolbar", []);  exportgraphics(fg, OUT_PNG, "Resolution", 200);

save(OUT_MAT, "J", "GF", "COND", "T_nodes", "T_fine", "T_233", "S233", "P", "N_FINE", "SPAN", "LOCAL_WIN", "CLIP", ...
    "Q1_REL", "Q2_BAND", "-v7.3");
fprintf("Saved: %s\nFigure: %s\n", OUT_MAT, OUT_PNG);

%% =========================================================================
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
% First-order gain from the weights and the chain's magnitude response relative to ideal (as v002).
f0 = f0(:);  nq = size(Wab, 3);  g = NaN(numel(f0), nq);
kUse = K(max(abs(Wab), [], [1 3]) > wMin);
for i = 1:numel(f0)
    if max(kUse) * f0(i) >= FS/2, continue, end
    fq = [f0(i); K(:) * f0(i)];  [Hv, Ha] = chainResponse_local(fam, fc, Tw, FS, fq);
    D = 1i * 2 * pi * fq;  Rv = abs(Hv ./ D);  Ra = abs(Ha ./ D.^2);
    rv = Rv(2:end) / Rv(1);  ra = Ra(2:end) / Ra(1);
    for q = 1:nq, g(i, q) = sum(Wab(1, :, q)' .* rv) + sum(Wab(2, :, q)' .* ra); end
end
end

function f = signedBins_local(M, FS)
f = (0:M-1)' * FS / M;  f(f > FS/2) = f(f > FS/2) - FS;
end

function b = emulateBeta_local(A_MM, ab, f0, bt, fam, fc, Tw, FS, NCYC, REG, zeroPhase)
% The pipeline's exact linear response applied to a whole-cycle record (as v002).
M = round(NCYC * FS / f0);  f0e = NCYC * FS / M;
[x, y] = generatePowerLawEllipse_v001(A_MM, A_MM / ab, f0e, FS, bt, 'M', M);
f = signedBins_local(M, FS);  [Hv, Ha] = chainResponse_local(fam, fc, Tw, FS, f);
if zeroPhase
    if fam ~= "BW", error("harmGain3:zpFamily", "%s", "Zero-phase variant is defined for BW only"); end
    Hv = Hv .* exp(1i * pi * f / FS);
end
nyq = abs(abs(f) - FS/2) < 1e-9;  Hv(nyq) = 0;  Ha(nyq) = 0;
X = fft([x(:) y(:)]);
b = betaKin_local(real(ifft(Hv .* X)), real(ifft(Ha .* X)), REG);
end

function [b, ef] = jobBeta_local(blk, r, ab, f0, bg, fam, fcv, twv, A_MM, FS, NCYC, CLIP, REG)
% One trajectory: the registered pipeline (clip 20 or edge-safe) or the corrected emulation.
b = NaN(1, numel(REG));  ef = false;
if blk == "pipe" || blk == "pipeSafe"
    [x, y] = generatePowerLawEllipse_v001(A_MM, A_MM / ab, f0, FS, bg, 'nCycles', NCYC);
    cl = CLIP;  if blk == "pipeSafe", cl = (sgFrame_local(twv(r), FS) - 1) / 2 + 1; end
    try
        b = pipelineBeta_local(x, y, fam(r), fcv(r), twv(r), FS, cl, REG);
    catch
        ef = true;   % counted and reported; never silently absorbed
    end
else
    b = emulateBeta_local(A_MM, ab, f0, bg, fam(r), fcv(r), twv(r), FS, NCYC, REG, fam(r) == "BW");
end
end

function b = pipelineBeta_local(x, y, fam, fc, tw, FS, clip, REG)
% #233's code path, with the clip as a parameter (as v002).
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
if isempty(g), error("harmGain3:sel", "%s", sprintf("No gains for %s, filter %d, a/b %.2f", blk, r, ab)); end
end

function g = pipeRef_local(GG, G233, flt, r, ab, abRef)
% Clip-20 pipeline gains: #233's own table at the reference geometry, v002's elsewhere.
if ab == abRef, g = gainRows_local(G233, flt.family(r), flt.fc(r), flt.Tw(r));
else, g = sel_local(GG, "pipe", r, ab); end
end

function gm = gainRows_local(G, fam, fc, Tw)
gm = sortrows(G(G.family == fam & isequaln_local(G.fc, fc) & isequaln_local(G.Tw, Tw), :), "f0");
end

function f = onset_local(GF, blk, c, col, ref)
g = sortrows(GF(GF.block == blk & GF.cond == c, :), "f0");
if isempty(g), error("harmGain3:onset", "%s", sprintf("No fine-grid gains for %s, condition %d", blk, c)); end
f = crossing_local(g.f0, g.(col), ref);
end

function checkSame_local(a, b, tol, tag)
% tol: scalar, or one value per element of a.
if numel(a) ~= numel(b) || ~isequal(isnan(a), isnan(b)) || any(abs(a(:) - b(:)) > tol(:))
    error("harmGain3:anchor", "%s", sprintf("%s: values differ (max %.3g)", tag, max(abs(a - b), [], "omitnan")));
end
end

function s = scale_local(row)
if row.family == "BW", s = 1 / row.fc; else, s = row.Tw; end
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

function b = corrBlk_local(flt, r)
% v002's corrected emulation block: half-sample delay removed for BW; SG has none.
b = "emul";  if flt.family(r) == "BW", b = "emulZP"; end
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
