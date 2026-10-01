%% checkBWFDHalfSampleNoise_v002.m
% Finishes #250. checkBWFDHalfSampleNoise_v001 found that realigning BWFD's velocity and acceleration
% (CEN, AVG) moves the recovered exponent by up to 0.013 (Block C) and that the shift looked
% deterministic: CEN's gain change under noise matched its zero-noise value. Two figures behind that
% reading were exploratory (a split of v001's saved results) and are not printed by v001:
%   (1) the noise-only part of each paired delta (noisy median paired delta minus the same cell's
%       zero-noise delta, from the sigma = 0 companion job);
%   (2) the beta_gen-band table (where in beta_gen the shift lives).
% This script reads v001's saved results, computes both, and gates them against predictions.
% It recomputes no trajectory, so it takes seconds.
% Predictions, fixed before the run (2026-10-01). They restate the exploratory read, so a HELD
% here replicates that read on the same data; it is not an independent test.
%  R1 noise-only part small: Block C, CEN, OLS and LMLS: max |noise-only| < NOISE_FRAC * MDC/2.77
%     (exploratory: <= 0.0004 OLS, <= 0.0007 LMLS). IRLS is reported, not gated (its failures and
%     non-additivity make zero-noise subtraction less clean, #249).
%  R2 the cells that reach MDC/2.77 in v001 (N1) have |noise-only| < NOISE_FRAC * MDC/2.77.
%  R3 beta_gen band, Block C, CEN, OLS: max |median paired delta| < 0.001 in the 0.27-0.40 band at
%     every f0 group, and >= MDC/2.77 in the 0.57-0.73 band at f0 = 2.5 Hz (exploratory: 0.0003 and
%     0.0006; 0.0129).
%  R4 deterministic gain: Block C, CEN, OLS and LMLS: max |dGain with noise - dGain at zero noise|
%     < NOISE_FRAC * MDC/2.77 (median-map gains from v001's MD table).
% Anchor (error if broken): the paired deltas rebuilt here from J equal v001's saved D exactly.
% Fail loud: a missing sigma = 0 companion is an error; non-finite deltas are counted and printed,
% not filled.
% Run from src/. Reads results/checkBWFDHalfSampleNoise_v001.mat; writes
% results/checkBWFDHalfSampleNoise_v002.mat

%% CONFIG
ROOT       = fileparts(fileparts(mfilename("fullpath")));
IN_MAT     = fullfile(ROOT, "results", "checkBWFDHalfSampleNoise_v001.mat");
OUT_MAT    = fullfile(ROOT, "results", "checkBWFDHalfSampleNoise_v002.mat");
NOISE_FRAC = 0.1;                                  % noise-only tolerance as a fraction of MDC/2.77
BAND_EDGES = [0 0.26 0.41 0.56 0.74];              % beta_gen grid is 0:1/30:0.7333
BAND_LAB   = ["0-0.23" "0.27-0.40" "0.43-0.53" "0.57-0.73"];
GRP_LAB    = ["f0 <= 1.5" "f0 = 2.5" "f0 = 4"];    % Block C f0 values: 0.25 0.5 1 1.5 2.5 4
BAND_LOW_MAX = 0.001;                              % R3, low band ceiling
VARIANTS   = ["REG" "CEN" "AVG"];
REG_NAMES  = ["OLS" "LMLS" "IRLS"];

%% Load
if ~isfile(IN_MAT), error("halfSampleNoise:in", "%s", "Missing input: " + IN_MAT); end
if ~isfolder(fileparts(OUT_MAT)), error("halfSampleNoise:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
S = load(IN_MAT);
for f = ["J" "D" "MD" "semAdequate"]
    if ~isfield(S, f), error("halfSampleNoise:field", "%s", "v001 results lack field " + f); end
end
J = S.J;  D = S.D;  MD = S.MD;  semAdequate = S.semAdequate;
thr = NOISE_FRAC * semAdequate;
fprintf("MDC/2.77 = %.5f; noise-only tolerance %.5f\n", semAdequate, thr);
col = @(v, q) (v - 1) * numel(REG_NAMES) + q;      % v001's betaRec column order: variant-major

%% Zero-noise companions
Z = J(J.sigma == 0 & ~J.lengthExcluded, :);
zk = arrayfun(@(i) key_local(Z.block(i), Z.dataset(i), Z.FS(i), Z.f0(i), Z.betaGen(i)), 1:height(Z))';
if numel(unique(zk)) ~= numel(zk), error("halfSampleNoise:dupCompanion", "Duplicate zero-noise companion keys."); end
zmap = containers.Map(cellstr(zk), num2cell(1:height(Z)));

%% Paired deltas, with and without noise (same loop order as v001's D)
Nz = J(J.sigma > 0 & ~J.lengthExcluded, :);
nR = height(Nz) * 2 * numel(REG_NAMES);
blk = strings(nR, 1);  dsn = strings(nR, 1);  vr = strings(nR, 1);  rg = strings(nR, 1);
FSv = zeros(nR, 1);  f0v = FSv;  bgv = FSv;  med = FSv;  d0 = FSv;  nP = FSv;
r = 0;
for k = 1:height(Nz)
    kk = key_local(Nz.block(k), Nz.dataset(k), Nz.FS(k), Nz.f0(k), Nz.betaGen(k));
    if ~isKey(zmap, char(kk)), error("halfSampleNoise:companion", "%s", "No sigma = 0 companion for " + kk); end
    b0 = Z.betaRec{zmap(char(kk))}(1, :);
    br = Nz.betaRec{k};
    for v = 2:numel(VARIANTS)
        for q = 1:numel(REG_NAMES)
            r = r + 1;
            d = br(:, col(v, q)) - br(:, col(1, q));
            blk(r) = Nz.block(k);  dsn(r) = Nz.dataset(k);  FSv(r) = Nz.FS(k);  f0v(r) = Nz.f0(k);
            bgv(r) = Nz.betaGen(k);  vr(r) = VARIANTS(v);  rg(r) = REG_NAMES(q);
            med(r) = median(d, "omitnan");  nP(r) = sum(isfinite(d));
            d0(r) = b0(col(v, q)) - b0(col(1, q));
        end
    end
end
T = table(blk, dsn, FSv, f0v, bgv, vr, rg, med, d0, nP, 'VariableNames', ...
    ["block" "dataset" "FS" "f0" "betaGen" "variant" "regression" "medDelta" "delta0" "nPaired"]);
T.noiseOnly = T.medDelta - T.delta0;

%% Anchor: rebuilt deltas equal v001's D
if height(D) ~= height(T), error("halfSampleNoise:anchor", "Row count %d against v001's D %d.", height(T), height(D)); end
if ~isequal(isnan(D.medDelta), isnan(T.medDelta)) || max(abs(D.medDelta - T.medDelta), [], "omitnan") > 1e-15 ...
        || ~isequal(D.variant, T.variant) || ~isequal(D.regression, T.regression)
    error("halfSampleNoise:anchor", "Rebuilt paired deltas differ from v001's D.");
end
fprintf("Anchor passed: %d paired deltas rebuilt from J equal v001's D exactly\n", height(T));
fprintf("Non-finite median paired deltas by regression: %s; zero-noise deltas: %s\n", ...
    strjoin(REG_NAMES + " " + arrayfun(@(q) sum(~isfinite(T.medDelta(T.regression == REG_NAMES(q)))), 1:3), ", "), ...
    strjoin(REG_NAMES + " " + arrayfun(@(q) sum(~isfinite(T.delta0(T.regression == REG_NAMES(q)))), 1:3), ", "));

%% (1) Noise-only part
for bl = ["C" "S"]
    Tb = T(T.block == bl, :);
    [g, gv, gr] = findgroups(Tb.variant, Tb.regression);
    O = table(gv, gr, splitapply(@(x) max(abs(x), [], "omitnan"), Tb.medDelta, g), ...
        splitapply(@(x) max(abs(x), [], "omitnan"), Tb.delta0, g), ...
        splitapply(@(x) max(abs(x), [], "omitnan"), Tb.noiseOnly, g), 'VariableNames', ...
        ["variant" "regression" "maxAbsMedDelta" "maxAbsZeroNoise" "maxAbsNoiseOnly"]);
    fprintf("\nBlock %s: largest |median paired delta|, its zero-noise part and its noise-only part:\n", bl);
    disp(O)
    if bl == "C", OC = O; else, OS = O; end
end

%% Cells reaching MDC/2.77 (v001's N1)
TC = T(T.block == "C", :);
hit = TC(TC.variant == "CEN" & abs(TC.medDelta) >= semAdequate, :);
fprintf("\nBlock C, CEN, cells with |median paired delta| >= MDC/2.77: %d\n", height(hit));
disp(hit(:, ["dataset" "f0" "betaGen" "regression" "medDelta" "delta0" "noiseOnly"]))

%% (2) beta_gen-band table
TB = TC(TC.variant == "CEN", :);
band = discretize(TB.betaGen, BAND_EDGES);
grp = NaN(height(TB), 1);  grp(TB.f0 <= 1.5) = 1;  grp(TB.f0 == 2.5) = 2;  grp(TB.f0 == 4) = 3;
if any(isnan(band)) || any(isnan(grp)), error("halfSampleNoise:band", "A Block C row falls outside the band or f0 groups."); end
[g, gr, gb, gg] = findgroups(TB.regression, band, grp);
BT = table(gr, BAND_LAB(gb)', GRP_LAB(gg)', splitapply(@(x) max(abs(x), [], "omitnan"), TB.medDelta, g), ...
    splitapply(@(x) max(abs(x), [], "omitnan"), TB.noiseOnly, g), splitapply(@numel, TB.medDelta, g), 'VariableNames', ...
    ["regression" "betaGenBand" "f0Group" "maxAbsMedDelta" "maxAbsNoiseOnly" "nCells"]);
BT = sortrows(BT, ["regression" "betaGenBand" "f0Group"]);
fprintf("\nBlock C, CEN against REG, by beta_gen band and f0 group:\n");
disp(BT)

%% (4) Gain change with and without noise (median maps)
keys = ["block" "dataset" "FS" "f0" "variant" "regression"];
MDn = MD(MD.sigma > 0, :);
MDz = [MD(MD.sigma == 0, cellstr(keys)), MD(MD.sigma == 0, "dGain")];
MDz = renamevars(MDz, "dGain", "dGain0");
GJ = innerjoin(MDn(:, [cellstr(keys), {'dGain'}]), MDz, 'Keys', cellstr(keys));
if height(GJ) ~= height(MDn), error("halfSampleNoise:joinGain", "Gain join lost rows: %d of %d.", height(GJ), height(MDn)); end
GJ.dNoise = GJ.dGain - GJ.dGain0;
GC = GJ(GJ.block == "C", :);
[g, gv, gr] = findgroups(GC.variant, GC.regression);
GS = table(gv, gr, splitapply(@(x) max(abs(x), [], "omitnan"), GC.dNoise, g), 'VariableNames', ["variant" "regression" "maxAbsGainChangeFromNoise"]);
fprintf("\nBlock C: |dGain with noise - dGain at zero noise|, largest per variant and regression:\n");
disp(GS)
ex = GC(GC.dataset == "Cook_CTRL" & GC.f0 == 2.5 & GC.variant == "CEN", ["regression" "dGain" "dGain0" "dNoise"]);
fprintf("Cook CTRL, 2.5 Hz, CEN:\n");  disp(ex)

%% Predictions
c1 = OC(OC.variant == "CEN" & ismember(OC.regression, ["OLS" "LMLS"]), :);
verdict_local("R1 noise-only part small", all(c1.maxAbsNoiseOnly < thr), sprintf("CEN Block C max |noise-only|: OLS %.5f, LMLS %.5f, IRLS %.5f (not gated); tolerance %.5f", ...
    c1.maxAbsNoiseOnly(c1.regression == "OLS"), c1.maxAbsNoiseOnly(c1.regression == "LMLS"), OC.maxAbsNoiseOnly(OC.variant == "CEN" & OC.regression == "IRLS"), thr));
if isempty(hit)
    fprintf("R2 cells at MDC/2.77: not evaluable (no cell reaches it)\n");
    R2 = NaN;
else
    R2 = all(abs(hit.noiseOnly) < thr);
    verdict_local("R2 cells at MDC/2.77", R2, sprintf("%d cells, max |noise-only| %.5f (tolerance %.5f)", height(hit), max(abs(hit.noiseOnly)), thr));
end
bo = BT(BT.regression == "OLS", :);
lowOK = all(bo.maxAbsMedDelta(bo.betaGenBand == "0.27-0.40") < BAND_LOW_MAX);
hiRow = bo(bo.betaGenBand == "0.57-0.73" & bo.f0Group == "f0 = 2.5", :);
R3 = lowOK && height(hiRow) == 1 && hiRow.maxAbsMedDelta >= semAdequate;
verdict_local("R3 beta_gen band", R3, sprintf("0.27-0.40 band max %.5f (ceiling %.3f); 0.57-0.73 band at 2.5 Hz %.5f (MDC/2.77 %.5f)", ...
    max(bo.maxAbsMedDelta(bo.betaGenBand == "0.27-0.40")), BAND_LOW_MAX, hiRow.maxAbsMedDelta, semAdequate));
c4 = GS(GS.variant == "CEN" & ismember(GS.regression, ["OLS" "LMLS"]), :);
verdict_local("R4 deterministic gain", all(c4.maxAbsGainChangeFromNoise < thr), sprintf("CEN Block C max gain change from noise: OLS %.5f, LMLS %.5f (tolerance %.5f)", ...
    c4.maxAbsGainChangeFromNoise(c4.regression == "OLS"), c4.maxAbsGainChangeFromNoise(c4.regression == "LMLS"), thr));
P = struct("R1", all(c1.maxAbsNoiseOnly < thr), "R2", R2, "R3", R3, "R4", all(c4.maxAbsGainChangeFromNoise < thr));

save(OUT_MAT, "T", "OC", "OS", "hit", "BT", "GS", "GJ", "P", "thr", "semAdequate", "-v7.3");
fprintf("Saved: %s\n", OUT_MAT);

%% =========================================================================
function k = key_local(blk, ds, fs, f0, bg)
k = string(blk) + "|" + string(ds) + "|" + string(sprintf("%d|%.4f|%.6f", fs, f0, bg));
end

function verdict_local(tag, held, detail)
s = "NOT HELD";  if held, s = "HELD"; end
fprintf("%s %s: %s\n", tag, s, detail);
end
