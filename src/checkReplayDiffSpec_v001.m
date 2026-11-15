%% checkReplayDiffSpec_v001  Deep Review GD17 item 8: size of the replay-vs-grid differences behind prediction P4 and the replay gate.
% Same inputs and join as analyseTempoFoldSweepE_v002 section 2 (results/tempoFoldSweepE_v001.mat block REPLAY against src/perCoordinateSEM_v2_001.mat
% at fs 120, VGF node E.VGF0), but keeps the node-level differences: median and maximum |engine - grid| in beta_rec, over (a) the noise nodes
% (alpha x sigma x pipeline x beta_gen; the "1,320 nodes" of the gate) and (b) the sigma = 0 nodes (P4 as written: "within 0.01 at every beta_gen node").
% Fails loud if the noise-node count is not 1,320 or the |z| > 3 count does not reproduce the reported 49.
% Run from src/. Writes results/checkReplayDiffSpec_v001.mat.   Fraser, D.S. (2026)  v001
ROOT = fileparts(fileparts(mfilename("fullpath")));
LABELS = ["BWFD-OLS" "BWFD-LMLS" "BWFD-IRLS" "SG-OLS" "SG-LMLS" "SG-IRLS"];  Z_GATE = 3;
fE = fullfile(ROOT, "results", "tempoFoldSweepE_v001.mat");  fS = fullfile(ROOT, "src", "perCoordinateSEM_v2_001.mat");
for f = [fE fS], if ~isfile(f), error("replayDiff:file", "%s", "FAILED PATH: " + f); end, end
E = load(fE);  R = E.R;  S = load(fS, "coordTable");  T = S.coordTable;
for v = ["betaGen" "VGF" "fs" "alpha" "sigma" "meanBetaRec" "sem" "nReps"], T.(v) = double(T.(v)); end
T = T(T.fs == E.FS & T.VGF == E.VGF0, :);
D = table();
for a = E.ALPHAS
    for s = [0 E.SIGMAS_MM]
        if s == 0 && a ~= E.ALPHAS(1), continue; end                         % sigma = 0 once, as in analyseTempoFoldSweepE_v002
        r = R(R.block == "REPLAY" & R.alpha == a & R.sigmaMM == s, :);
        for p = 1:numel(LABELS)
            gp = T(T.alpha == a & T.sigma == s & T.pipeline == LABELS(p), ["betaGen" "meanBetaRec" "sem" "nReps"]);
            rp = table(r.betaGen, r.betaMean(:, p), r.betaSD(:, p), r.nReps, 'VariableNames', ["betaGen" "eng" "engSD" "engN"]);
            m = innerjoin(rp, gp, "Keys", "betaGen");  m = m(isfinite(m.eng) & isfinite(m.meanBetaRec), :);
            v1 = m.engSD.^2 ./ m.engN;  v2 = m.sem.^2 ./ m.nReps;  z = (m.eng - m.meanBetaRec) ./ sqrt(v1 + v2);
            D = [D; table(repmat(s, height(m), 1), repmat(LABELS(p), height(m), 1), abs(m.eng - m.meanBetaRec), z, 'VariableNames', ["sigma" "pipeline" "absDiff" "z"])]; %#ok<AGROW>
    end, end
end
N = D(D.sigma > 0, :);  Z = D(D.sigma == 0, :);
R8 = struct();
R8.nNoise = height(N);  R8.noiseMedian = median(N.absDiff);  R8.noiseMax = max(N.absDiff);  R8.noiseOver01 = sum(N.absDiff > 0.01);  R8.noiseZover = sum(abs(N.z) > Z_GATE);
R8.nZero = height(Z);   R8.zeroMedian = median(Z.absDiff);   R8.zeroMax = max(Z.absDiff);   R8.zeroOver01 = sum(Z.absDiff > 0.01);
% (c) Block B2 as P4 is specified (spec 3: replay at the grid's own tempo path against the sigma = 0 grid), as in analyseTempoFoldSweep_v001 section 2
fA = fullfile(ROOT, "results", "tempoFoldSweep_v001.mat");  fG = fullfile(ROOT, "results", "gridTempoAtFold_v004.mat");
for f = [fA fG], if ~isfile(f), error("replayDiff:file", "%s", "FAILED PATH: " + f); end, end
RB = load(fA).R;  TG = load(fG).out.perCoordinate;  B = RB(RB.block == "B2", :);  FS_REF = 120;  B2 = table();
for p = 1:numel(LABELS)
    gp = TG(TG.sampling_rate == FS_REF & TG.pipeline == LABELS(p), ["generated_beta" "vgf_value" "betaRec"]);
    rp = table(B.betaGen, B.vgf, B.betaMean(:, p), 'VariableNames', ["generated_beta" "vgf_value" "eng"]);
    m = innerjoin(rp, gp, "Keys", ["generated_beta" "vgf_value"]);  B2 = [B2; table(repmat(LABELS(p), height(m), 1), abs(m.eng - m.betaRec), 'VariableNames', ["pipeline" "absDiff"])]; %#ok<AGROW>
end
B2 = B2(isfinite(B2.absDiff), :);
R8.nB2 = height(B2);  R8.b2Median = median(B2.absDiff);  R8.b2Max = max(B2.absDiff);  R8.b2Over01 = sum(B2.absDiff > 0.01);
R8.b2MaxPipeline = B2.pipeline(find(B2.absDiff == max(B2.absDiff), 1));
R8.b2PipeMedianMax = max(groupsummary(B2, "pipeline", "median", "absDiff").median_absDiff);
% (d) the cell-mean OLS refit of the registered VGF formula against the fitted model (checkLMMTempoConfound_v002 gate, M0)
fL = fullfile(ROOT, "results", "lmmTempoConfound_v002.mat");  if ~isfile(fL), error("replayDiff:file", "%s", "FAILED PATH: " + fL); end
LO = load(fL).out;  G = LO.gate;
assert(height(G) == 192 && all(G.pass), "replayDiff:gate", "The cell-mean gate does not hold for all 192 coefficients");
R8.gateN = height(G);  R8.gateMaxSE = max(G.zd);  R8.gateMaxAbs = max(abs(G.d));  R8.gateR = LO.gateR;
assert(R8.nNoise == 1320, "replayDiff:n", "Noise-node count is %d, not 1,320", R8.nNoise);
assert(R8.noiseZover == 49, "replayDiff:z", "|z| > 3 count is %d, not the reported 49", R8.noiseZover);
disp(R8);  save(fullfile(ROOT, "results", "checkReplayDiffSpec_v001.mat"), "-struct", "R8");
