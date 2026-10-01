%% checkIRLSGap_v001.m
% Do the non-finite IRLS fits in predictHarmonicGain_v003 move the IRLS gain onsets?
% v003 saved 99 non-finite IRLS betas (limitBreak = 0 fails loud; none for OLS, no error flags),
% in conditions 4, 5, 6, 9 and 10. Each gain is a straight-line slope over the 10 beta_gen points
% of a curve (0.2-0.5, polyfit on the finite points, at least 3 needed), so a gap removes real
% information, and the gain is silently a fit to fewer points. #249 flags the IRLS onsets as resting
% on curves with gaps; this script measures how much that can matter.
% Three views of each gapped condition's IRLS onset (gain 0.9, log-linear crossing as v003):
%   as is        v003's gains, fitted to the finite points;
%   adjusted     proxy: each gapped curve's IRLS gain plus (OLS gain on all points minus OLS gain on
%                the same finite subset). OLS has no gaps, so this estimates what dropping those
%                beta_gen points does to a gain. A proxy, not an imputation: nothing is written back;
%   complete     gapped curves left out, the crossing interpolated between complete curves.
% Context: single-point leave-one-out influence on complete IRLS curves (largest change in gain when
% one beta_gen point is dropped), per condition.
% Prediction, fixed before the run (2026-10-01):
%  G1 for every (block, condition) with a gapped IRLS curve, the adjusted and the complete-curves
%     onsets each differ from the as-is onset by less than TOL_ONSET (0.5%, the Q1 theory tolerance
%     of #249). NaN onsets are counted as failures of G1, not skipped.
% Anchor (error if broken): gains rebuilt here from J equal v003's saved GF (OLS and IRLS, both
% windows), and v003's LOCAL_WIN equals this script's.
% Fail loud: the script prints how many gains are NaN (local window) and onsets that cannot be
% computed; nothing is filled.
% Run from src/. Reads results/predictHarmonicGain_v003.mat and results/filterCutoffCollapse_v001.mat;
% writes results/checkIRLSGap_v001.mat

%% CONFIG
ROOT      = fileparts(fileparts(mfilename("fullpath")));
IN_V003   = fullfile(ROOT, "results", "predictHarmonicGain_v003.mat");
IN_233    = fullfile(ROOT, "results", "filterCutoffCollapse_v001.mat");
OUT_MAT   = fullfile(ROOT, "results", "checkIRLSGap_v001.mat");
LOCAL_WIN = [0.29 0.37];                           % as predictHarmonicGain_v003 L41
TOL_ONSET = 0.005;
TOL_GAIN  = 1e-12;
NEAR      = 0.10;                                  % a gapped curve is "near" the crossing within +-10% in f0

%% Load
for f = [IN_V003 IN_233]
    if ~isfile(f), error("irlsGap:in", "%s", "Missing input: " + f); end
end
if ~isfolder(fileparts(OUT_MAT)), error("irlsGap:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
S = load(IN_V003, "J", "GF", "LOCAL_WIN");
R233 = load(IN_233, "GAIN_WIN", "GAIN_REF");
J = S.J;  GF = S.GF;  GAIN_WIN = R233.GAIN_WIN;  GAIN_REF = R233.GAIN_REF;
if ~isequal(S.LOCAL_WIN, LOCAL_WIN), error("irlsGap:localWin", "v003's LOCAL_WIN differs from this script's."); end
fprintf("Gain window %s, local window %s, reference gain %.2f; %d curves, %d non-finite IRLS betas (OLS %d)\n", ...
    mat2str(GAIN_WIN), mat2str(LOCAL_WIN), GAIN_REF, height(GF), sum(~isfinite(J.bIRLS)), sum(~isfinite(J.bOLS)));

%% Per-curve gains, rebuilt, with the OLS-on-the-same-subset companion
[grp, K] = findgroups(J(:, ["block" "cond" "filt" "ab" "f0"]));
nC = height(K);
gO = NaN(nC, 2);  gI = NaN(nC, 2);  gOsub = NaN(nC, 2);  nGap = zeros(nC, 2);  loo = NaN(nC, 1);
WINS = {GAIN_WIN, LOCAL_WIN};
for c = 1:nC
    m = grp == c;  b = J.betaGen(m);  oy = J.bOLS(m);  iy = J.bIRLS(m);  fi = isfinite(iy);
    for w = 1:2
        inW = b >= WINS{w}(1) & b <= WINS{w}(2);
        nGap(c, w) = sum(inW & ~fi);
        gO(c, w) = gainAt_local(b, oy, WINS{w});
        gI(c, w) = gainAt_local(b, iy, WINS{w});
        gOsub(c, w) = gainAt_local(b(fi), oy(fi), WINS{w});
    end
    if all(fi)                                      % single-point influence on a complete curve (global window)
        inW = find(b >= GAIN_WIN(1) & b <= GAIN_WIN(2));  d = zeros(numel(inW), 1);
        for q = 1:numel(inW)
            keep = true(size(b));  keep(inW(q)) = false;
            d(q) = gainAt_local(b(keep), iy(keep), GAIN_WIN) - gI(c, 1);
        end
        loo(c) = max(abs(d));
    end
end

%% Anchor: rebuilt gains equal v003's GF
if ~isequal(K.block, GF.block) || ~isequal(K.cond, GF.cond) || ~isequal(K.f0, GF.f0)
    error("irlsGap:anchor", "Curve keys differ from v003's GF.");
end
chk = {gO(:, 1), GF.gainOLS, "OLS"; gO(:, 2), GF.gainOLSloc, "OLS local"; gI(:, 1), GF.gainIRLS, "IRLS"; gI(:, 2), GF.gainIRLSloc, "IRLS local"};
for q = 1:size(chk, 1)
    if ~isequal(isnan(chk{q, 1}), isnan(chk{q, 2})) || max(abs(chk{q, 1} - chk{q, 2}), [], "omitnan") > TOL_GAIN
        error("irlsGap:anchor", "%s", "Rebuilt " + chk{q, 3} + " gains differ from v003's GF.");
    end
end
fprintf("Anchor passed: %d curves, OLS and IRLS gains in both windows equal v003's GF\n", nC);
fprintf("Local-window gains that are NaN: OLS %d, IRLS %d (these curves have fewer than 3 finite points in 0.29-0.37)\n", ...
    sum(isnan(gO(:, 2))), sum(isnan(gI(:, 2))));

%% Onsets under the three views, per (block, condition) with a gapped IRLS curve
useLoc = K.block == "emulC";                        % v003 L181-183: emulC onsets use the local window
w = 1 + useLoc;
idxSel = sub2ind(size(gI), (1:nC)', w);
gSel = gI(idxSel);  adj = gSel + (gO(idxSel) - gOsub(idxSel));  ng = nGap(idxSel);
adj(ng == 0) = gSel(ng == 0);                       % complete curves unchanged
BC = unique(K(:, ["block" "cond"]), "rows");
rows = table();
for i = 1:height(BC)
    ix = find(K.block == BC.block(i) & K.cond == BC.cond(i));
    [~, o] = sort(K.f0(ix));  ix = ix(o);
    if ~any(ng(ix) > 0), continue, end
    f0 = K.f0(ix);  gs = gSel(ix);  ga = adj(ix);  cmp = ng(ix) == 0;
    on0 = crossing_local(f0, gs, GAIN_REF);
    onA = crossing_local(f0, ga, GAIN_REF);
    onC = crossing_local(f0(cmp), gs(cmp), GAIN_REF);
    near = NaN;  if isfinite(on0), near = sum(ng(ix) > 0 & abs(log(f0 ./ on0)) < log(1 + NEAR)); end
    rows = [rows; table(BC.block(i), BC.cond(i), K.filt(ix(1)), K.ab(ix(1)), numel(ix), sum(ng(ix) > 0), sum(ng(ix)), near, on0, onA, onC, ...
        (onA - on0) / on0, (onC - on0) / on0, max(loo(ix), [], "omitnan"), max(abs(ga - gs), [], "omitnan"), ...
        'VariableNames', ["block" "cond" "filt" "ab" "nCurves" "nGapped" "nGapPoints" "nGappedNearOnset" "onsetAsIs" "onsetAdjusted" ...
        "onsetComplete" "relAdjusted" "relComplete" "maxLOOGain" "maxProxyGainShift"])]; %#ok<AGROW>
end
if isempty(rows), error("irlsGap:none", "No gapped IRLS curve found; v003's saved results changed?"); end
fprintf("\nIRLS onset (Hz) for each gapped condition: as is, adjusted (OLS proxy), complete curves only:\n");
disp(rows)

%% Prediction
bad = ~(abs(rows.relAdjusted) < TOL_ONSET & abs(rows.relComplete) < TOL_ONSET);   % NaN counts as failing
if ~any(bad)
    verdict_local("G1 IRLS onsets robust to the gaps", true, sprintf("%d gapped conditions; largest |relative change| %.2f%% (adjusted), %.2f%% (complete); tolerance %.1f%%", ...
        height(rows), 100 * max(abs(rows.relAdjusted)), 100 * max(abs(rows.relComplete)), 100 * TOL_ONSET));
else
    verdict_local("G1 IRLS onsets robust to the gaps", false, sprintf("%d of %d gapped conditions exceed %.1f%% or cannot be computed: %s", ...
        sum(bad), height(rows), 100 * TOL_ONSET, strjoin(rows.block(bad) + " c" + rows.cond(bad), ", ")));
end
fprintf("Largest single-point leave-one-out change in an IRLS gain on a complete curve, per gapped condition: %.4f to %.4f\n", ...
    min(rows.maxLOOGain), max(rows.maxLOOGain));

save(OUT_MAT, "rows", "K", "gI", "gO", "gOsub", "nGap", "loo", "TOL_ONSET", "GAIN_WIN", "LOCAL_WIN", "GAIN_REF", "-v7.3");
fprintf("Saved: %s\n", OUT_MAT);

%% =========================================================================
function g = gainAt_local(betaGen, bRec, win)          % as predictHarmonicGain_v003 L350-353
w = betaGen >= win(1) & betaGen <= win(2) & isfinite(bRec);
g = NaN;  if sum(w) >= 3, p = polyfit(betaGen(w), bRec(w), 1); g = p(1); end
end

function f = crossing_local(f0, g, ref)                % as predictHarmonicGain_v003 L406-412
f = NaN;
i = find(g < ref, 1);
if isempty(i) || i == 1, return, end
t = (ref - g(i-1)) / (g(i) - g(i-1));
f = exp(log(f0(i-1)) + t * (log(f0(i)) - log(f0(i-1))));
end

function verdict_local(tag, held, detail)
s = "NOT HELD";  if held, s = "HELD"; end
fprintf("%s %s: %s\n", tag, s, detail);
end
