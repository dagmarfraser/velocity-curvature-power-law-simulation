function [predBetaRec, predDeltaBeta, info] = predictedDeltaBeta_v001( ...
    alpha, sigmaMM, fs, betaTrue, vgfTrue, lookupGrid, gridSpec)
% predictedDeltaBeta_v001  Look up the predicted constellation at one
% empirical trial's full 6D coordinate.
%
% Snaps (alpha, sigmaMM, fs, betaTrue, vgfTrue) to nearest grid points in
% the simulation parameter space and reads the predicted beta_rec for each
% of the 6 canonical pipelines, then computes the 15 pairwise predicted
% Delta beta values in canonical order (i<j enumeration over 6 pipelines).
%
% No marginalisation: each empirical trial is matched to its single
% nearest-neighbour simulation cell.
%
% USAGE:
%   [predBetaRec, predDeltaBeta, info] = predictedDeltaBeta_v001( ...
%       alpha, sigmaMM, fs, betaTrue, vgfTrue, lookupGrid, gridSpec)
%
% INPUTS:
%   alpha       scalar empirical noise colour exponent (position space)
%   sigmaMM     scalar empirical noise magnitude (mm)
%   fs          scalar empirical sampling rate (Hz)
%   betaTrue    scalar best estimate of true beta for this trial
%   vgfTrue     scalar best estimate of true VGF for this trial
%   lookupGrid  6D double from buildPredictionLookup_v001
%               [nAlpha x nSigma x nFs x nBetaGen x nVGF x 6]
%   gridSpec    struct with .alphaGrid, .sigmaGrid, .fsGrid,
%               .betaGenGrid, .vgfGrid, .pipeOrder
%
% OUTPUTS:
%   predBetaRec    [6x1] predicted beta_rec per pipeline (canonical order)
%   predDeltaBeta  [15x1] predicted pairwise differences (i<j enumeration)
%   info           struct: .alphaSnap .sigmaSnap .fsSnap .betaGenSnap .vgfSnap
%                          .iA .iS .iF .iB .iV .warnings (string array)
%
% PREREG REF: Section 7.4 (constellation prediction generation)
% Fraser, D.S. (2026)

arguments
    alpha       (1,1) double
    sigmaMM     (1,1) double
    fs          (1,1) double
    betaTrue    (1,1) double
    vgfTrue     (1,1) double
    lookupGrid  double
    gridSpec    struct
end

warns = strings(0, 1);

[aSnap, iA, wA] = snapAndCheck(alpha,    gridSpec.alphaGrid,   "alpha");
[sSnap, iS, wS] = snapAndCheck(sigmaMM,  gridSpec.sigmaGrid,   "sigma_mm");
[fSnap, iF, wF] = snapAndCheck(fs,       gridSpec.fsGrid,      "fs");
[bSnap, iB, wB] = snapAndCheck(betaTrue, gridSpec.betaGenGrid, "betaTrue");
[vSnap, iV, wV] = snapAndCheck(vgfTrue,  gridSpec.vgfGrid,     "vgfTrue");

warns = [warns; wA; wS; wF; wB; wV];

if any(isnan([iA iS iF iB iV]))
    predBetaRec   = NaN(6, 1);
    predDeltaBeta = NaN(15, 1);
    info = struct("alphaSnap", aSnap, "sigmaSnap", sSnap, "fsSnap", fSnap, ...
                  "betaGenSnap", bSnap, "vgfSnap", vSnap, ...
                  "iA", iA, "iS", iS, "iF", iF, "iB", iB, "iV", iV, ...
                  "warnings", warns);
    return;
end

predBetaRec = squeeze(lookupGrid(iA, iS, iF, iB, iV, :));
if numel(predBetaRec) ~= 6
    error("predictedDeltaBeta:Shape", ...
        "Expected 6 pipelines from lookup, got %d.", numel(predBetaRec));
end

if any(isnan(predBetaRec))
    nanIdx = find(isnan(predBetaRec));
    warns(end+1, 1) = sprintf( ...
        "predBetaRec NaN for pipelines [%s] at (a=%.2f, s=%.3g, fs=%d, b=%.3f, v=%.1f).", ...
        strjoin(gridSpec.pipeOrder(nanIdx), ","), aSnap, sSnap, fSnap, bSnap, vSnap);
end

% 15 pairwise diffs (i<j over 6 pipelines)
predDeltaBeta = NaN(15, 1);
k = 0;
for p = 1:6
    for q = (p+1):6
        k = k + 1;
        predDeltaBeta(k) = predBetaRec(p) - predBetaRec(q);
    end
end

info = struct( ...
    "alphaSnap",   aSnap,   "sigmaSnap",   sSnap, ...
    "fsSnap",      fSnap,   "betaGenSnap", bSnap, ...
    "vgfSnap",     vSnap, ...
    "iA", iA, "iS", iS, "iF", iF, "iB", iB, "iV", iV, ...
    "warnings", warns);
end

% --------------------------------------------------------------------------
function [snap, idx, warns] = snapAndCheck(val, grid, label)
% Snap val to nearest grid value. Returns (snap, linearIndex, warnings).
% Out-of-range inputs are clamped to the nearest endpoint and a warning is
% returned (visible-fallback per fail-loud convention).
    warns = strings(0, 1);
    if isnan(val)
        snap = NaN; idx = NaN;
        warns(end+1, 1) = sprintf("%s is NaN.", label);
        return;
    end
    [delta, idx] = min(abs(grid - val));
    snap = grid(idx);
    if val < min(grid) || val > max(grid)
        warns(end+1, 1) = sprintf( ...
            "%s=%.4g outside grid [%.4g, %.4g]; clamped to %.4g (delta=%.4g).", ...
            label, val, min(grid), max(grid), snap, delta);
    end
end
