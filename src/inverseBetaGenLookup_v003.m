function out = inverseBetaGenLookup_v003(p, alphaT, sigmaT, fsT, f0K, betaRec, lookupBundle)
% inverseBetaGenLookup_v003  PCHIP inversion using f0-indexed lookup.
%
% Replaces inverseBetaGenLookup_v002 (VGF-indexed with per-candidate G_sim
% recomputation, 2026-05-18). The VGF axis has been replaced by f0 in
% predictionLookup_v2_002, so:
%   - No gSimBundle argument required.
%   - f0K is snapped once to the nearest f0Grid node.
%   - The composite beta_rec curve is read directly as
%       lookupGrid[aIdx, sIdx, fIdx, :, f0Idx, p]
%     where VGF = f0 x G_sim(beta_gen) is already encoded in the grid by
%     construction (buildPredictionLookup_v2_002).
%   - Self-consistency is guaranteed structurally, not numerically.
%
% INPUTS:
%   p            pipeline index (integer, 1-6)
%   alphaT       empirical noise alpha (snapped to lookupBundle.alphaGrid)
%   sigmaT       empirical noise sigma mm (snapped to lookupBundle.sigmaGrid)
%   fsT          empirical sampling rate Hz (snapped to lookupBundle.fsGrid)
%   f0K          empirical fundamental frequency Hz (snapped to f0Grid)
%   betaRec      observed beta_rec for pipeline p on this trial
%   lookupBundle struct from predictionLookup_v2_002.mat:
%                  .lookupGrid [nA x nS x nF x nB x nF0 x 6]
%                  .alphaGrid, .sigmaGrid, .fsGrid, .betaGenGrid, .f0Grid
%
% OUTPUT struct (same field names as v001/v002 for drop-in compatibility):
%   .betaGenHat          double   NaN if not invertible
%   .status              string   "INVERTIBLE"|"NO_RISE_AT_COORD"|
%                                 "BELOW_RISE"|"ABOVE_RISE"
%   .branch              string   "rise" | "none"
%   .meanSlope           double   NaN if not invertible
%   .descInvert          logical  false
%   .betaRecInDescImage  logical  false
%
% PREREG REF:
%   claude.md Session 18 (VGF grid coverage failure; f0-axis fix).
%   Caller: constellationCCC_v005.
%
% Fraser, D.S. (2026)

arguments
    p            (1,1) double {mustBeInteger, mustBePositive}
    alphaT       (1,1) double
    sigmaT       (1,1) double
    fsT          (1,1) double
    f0K          (1,1) double
    betaRec      (1,1) double
    lookupBundle struct
end

%% Initialise output
out.betaGenHat         = NaN;
out.status             = "NO_RISE_AT_COORD";
out.branch             = "none";
out.meanSlope          = NaN;
out.descInvert         = false;
out.betaRecInDescImage = false;

%% Guards
if isnan(betaRec) || isnan(f0K) || f0K <= 0
    return
end

%% Validate f0Grid present
if ~isfield(lookupBundle, "f0Grid")
    error("inverseBetaGenLookup_v003:NoF0Grid", "%s", ...
        "lookupBundle has no 'f0Grid'. Load predictionLookup_v2_002.mat, not v2_001.");
end

%% Snap all coordinates
[~, aIdx]  = min(abs(lookupBundle.alphaGrid   - alphaT));
[~, sIdx]  = min(abs(lookupBundle.sigmaGrid   - sigmaT));
[~, fIdx]  = min(abs(lookupBundle.fsGrid      - fsT));
[~, f0Idx] = min(abs(log(lookupBundle.f0Grid) - log(f0K)));   % log-distance snap

%% Extract composite beta_rec curve across betaGen at snapped coords
betaGenGrid = double(lookupBundle.betaGenGrid(:)');
nB = numel(betaGenGrid);

betaRecComposite = double(squeeze( ...
    lookupBundle.lookupGrid(aIdx, sIdx, fIdx, :, f0Idx, p)))';   % 1 x nB

%% Find rising segment
valid = find(isfinite(betaRecComposite));
if numel(valid) < 3
    return
end

bG = betaGenGrid(valid);
bR = betaRecComposite(valid);

[~, peakLocal] = max(bR);
if peakLocal < 2
    return
end

segBG = bG(1:peakLocal);
segBR = bR(1:peakLocal);

out.meanSlope = (segBR(end) - segBR(1)) / max(segBG(end) - segBG(1), 1e-8);

%% Range check
if betaRec < segBR(1)
    out.status = "BELOW_RISE";
    return
end
if betaRec > segBR(end)
    out.status = "ABOVE_RISE";
    return
end

%% PCHIP inversion
[segBR_u, ia] = unique(segBR, "stable");
segBG_u = segBG(ia);
if numel(segBR_u) < 2
    return
end

betaGenHat = interp1(segBR_u, segBG_u, betaRec, "pchip", NaN);

if isnan(betaGenHat)
    error("inverseBetaGenLookup_v003:InterpFailed", "%s", ...
        sprintf("interp1 returned NaN for betaRec=%.4f in [%.4f,%.4f] p=%d f0=%.3f.", ...
        betaRec, segBR(1), segBR(end), p, f0K));
end

out.betaGenHat = betaGenHat;
out.status     = "INVERTIBLE";
out.branch     = "rise";
end
