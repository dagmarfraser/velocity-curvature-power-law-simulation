function out = checkTrialInvertibility_v001(alphaT, sigmaT, fsT, vgfT, betaRecPerPipe, monoBundle)
% checkTrialInvertibility_v001  Per-trial branch membership for one trial.
%
% Snaps (alpha, sigma, fs, VGF) to the simulation grid, then asks for each
% of the 6 pipelines: does this trial's recovered beta_rec sit inside the
% rising-branch image at this coord? The descending-branch image? Both?
% Neither? Returns a flat struct with per-pipeline logical arrays plus the
% snapped grid indices and any warning text accrued during snap.
%
% INPUTS:
%   alphaT, sigmaT, fsT, vgfT  scalar trial coordinates in canonical units
%                              (sigma in mm, fs in Hz, VGF in mm-equivalent
%                              after sigmaToMM and beta-aware unit conversion)
%   betaRecPerPipe             [1 x 6] beta_rec per pipeline, ordered as
%                              monoBundle.pipeOrder
%                              ["BWFD-OLS","BWFD-LMLS","BWFD-IRLS",
%                               "SG-OLS","SG-LMLS","SG-IRLS"]
%   monoBundle                 struct from monotonicSegments_v001.mat
%                              required fields: alphaGrid, sigmaGrid, fsGrid,
%                              VGFGrid, betaSeg (with .rise, .desc, .ambiguous),
%                              pipeOrder
%
% OUTPUT struct fields:
%   .inRise        [1 x 6] logical, true if beta_rec in rising image AND coord rise.invertible
%   .inDesc        [1 x 6] logical, mirror for descending
%   .inAmbiguous   [1 x 6] logical, inRise & inDesc & coord ambiguous
%   .inNeither     [1 x 6] logical, ~inRise & ~inDesc
%   .aIdx, .sIdx, .fIdx, .vIdx   snapped grid indices (scalars)
%   .snapWarnings                string scalar (joined "|" if multiple)
%   .anyBetaNaN                  logical, true if any betaRecPerPipe entry NaN
%
% Fraser, D.S. (2026)

nP = numel(monoBundle.pipeOrder);

out.inRise       = false(1, nP);
out.inDesc       = false(1, nP);
out.inAmbiguous  = false(1, nP);
out.inNeither    = false(1, nP);
out.aIdx         = NaN;
out.sIdx         = NaN;
out.fIdx         = NaN;
out.vIdx         = NaN;
out.snapWarnings = "";
out.anyBetaNaN   = any(isnan(betaRecPerPipe));

% --- Coord NaN guard ---
coords = [alphaT, sigmaT, fsT, vgfT];
if any(isnan(coords))
    out.snapWarnings = "Trial coord NaN (alpha|sigma|fs|VGF).";
    return
end

% --- Snap each axis ---
warns = strings(0, 1);
[aIdx, w] = nearestIdx(alphaT, monoBundle.alphaGrid, "alpha");
if ~isempty(w), warns(end+1, 1) = w; end
[sIdx, w] = nearestIdx(sigmaT, monoBundle.sigmaGrid, "sigma");
if ~isempty(w), warns(end+1, 1) = w; end
[fIdx, w] = nearestIdx(fsT,    monoBundle.fsGrid,    "fs");
if ~isempty(w), warns(end+1, 1) = w; end
[vIdx, w] = nearestIdx(vgfT,   monoBundle.VGFGrid,   "VGF");
if ~isempty(w), warns(end+1, 1) = w; end

out.aIdx = aIdx; out.sIdx = sIdx; out.fIdx = fIdx; out.vIdx = vIdx;
if ~isempty(warns), out.snapWarnings = strjoin(warns, " | "); end

% --- Per-pipeline branch test ---
% betaSeg.rise.segImage indexed [P, A, S, F, V, 2]; (:,:,:,:,:,1)=y_min, 2=y_max
for p = 1:nP
    riseInv = monoBundle.betaSeg.rise.invertible(p, aIdx, sIdx, fIdx, vIdx);
    descInv = monoBundle.betaSeg.desc.invertible(p, aIdx, sIdx, fIdx, vIdx);
    ambig   = monoBundle.betaSeg.ambiguous(p, aIdx, sIdx, fIdx, vIdx);
    bRec    = betaRecPerPipe(p);

    if isnan(bRec)
        % leave all flags false; flagged by anyBetaNaN at trial level
        continue
    end

    if riseInv
        riseLo = monoBundle.betaSeg.rise.segImage(p, aIdx, sIdx, fIdx, vIdx, 1);
        riseHi = monoBundle.betaSeg.rise.segImage(p, aIdx, sIdx, fIdx, vIdx, 2);
        out.inRise(p) = bRec >= riseLo && bRec <= riseHi;
    end

    if descInv
        descLo = monoBundle.betaSeg.desc.segImage(p, aIdx, sIdx, fIdx, vIdx, 1);
        descHi = monoBundle.betaSeg.desc.segImage(p, aIdx, sIdx, fIdx, vIdx, 2);
        out.inDesc(p) = bRec >= descLo && bRec <= descHi;
    end

    out.inAmbiguous(p) = out.inRise(p) && out.inDesc(p) && ambig;
    out.inNeither(p)   = ~out.inRise(p) && ~out.inDesc(p);
end
end


% =============================================================
function [idx, warnStr] = nearestIdx(v, grid, name)
    % Nearest-grid snap with out-of-range warning. Grid assumed sorted.
    [d, idx] = min(abs(grid - v));
    warnStr = "";
    if v < grid(1) || v > grid(end)
        warnStr = sprintf("%s=%.4g outside grid [%.4g, %.4g]; snapped to %.4g (d=%.4g)", ...
            name, v, grid(1), grid(end), grid(idx), d);
    end
end
