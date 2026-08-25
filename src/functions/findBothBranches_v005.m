function both = findBothBranches_v005(xGrid, yRaw, sigmaY, params)
%FINDBOTHBRANCHES_V005 Classify a 1D forward-map curve into ALL invertible
% monotonic branches -- grid-invariant segment width and smoothing (v004)
% PLUS a noise-calibrated slope threshold (new in v005), replacing the
% flat params.slopeTol constant used in v001-v004.
%
% RATIONALE (2026-08-11, docs/TODO_NoiseCalibratedInvertibilityGate_v001.md):
% v002/v003/v004 all used a single flat slopeTol=0.05 everywhere,
% regardless of how noisy the curve actually is at a given point. This
% project already has a noise-calibrated adequacy threshold sitting
% unused for this purpose: prereg_v101 Section 3.2's SEM<0.011 criterion
% for a recovered beta_gen to be adequately measured. Error propagation
% through a locally-linear inverse gives sigma_x ~= sigma_y/|slope|;
% requiring sigma_x < 0.011 rearranges to |slope| > sigma_y/0.011. This
% makes the bar a segment must clear ADAPT to the actual local noise
% level (estimated directly from the N_REPS replicate draws already
% generated at each beta_gen point), rather than a single constant
% applied everywhere regardless of context.
%
% API CHANGE vs v001-v004: sigmaY (a vector, one value per xGrid point --
% the noise SD in yRaw at that point, typically std across N_REPS
% replicate draws) is now a REQUIRED input. params.slopeTol is no longer
% used for this purpose; params.mdcThreshold (default 0.011, the
% project's own existing MDC/2.77 adequacy criterion) replaces it as the
% tunable constant.
%
% For the finite-difference slope between adjacent points i, i+1, the
% effective local noise used is mean([sigmaY(i), sigmaY(i+1)]) -- the
% average of the two points' own noise estimates, since the slope itself
% is a comparison of exactly those two values.
%
% v004 is left completely unmodified (new file, matching this project's
% version-lineage-immutability convention throughout this thread).
%
% HONEST EXPECTATION, STATED BEFORE TESTING: Finding #151 found the
% residual sensitivity remaining after v004's fix (16/36 cells, 44.4%) is
% UNIFORMLY one direction ([0,1,1] -- coarser grids always show MORE
% coverage), not a mixed pattern. That consistency is more suggestive of
% structural aliasing (the coarse grid's own sample points do not include
% where a real reversal occurs -- no threshold recalibration can reveal a
% feature that was never sampled) than of noise-threshold miscalibration.
% This version is built and tested regardless, since the alternative is
% assuming rather than checking -- but do not be surprised if the 44.4%
% figure does not drop substantially further.
%
% See also: findBothBranches_v004 (base this extends; flat slopeTol still
% used there)
%
% Fraser, D.S. (2026)

both.rise = findMonotonicRuns_v005(xGrid, yRaw, sigmaY, params, "rising");
both.desc = findMonotonicRuns_v005(xGrid, yRaw, sigmaY, params, "descending");

end

% --------------------------------------------------------------------------
function segs = findMonotonicRuns_v005(xGrid, yRaw, sigmaY, params, direction)
    segTemplate = struct("x", [], "y", [], "n", 0, "meanSlope", NaN, "invertible", false);
    segs = segTemplate([]);

    good = ~isnan(yRaw) & ~isnan(sigmaY);
    if nnz(good) < 2, return, end
    x  = xGrid(good);
    y  = yRaw(good);
    sy = sigmaY(good);

    if numel(x) >= 2
        gridSpacing = median(diff(sort(x)));
    else
        gridSpacing = Inf;
    end
    if isfield(params, 'smoothWidth') && isfinite(gridSpacing) && gridSpacing > 0
        windowPts = max(1, round(params.smoothWidth / gridSpacing));
    else
        windowPts = 1;
    end

    if windowPts > 1 && numel(y) >= windowPts
        ys = movmean(y, windowPts);
    else
        ys = y;
    end

    if isfield(params, 'mdcThreshold')
        mdcThreshold = params.mdcThreshold;
    else
        mdcThreshold = 0.011;  % project's own MDC/2.77 adequacy criterion
    end

    slope        = diff(ys) ./ diff(x);
    sigmaEff     = (sy(1:end-1) + sy(2:end)) / 2;   % avg local noise per finite-difference step
    slopeThresh  = sigmaEff / mdcThreshold;         % adaptive, replaces flat slopeTol

    switch direction
        case "rising",     isValid = slope >  slopeThresh;
        case "descending", isValid = slope < -slopeThresh;
        otherwise
            error("findMonotonicRuns_v005:BadDirection", "%s", ...
                "direction must be 'rising' or 'descending', got: " + direction);
    end

    [starts, ends] = allTrueRuns(isValid);

    for k = 1:numel(starts)
        runStart = starts(k); runEnd = ends(k);

        if (x(runEnd+1) - x(runStart)) < params.minSegWidth, continue, end

        sx = x(runStart:runEnd+1);
        syv = y(runStart:runEnd+1);

        if direction == "rising"
            bad = find(diff(syv) <= 0, 1, "first");
        else
            bad = find(diff(syv) >= 0, 1, "first");
        end
        if ~isempty(bad)
            sx  = sx(1:bad);
            syv = syv(1:bad);
        end

        if numel(syv) < 2 || (max(sx) - min(sx)) < params.minSegWidth
            continue
        end

        if direction == "descending"
            sx  = flipud(sx(:));
            syv = flipud(syv(:));
        else
            sx  = sx(:);
            syv = syv(:);
        end

        meanSlope = mean(diff(syv) ./ diff(sx));
        % Segment-level check against the MEAN effective threshold across
        % the segment's own points (consistent with the per-step gating
        % already applied above; re-checked here since trimming can
        % change which points are included).
        segIdxInGood = ismember(x, sx);
        segSigma = mean(sy(segIdxInGood));
        segThresh = segSigma / mdcThreshold;
        if abs(meanSlope) <= segThresh, continue, end

        seg = struct("x", sx, "y", syv, "n", numel(syv), ...
            "meanSlope", meanSlope, "invertible", true);
        segs = [segs, seg]; %#ok<AGROW>
    end
end

% --------------------------------------------------------------------------
function [starts, ends] = allTrueRuns(v)
    d      = diff([false; v(:); false]);
    starts = find(d ==  1);
    ends   = find(d == -1) - 1;
end
