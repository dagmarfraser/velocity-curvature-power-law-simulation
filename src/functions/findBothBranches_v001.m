function both = findBothBranches_v001(xGrid, yRaw, params)
%FINDBOTHBRANCHES_V001 Classify a 1D forward-map curve into invertible
% monotonic branches (rising and/or descending).
%
% Faithful extraction of the local functions findBothBranches/
% findMonotonicRun/longestTrueRun from buildMonotonicSegments_v2_002.m,
% where this algorithm was applied to pre-built v058 SIMULATION GRID
% cells (one call per (pipeline, alpha, sigma, fs, VGF) grid coordinate,
% sweeping the beta_gen axis). Algorithm is UNCHANGED here -- only the
% packaging (private local function -> standalone reusable file) differs.
%
% Extracted 2026-08-09 specifically so this same algorithm can be applied
% to a PER-TRIAL forward-map curve (this trial's own exact alpha/sigma/
% geometry, no grid snapping) rather than only to the pre-built grid --
% see runLoopClosureFftnoise_v009's invertBeta_local. Before this
% extraction, the only place this algorithm ran was on the grid, via
% checkInvertibilityForEmpirical_v2_002's nearest-neighbour snap, which
% is a materially different (and, per the 2026-08-09 Core Function
% correction, superseded) methodology from genuine per-trial loop closure.
%
% buildMonotonicSegments_v2_002.m's own internal copy is deliberately left
% untouched, not refactored to call this file -- its existing
% monotonicSegments_v2_002.mat output must remain reproducible from a
% script that doesn't depend on files outside its own version lineage,
% matching this project's established pattern (see e.g. estimateIRASA.m's
% "Bundled here, not called externally" rationale).
%
% SYNTAX:
%   both = findBothBranches_v001(xGrid, yRaw, params)
%
% INPUTS:
%   xGrid  - Generator-axis grid values (e.g. betaGenVec), real vector.
%            Must be sorted ascending.
%   yRaw   - Recovered-value curve at each xGrid point (e.g. median
%            beta_rec across reps), same length as xGrid. NaN entries are
%            dropped before analysis (caller should pre-filter if a
%            different NaN-handling policy is wanted).
%   params - struct with fields:
%              .slopeTol     - minimum |slope| (dy/dx) to count as part of
%                               a monotonic run (default context:
%                               buildMonotonicSegments_v2_002 uses 0.05)
%              .smoothWindow - movmean window applied to yRaw before slope
%                               thresholding, purely for run-detection
%                               (the returned seg.y is NOT smoothed -- see
%                               NOTE below). buildMonotonicSegments_v2_002
%                               uses 3.
%              .minSegLength - minimum number of points for a segment to
%                               be considered invertible. Also used as the
%                               minimum finite-point floor by callers.
%                               buildMonotonicSegments_v2_002 uses 3.
%
% OUTPUT:
%   both.rise, both.desc - each a struct:
%     .x          - generator-axis values spanning the monotonic run
%                    (ascending order in both cases; for the descending
%                    branch this means the branch is stored REVERSED
%                    relative to xGrid, so that .y below is increasing)
%     .y          - recovered-value curve over .x. GUARANTEED strictly
%                    monotonically increasing (by construction: trimmed to
%                    a strictly-monotonic prefix on the RAW, unsmoothed y,
%                    and the descending branch is flipped) -- callers can
%                    rely on this for interp1(seg.y, seg.x, ...) inversion
%                    without re-checking.
%     .n          - number of points in the segment
%     .meanSlope  - mean(diff(y)./diff(x)) over the segment
%     .invertible - true iff n >= params.minSegLength AND
%                    abs(meanSlope) > params.slopeTol
%
% NOTE ON SMOOTHING: smoothWindow affects only which samples are classed
% as "rising enough" / "descending enough" to seed the longest-run search;
% the segment's own .y values returned are the RAW (unsmoothed) yRaw at
% those points, so seg.y reflects the actual recovered-value curve, not a
% smoothed proxy of it.
%
% See also: findMonotonicRun (local), buildMonotonicSegments_v2_002
%
% Fraser, D.S. (2026)

both.rise = findMonotonicRun(xGrid, yRaw, params, "rising");
both.desc = findMonotonicRun(xGrid, yRaw, params, "descending");

end

% --------------------------------------------------------------------------
function seg = findMonotonicRun(xGrid, yRaw, params, direction)
    seg = struct("x", [], "y", [], "n", 0, "meanSlope", NaN, "invertible", false);

    good = ~isnan(yRaw);
    if nnz(good) < params.minSegLength, return, end
    x = xGrid(good);
    y = yRaw(good);

    if params.smoothWindow > 1 && numel(y) >= params.smoothWindow
        ys = movmean(y, params.smoothWindow);
    else
        ys = y;
    end

    slope = diff(ys) ./ diff(x);
    switch direction
        case "rising",     isValid = slope >  params.slopeTol;
        case "descending", isValid = slope < -params.slopeTol;
        otherwise
            error("findMonotonicRun:BadDirection", "%s", ...
                "direction must be 'rising' or 'descending', got: " + direction);
    end

    [runStart, runEnd, runLen] = longestTrueRun(isValid);
    if isempty(runStart) || runLen < params.minSegLength - 1, return, end

    sx = x(runStart:runEnd+1);
    sy = y(runStart:runEnd+1);

    % Trim to strictly-monotonic prefix on raw y
    if direction == "rising"
        bad = find(diff(sy) <= 0, 1, "first");
    else
        bad = find(diff(sy) >= 0, 1, "first");
    end
    if ~isempty(bad)
        sx = sx(1:bad);
        sy = sy(1:bad);
    end
    if numel(sy) < params.minSegLength, return, end

    % Flip descending branch so segY is monotonically increasing for interp1
    if direction == "descending"
        sx = flipud(sx(:));
        sy = flipud(sy(:));
    else
        sx = sx(:);
        sy = sy(:);
    end

    seg.x          = sx;
    seg.y          = sy;
    seg.n          = numel(sy);
    seg.meanSlope  = mean(diff(sy) ./ diff(sx));
    seg.invertible = seg.n >= params.minSegLength && abs(seg.meanSlope) > params.slopeTol;
end

% --------------------------------------------------------------------------
function [runStart, runEnd, runLen] = longestTrueRun(v)
    runStart = []; runEnd = []; runLen = 0;
    d      = diff([false; v(:); false]);
    starts = find(d ==  1);
    ends   = find(d == -1) - 1;
    if isempty(starts), return, end
    lens = ends - starts + 1;
    [runLen, best] = max(lens);
    runStart = starts(best);
    runEnd   = ends(best);
end
