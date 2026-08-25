function both = findBothBranches_v004(xGrid, yRaw, params)
%FINDBOTHBRANCHES_V004 Classify a 1D forward-map curve into ALL invertible
% monotonic branches -- grid-invariant minimum segment size (from v003)
% PLUS grid-invariant smoothing window (new in v004).
%
% HYPOTHESIS UNDER TEST (2026-08-11, following on from Finding #149/#150):
% params.smoothWindow (movmean, used for classification only, not the
% returned segment's raw y-values) was, like v002/v003's minSegLength, a
% grid POINT COUNT, not a beta_gen WIDTH. At N_BETA=25 (canonical
% default), smoothWindow=3 spans ~0.063 in beta_gen; at N_BETA=9 (one of
% the resolutions tested in Finding #149), the same "3 points" spans
% nearly a third of the whole sweep. This is a more direct match to
% Finding #149's identified mechanism (movmean smoothing over a
% disproportionate fraction of a shorter curve, hiding a genuine local
% reversal -- aliasing) than a slope-threshold recalibration would be.
%
% v004 replaces params.smoothWindow (point count) with params.smoothWidth
% (beta_gen units), converted internally to an integer window size at
% whatever grid spacing the caller's xGrid actually has:
%   windowPts = max(1, round(smoothWidth / median(diff(xGrid))))
% Default smoothWidth=0.0625, matching what smoothWindow=3 meant at the
% canonical N_BETA=25 default (continuity with existing behaviour at that
% specific resolution, same convention as v003's minSegWidth default).
%
% params.slopeTol is UNCHANGED (still a flat constant) -- this version
% deliberately isolates the smoothing-window fix on its own, to measure
% its individual contribution against Finding #149's 36-cell benchmark
% before deciding whether the noise-calibrated slopeTol redesign
% (docs/TODO_NoiseCalibratedInvertibilityGate_v001.md) is still needed on
% top.
%
% v003 is left completely unmodified (new file, matching this project's
% version-lineage-immutability convention throughout this thread).
%
% See also: findBothBranches_v003 (minSegWidth only; smoothWindow still
% point-count-based there)
%
% Fraser, D.S. (2026)

both.rise = findMonotonicRuns_v004(xGrid, yRaw, params, "rising");
both.desc = findMonotonicRuns_v004(xGrid, yRaw, params, "descending");

end

% --------------------------------------------------------------------------
function segs = findMonotonicRuns_v004(xGrid, yRaw, params, direction)
    segTemplate = struct("x", [], "y", [], "n", 0, "meanSlope", NaN, "invertible", false);
    segs = segTemplate([]);

    good = ~isnan(yRaw);
    if nnz(good) < 2, return, end
    x = xGrid(good);
    y = yRaw(good);

    % Grid-invariant smoothing window: convert smoothWidth (beta_gen units)
    % to a point count at THIS xGrid's own spacing, rather than trusting a
    % fixed point count across arbitrary N_BETA.
    if numel(x) >= 2
        gridSpacing = median(diff(sort(x)));
    else
        gridSpacing = Inf;
    end
    if isfield(params, 'smoothWidth') && isfinite(gridSpacing) && gridSpacing > 0
        windowPts = max(1, round(params.smoothWidth / gridSpacing));
    else
        windowPts = params.smoothWindow;  % fallback: behave like v002/v003
    end

    if windowPts > 1 && numel(y) >= windowPts
        ys = movmean(y, windowPts);
    else
        ys = y;
    end

    % Slope threshold UNCHANGED from v002/v003 -- deliberately isolating
    % the smoothing-window fix for this test.
    slope = diff(ys) ./ diff(x);
    switch direction
        case "rising",     isValid = slope >  params.slopeTol;
        case "descending", isValid = slope < -params.slopeTol;
        otherwise
            error("findMonotonicRuns_v004:BadDirection", "%s", ...
                "direction must be 'rising' or 'descending', got: " + direction);
    end

    [starts, ends] = allTrueRuns(isValid);

    for k = 1:numel(starts)
        runStart = starts(k); runEnd = ends(k);

        if (x(runEnd+1) - x(runStart)) < params.minSegWidth, continue, end

        sx = x(runStart:runEnd+1);
        sy = y(runStart:runEnd+1);

        if direction == "rising"
            bad = find(diff(sy) <= 0, 1, "first");
        else
            bad = find(diff(sy) >= 0, 1, "first");
        end
        if ~isempty(bad)
            sx = sx(1:bad);
            sy = sy(1:bad);
        end

        if numel(sy) < 2 || (max(sx) - min(sx)) < params.minSegWidth
            continue
        end

        if direction == "descending"
            sx = flipud(sx(:));
            sy = flipud(sy(:));
        else
            sx = sx(:);
            sy = sy(:);
        end

        meanSlope = mean(diff(sy) ./ diff(sx));
        if abs(meanSlope) <= params.slopeTol, continue, end

        seg = struct("x", sx, "y", sy, "n", numel(sy), ...
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
