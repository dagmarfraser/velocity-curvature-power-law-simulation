function both = findBothBranches_v003(xGrid, yRaw, params)
%FINDBOTHBRANCHES_V003 Classify a 1D forward-map curve into ALL invertible
% monotonic branches (rising and/or descending) -- grid-invariant minimum
% segment size (WIDTH in x-units, not grid-point count).
%
% FIX vs v002 (2026-08-11, Finding #149): v002's minimum qualifying
% segment was defined as a POINT COUNT (params.minSegLength, e.g. 3
% points), which is not a property of the curve -- it is a property of
% how many grid points the caller happened to sample. At N_BETA=25 (grid
% spacing ~0.031), 3 points spans ~0.063 in beta_gen; at N_BETA=50 the
% same "3 points" spans half that physical width. Confirmed on real data
% (Finding #149): resampling the SAME captured curve at different grid
% densities flipped classification in 72.2% of tested cells (26/36).
%
% v003 replaces the point-count floor with a WIDTH floor
% (params.minSegWidth, in beta_gen units): a segment must span at least
% this much of the x-axis after trimming, regardless of how many grid
% points that happens to take. minSegWidth=0.0625 default matches what
% v001/v002's minSegLength=3 meant at the canonical N_BETA=25 default (3
% points = 2 grid intervals = 2*0.75/24 = 0.0625), preserving continuity
% with existing results AT that specific resolution while making the
% criterion explicit and grid-invariant going forward.
%
% IMPORTANT, STATED PLAINLY: this fixes only the SPATIAL-EXTENT half of
% the grid-dependence problem. It does NOT fix the slope/smoothing step
% (params.slopeTol, params.smoothWindow) where Finding #149's own
% diagnosis identified the DOMINANT mechanism -- aliasing, where a
% coarser grid's movmean smooths past a genuine local reversal that a
% finer grid correctly detects. A segment that only ever forms because of
% that aliasing will still form under v003 exactly as it did under v002;
% v003 only changes whether a segment, once formed, is long ENOUGH to
% count. Do not expect this alone to resolve most of the 72.2% figure --
% see Finding #149's follow-up empirical check for the measured effect.
%
% SYNTAX:
%   both = findBothBranches_v003(xGrid, yRaw, params)
%
% INPUTS:
%   xGrid, yRaw     -- identical to v002
%   params.slopeTol, params.smoothWindow -- identical to v002, UNCHANGED
%   params.minSegWidth -- NEW, replaces params.minSegLength. Minimum
%                         x-axis span (same units as xGrid) a trimmed
%                         segment must cover to qualify. Default 0.0625.
%
% OUTPUT: identical shape to v002 (both.rise/both.desc struct arrays).
%
% v002 is left completely unmodified -- this is a new file, matching the
% established version-lineage-immutability convention (v001 was left
% unmodified when v002 was built, in turn).
%
% See also: findBothBranches_v002 (superseded for the spatial-extent
% criterion; slope/smoothing logic unchanged and reused here verbatim)
%
% Fraser, D.S. (2026)

both.rise = findMonotonicRuns_v003(xGrid, yRaw, params, "rising");
both.desc = findMonotonicRuns_v003(xGrid, yRaw, params, "descending");

end

% --------------------------------------------------------------------------
function segs = findMonotonicRuns_v003(xGrid, yRaw, params, direction)
    segTemplate = struct("x", [], "y", [], "n", 0, "meanSlope", NaN, "invertible", false);
    segs = segTemplate([]);  % 0x0 struct array with the right fieldnames

    good = ~isnan(yRaw);
    if nnz(good) < 2, return, end  % bare numerical minimum for any slope
    x = xGrid(good);
    y = yRaw(good);

    if params.smoothWindow > 1 && numel(y) >= params.smoothWindow
        ys = movmean(y, params.smoothWindow);
    else
        ys = y;
    end

    % Slope/smoothing logic UNCHANGED from v002 -- this is deliberately
    % not where this version's fix operates. See header note above.
    slope = diff(ys) ./ diff(x);
    switch direction
        case "rising",     isValid = slope >  params.slopeTol;
        case "descending", isValid = slope < -params.slopeTol;
        otherwise
            error("findMonotonicRuns_v003:BadDirection", "%s", ...
                "direction must be 'rising' or 'descending', got: " + direction);
    end

    [starts, ends] = allTrueRuns(isValid);

    for k = 1:numel(starts)
        runStart = starts(k); runEnd = ends(k);

        % Early-exit on the UN-TRIMMED candidate span: trimming can only
        % shrink a segment, never grow it, so if even the untrimmed span
        % fails the width floor, the trimmed one will too.
        if (x(runEnd+1) - x(runStart)) < params.minSegWidth, continue, end

        sx = x(runStart:runEnd+1);
        sy = y(runStart:runEnd+1);

        % Trim to strictly-monotonic prefix on raw y (identical to v002)
        if direction == "rising"
            bad = find(diff(sy) <= 0, 1, "first");
        else
            bad = find(diff(sy) >= 0, 1, "first");
        end
        if ~isempty(bad)
            sx = sx(1:bad);
            sy = sy(1:bad);
        end

        % The actual (post-trim) width floor -- this is the real fix vs
        % v002's numel(sy) < params.minSegLength point-count check.
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
% Identical to v002's allTrueRuns.
    d      = diff([false; v(:); false]);
    starts = find(d ==  1);
    ends   = find(d == -1) - 1;
end
