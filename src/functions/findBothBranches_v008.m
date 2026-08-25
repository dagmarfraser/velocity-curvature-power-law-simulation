function both = findBothBranches_v008(xGrid, yRaw, sigmaY, params)
%FINDBOTHBRANCHES_V008 Classify a 1D forward-map curve into ALL invertible
% monotonic branches -- v004's detection (width-based smoothing/segment
% floor, flat slopeTol) UNCHANGED, but the strict-monotonicity TRIM step
% is now noise-tolerant instead of zero-tolerance.
%
% THE MECHANISM THIS FIXES (Finding #153, 2026-08-11): v001-v007 all trim
% a candidate run to its strictly-monotonic prefix by hard-stopping the
% instant the RAW curve decreases AT ALL (diff(sy)<=0 for rising),
% however tiny. Concrete real-data example (trial F309_V2_S1_T031_sh3,
% BWFD-OLS): a genuine but tiny -0.0051 decrease truncated an otherwise
% robust 0.66-span rising segment, excluding betaObs by a margin of just
% 0.0036 -- almost certainly within noise. At coarser grid resolution the
% sampled points simply skipped over that exact location, so the segment
% ran uninterrupted and covered betaObs easily. This zero-tolerance trim,
% not slope-detection miscalibration (which v005/v006/v007 tried and
% failed to fix, see Finding #152/#153), is very likely the dominant
% remaining mechanism after v004's width fixes.
%
% v008's TRIM criterion: a decrease only ends the segment if it exceeds a
% noise-calibrated tolerance, not if it is merely negative.
%   rising:     stop when  diff(sy) <= -k * sigmaTolLocal
%   descending: stop when  diff(sy) >=  k * sigmaTolLocal
% where sigmaTolLocal is the local noise scale (mean of sigmaY at the two
% points spanning the step) and k (params.trimToleranceK, default 1) is a
% tunable multiplier. This compares a RAW VALUE difference to a RAW noise
% scale -- unlike v006/v007's slope-significance formula, there is no
% division by dx here, so this does NOT inherit the 1/dx structural bias
% that made v007 worse than the original problem (Finding #153).
%
% DETECTION (which points get grouped into a candidate run in the first
% place) is UNCHANGED from v004: width-based smoothing/segment floor,
% flat params.slopeTol. This version is deliberately scoped to fix ONLY
% the trim step, isolating its own contribution the same way v004's
% smoothing fix and v006's detection fix were each tested in isolation
% before being combined.
%
% API: sigmaY (vector, one value per xGrid point) is a REQUIRED input,
% same as v005-v007 -- used here only for trim tolerance, not detection.
%
% v004 (and v005/v006/v007) are left completely unmodified.
%
% PRODUCTION STATUS (updated 2026-08-11 as later findings landed, not
% rewritten to hide the sequence): built and validated on real data
% (Finding #156, 0/36 flip on the grid-density-sensitivity benchmark).
% Also tested inside a per-replicate consistency-voting wrapper
% (checkMonotonicityConsistency_v002) as part of validating the fix
% itself -- but Finding #158 found that wrapper is NOT required: this
% function applied once, directly to the median curve, at the SAME
% compute cost as every prior version, already achieves the full result.
% Independently cross-validated against R's MonotonicityTest package
% (Hall & Heckman 2000) on synthetic ground truth, Finding #157. This is
% the version wired into production: runLoopClosureFftnoise_v012.m's
% invertBeta_local calls this function directly, once per trial per
% pipeline, exactly the way v002/v003/v004 were called before it.
%
% See also: findBothBranches_v004 (the detection logic this reuses
% unchanged); runLoopClosureFftnoise_v012 (the production caller);
% checkMonotonicityConsistency_v002 (the replicate-voting wrapper this
% was validated inside, kept as a diagnostic/robustness cross-check but
% not required in the production path); docs/FINDINGS_REFERENCE.md
% Findings #153/#156/#157/#158 (the full diagnosis, fix, and validation)
%
% Fraser, D.S. (2026)

both.rise = findMonotonicRuns_v008(xGrid, yRaw, sigmaY, params, "rising");
both.desc = findMonotonicRuns_v008(xGrid, yRaw, sigmaY, params, "descending");

end

% --------------------------------------------------------------------------
function segs = findMonotonicRuns_v008(xGrid, yRaw, sigmaY, params, direction)
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

    if isfield(params, 'trimToleranceK'), trimK = params.trimToleranceK; else, trimK = 1; end

    % Detection: UNCHANGED from v004 -- flat slopeTol on the smoothed curve.
    slope = diff(ys) ./ diff(x);
    switch direction
        case "rising",     isValid = slope >  params.slopeTol;
        case "descending", isValid = slope < -params.slopeTol;
        otherwise
            error("findMonotonicRuns_v008:BadDirection", "%s", ...
                "direction must be 'rising' or 'descending', got: " + direction);
    end

    [starts, ends] = allTrueRuns(isValid);

    for k = 1:numel(starts)
        runStart = starts(k); runEnd = ends(k);

        if (x(runEnd+1) - x(runStart)) < params.minSegWidth, continue, end

        sx  = x(runStart:runEnd+1);
        syv = y(runStart:runEnd+1);
        ssy = sy(runStart:runEnd+1);

        % NOISE-TOLERANT TRIM (the fix): a decrease only ends the segment
        % if it exceeds a noise-calibrated tolerance, not if it is merely
        % negative.
        dY = diff(syv);
        sigmaTolLocal = (ssy(1:end-1) + ssy(2:end)) / 2;
        if direction == "rising"
            bad = find(dY <= -trimK * sigmaTolLocal, 1, "first");
        else
            bad = find(dY >=  trimK * sigmaTolLocal, 1, "first");
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
        if abs(meanSlope) <= params.slopeTol, continue, end

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
