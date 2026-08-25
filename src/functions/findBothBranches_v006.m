function both = findBothBranches_v006(xGrid, yRaw, sigmaY, params)
%FINDBOTHBRANCHES_V006 Classify a 1D forward-map curve into ALL invertible
% monotonic branches -- grid-invariant segment width/smoothing (v004) PLUS
% a statistically-motivated slope-significance threshold for DETECTION
% (corrects v005's design error), with the project's MDC/SEM=0.011
% precision criterion kept as a SEPARATE, downstream annotation rather
% than the detection bar itself.
%
% WHY v005 WAS WRONG, NOT JUST INSUFFICIENT (2026-08-11): v005 used
% |slope| > sigma_y/0.011 as the criterion for whether a segment exists
% at all. Tested against real data: ALL 36 cells returned covered=[0,0,0]
% at every grid density -- not a fix, a degenerate collapse. Diagnosis:
% sigma_y/0.011 conflates two different questions. "Is this apparent
% trend distinguishable from noise" (detection) is a statistical
% significance question, answered by comparing the observed slope to the
% ESTIMATION UNCERTAINTY of that slope. "Would inverting through this
% trend meet clinical precision" (the project's actual MDC/SEM=0.011
% adequacy criterion, prereg_v101 Section 3.2) is a different, downstream
% question about how good the resulting estimate is, not whether a trend
% exists in the first place. Checked directly: at this project's real
% N_REPS=8 noise scale, sigma_y/0.011 gives implied thresholds (median
% ~4, up to ~10) far exceeding actual observed slopes (-2.7 to 4.5) --
% essentially unclearable, hence the collapse.
%
% v006's DETECTION criterion: for the finite-difference slope between
% adjacent points i, i+1, the slope's own ESTIMATION uncertainty is
%   sigma_slope = sqrt(sigmaY(i)^2 + sigmaY(i+1)^2) / dx
% (independent-noise error propagation through the finite-difference
% formula itself). Require |slope| > z * sigma_slope (params.zThreshold,
% default 2 -- roughly a 95% one-sided significance bar that the true
% slope is nonzero, not tied to any clinical constant). This directly
% targets Finding #149's aliasing mechanism: at coarser grids, dx grows,
% so sigma_slope shrinks relative to a given sigma_y, making a spurious
% coarse-grid slope MORE likely to clear a flat bar and LESS likely to
% clear this one at the same sigma_y -- the opposite of what was
% happening under the old flat slopeTol.
%
% v006's PRECISION annotation (separate from detection, not gating it):
% each returned segment also reports .meetsClinicalPrecision, a logical
% -- whether THIS segment's own local slope, if used to invert an
% observed value, would give sigma_x = sigma_y_local/|slope| < the
% project's own MDC/2.77 = 0.011 threshold (params.mdcThreshold). This
% keeps the clinical-adequacy question available downstream (e.g. to
% flag a betaGenStar as detected-but-imprecise) without corrupting
% detection itself.
%
% v005 is left completely unmodified (new file, matching this project's
% version-lineage-immutability convention throughout this thread) --
% v005's design error is documented here, not silently fixed in place.
%
% See also: findBothBranches_v004 (base this extends), findBothBranches_v005
% (the design error this corrects -- read its own docstring for the
% detail, not repeated here)
%
% Fraser, D.S. (2026)

both.rise = findMonotonicRuns_v006(xGrid, yRaw, sigmaY, params, "rising");
both.desc = findMonotonicRuns_v006(xGrid, yRaw, sigmaY, params, "descending");

end

% --------------------------------------------------------------------------
function segs = findMonotonicRuns_v006(xGrid, yRaw, sigmaY, params, direction)
    segTemplate = struct("x", [], "y", [], "n", 0, "meanSlope", NaN, ...
        "invertible", false, "meetsClinicalPrecision", false);
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

    if isfield(params, 'zThreshold'),   zThreshold   = params.zThreshold;   else, zThreshold   = 2;     end
    if isfield(params, 'mdcThreshold'), mdcThreshold = params.mdcThreshold; else, mdcThreshold = 0.011; end

    slope       = diff(ys) ./ diff(x);
    dx          = diff(x);
    sigmaSlope  = sqrt(sy(1:end-1).^2 + sy(2:end).^2) ./ dx;  % estimation uncertainty of the slope itself
    detectThresh = zThreshold * sigmaSlope;

    switch direction
        case "rising",     isValid = slope >  detectThresh;
        case "descending", isValid = slope < -detectThresh;
        otherwise
            error("findMonotonicRuns_v006:BadDirection", "%s", ...
                "direction must be 'rising' or 'descending', got: " + direction);
    end

    [starts, ends] = allTrueRuns(isValid);

    for k = 1:numel(starts)
        runStart = starts(k); runEnd = ends(k);

        if (x(runEnd+1) - x(runStart)) < params.minSegWidth, continue, end

        sx  = x(runStart:runEnd+1);
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

        % Re-check significance at the segment level (trimming can change
        % which points are included -- re-derive sigma_slope over the
        % final segment's own span, not just re-use the per-step gate).
        segIdxInGood = ismember(x, sx);
        segSy = sy(segIdxInGood);
        segSpan = max(sx) - min(sx);
        segSigmaSlope = sqrt(sum(segSy.^2)) / segSpan / sqrt(max(numel(segSy)-1,1));
        if abs(meanSlope) <= zThreshold * segSigmaSlope, continue, end

        % Precision annotation, separate from detection.
        segSigmaLocal = mean(segSy);
        meetsPrecision = (segSigmaLocal / abs(meanSlope)) < mdcThreshold;

        seg = struct("x", sx, "y", syv, "n", numel(syv), ...
            "meanSlope", meanSlope, "invertible", true, ...
            "meetsClinicalPrecision", meetsPrecision);
        segs = [segs, seg]; %#ok<AGROW>
    end
end

% --------------------------------------------------------------------------
function [starts, ends] = allTrueRuns(v)
    d      = diff([false; v(:); false]);
    starts = find(d ==  1);
    ends   = find(d == -1) - 1;
end
