function both = findBothBranches_v007(xGrid, yRaw, sigmaY, nReps, params)
%FINDBOTHBRANCHES_V007 Classify a 1D forward-map curve into ALL invertible
% monotonic branches -- corrects a real error in v005/v006's slope
% estimation-uncertainty formula (used raw per-draw SD where the
% standard error of the MEDIAN estimate was required).
%
% THE BUG, FOUND BEFORE SHIPPING v006 AS FINAL (2026-08-11): v005/v006
% both computed sigma_slope = sqrt(sigmaY(i)^2+sigmaY(i+1)^2)/dx, treating
% sigmaY (the raw SD across N_REPS individual surrogate draws at a point)
% as if it were the standard error of medC(i) (the MEDIAN of those N_REPS
% draws, which is what the curve actually reports and what the slope is
% computed from). It is not -- medC(i)'s own uncertainty is smaller than
% the raw per-draw spread by a factor of approximately
% 1.2533/sqrt(N_REPS) (1.2533 = sqrt(pi/2), the median's efficiency loss
% vs the mean under approximate normality). Checked directly at
% N_REPS=8 (this thread's smoke-config captures): the correction factor
% is 0.443, meaning v006's detection threshold was roughly 2.3x too
% strict -- almost certainly the dominant reason coverage sat at only
% 19.4% rather than higher.
%
% v007's DETECTION criterion, corrected:
%   SE_median(i) = 1.2533 * sigmaY(i) / sqrt(nReps)
%   sigma_slope   = sqrt(SE_median(i)^2 + SE_median(i+1)^2) / dx
%   |slope| > z * sigma_slope   (z=2, unchanged from v006)
%
% API CHANGE vs v006: nReps (a scalar, the N_REPS used to generate the
% surrogate draws sigmaY was computed from) is now a REQUIRED input --
% the correction factor cannot be applied without knowing it.
%
% v006 is left completely unmodified (new file, matching this project's
% version-lineage-immutability convention throughout this thread) -- the
% error is documented here and in Finding #153, not silently patched in
% place.
%
% See also: findBothBranches_v006 (the version this corrects -- read its
% own docstring for the full detection/precision-annotation design, all
% unchanged here except the SE formula itself)
%
% Fraser, D.S. (2026)

both.rise = findMonotonicRuns_v007(xGrid, yRaw, sigmaY, nReps, params, "rising");
both.desc = findMonotonicRuns_v007(xGrid, yRaw, sigmaY, nReps, params, "descending");

end

% --------------------------------------------------------------------------
function segs = findMonotonicRuns_v007(xGrid, yRaw, sigmaY, nReps, params, direction)
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

    MEDIAN_SE_FACTOR = 1.2533;  % sqrt(pi/2)
    seMedian = MEDIAN_SE_FACTOR * sy / sqrt(nReps);   % SE of medC(i), not raw per-draw SD

    slope        = diff(ys) ./ diff(x);
    dx           = diff(x);
    sigmaSlope   = sqrt(seMedian(1:end-1).^2 + seMedian(2:end).^2) ./ dx;
    detectThresh = zThreshold * sigmaSlope;

    switch direction
        case "rising",     isValid = slope >  detectThresh;
        case "descending", isValid = slope < -detectThresh;
        otherwise
            error("findMonotonicRuns_v007:BadDirection", "%s", ...
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

        segIdxInGood  = ismember(x, sx);
        segSeMedian   = seMedian(segIdxInGood);
        segSpan       = max(sx) - min(sx);
        segSigmaSlope = sqrt(sum(segSeMedian.^2)) / segSpan / sqrt(max(numel(segSeMedian)-1,1));
        if abs(meanSlope) <= zThreshold * segSigmaSlope, continue, end

        segSigmaLocal  = mean(sy(segIdxInGood));  % raw noise, for the CLINICAL precision check
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
