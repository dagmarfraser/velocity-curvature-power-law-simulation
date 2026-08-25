function both = findBothBranches_v002(xGrid, yRaw, params)
%FINDBOTHBRANCHES_V002 Classify a 1D forward-map curve into ALL invertible
% monotonic branches (rising and/or descending) -- not just the longest.
%
% FIX vs v001 (2026-08-11, Finding #146): v001's findMonotonicRun kept
% only the single LONGEST qualifying run per direction (via
% longestTrueRun's max(lens)), discarding every other run entirely. For
% the actual use case -- given a trial's own observed beta, is THIS
% VALUE locally recoverable -- that is the wrong question. A curve can
% have a long qualifying run elsewhere that does not contain the
% observed value, and a shorter, entirely valid qualifying run nearby
% that does. Confirmed on real data, not just argued: checked 8 real
% Fraser trials x 6 pipelines (48 cells); 6 cells had >1 qualifying run
% in some direction, and in 4 of those v001's longest-only design
% produced a false-negative "neither" where an alternate run genuinely
% covered the observed value (e.g. trial F320_V2_S3_T019_sh3, SG-IRLS:
% betaObs=0.2331 fell outside the longest run's y=[0.062,0.222] but
% inside a second, shorter, 3-point run's y=[0.158,0.254] that v001
% never inspected). v001 is left unmodified (matches
% buildMonotonicSegments_v2_002's own internal-copy-left-untouched
% precedent) -- this is a new file, not a patch.
%
% v001 remains valid and is NOT deprecated for any use where "does this
% curve have a usably long monotonic run at all" is the actual question
% (e.g. general curve characterisation). It is specifically wrong for
% per-trial local-invertibility gating, which is the only thing this
% project uses it for (runLoopClosureFftnoise_v010's invertBeta_local,
% superseding v009's use of v001 for this purpose).
%
% SYNTAX:
%   both = findBothBranches_v002(xGrid, yRaw, params)
%
% INPUTS:  identical to v001 -- xGrid, yRaw, params.slopeTol/
%          smoothWindow/minSegLength.
%
% OUTPUT:
%   both.rise, both.desc -- each a STRUCT ARRAY (zero or more elements,
%   vs v001's single struct), one element per qualifying run, each with
%   the same fields as v001's single segment (.x, .y, .n, .meanSlope,
%   .invertible -- all still true by construction for elements present
%   in the array; a direction with no qualifying runs returns a 0x0
%   struct array with the same fieldnames, not an empty invertible-false
%   placeholder like v001 used). Callers must iterate, not assume a
%   single segment -- see runLoopClosureFftnoise_v010's invertBeta_local
%   for the reference caller pattern (checks membership across every
%   element of both.rise and both.desc, flags "ambiguous" if the
%   observed value falls inside more than one across the combined set,
%   not just "falls in both rise AND desc" as v001's caller did).
%
% See also: findBothBranches_v001 (superseded for this project's own
% per-trial gating use case, left unmodified), findMonotonicRun_v002 (local)
%
% Fraser, D.S. (2026)

both.rise = findMonotonicRuns_v002(xGrid, yRaw, params, "rising");
both.desc = findMonotonicRuns_v002(xGrid, yRaw, params, "descending");

end

% --------------------------------------------------------------------------
function segs = findMonotonicRuns_v002(xGrid, yRaw, params, direction)
    segTemplate = struct("x", [], "y", [], "n", 0, "meanSlope", NaN, "invertible", false);
    segs = segTemplate([]);  % 0x0 struct array with the right fieldnames

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
            error("findMonotonicRuns_v002:BadDirection", "%s", ...
                "direction must be 'rising' or 'descending', got: " + direction);
    end

    [starts, ends] = allTrueRuns(isValid);

    for k = 1:numel(starts)
        runStart = starts(k); runEnd = ends(k);
        if (runEnd - runStart + 1) < params.minSegLength - 1, continue, end

        sx = x(runStart:runEnd+1);
        sy = y(runStart:runEnd+1);

        % Trim to strictly-monotonic prefix on raw y (identical to v001)
        if direction == "rising"
            bad = find(diff(sy) <= 0, 1, "first");
        else
            bad = find(diff(sy) >= 0, 1, "first");
        end
        if ~isempty(bad)
            sx = sx(1:bad);
            sy = sy(1:bad);
        end
        if numel(sy) < params.minSegLength, continue, end

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
% Generalises v001's longestTrueRun: returns EVERY contiguous true run's
% [start,end] indices, not just the longest.
    d      = diff([false; v(:); false]);
    starts = find(d ==  1);
    ends   = find(d == -1) - 1;
end
