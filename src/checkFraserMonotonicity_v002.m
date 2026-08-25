function out = checkFraserMonotonicity_v002(opts)
% checkFraserMonotonicity_v002  R6 (local monotonicity margin) + R7
% (surrogate alpha-gap check) for Fraser's ALL-TRIALS loop-closure result
% (N=2829), superseding checkFraserMonotonicity_v001.m's per_subject-based
% result (N=48) now that the definitive all-trials tier exists (matching
% Cook CTRL/ASD, Hickman PLAC/HALO's own TRIAL_SELECTION='all' basis).
%
% v001 is left untouched on disk; this is a new file, not an edit, per
% project versioning convention. R5 (SEM lookup) is NOT re-run here: it
% depends only on Fraser's noise centroid (alpha=4.289, sigma=2.009mm,
% fs=240Hz), computed once at import from noiseCharacterisation_fraser.mat
% across the full imported trial set -- independent of which subset the
% downstream loop-closure regression is run on, so it is unaffected by the
% per_subject -> all switch and needs no re-check.
%
% Methodology (R6) identical to v001: same indirect-slope and E5-style
% boundary formulas as loopClosureSlopeMargin_v001.m, verified there against
% Cook CTRL/Pilot in Session 85. Methodology (R7) identical to v001's R7
% section: mean |surrogateAlphaGapMaj/Min| read directly from the results
% struct, compared against the README_ShapedXu_v001.md Sec 6.1 gate (0.5)
% and the in-runner range documented for the other five datasets (Finding
% #75: 0.199-0.325).
%
% OUTPUT: console tables + struct `out` + checkFraserMonotonicity_v002.mat.
%
% USAGE: checkFraserMonotonicity_v002()
%
% Fraser, D.S. (2026)  v002

    arguments
        opts.NBoot    (1,1) double {mustBeInteger, mustBePositive} = 2000
        opts.NInner   (1,1) double {mustBeInteger, mustBePositive} = 200
        opts.RngSeed  (1,1) double {mustBeInteger}                 = 20260807
        opts.PrimaryPP(1,1) double {mustBeInRange(opts.PrimaryPP, 1, 6)} = 4
    end

    PP_NAMES     = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];
    PP_PRIMARY   = opts.PrimaryPP;
    N_PIPELINES  = 6;
    Z90          = 1.644853626951472;
    PUREXU_BOUND = 0.30;
    ONE_THIRD    = 1/3;
    LABEL        = "Fraser (all-trials)";

    srcDir  = fileparts(mfilename("fullpath"));
    if isempty(srcDir), srcDir = pwd; end
    matFile = fullfile(srcDir, "loopClosureResults_Fraser_all_shaped_xu_v008.mat");
    outFile = fullfile(srcDir, "checkFraserMonotonicity_v002.mat");

    rng(opts.RngSeed, "twister");

    if ~isfile(matFile)
        error("checkFraserMonotonicity_v002:notFound", "%s", ...
            sprintf("FAILED PATH: %s not found.", matFile));
    end
    D = load(matFile, "results", "pipelineLabels");
    r = D.results;
    N = numel(r);

    if ~isequal(string(D.pipelineLabels), PP_NAMES)
        error("checkFraserMonotonicity_v002:pipelineOrder", "%s", ...
            "pipelineLabels in the mat do not match the assumed PP_NAMES order -- check before trusting PrimaryPP index.");
    end

    %% --- R6: per-trial arrays --------------------------------------------
    bgStarMed = nan(N,1);
    bgStarPP  = nan(N, N_PIPELINES);
    ciHiMat   = nan(N, N_PIPELINES);
    ciLoMat   = nan(N, N_PIPELINES);
    sliceP5   = nan(N, N_PIPELINES);
    sliceP95  = nan(N, N_PIPELINES);
    gapMaj    = nan(N,1);
    gapMin    = nan(N,1);

    for ti = 1:N
        bgStarMed(ti)  = r(ti).betaGenStarMed;
        bgStarPP(ti,:) = r(ti).betaGenStar(:)';
        ciHiMat(ti,:)  = r(ti).ciHi(:)';
        ciLoMat(ti,:)  = r(ti).ciLo(:)';
        if isfield(r(ti),'surrogateAlphaGapMaj'), gapMaj(ti) = r(ti).surrogateAlphaGapMaj; end
        if isfield(r(ti),'surrogateAlphaGapMin'), gapMin(ti) = r(ti).surrogateAlphaGapMin; end

        sl = r(ti).betaRecSlice;
        if isnumeric(sl) && size(sl,2) == N_PIPELINES && any(isfinite(sl(:)))
            for pp = 1:N_PIPELINES
                col = sl(:, pp);
                if sum(isfinite(col)) >= 2
                    q = prctile(col(isfinite(col)), [5 95]);
                    sliceP5(ti, pp)  = q(1);
                    sliceP95(ti, pp) = q(2);
                end
            end
        end
    end

    %% --- R6: indirect slope -----------------------------------------------
    ciWidth  = ciLoMat - ciHiMat;
    vSpread  = sliceP95 - sliceP5;
    slopeInd = nan(N, N_PIPELINES);
    okSlope  = isfinite(ciWidth) & isfinite(vSpread) & (abs(ciWidth) > eps);
    slopeInd(okSlope) = vSpread(okSlope) ./ ciWidth(okSlope);

    medBgStar = median(bgStarMed, "omitnan");

    %% --- R6: forward-map-aware bootstrap -----------------------------------
    sdInvPP  = ciWidth / (2*Z90);
    sdInvMed = innerMedianSD_local(bgStarPP, sdInvPP, opts.NInner);

    bootPrim = bootstrapBeta_local(bgStarPP(:,PP_PRIMARY), sdInvPP(:,PP_PRIMARY), opts.NBoot, ONE_THIRD);
    bootMed  = bootstrapBeta_local(bgStarMed, sdInvMed, opts.NBoot, ONE_THIRD);

    %% --- R6: boundary metrics ----------------------------------------------
    finP    = isfinite(bgStarMed) & isfinite(ciHiMat(:,PP_PRIMARY));
    shapedLB = median(ciHiMat(finP, PP_PRIMARY), "omitnan");
    marginPureXu = medBgStar - PUREXU_BOUND;
    marginShaped = medBgStar - shapedLB;
    pctMonoPrim  = 100 * mean(slopeInd(isfinite(slopeInd(:,PP_PRIMARY)), PP_PRIMARY) > 0);
    nUndefPrim   = sum(~isfinite(slopeInd(:,PP_PRIMARY)) & finP);

    isGood = (marginShaped >= 0.03) && (median(slopeInd(:,PP_PRIMARY),"omitnan") > 0) && (pctMonoPrim >= 90);
    verdict = "BORDERLINE-GOOD (downgrade)";
    if isGood, verdict = "GOOD"; end

    %% --- R7: surrogate alpha-gap, recomputed from source -------------------
    meanGapMaj = mean(gapMaj, 'omitnan');
    meanGapMin = mean(gapMin, 'omitnan');
    nMaj = sum(isfinite(gapMaj));
    nMin = sum(isfinite(gapMin));
    STANDALONE_GATE = 0.5;

    %% --- Pack ---------------------------------------------------------------
    out.label = LABEL; out.N_total = N; out.medBgStar = medBgStar;
    out.slopeInd_med  = median(slopeInd, 1, "omitnan");
    out.slopeInd_prim = median(slopeInd(:,PP_PRIMARY), "omitnan");
    out.pctMono_prim  = pctMonoPrim;
    out.nUndef_prim   = nUndefPrim;
    out.boot_primary  = bootPrim;
    out.boot_median   = bootMed;
    out.shapedLB      = shapedLB;
    out.marginPureXu  = marginPureXu;
    out.marginShaped  = marginShaped;
    out.e5_verdict    = verdict;
    out.r7_gapMaj     = meanGapMaj;
    out.r7_gapMin     = meanGapMin;

    %% --- Print: R6 -------------------------------------------------------
    fprintf("\n%s\n", repmat("=",1,78));
    fprintf("FRASER (ALL-TRIALS, N=%d): LOCAL MONOTONICITY MARGIN at beta_gen*  (SG-LMLS primary, pp=%d)\n", N, PP_PRIMARY);
    fprintf("%s\n", repmat("=",1,78));
    fprintf("%-20s  %8s  %9s  %7s  %7s\n", "Dataset", "bg*_med", "slope_ind", "%mono", "nUndef");
    fprintf("%s\n", repmat("-",1,78));
    fprintf("%-20s  %8.4f  %9.3f  %6.1f%%  %7d\n", LABEL, medBgStar, out.slopeInd_prim, pctMonoPrim, nUndefPrim);

    fprintf("\n--- Indirect slope (median over trials), all six pipelines ---\n");
    fprintf("%-20s", "Dataset");
    for pp = 1:N_PIPELINES, fprintf("  %9s", PP_NAMES(pp)); end
    fprintf("\n%-20s", LABEL);
    for pp = 1:N_PIPELINES, fprintf("  %9.3f", out.slopeInd_med(pp)); end
    fprintf("\n");

    fprintf("\n%s\n", repmat("=",1,78));
    fprintf("FRASER (ALL-TRIALS): FORWARD-MAP-AWARE BOOTSTRAP CI on beta_gen*  (NBoot=%d)\n", opts.NBoot);
    fprintf("%s\n", repmat("=",1,78));
    fprintf("%-24s  %5s  %8s  %20s  %20s  %8s\n", ...
        "Dataset [estimator]", "N", "median", "FM-aware 95% CI", "trial-only 95% CI", "p(>=1/3)");
    fprintf("%s\n", repmat("-", 1, 92));
    printBootRow_local(LABEL + " [SG-LMLS]", bootPrim);
    printBootRow_local(LABEL + " [median]",  bootMed);

    fprintf("\n%s\n", repmat("=",1,78));
    fprintf("FRASER (ALL-TRIALS): BOUNDARY DECISION (E5-style, cf. Finding #88)\n");
    fprintf("%s\n", repmat("=",1,78));
    fprintf("%-20s  bg*=%.4f  shaped_xuLB=%.4f  d(0.30)=%+.4f  d(LB)=%+.4f  slope=%+.3f\n", ...
        LABEL, medBgStar, shapedLB, marginPureXu, marginShaped, out.slopeInd_prim);
    fprintf("\nVERDICT: %s   (rule: GOOD if d(LB) >= 0.03 AND slope > 0 AND %%mono >= 90)\n", verdict);
    fprintf("\nFor comparison, per_subject result (checkFraserMonotonicity_v001.m, N=48):\n");
    fprintf("  bg*=0.3283, slope=+0.313, %%mono=100.0%%, d(LB)=+0.1038, verdict GOOD\n");

    %% --- Print: R7 ---------------------------------------------------------
    fprintf("\n%s\n", repmat("=",1,78));
    fprintf("FRASER (ALL-TRIALS): R7 SURROGATE ALPHA-GAP CHECK\n");
    fprintf("%s\n", repmat("=",1,78));
    fprintf("  N trials = %d\n", N);
    fprintf("  mean |gap| major axis = %.4f  (N=%d finite)\n", meanGapMaj, nMaj);
    fprintf("  mean |gap| minor axis = %.4f  (N=%d finite)\n", meanGapMin, nMin);
    fprintf("  Standalone gate PASS threshold (README_ShapedXu_v001.md Sec 6.1): %.3f\n", STANDALONE_GATE);
    fprintf("  In-runner range, other datasets (Finding #75): 0.199-0.325\n");
    fprintf("  Per_subject result (v001, N=48): 0.297 / 0.344\n");
    if meanGapMaj < STANDALONE_GATE && meanGapMin < STANDALONE_GATE
        fprintf("  Verdict: comfortably below the %.1f gate threshold.\n", STANDALONE_GATE);
    else
        fprintf("  WARNING: exceeds the standalone gate threshold. Flag before citing loopCCC.\n");
    end

    save(outFile, "out");
    fprintf("\nSaved: %s\n", outFile);
end

%% ======================== local functions ============================

function sdMed = innerMedianSD_local(bgStarPP, sdInvPP, NInner)
    N = size(bgStarPP, 1);
    sdMed = zeros(N, 1);
    for ti = 1:N
        pts = bgStarPP(ti, :);
        sds = sdInvPP(ti, :);
        conv = isfinite(pts);
        if ~any(conv), sdMed(ti) = 0; continue; end
        p = pts(conv);
        d = sds(conv);
        d(~isfinite(d) | d < 0) = 0;
        if all(d == 0), sdMed(ti) = 0; continue; end
        draws = p(:)' + d(:)' .* randn(NInner, numel(p));
        meds  = median(draws, 2);
        sdMed(ti) = std(meds, 0);
    end
end

function b = bootstrapBeta_local(point, sdInv, NBoot, refVal)
    fin   = isfinite(point);
    p     = point(fin);
    s     = sdInv(fin);
    s(~isfinite(s) | s < 0) = 0;
    nUndefCI = sum(fin & (~isfinite(sdInv) | sdInv < 0));
    N     = numel(p);

    b = struct("N", N, "nUndefCI", nUndefCI, "median", NaN, ...
        "ciFM", [NaN NaN], "ciTrial", [NaN NaN], "pGE", NaN, "refVal", refVal);
    if N < 3, return; end

    b.median = median(p);
    medFM    = nan(NBoot, 1);
    medTrial = nan(NBoot, 1);
    for bi = 1:NBoot
        idx = randi(N, N, 1);
        pr  = p(idx);
        medTrial(bi) = median(pr);
        medFM(bi)    = median(pr + s(idx) .* randn(N, 1));
    end
    b.ciFM    = prctile(medFM,    [2.5 97.5]);
    b.ciTrial = prctile(medTrial, [2.5 97.5]);
    b.pGE     = mean(medFM >= refVal);
end

function printBootRow_local(name, b)
    fprintf("%-24s  %5d  %8.4f  [%7.4f, %7.4f]  [%7.4f, %7.4f]  %8.4f\n", ...
        name, b.N, b.median, b.ciFM(1), b.ciFM(2), b.ciTrial(1), b.ciTrial(2), b.pGE);
end
