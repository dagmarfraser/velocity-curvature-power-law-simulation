function out = compareHickmanPlacHalo_v002(opts)
% compareHickmanPlacHalo_v002  Hickman PLAC vs HALO, properly powered AND
% properly paired.
%
% NOT a parameter change on compareCookCtrlAsd_v001.m -- that script's joint
% bootstrap assumes two independent groups, which is correct for Cook
% CTRL/ASD (different people) but WRONG here: PLAC and HALO are drawn from
% loopClosureResults_Hickman_PLAC/HALO_all_shaped_xu_v007.mat, and 28 of
% their 35 subjects each are the SAME people under two drug conditions (a
% crossover design, confirmed directly by subjectID intersection: 28
% shared, 7 PLAC-only, 7 HALO-only). Finding #30's own original check
% already used a paired t-test on the 28 shared subjects for exactly this
% reason. This reruns that same paired comparison with the cluster-
% bootstrap machinery built today, rather than a flat paired t-test.
%
% Design: a ONE-sample bootstrap of 28 within-subject paired differences,
% not a two-sample joint bootstrap. For each of the 28 shared subjects:
% resample their own PLAC trials, resample their own HALO trials
% (independently of each other -- different sessions), apply FM-aware
% perturbation to both, take each condition's median, then that subject's
% difference for this iteration. Resample which subjects (with replacement)
% at the outer level, take the across-subject median of differences. This
% is structurally identical to loopClosureVarDecomp_v008/v009's
% clusterEstimate_local applied to a single dataset, just with two trial-
% index-sets per subject instead of one -- t(M-1) with M=28 (df=27) is the
% correct reference, not the two-sample Welch correction
% compareCookCtrlAsd_v001.m needed.
%
% The 7 PLAC-only and 7 HALO-only subjects cannot contribute to a paired
% difference at all and are excluded from the primary analysis (matching
% Finding #30's own precedent) -- reported, not silently dropped, but not
% folded into a more complex hybrid paired/unpaired model without being
% asked for one.
%
% Also reports the NAIVE independent-groups estimate (all 35 vs all 35,
% ignoring the shared-subject structure entirely) alongside the correct
% paired estimate, so the size of the correction is visible, not silent --
% same convention as the flat-vs-cluster SE tables elsewhere this session.
%
% Same two-CI-level structure as compareCookCtrlAsd_v001.m: 95% CI tests
% whether the difference excludes zero; 90% CI (correct TOST construction)
% tests whether it fits entirely within +/-MDC.
%
%
% v002 (2026-09-19, Session 113): reruns Finding #139's comparison on the
% gated, N_REPS=200 corpora, replacing legacy v007. Default is v015, the
% paper's N_REPS=200 basis from Session 113 (decided before these
% contrasts were run); v013/v014 ("confirmatory") is kept as a replicate.
% Three changes only; the bootstrap logic is untouched:
%   1. loadResults_local loads v015 (default) or v013/v014 (exactly one).
%   2. CI width is ciHi - ciLo. Post-v009 corpora are ascending (Finding
%      #196); v001's ciLo - ciHi would be negative there, and the downstream
%      clamp would silently zero every FM-aware perturbation. The ascending
%      convention is now ASSERTED, not assumed.
%   3. The constellation-median point estimate(s) are checked against
%      cluster_median.pointEst in loopClosureVarDecomp_v012_gated_all6_
%      <corpus>.mat (v015 or confirmatory), so this uses the same
%      estimator as the pool built on the same corpus.
% Context: Finding #215 (per-trial correction inadequate; this asks whether
% the clinical contrast is resolvable in aggregate).
% opts.Corpus: "v015" (default, paper basis) or "confirmatory"
%   (v013/v014 replicate; output mat suffixed _confirmatory). Finding #204
%   drift means the two differ slightly; each is checked against its own
%   pool mat.
%
% Fraser, D.S. (2026)  v002

    arguments
        opts.NBoot     (1,1) double {mustBeInteger, mustBePositive} = 2000
        opts.NInner    (1,1) double {mustBeInteger, mustBePositive} = 200
        opts.RngSeed   (1,1) double {mustBeInteger}                 = 20260619
        opts.PrimaryPP (1,1) double {mustBeInRange(opts.PrimaryPP, 1, 6)} = 4
        opts.MDC       (1,1) double {mustBePositive}                = 0.03
        opts.SaveMat   (1,1) logical = true
        opts.Corpus    (1,1) string {mustBeMember(opts.Corpus, ["confirmatory","v015"])} = "v015"
    end

    N_PIPELINES = 6;
    Z90 = 1.644853626951472;
    PP_PRIMARY = opts.PrimaryPP;
    MDC = opts.MDC;
    NBoot = opts.NBoot;

    srcDir = fileparts(mfilename("fullpath"));
    rng(opts.RngSeed, "twister");

    rP = loadResults_local(srcDir, opts.Corpus, "Hickman_PLAC");
    rH = loadResults_local(srcDir, opts.Corpus, "Hickman_HALO");

    [bgMedP, bgPPP, sdInvMedP, sdInvPPP, subjP] = extractDataset_local(rP, N_PIPELINES, Z90, opts.NInner);
    [bgMedH, bgPPH, sdInvMedH, sdInvPPH, subjH] = extractDataset_local(rH, N_PIPELINES, Z90, opts.NInner);

    uP = unique(subjP); uH = unique(subjH);
    shared  = intersect(uP, uH);
    placOnly = setdiff(uP, uH);
    haloOnly = setdiff(uH, uP);
    fprintf("Hickman PLAC: %d trials, %d subjects\n", numel(bgMedP), numel(uP));
    fprintf("Hickman HALO: %d trials, %d subjects\n", numel(bgMedH), numel(uH));
    fprintf("Shared (paired analysis uses these): %d\n", numel(shared));
    fprintf("PLAC-only (excluded from paired analysis): %d\n", numel(placOnly));
    fprintf("HALO-only (excluded from paired analysis): %d\n\n", numel(haloOnly));

    out = struct();
    out.median = struct();
    out.median.paired = pairedOneEstimator_local(bgMedP, sdInvMedP, subjP, bgMedH, sdInvMedH, subjH, shared, NBoot, opts.NInner, MDC, "constellation median, PAIRED (correct)");
    out.median.naiveIndependent = compareIndependent_local(bgMedP, sdInvMedP, subjP, bgMedH, sdInvMedH, subjH, NBoot, MDC, "constellation median, naive independent-groups (WRONG -- for comparison only)");

    out.primary = struct();
    out.primary.paired = pairedOneEstimator_local(bgPPP(:,PP_PRIMARY), sdInvPPP(:,PP_PRIMARY), subjP, ...
                                                    bgPPH(:,PP_PRIMARY), sdInvPPH(:,PP_PRIMARY), subjH, shared, NBoot, opts.NInner, MDC, "SG-LMLS, PAIRED (correct)");
    out.primary.naiveIndependent = compareIndependent_local(bgPPP(:,PP_PRIMARY), sdInvPPP(:,PP_PRIMARY), subjP, ...
                                                              bgPPH(:,PP_PRIMARY), sdInvPPH(:,PP_PRIMARY), subjH, NBoot, MDC, "SG-LMLS, naive independent-groups (WRONG -- for comparison only)");

    % All-subject arm medians (naive block) are the pool's estimator; the paired block uses shared subjects only.
    checkPoint_local(srcDir, opts.Corpus, "Hickman PLAC", out.median.naiveIndependent.pointEstPlac, "compareHickmanPlacHalo");
    checkPoint_local(srcDir, opts.Corpus, "Hickman HALO", out.median.naiveIndependent.pointEstHalo, "compareHickmanPlacHalo");

    out.sharedSubjects = shared; out.placOnlySubjects = placOnly; out.haloOnlySubjects = haloOnly;
    out.config = struct("Corpus", opts.Corpus, "NBoot", NBoot, "NInner", opts.NInner, "RngSeed", opts.RngSeed, ...
        "PP_PRIMARY", PP_PRIMARY, "MDC", MDC);
    out.runDate = string(datetime("now"));

    if opts.SaveMat
        outFile = fullfile(srcDir, "compareHickmanPlacHalo_v002" + pick_local(opts.Corpus == "confirmatory", "_confirmatory", "") + ".mat");
        save(outFile, "out", "-v7.3");
        fprintf("\nSaved: %s\n", outFile);
    end
end

%% ======================== per-dataset extraction (reused logic) ============

function [bgMed, bgPP, sdInvMed, sdInvPP, subjID] = extractDataset_local(r, N_PIPELINES, Z90, NInner)
    N = numel(r);
    bgPP  = nan(N, N_PIPELINES);
    ciHiMat = nan(N, N_PIPELINES);
    ciLoMat = nan(N, N_PIPELINES);
    subjID  = strings(N,1);
    for ti = 1:N
        bgPP(ti,:)    = getRowVec_local(r(ti), "betaGenStar", N_PIPELINES);
        ciHiMat(ti,:) = getRowVec_local(r(ti), "ciHi",        N_PIPELINES);
        ciLoMat(ti,:) = getRowVec_local(r(ti), "ciLo",        N_PIPELINES);
        subjID(ti)    = getSubjID_local(r(ti));
    end
    if any(subjID == "")
        error("compareHickmanPlacHalo:missingSubjectID", "%s", "Missing subjectID for one or more trials.");
    end
    m = isfinite(ciLoMat) & isfinite(ciHiMat);
    if any(ciLoMat(m) > ciHiMat(m))
        error("compareHickmanPlacHalo:ciConvention", "%s", "ciLo > ciHi found: expected the ascending post-v009 convention (Finding #196).");
    end
    ciWidth = ciHiMat - ciLoMat;
    sdInvPP = ciWidth / (2*Z90);
    bgMed = median(bgPP, 2, "omitnan");
    sdNaive = innerMedianSD_local(bgPP, sdInvPP, NInner);
    sdInvMed = floorMedianSD_local(bgPP, sdNaive);
end

function sdMed = innerMedianSD_local(bgStarPP, sdInvPP, NInner)
    N = size(bgStarPP, 1);
    sdMed = zeros(N, 1);
    for ti = 1:N
        pts = bgStarPP(ti, :); sds = sdInvPP(ti, :);
        conv = isfinite(pts);
        if ~any(conv), sdMed(ti) = 0; continue; end
        p = pts(conv); d = sds(conv);
        d(~isfinite(d) | d < 0) = 0;
        if all(d == 0), sdMed(ti) = 0; continue; end
        draws = p(:)' + d(:)' .* randn(NInner, numel(p));
        sdMed(ti) = std(median(draws, 2), 0);
    end
end

function sdMed = floorMedianSD_local(bgStarPP, sdNaive)
    MEDIAN_SE_CONST = 1.2533141373155003;
    nConv  = sum(isfinite(bgStarPP), 2);
    empSD  = std(bgStarPP, 0, 2, "omitnan");
    empSD(nConv < 2) = 0;
    floorSD = empSD ./ sqrt(max(nConv,1)) * MEDIAN_SE_CONST;
    floorSD(nConv < 2) = 0;
    sdMed = max(sdNaive, floorSD);
end

%% ======================== PAIRED one-sample cluster bootstrap (correct) ====

function res = pairedOneEstimator_local(pointP, sdInvP, subjP, pointH, sdInvH, subjH, sharedSubj, NBoot, NInner, MDC, label)
    finP = isfinite(pointP); pP = pointP(finP); sP = sdInvP(finP); sP(~isfinite(sP)|sP<0)=0; sidP = subjP(finP);
    finH = isfinite(pointH); pH = pointH(finH); sH = sdInvH(finH); sH(~isfinite(sH)|sH<0)=0; sidH = subjH(finH);

    M = numel(sharedSubj);
    idxP = cell(M,1); idxH = cell(M,1);
    for i = 1:M
        idxP{i} = find(sidP == sharedSubj(i));
        idxH{i} = find(sidH == sharedSubj(i));
    end
    if any(cellfun(@isempty, idxP)) || any(cellfun(@isempty, idxH))
        error("compareHickmanPlacHalo:missingSharedTrials", "%s", ...
            "A shared subject has zero finite trials in one condition -- cannot pair.");
    end

    subjP_pt = nan(M,1); subjH_pt = nan(M,1);
    for i = 1:M
        subjP_pt(i) = median(pP(idxP{i}), "omitnan");
        subjH_pt(i) = median(pH(idxH{i}), "omitnan");
    end
    subjDiff = subjP_pt - subjH_pt;
    pointEstDiff = median(subjDiff, "omitnan");

    diffBoot = nan(NBoot,1);
    for bi = 1:NBoot
        subjDraw = randi(M, M, 1);
        dBoot = nan(M,1);
        for i = 1:M
            src = subjDraw(i);
            iiP = idxP{src}; nP = numel(iiP); ddP = iiP(randi(nP,nP,1));
            iiH = idxH{src}; nH = numel(iiH); ddH = iiH(randi(nH,nH,1));
            mP = median(pP(ddP) + sP(ddP).*randn(nP,1), "omitnan");
            mH = median(pH(ddH) + sH(ddH).*randn(nH,1), "omitnan");
            dBoot(i) = mP - mH;
        end
        diffBoot(bi) = median(dBoot, "omitnan");
    end

    seDiff = std(diffBoot, 0, "omitnan");
    ci95pct = prctile(diffBoot, [2.5 97.5]);
    ci90pct = prctile(diffBoot, [5 95]);

    df = M - 1;
    tcrit95 = tinv(0.975, df);
    tcrit90 = tinv(0.95,  df);
    ci95t = [pointEstDiff - tcrit95*seDiff, pointEstDiff + tcrit95*seDiff];
    ci90t = [pointEstDiff - tcrit90*seDiff, pointEstDiff + tcrit90*seDiff];

    diffExcludesZero = (ci95t(1) > 0) || (ci95t(2) < 0);
    equivalentWithinMDC = (ci90t(1) > -MDC) && (ci90t(2) < MDC);

    if diffExcludesZero
        verdict = "REAL DIFFERENCE (95% CI excludes 0)";
    elseif equivalentWithinMDC
        verdict = "EQUIVALENT within MDC (90% TOST CI inside +/-MDC)";
    else
        verdict = "INCONCLUSIVE -- neither excludes 0 nor fits within +/-MDC";
    end

    res = struct("label", label, "M", M, "pointEstPlac", median(subjP_pt,"omitnan"), ...
        "pointEstHalo", median(subjH_pt,"omitnan"), "diffPoint", pointEstDiff, ...
        "seDiff", seDiff, "df", df, ...
        "ci95_percentile", ci95pct, "ci90_percentile", ci90pct, ...
        "ci95_t", ci95t, "ci90_t", ci90t, ...
        "diffExcludesZero", diffExcludesZero, "equivalentWithinMDC", equivalentWithinMDC, ...
        "verdict", verdict);

    fprintf("--- %s ---\n", label);
    fprintf("  M=%d shared subjects (paired). PLAC point=%.4f, HALO point=%.4f\n", M, res.pointEstPlac, res.pointEstHalo);
    fprintf("  Paired difference (PLAC-HALO): %.4f, bootstrap SE=%.4f (df=%d)\n", pointEstDiff, seDiff, df);
    fprintf("  95%% CI (t): [%.4f, %.4f]  -- excludes 0? %d\n", ci95t(1), ci95t(2), diffExcludesZero);
    fprintf("  90%% CI (t), TOST vs +/-%.3f: [%.4f, %.4f]  -- inside band? %d\n", MDC, ci90t(1), ci90t(2), equivalentWithinMDC);
    fprintf("  VERDICT: %s\n\n", verdict);
end

%% ======================== NAIVE independent-groups (for comparison only) ==

function res = compareIndependent_local(pointP, sdInvP, subjP, pointH, sdInvH, subjH, NBoot, MDC, label)
% Deliberately ignores the shared-subject structure -- treats all 35 PLAC
% and all 35 HALO as if independent. Reported ONLY to show the size of the
% correction the paired analysis makes; not to be cited as the result.
    finP = isfinite(pointP); pP = pointP(finP); sP = sdInvP(finP); sP(~isfinite(sP)|sP<0)=0; sidP = subjP(finP);
    finH = isfinite(pointH); pH = pointH(finH); sH = sdInvH(finH); sH(~isfinite(sH)|sH<0)=0; sidH = subjH(finH);

    uP = unique(sidP); Mp = numel(uP);
    uH = unique(sidH); Mh = numel(uH);
    idxP = cell(Mp,1); for i=1:Mp, idxP{i} = find(sidP==uP(i)); end
    idxH = cell(Mh,1); for i=1:Mh, idxH{i} = find(sidH==uH(i)); end

    subjMedP = nan(Mp,1); for i=1:Mp, subjMedP(i) = median(pP(idxP{i}), "omitnan"); end
    subjMedH = nan(Mh,1); for i=1:Mh, subjMedH(i) = median(pH(idxH{i}), "omitnan"); end
    pointEstP = median(subjMedP, "omitnan");
    pointEstH = median(subjMedH, "omitnan");
    diffPoint = pointEstP - pointEstH;

    diffBoot = nan(NBoot,1); pBoot = nan(NBoot,1); hBoot = nan(NBoot,1);
    for bi = 1:NBoot
        drawP = randi(Mp, Mp, 1); mP = nan(Mp,1);
        for i = 1:Mp
            ii = idxP{drawP(i)}; nI = numel(ii); dd = ii(randi(nI,nI,1));
            mP(i) = median(pP(dd) + sP(dd).*randn(nI,1), "omitnan");
        end
        drawH = randi(Mh, Mh, 1); mH = nan(Mh,1);
        for i = 1:Mh
            ii = idxH{drawH(i)}; nI = numel(ii); dd = ii(randi(nI,nI,1));
            mH(i) = median(pH(dd) + sH(dd).*randn(nI,1), "omitnan");
        end
        pBoot(bi) = median(mP, "omitnan"); hBoot(bi) = median(mH, "omitnan");
        diffBoot(bi) = pBoot(bi) - hBoot(bi);
    end

    seDiff = std(diffBoot, 0, "omitnan");
    sePBoot = std(pBoot, 0, "omitnan"); sHBoot = std(hBoot, 0, "omitnan");
    dfWelch = (sePBoot^2 + sHBoot^2)^2 / ( (sePBoot^4)/(Mp-1) + (sHBoot^4)/(Mh-1) );
    tcrit95 = tinv(0.975, dfWelch); tcrit90 = tinv(0.95, dfWelch);
    ci95t = [diffPoint - tcrit95*seDiff, diffPoint + tcrit95*seDiff];
    ci90t = [diffPoint - tcrit90*seDiff, diffPoint + tcrit90*seDiff];
    diffExcludesZero = (ci95t(1) > 0) || (ci95t(2) < 0);
    equivalentWithinMDC = (ci90t(1) > -MDC) && (ci90t(2) < MDC);

    res = struct("label", label, "Mp", Mp, "Mh", Mh, "pointEstPlac", pointEstP, "pointEstHalo", pointEstH, "diffPoint", diffPoint, ...
        "seDiff", seDiff, "dfWelch", dfWelch, "ci95_t", ci95t, "ci90_t", ci90t, ...
        "diffExcludesZero", diffExcludesZero, "equivalentWithinMDC", equivalentWithinMDC);

    fprintf("--- %s ---\n", label);
    fprintf("  (all %d PLAC vs all %d HALO, shared-subject structure IGNORED)\n", Mp, Mh);
    fprintf("  Difference: %.4f, bootstrap SE=%.4f  95%% CI [%.4f, %.4f]  90%% CI [%.4f, %.4f]\n\n", ...
        diffPoint, seDiff, ci95t(1), ci95t(2), ci90t(1), ci90t(2));
end

%% ======================== shared IO helpers =================================

function r = loadResults_local(srcDir, corpus, tag)
    if corpus == "v015", cand = "v015"; else, cand = ["v013","v014"]; end
    f = arrayfun(@(v) fullfile(srcDir, sprintf("loopClosureResults_%s_all_shaped_xu_%s.mat", tag, v)), cand);
    hit = isfile(f);
    if nnz(hit) ~= 1
        error("compareHickmanPlacHalo:notFound", "%s", sprintf("FAILED PATH: %s -- %d of {%s} present, expected exactly 1.", tag, nnz(hit), strjoin(cand, ",")));
    end
    matFile = f(hit);
    fprintf("%s: %s\n", tag, matFile);
    D = load(matFile, "results");
    r = D.results(:);
end

function s = getSubjID_local(r)
    s = "";
    if isfield(r, "subjectID") && ~isempty(r.subjectID)
        s = string(r.subjectID);
    end
end

function row = getRowVec_local(s, fld, nExpected)
    row = nan(1, nExpected);
    if isfield(s, fld)
        x = s.(fld);
        if isnumeric(x) && numel(x) == nExpected
            row = double(x(:)');
        end
    end
end

function checkPoint_local(srcDir, corpus, label, pointEst, errid)
% Same estimator as the reported pool? Compare with the confirmatory mat.
    f = fullfile(srcDir, "loopClosureVarDecomp_v012_gated_all6_" + corpus + ".mat");
    if ~isfile(f), error(errid + ":refMissing", "%s", "FAILED PATH: " + f); end
    O = load(f, "out"); O = O.out;
    k = find(string({O.label}) == label);
    if numel(k) ~= 1, error(errid + ":refLabel", "%s", "No unique '" + label + "' in " + corpus + " pool mat."); end
    ref = O(k).cluster_median.pointEst;
    fprintf("  Reproduction check, %s: this %.6f vs pool %.6f (diff %.2g)\n", label, pointEst, ref, pointEst - ref);
    if abs(pointEst - ref) > 1e-9
        warning(errid + ":refMismatch", "%s", sprintf("%s point estimate differs from the pool by %.3g: estimator or trial set is not identical.", label, pointEst - ref));
    end
end

function out = pick_local(cond, a, b)
    if cond, out = a; else, out = b; end
end
