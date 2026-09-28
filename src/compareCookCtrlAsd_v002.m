function out = compareCookCtrlAsd_v002(opts)
% compareCookCtrlAsd_v002  Cook CTRL vs Cook ASD, properly powered.
%
% Motivation (2026-08-08): Finding #114's original comparison (t(37)=0.264,
% p=0.793, flat between-subject t-test) was explicitly flagged there as
% underpowered even for a proper TOST equivalence test, and as an absence-
% of-evidence result, not evidence of absence. That comparison also predates
% the two-stage cluster bootstrap entirely (Finding #136). This reruns the
% comparison with the machinery actually built today: subject-median point
% estimates and a JOINT two-stage cluster bootstrap of the difference
% (resamples both groups' subjects and within-subject trials independently
% per iteration, with FM-aware perturbation -- byte-identical logic to
% loopClosureVarDecomp_v008/v009's clusterEstimate_local, reused rather than
% reimplemented), rather than combining two marginal SEs by formula.
%
% Two genuinely different tests, not one CI level reused twice:
%   - 95% CI on the difference: does it exclude 0? (a real difference)
%   - 90% CI on the difference: does it fall entirely within +/-MDC?
%     (TOST equivalence -- the textbook-correct construction: two one-sided
%     alpha=0.05 tests intersect to a (1-2*0.05)=90% two-sided CI, not 95%.
%     Using 95% for both would be a stricter-than-standard equivalence bound,
%     not simply "more conservative" -- it changes what's being claimed.)
% A result excluding neither is a real, reportable outcome: genuinely
% inconclusive even under the more rigorous machinery, not a failure of the
% script. This is what "specificity/selectivity" actually asks: can the new
% tools resolve what the old flat test structurally could not?
%
% Cook CTRL and Cook ASD are drawn from `loopClosureResults_Cook_CTRL_all_
% shaped_xu_v007.mat` / `..._Cook_ASD_..._v007.mat` (unchanged, proven files
% from throughout today's session).
%
% OUTPUT
%   out.median / out.primary : each with pointEst_ctrl, pointEst_asd,
%     diffPoint, seDiff, ci95 (excludes-zero test), ci90 (TOST-equivalence
%     test), diffExcludesZero, equivalentWithinMDC, verdict (string)
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

    rCTRL = loadResults_local(srcDir, opts.Corpus, "Cook_CTRL");
    rASD  = loadResults_local(srcDir, opts.Corpus, "Cook_ASD");

    [bgMedC, bgPPC, sdInvMedC, sdInvPPC, subjC] = extractDataset_local(rCTRL, N_PIPELINES, Z90, opts.NInner);
    [bgMedA, bgPPA, sdInvMedA, sdInvPPA, subjA] = extractDataset_local(rASD,  N_PIPELINES, Z90, opts.NInner);

    fprintf("Cook CTRL: %d trials, %d subjects\n", numel(bgMedC), numel(unique(subjC)));
    fprintf("Cook ASD:  %d trials, %d subjects\n\n", numel(bgMedA), numel(unique(subjA)));

    out = struct();
    out.median  = compareOneEstimator_local(bgMedC, sdInvMedC, subjC, bgMedA, sdInvMedA, subjA, NBoot, MDC, "constellation median");
    out.primary = compareOneEstimator_local(bgPPC(:,PP_PRIMARY), sdInvPPC(:,PP_PRIMARY), subjC, ...
                                             bgPPA(:,PP_PRIMARY), sdInvPPA(:,PP_PRIMARY), subjA, NBoot, MDC, "SG-LMLS (sensitivity)");

    checkPoint_local(srcDir, opts.Corpus, "Cook CTRL", out.median.pointEstCtrl, "compareCookCtrlAsd");
    checkPoint_local(srcDir, opts.Corpus, "Cook ASD",  out.median.pointEstAsd,  "compareCookCtrlAsd");

    out.config = struct("Corpus", opts.Corpus, "NBoot", NBoot, "NInner", opts.NInner, "RngSeed", opts.RngSeed, ...
        "PP_PRIMARY", PP_PRIMARY, "MDC", MDC);
    out.runDate = string(datetime("now"));

    if opts.SaveMat
        outFile = fullfile(srcDir, "compareCookCtrlAsd_v002" + pick_local(opts.Corpus == "confirmatory", "_confirmatory", "") + ".mat");
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
        error("compareCookCtrlAsd:missingSubjectID", "%s", "Missing subjectID for one or more trials.");
    end
    m = isfinite(ciLoMat) & isfinite(ciHiMat);
    if any(ciLoMat(m) > ciHiMat(m))
        error("compareCookCtrlAsd:ciConvention", "%s", "ciLo > ciHi found: expected the ascending post-v009 convention (Finding #196).");
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

%% ======================== joint two-sample cluster bootstrap ===============

function res = compareOneEstimator_local(pointC, sdInvC, subjC, pointA, sdInvA, subjA, NBoot, MDC, label)
    finC = isfinite(pointC); pC = pointC(finC); sC = sdInvC(finC); sC(~isfinite(sC)|sC<0)=0; sidC = subjC(finC);
    finA = isfinite(pointA); pA = pointA(finA); sA = sdInvA(finA); sA(~isfinite(sA)|sA<0)=0; sidA = subjA(finA);

    uC = unique(sidC); Mc = numel(uC);
    uA = unique(sidA); Ma = numel(uA);
    idxC = cell(Mc,1); for i=1:Mc, idxC{i} = find(sidC==uC(i)); end
    idxA = cell(Ma,1); for i=1:Ma, idxA{i} = find(sidA==uA(i)); end

    subjMedC = nan(Mc,1); for i=1:Mc, subjMedC(i) = median(pC(idxC{i}), "omitnan"); end
    subjMedA = nan(Ma,1); for i=1:Ma, subjMedA(i) = median(pA(idxA{i}), "omitnan"); end
    pointEstC = median(subjMedC, "omitnan");
    pointEstA = median(subjMedA, "omitnan");
    diffPoint = pointEstC - pointEstA;

    diffBoot = nan(NBoot,1);
    ctrlBoot = nan(NBoot,1);   % tracked separately (not just the difference) so the
    asdBoot  = nan(NBoot,1);   % true Welch-Satterthwaite df can be computed below,
                               % not approximated -- CTRL and ASD are resampled fully
                               % independently each iteration, so var(diff) = var(ctrl)
                               % + var(asd) exactly, letting both marginal SEs be
                               % recovered from this single joint bootstrap at no extra cost.
    for bi = 1:NBoot
        drawC = randi(Mc, Mc, 1);
        mC = nan(Mc,1);
        for i = 1:Mc
            ii = idxC{drawC(i)}; nI = numel(ii); dd = ii(randi(nI,nI,1));
            mC(i) = median(pC(dd) + sC(dd).*randn(nI,1), "omitnan");
        end
        drawA = randi(Ma, Ma, 1);
        mA = nan(Ma,1);
        for i = 1:Ma
            ii = idxA{drawA(i)}; nI = numel(ii); dd = ii(randi(nI,nI,1));
            mA(i) = median(pA(dd) + sA(dd).*randn(nI,1), "omitnan");
        end
        ctrlBoot(bi) = median(mC, "omitnan");
        asdBoot(bi)  = median(mA, "omitnan");
        diffBoot(bi) = ctrlBoot(bi) - asdBoot(bi);
    end

    seDiff = std(diffBoot, 0, "omitnan");
    seCtrl = std(ctrlBoot, 0, "omitnan");
    seAsd  = std(asdBoot,  0, "omitnan");
    ci95pct = prctile(diffBoot, [2.5 97.5]);
    ci90pct = prctile(diffBoot, [5 95]);

    % True Welch-Satterthwaite df, not an approximation: both marginal SEs are
    % directly available above (independent resampling per iteration), so the
    % standard two-sample unequal-variance df formula applies exactly as it
    % would for any two-sample t-test with unequal n and unequal variance.
    dfWelch = (seCtrl^2 + seAsd^2)^2 / ( (seCtrl^4)/(Mc-1) + (seAsd^4)/(Ma-1) );
    varCheck = abs(seDiff^2 - (seCtrl^2 + seAsd^2)) / (seCtrl^2 + seAsd^2);
    if varCheck > 0.05
        warning("compareCookCtrlAsd:varianceMismatch", "%s", sprintf( ...
            "seDiff^2 (%.6f) and seCtrl^2+seAsd^2 (%.6f) disagree by %.1f%% -- " + ...
            "independence assumption may not hold as expected; Welch df below may not be reliable.", ...
            seDiff^2, seCtrl^2+seAsd^2, 100*varCheck));
    end

    tcrit95 = tinv(0.975, dfWelch);
    tcrit90 = tinv(0.95,  dfWelch);
    ci95t = [diffPoint - tcrit95*seDiff, diffPoint + tcrit95*seDiff];
    ci90t = [diffPoint - tcrit90*seDiff, diffPoint + tcrit90*seDiff];

    diffExcludesZero = (ci95t(1) > 0) || (ci95t(2) < 0);
    equivalentWithinMDC = (ci90t(1) > -MDC) && (ci90t(2) < MDC);

    if diffExcludesZero
        verdict = "REAL DIFFERENCE (95% CI excludes 0)";
    elseif equivalentWithinMDC
        verdict = "EQUIVALENT within MDC (90% TOST CI inside +/-MDC)";
    else
        verdict = "INCONCLUSIVE -- neither excludes 0 nor fits within +/-MDC";
    end

    res = struct("label", label, "Mc", Mc, "Ma", Ma, ...
        "pointEstCtrl", pointEstC, "pointEstAsd", pointEstA, "diffPoint", diffPoint, ...
        "seDiff", seDiff, "seCtrl", seCtrl, "seAsd", seAsd, "varianceCheckPctDiff", 100*varCheck, ...
        "dfWelch", dfWelch, ...
        "ci95_percentile", ci95pct, "ci90_percentile", ci90pct, ...
        "ci95_t", ci95t, "ci90_t", ci90t, ...
        "diffExcludesZero", diffExcludesZero, "equivalentWithinMDC", equivalentWithinMDC, ...
        "verdict", verdict);

    fprintf("--- %s ---\n", label);
    fprintf("  Cook CTRL: M=%d subjects, point=%.4f, bootstrap SE=%.4f\n", Mc, pointEstC, seCtrl);
    fprintf("  Cook ASD:  M=%d subjects, point=%.4f, bootstrap SE=%.4f\n", Ma, pointEstA, seAsd);
    fprintf("  Difference (CTRL-ASD): %.4f, bootstrap SE=%.4f (Welch df=%.1f, variance-additivity check: %.1f%% off)\n", ...
        diffPoint, seDiff, dfWelch, 100*varCheck);
    fprintf("  95%% CI (t): [%.4f, %.4f]  -- excludes 0? %d\n", ci95t(1), ci95t(2), diffExcludesZero);
    fprintf("  90%% CI (t), TOST vs +/-%.3f: [%.4f, %.4f]  -- inside band? %d\n", MDC, ci90t(1), ci90t(2), equivalentWithinMDC);
    fprintf("  VERDICT: %s\n\n", verdict);
end

%% ======================== shared IO helpers =================================

function r = loadResults_local(srcDir, corpus, tag)
    if corpus == "v015", cand = "v015"; else, cand = ["v013","v014"]; end
    f = arrayfun(@(v) fullfile(srcDir, sprintf("loopClosureResults_%s_all_shaped_xu_%s.mat", tag, v)), cand);
    hit = isfile(f);
    if nnz(hit) ~= 1
        error("compareCookCtrlAsd:notFound", "%s", sprintf("FAILED PATH: %s -- %d of {%s} present, expected exactly 1.", tag, nnz(hit), strjoin(cand, ",")));
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
