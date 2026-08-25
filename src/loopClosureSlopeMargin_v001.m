function out = loopClosureSlopeMargin_v001(opts)
% loopClosureSlopeMargin_v001  E1 + E5 claim-hardening: local monotonicity margin at beta_gen*.
%
% Hardens the central claim ("canonical beta ~ 1/3 is a pipeline artefact; true beta ~ 0.30")
% against the circularity charge a tough reviewer will swing at hardest (TODO Section E, E1).
% Extends D16-6 (loopClosureMonotoneZone_v001 / Finding #87).
%
% Loads loopClosureResults_<tag>_all_shaped_xu_v007.mat for the five clinical datasets
% (Cook CTRL/ASD, Hickman PLAC/HALO, Pilot) and, per dataset, reports:
%
%   (1) MONOTONICITY MARGIN at beta_gen*  -- two independent estimators:
%       (a) Indirect (loop-closure-honest, PRIMARY).  The full forward-map cube is discarded
%           by v007 (R1), but the inversion CI geometry recovers the local slope of the actual
%           shaped_xu forward map at beta_gen*:
%                 slope_pp = (p95 - p5 of betaRecSlice(:,pp)) / (ciLo(pp) - ciHi(pp))
%           because the inversion maps the vertical beta_rec band [p5,p95] at the operating
%           point onto the horizontal beta_gen* band [ciHi,ciLo] via the inverse slope.
%           NB ciHi is the LOWER beta_gen* bound and ciLo the UPPER (naming inverted in the
%           mat: they are named for the beta_rec percentile, not the resulting beta_gen bound;
%           see invertBeta_local in runLoopClosureFftnoise_v007 and loopClosureMonotoneZone_v001).
%           Positive slope = locally monotone = inversion is locally well-posed = NOT circular.
%           This estimator is evaluated at each trial's TRUE (alpha, sigma) with the REAL
%           surrogate, and is independent of any simulation grid.
%       (b) Cross-check (generic forward map).  Finite-difference d(beta_rec)/d(beta_gen) at
%           the dataset operating point on perCoordinateSEM_v2_001 (reuses the machinery of
%           plotForwardMapAtCentroid_v001). Evaluated at the LC residual alpha (Finding #83,
%           canonical) and, as a sensitivity, the noiseChar alpha. Degraded-mode (banner, not
%           crash) if perCoordinateSEM is unavailable or its grid clips the high-alpha datasets.
%
%   (2) FORWARD-MAP-AWARE BOOTSTRAP CI on beta_gen*  -- carries forward-map (inversion)
%       uncertainty, not just trial scatter (E1).  Double bootstrap: resample trials, and
%       within each draw perturb each trial's beta_gen* by its per-trial inversion SD derived
%       from the 90% inversion CI. Reported for SG-LMLS (primary) and for betaGenStarMed
%       (the constellation-centre headline). A trial-scatter-only bootstrap is reported
%       alongside so the forward-map widening is explicit. A bootstrap overlap test against
%       1/3 is included (light precursor to E3).
%
%   (3) E5 Cook CTRL boundary decision.  Distance from beta_gen* to (i) the external 0.30
%       pure-Xu conservative bound and (ii) the data-driven shaped_xu lower bound
%       = median(ciHi_SG-LMLS) (this is D16-6's primary_margin). Plus per-trial slope
%       positivity at the boundary. Applies a Good / borderline-Good rule (Finding #87
%       upgraded Cook to Tier 2b; the rule re-derives that from the metrics).
%
% OUTPUT
%   Console tables + struct array `out` + src/loopClosureSlopeMargin_v001.mat
%
% USAGE
%   loopClosureSlopeMargin_v001()
%   loopClosureSlopeMargin_v001(NBoot=5000, RngSeed=42)
%
% Fraser, D.S. (2026)  v001

    arguments
        opts.NBoot    (1,1) double {mustBeInteger, mustBePositive} = 2000
        opts.NInner   (1,1) double {mustBeInteger, mustBePositive} = 200
        opts.RngSeed  (1,1) double {mustBeInteger}                 = 20260619
        opts.PrimaryPP(1,1) double {mustBeInRange(opts.PrimaryPP, 1, 6)} = 4
    end

    %% --- Config ----------------------------------------------------------
    PP_NAMES_V007 = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];
    PP_PRIMARY    = opts.PrimaryPP;     % 4 = SG-LMLS in the v007 loop-closure order
    N_PIPELINES   = 6;
    Z90           = 1.644853626951472; % 90% two-sided normal quantile
    PUREXU_BOUND  = 0.30;              % external pure-Xu conservative lower bound (E5)
    ONE_THIRD     = 1/3;

    % tag, label, fs (fs mirrors FS_DS in runLoopClosureFftnoise_v007)
    DATASETS = { ...
        "Cook_CTRL",    "Cook CTRL",    133; ...
        "Cook_ASD",     "Cook ASD",     133; ...
        "Hickman_PLAC", "Hickman PLAC", 133; ...
        "Hickman_HALO", "Hickman HALO", 133; ...
        "Pilot",        "Pilot",        240 };
    nDS = size(DATASETS, 1);

    srcDir  = fileparts(mfilename("fullpath"));
    outFile = fullfile(srcDir, "loopClosureSlopeMargin_v001.mat");

    rng(opts.RngSeed, "twister");

    %% --- Load the generic forward map once (degraded-mode on failure) -----
    [coordTable, fwdMapOK, fwdMapMsg] = loadForwardMap_local(srcDir);
    if ~fwdMapOK
        printBanner_local("FORWARD-MAP CROSS-CHECK UNAVAILABLE (degraded mode)", fwdMapMsg);
    end

    %% --- Per-dataset computation -----------------------------------------
    out(nDS) = struct();
    for k = 1:nDS
        tag   = DATASETS{k,1};
        label = DATASETS{k,2};
        fsDS  = DATASETS{k,3};

        matFile = fullfile(srcDir, sprintf("loopClosureResults_%s_all_shaped_xu_v007.mat", tag));
        if ~isfile(matFile)
            error("loopClosureSlopeMargin_v001:notFound", "%s", ...
                sprintf("Loop-closure mat not found (FAILED PATH): %s", matFile));
        end
        D = load(matFile, "results");
        r = D.results;
        N = numel(r);

        % --- Per-trial arrays -------------------------------------------
        bgStarMed = nan(N, 1);
        bgStarPP  = nan(N, N_PIPELINES);
        ciHiMat   = nan(N, N_PIPELINES);   % lower beta_gen* bounds
        ciLoMat   = nan(N, N_PIPELINES);   % upper beta_gen* bounds
        sliceP5   = nan(N, N_PIPELINES);   % p5  of betaRecSlice column
        sliceP95  = nan(N, N_PIPELINES);   % p95 of betaRecSlice column
        alphaLC   = nan(N, 1);             % LC residual alpha (axis mean)
        alphaNC   = nan(N, 1);             % noiseChar alpha
        sigmaMM   = nan(N, 1);

        for ti = 1:N
            bgStarMed(ti) = getfieldScalar_local(r(ti), "betaGenStarMed");
            bgStarPP(ti,:) = getRowVec_local(r(ti), "betaGenStar", N_PIPELINES);
            ciHiMat(ti,:)  = getRowVec_local(r(ti), "ciHi",        N_PIPELINES);
            ciLoMat(ti,:)  = getRowVec_local(r(ti), "ciLo",        N_PIPELINES);

            sl = r(ti).betaRecSlice;       % N_REPS x 6
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

            aMaj = getfieldScalar_local(r(ti), "alphaMaj");
            aMin = getfieldScalar_local(r(ti), "alphaMin");
            alphaLC(ti) = mean([aMaj, aMin], "omitnan");
            alphaNC(ti) = getfieldScalar_local(r(ti), "alphaIRA");
            sigmaMM(ti) = getfieldScalar_local(r(ti), "sigmaMM");
        end

        % --- (1a) Indirect slope: slope_pp = vSpread / ciWidth ----------
        ciWidth   = ciLoMat - ciHiMat;                 % upper - lower (>0 if monotone)
        vSpread   = sliceP95 - sliceP5;                % >= 0
        slopeInd  = nan(N, N_PIPELINES);
        okSlope   = isfinite(ciWidth) & isfinite(vSpread) & (abs(ciWidth) > eps);
        slopeInd(okSlope) = vSpread(okSlope) ./ ciWidth(okSlope);

        % --- (1b) Generic-map finite-difference slope at operating point -
        medAlphaLC = median(alphaLC, "omitnan");
        medAlphaNC = median(alphaNC, "omitnan");
        medSigma   = median(sigmaMM, "omitnan");
        medBgStar  = median(bgStarMed, "omitnan");
        fwdSlope   = struct("ok", false);
        if fwdMapOK
            fwdSlope = fwdMapSlope_local(coordTable, medAlphaLC, medAlphaNC, medSigma, ...
                fsDS, medBgStar, PP_NAMES_V007);
        end

        % --- (2) Forward-map-aware bootstrap CIs ------------------------
        % Per-trial inversion SD from the 90% CI (width / (2*z90)).
        sdInvPP = ciWidth / (2*Z90);                   % N x 6; may be <0/NaN -> cleaned in boot
        % Per-trial inversion SD on the median-of-pipelines (inner MC).
        sdInvMed = innerMedianSD_local(bgStarPP, sdInvPP, opts.NInner);

        bootPrim = bootstrapBeta_local(bgStarPP(:,PP_PRIMARY), sdInvPP(:,PP_PRIMARY), ...
            opts.NBoot, ONE_THIRD);
        bootMed  = bootstrapBeta_local(bgStarMed, sdInvMed, opts.NBoot, ONE_THIRD);

        % --- (3 / E5) Boundary metrics (data-driven shaped_xu LB) -------
        finP    = isfinite(bgStarMed) & isfinite(ciHiMat(:,PP_PRIMARY));
        shapedLB = median(ciHiMat(finP, PP_PRIMARY), "omitnan");   % = D16-6 primary ciHi median
        marginPureXu = medBgStar - PUREXU_BOUND;
        marginShaped = medBgStar - shapedLB;                       % = D16-6 primary_margin
        pctMonoPrim  = 100 * mean(slopeInd(isfinite(slopeInd(:,PP_PRIMARY)), PP_PRIMARY) > 0);
        nUndefPrim   = sum(~isfinite(slopeInd(:,PP_PRIMARY)) & finP);

        % --- Pack -------------------------------------------------------
        out(k).label        = label;
        out(k).tag          = tag;
        out(k).N_total      = N;
        out(k).fs           = fsDS;
        out(k).medAlphaLC   = medAlphaLC;
        out(k).medAlphaNC   = medAlphaNC;
        out(k).medSigma     = medSigma;
        out(k).medBgStar    = medBgStar;

        out(k).slopeInd_med  = median(slopeInd, 1, "omitnan");      % 1 x 6
        out(k).slopeInd_prim = median(slopeInd(:,PP_PRIMARY), "omitnan");
        out(k).pctMono_prim  = pctMonoPrim;
        out(k).nUndef_prim   = nUndefPrim;
        out(k).fwdSlope      = fwdSlope;

        out(k).boot_primary  = bootPrim;
        out(k).boot_median   = bootMed;

        out(k).shapedLB      = shapedLB;
        out(k).marginPureXu  = marginPureXu;
        out(k).marginShaped  = marginShaped;
    end

    %% --- E5 Cook CTRL verdict (apply rule) -------------------------------
    cookIdx = find([out.label] == "Cook CTRL", 1);
    if ~isempty(cookIdx)
        c = out(cookIdx);
        % Good if: above shaped_xu LB with margin >= MDC (0.03), slope positive (median>0),
        % and per-trial monotone fraction >= 90%. Else borderline-Good.
        isGood = (c.marginShaped >= 0.03) && (c.slopeInd_prim > 0) && (c.pctMono_prim >= 90);
        if isGood
            cookVerdict = "GOOD (Tier 2b retained)";
        else
            cookVerdict = "BORDERLINE-GOOD (downgrade)";
        end
        out(cookIdx).e5_verdict = cookVerdict;
    end

    %% --- Print: monotonicity margin --------------------------------------
    fprintf("\n%s\n", repmat("=", 1, 78));
    fprintf("E1  LOCAL MONOTONICITY MARGIN at beta_gen*  (SG-LMLS primary, pp=%d)\n", PP_PRIMARY);
    fprintf("%s\n", repmat("=", 1, 78));
    fprintf("Indirect slope = (p95-p5 betaRecSlice) / (ciLo-ciHi);  +ve => locally monotone.\n");
    fprintf("Evaluated at each trial's true (alpha,sigma) with the real shaped_xu surrogate.\n\n");
    fprintf("%-14s  %8s  %9s  %7s  %7s  %9s  %9s\n", ...
        "Dataset", "bg*_med", "slope_ind", "%mono", "nUndef", "fwdSlpLC", "fwdSlpNC");
    fprintf("%s\n", repmat("-", 1, 78));
    for k = 1:nDS
        [sLC, sNC] = fwdSlopePrimary_local(out(k).fwdSlope, "SG-LMLS");
        fprintf("%-14s  %8.4f  %9.3f  %6.1f%%  %7d  %9s  %9s\n", ...
            out(k).label, out(k).medBgStar, out(k).slopeInd_prim, ...
            out(k).pctMono_prim, out(k).nUndef_prim, fmtSlope_local(sLC), fmtSlope_local(sNC));
    end
    fprintf("\nfwdSlpLC / fwdSlpNC = generic-map finite-diff slope at LC / noiseChar alpha");
    fprintf(" ('--' = unavailable/clipped).\n");

    %% --- Print: all-pipeline indirect slope ------------------------------
    fprintf("\n--- Indirect slope (median over trials), all six pipelines ---\n");
    fprintf("%-14s", "Dataset");
    for pp = 1:N_PIPELINES, fprintf("  %9s", PP_NAMES_V007(pp)); end
    fprintf("\n%s\n", repmat("-", 1, 14 + N_PIPELINES*11));
    for k = 1:nDS
        fprintf("%-14s", out(k).label);
        for pp = 1:N_PIPELINES, fprintf("  %9.3f", out(k).slopeInd_med(pp)); end
        fprintf("\n");
    end

    %% --- Print: forward-map-aware bootstrap CI ---------------------------
    fprintf("\n%s\n", repmat("=", 1, 78));
    fprintf("E1  FORWARD-MAP-AWARE BOOTSTRAP CI on beta_gen*  (NBoot=%d)\n", opts.NBoot);
    fprintf("%s\n", repmat("=", 1, 78));
    fprintf("FM-aware carries per-trial inversion SD (from 90%% CI); trial-only does not.\n");
    fprintf("p(>=1/3) = bootstrap fraction of medians at/above 1/3 (overlap test, precursor to E3).\n\n");
    printBootHeader_local();
    for k = 1:nDS
        printBootRow_local(out(k).label + " [SG-LMLS]", out(k).boot_primary);
        printBootRow_local(out(k).label + " [median]",  out(k).boot_median);
    end

    %% --- Print: E5 Cook CTRL boundary ------------------------------------
    fprintf("\n%s\n", repmat("=", 1, 78));
    fprintf("E5  COOK CTRL BOUNDARY DECISION\n");
    fprintf("%s\n", repmat("=", 1, 78));
    for k = 1:nDS
        fprintf("%-14s  bg*=%.4f  shaped_xuLB=%.4f  d(0.30)=%+.4f  d(LB)=%+.4f  slope=%+.3f\n", ...
            out(k).label, out(k).medBgStar, out(k).shapedLB, ...
            out(k).marginPureXu, out(k).marginShaped, out(k).slopeInd_prim);
    end
    if ~isempty(cookIdx)
        fprintf("\nd(0.30) = bg* - 0.30 (external pure-Xu bound; informational, Finding #87).\n");
        fprintf("d(LB)   = bg* - median(ciHi_SG-LMLS) = data-driven shaped_xu lower bound\n");
        fprintf("          (this is D16-6 primary_margin; the operative bound per Finding #87).\n");
        fprintf("\nCOOK CTRL VERDICT: %s\n", out(cookIdx).e5_verdict);
        fprintf("  rule: GOOD if d(LB) >= 0.03 AND slope > 0 AND %%mono >= 90.\n");
    end

    %% --- Save ------------------------------------------------------------
    meta = struct("NBoot", opts.NBoot, "NInner", opts.NInner, "RngSeed", opts.RngSeed, ...
        "PP_PRIMARY", PP_PRIMARY, "pipelineNames", PP_NAMES_V007, ...
        "fwdMapAvailable", fwdMapOK, "fwdMapMsg", string(fwdMapMsg), ...
        "generated", string(datetime("now")));
    save(outFile, "out", "meta", "-v7.3");
    fprintf("\nSaved: %s\n", outFile);
end

%% ========================  local functions  ============================

function v = getfieldScalar_local(s, fld)
% Return a finite scalar field or NaN.
    v = NaN;
    if isfield(s, fld)
        x = s.(fld);
        if isnumeric(x) && isscalar(x), v = double(x); end
    end
end

function row = getRowVec_local(s, fld, nExpected)
% Return a 1 x nExpected row from a field of length nExpected, else NaNs.
    row = nan(1, nExpected);
    if isfield(s, fld)
        x = s.(fld);
        if isnumeric(x) && numel(x) == nExpected
            row = double(x(:)');
        end
    end
end

function [T, ok, msg] = loadForwardMap_local(srcDir)
% Load perCoordinateSEM_v2_001 coordTable; return ok=false with a reason if unavailable.
    T = table(); ok = false; msg = "";
    f = fullfile(srcDir, "perCoordinateSEM_v2_001.mat");
    if ~isfile(f)
        msg = sprintf("perCoordinateSEM_v2_001.mat not found in %s", srcDir);
        return
    end
    try
        S = load(f, "coordTable");
    catch ME
        msg = sprintf("load failed: %s", ME.message);
        return
    end
    if ~isfield(S, "coordTable")
        msg = "coordTable variable absent from perCoordinateSEM_v2_001.mat";
        return
    end
    T = S.coordTable;
    needed = ["alpha","sigma","fs","betaGen","pipeline","meanBetaRec"];
    have   = string(T.Properties.VariableNames);
    miss   = needed(~ismember(needed, have));
    if ~isempty(miss)
        msg = sprintf("coordTable missing columns: %s", strjoin(miss, ", "));
        T = table(); return
    end
    ok = true;
end

function fs = fwdMapSlope_local(T, alphaLC, alphaNC, sigmaMM, fsTgt, bgStar, ppNames)
% Generic-map finite-difference slope d(beta_rec)/d(beta_gen) at beta_gen=bgStar, per pipeline,
% at the LC alpha and the noiseChar alpha. Snaps (alpha,sigma,fs) to grid and warns on clip.
    fs = struct("ok", true);
    aGrid = sort(unique(T.alpha));
    sGrid = sort(unique(T.sigma));
    fGrid = sort(unique(T.fs));

    [snapA_LC, clipA_LC] = snapWarn_local(aGrid, alphaLC);
    [snapA_NC, clipA_NC] = snapWarn_local(aGrid, alphaNC);
    [snapS, ~]           = snapWarn_local(sGrid, sigmaMM);
    [snapF, ~]           = snapWarn_local(fGrid, fsTgt);

    fs.snapAlphaLC = snapA_LC; fs.clipAlphaLC = clipA_LC;
    fs.snapAlphaNC = snapA_NC; fs.clipAlphaNC = clipA_NC;
    fs.snapSigma   = snapS;    fs.snapFs = snapF;
    fs.slopeLC = struct(); fs.slopeNC = struct();

    if clipA_LC
        fprintf("  [fwdmap] LC alpha=%.2f snapped to grid alpha=%.2f (CLIP by %.2f; slope grid-limited)\n", ...
            alphaLC, snapA_LC, abs(alphaLC - snapA_LC));
    end

    for pi = 1:numel(ppNames)
        nm = ppNames(pi);
        fs.slopeLC.(matlab.lang.makeValidName(nm)) = ...
            curveSlope_local(T, snapA_LC, snapS, snapF, nm, bgStar);
        fs.slopeNC.(matlab.lang.makeValidName(nm)) = ...
            curveSlope_local(T, snapA_NC, snapS, snapF, nm, bgStar);
    end
end

function [snapV, clipped] = snapWarn_local(grid, tgt)
% Nearest grid node; clipped=true if the snap distance exceeds 1.5x the local grid spacing.
    [~, idx] = min(abs(grid - tgt));
    snapV = grid(idx);
    if numel(grid) >= 2
        dns = diff(grid);
        localSpacing = median(dns);
        clipped = abs(tgt - snapV) > 1.5*localSpacing;
    else
        clipped = false;
    end
end

function s = curveSlope_local(T, snapA, snapS, snapF, pipeName, bgStar)
% Marginalise VGF (mean meanBetaRec), build beta_gen->beta_rec curve, local-linear slope at bgStar.
    sub = T(T.alpha == snapA & T.sigma == snapS & T.fs == snapF & T.pipeline == pipeName, :);
    if isempty(sub), s = NaN; return; end
    bg = sort(unique(sub.betaGen));
    br = nan(numel(bg), 1);
    for bi = 1:numel(bg)
        rr = sub(sub.betaGen == bg(bi), :);
        br(bi) = mean(rr.meanBetaRec, "omitnan");
    end
    fin = isfinite(bg) & isfinite(br);
    bg = bg(fin); br = br(fin);
    if numel(bg) < 3, s = NaN; return; end
    % Local linear fit over the three grid nodes nearest bgStar.
    [~, ord] = sort(abs(bg - bgStar));
    sel = sort(ord(1:min(3, numel(ord))));
    P = polyfit(bg(sel), br(sel), 1);
    s = P(1);
end

function [sLC, sNC] = fwdSlopePrimary_local(fwdSlope, pipeName)
% Pull the primary pipeline's LC/noiseChar finite-diff slope from the fwdSlope struct.
    sLC = NaN; sNC = NaN;
    if isfield(fwdSlope, "ok") && fwdSlope.ok
        key = matlab.lang.makeValidName(pipeName);
        if isfield(fwdSlope.slopeLC, key), sLC = fwdSlope.slopeLC.(key); end
        if isfield(fwdSlope.slopeNC, key), sNC = fwdSlope.slopeNC.(key); end
    end
end

function sdMed = innerMedianSD_local(bgStarPP, sdInvPP, NInner)
% Per-trial inversion SD on the median-of-converging-pipelines beta_gen*, by inner Monte-Carlo.
% Each converging pipeline's beta_gen* is drawn from N(point, inversionSD); the median is
% re-taken NInner times; the SD across draws is that trial's forward-map uncertainty on the median.
    N = size(bgStarPP, 1);
    sdMed = zeros(N, 1);
    for ti = 1:N
        pts = bgStarPP(ti, :);
        sds = sdInvPP(ti, :);
        conv = isfinite(pts);
        if ~any(conv), sdMed(ti) = 0; continue; end
        p = pts(conv);
        d = sds(conv);
        d(~isfinite(d) | d < 0) = 0;        % undefined CI -> point mass (counted elsewhere)
        if all(d == 0), sdMed(ti) = 0; continue; end
        draws = p(:)' + d(:)' .* randn(NInner, numel(p));
        meds  = median(draws, 2);
        sdMed(ti) = std(meds, 0);
    end
end

function b = bootstrapBeta_local(point, sdInv, NBoot, refVal)
% Double bootstrap on the dataset median of `point`.
%   FM-aware : resample trials, perturb each by N(0,sdInv), take median.
%   trial-only: resample trials, take median (no perturbation).
% refVal: value to test overlap against (1/3). Returns CIs and bootstrap p(>=refVal).
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
    b.pGE     = mean(medFM >= refVal);   % overlap test against 1/3 (FM-aware distribution)
end

function printBootHeader_local()
    fprintf("%-22s  %5s  %8s  %20s  %20s  %8s\n", ...
        "Dataset [estimator]", "N", "median", "FM-aware 95% CI", "trial-only 95% CI", "p(>=1/3)");
    fprintf("%s\n", repmat("-", 1, 92));
end

function printBootRow_local(name, b)
    fprintf("%-22s  %5d  %8.4f  [%7.4f, %7.4f]  [%7.4f, %7.4f]  %8.4f\n", ...
        name, b.N, b.median, b.ciFM(1), b.ciFM(2), b.ciTrial(1), b.ciTrial(2), b.pGE);
end

function str = fmtSlope_local(v)
% Format a finite-diff slope, or '--' if unavailable.
    if isfinite(v), str = sprintf("%+.3f", v); else, str = "--"; end
end

function printBanner_local(title, body)
    fprintf("\n%s\n", repmat("!", 1, 78));
    fprintf("!! %s\n", title);
    if strlength(string(body)) > 0
        fprintf("!! %s\n", body);
    end
    fprintf("%s\n", repmat("!", 1, 78));
end
