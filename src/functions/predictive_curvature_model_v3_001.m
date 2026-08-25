function out = predictive_curvature_model_v3_001(x, y, fs, params)
% PREDICTIVE_CURVATURE_MODEL_V3_001  Vestibular pq noise model with Menger curvature.
%
% V3 architecture (2026-06-08). Key changes from v2_005:
%   PRIMITIVES: a_c = v^2 * kappa_M (Type II, centripetal) and
%               j_path = v^2 * slope(kappa_M over preview) (Type I, prospective jerk).
%               kappa_M from Menger positional triplets (zero differentiations).
%               v from differentiateKinematicsEBR filterType 6 (one differentiation).
%               Total differentiation budget: 1 (vs 2 in v2_005).
%
%   SPEED:      speedPred = vmax * kappaEff^(-beta), beta=1/3 default.
%               Pipeline-honest: 1/3 reflects pipeline compression of biological
%               beta_gen* ~0.4-0.5 by SG differentiation. See TODO_VestibularModel_v001
%               for the 1/3 vs 0.5 architectural rationale and upgrade path to v3_002
%               (Wiener filter, beta emergent) and v3_003 (PLN2020, zero differentiations).
%
%   RELIABILITY: 1 / (1 + noiseGain * acNorm^p * abs(jNorm)^q)
%               p weights Type II (centripetal acc channel, Jones et al. 2002).
%               q weights Type I (jerk channel, Manzari 2026; Bockisch & Haslwanter 2007).
%               DA prediction:  q_PLAC > q_HALO (reduced DA -> less anticipatory braking).
%               ASD prediction: q_ASD  < q_CTRL (degraded Type I, Mansour 2021).
%
%   HORIZON:    kappaEff integrates kappa_M over [KMemory past, KPreview future] with
%               exponential weighting. Future window uses Menger directly (exact; no
%               sensory/imagination fusion needed for k <= KPreview).
%               j_path is the exponentially-weighted average slope of kappa_M over the
%               preview window, dual-use: (1) Type I reliability channel;
%               (2) forward model extrapolation beyond KPreview (reserved, v3_002).
%
%   ADAPT RATE: drives EMA of total pq noise (acNorm^p * abs(jNorm)^q), expanding
%               the horizon when vestibular demand is sustained.
%
%   REMOVED vs v2_005: goalDecay/goalDecayMode (D2 null, Finding #60),
%               volatilityGain/omegaEMA (D15, subsumed into j channel),
%               sigmaEMA/reliabilityGain (replaced by pq model),
%               sensoryDecay/forwardModelDecay (no channel fusion in preview window),
%               ax/ay usage (Menger replaces acceleration-based curvature).
%
% INPUTS
%   x, y    Position vectors (mm recommended; any consistent unit). Column or row.
%   fs      Sampling rate in Hz (> 0).
%   params  Name-value pairs (see arguments block).
%
% OUTPUT  (single struct; 4-output cap)
%   out.speedPred    N-by-1 predicted speed (NaN at unmodellable edges).
%   out.kappaEff     N-by-1 horizon-integrated Menger curvature.
%   out.kappaM       N-by-1 Menger curvature (NaN at edges of width kMenger).
%   out.acEff        N-by-1 centripetal demand: speedNow^2 * kappaEff.
%   out.acRaw        N-by-1 instantaneous centripetal demand: speedNow^2 * kappaM.
%   out.jPath        N-by-1 prospective jerk: v^2 * weighted-avg slope(kappa_M).
%   out.acNorm       N-by-1 soft-normalised centripetal demand, in (0,1).
%   out.jNorm        N-by-1 soft-normalised jerk demand, in (0,1).
%   out.reliability  N-by-1 pq noise-derived speed gain, in (0,1].
%   out.lambdaFuture N-by-1 adaptive future decay rate (1/sample).
%   out.lambdaPast   N-by-1 adaptive past decay rate (1/sample).
%   out.alphaFuture  N-by-1 adaptive horizon state.
%   out.fit          Struct: betaRecovered, yGain, r2, nUsed, method, recoveryHalfLife.
%   out.diff         Struct: differentiation method and diagnostics.
%   out.params       Resolved parameter struct.
%
% FAIL-LOUD CONTRACT (project policy)
%   - Differentiation errors: hard-fail unless allowFallbackDiff=true (disclosed banner).
%   - Menger edges: NaN-padded, never zero-filled.
%   - kappaFloor hit rate > 5%: warning, not silent.
%   - Degenerate acScale (all-zero centripetal demand): hard error.
%   - jScale fallback (cannot estimate from proxy): disclosed warning.
%
% NOTE ON BETA RECOVERY
%   out.fit.betaRecovered is a self-consistency check only. With beta=1/3 input and
%   noiseGain>0, the recovered value will deviate due to reliability modulation.
%   Set noiseGain=0 to isolate pure power-law recovery. Any empirical claim about
%   the biological beta requires v058 loop-closure validation (Finding #53).
%
% See also:
%   predictive_curvature_model_v2_005.m  (predecessor)
%   TODO_VestibularModel_v001.md         (architectural decisions, upgrade path)
%   README_VestibularNoise_v001.md       (vestibular pq noise model)
%   differentiateKinematicsEBR.m         (Fraser et al. 2025)
%   regressDataEBR.m                     (Fraser et al. 2025)
%
% D.S. Fraser, University of Birmingham, 2026-06-08

    arguments
        x  double
        y  double
        fs (1,1) double {mustBePositive}
        params.KPreview      (1,1) double {mustBeInteger, mustBePositive}          = 48
        params.KMemory       (1,1) double {mustBeInteger, mustBePositive}          = 48
        params.beta          (1,1) double {mustBeNonnegative}                      = 1/3
        params.vmax          (1,1) double {mustBePositive}                         = 1.0
        params.lambda0       (1,1) double {mustBeNonnegative}                      = 0.035
        params.adaptRate     (1,1) double {mustBeInRange(params.adaptRate, 0, 1)}  = 0.05
        params.alpha0        (1,1) double {mustBeNonnegative}                      = 0.0
        params.memoryGain    (1,1) double {mustBeNonnegative}                      = 0.0
        params.p             (1,1) double {mustBeNonnegative}                      = 2.0
        params.q             (1,1) double {mustBeNonnegative}                      = 1.0
        params.noiseGain     (1,1) double {mustBeNonnegative}                      = 1.0
        params.kappaFloor    (1,1) double {mustBePositive}                         = 1e-4
        params.acFloor       (1,1) double {mustBePositive}                         = 1e-6
        params.kMenger       (1,1) double {mustBeInteger, mustBePositive}          = 1
        params.allowFallbackDiff (1,1) logical                                     = false
        params.filterType    (1,1) double {mustBeInteger, mustBePositive}          = 6
        params.filterParams  (1,:) double                                          = [4, 17]
        params.regressType   (1,1) double {mustBeInteger, mustBePositive}          = 3
    end

    %% ── Input normalisation and validation ───────────────────────────────
    x = x(:);
    y = y(:);
    if numel(x) ~= numel(y)
        error("PCMv3:lengthMismatch", "%s", "x and y must have the same length.");
    end
    N      = numel(x);
    kLag   = params.kMenger;
    minLen = params.KMemory + params.KPreview + 1 + 2*kLag;
    if N <= minLen
        error("PCMv3:tooShort", "%s", ...
            sprintf("Trajectory too short: N=%d, minimum=%d " + ...
            "(KMemory+KPreview+1+2*kMenger).", N, minLen + 1));
    end
    if any(~isfinite(x)) || any(~isfinite(y))
        error("PCMv3:nonFinitePosition", "%s", ...
            "x and y must be finite (no NaN/Inf).");
    end

    epsP = 1e-9;

    %% ── Velocity via EBR ─────────────────────────────────────────────────
    % ax, ay are returned by computeKinematics but not used in v3_001.
    % Menger curvature replaces acceleration-based curvature entirely.
    [vx, vy, ~, ~, diffInfo] = computeKinematics(x, y, fs, ...
        params.filterType, params.filterParams, params.allowFallbackDiff);
    speedNow = hypot(vx, vy);                   % N-by-1

    %% ── Menger curvature (zero additional differentiations) ──────────────
    % kappaM_raw: length N-2*kLag, valid positions kLag+1 : N-kLag.
    % Padded to N with NaN so direct indexing kappaM(t) works in main loop.
    kappaM_raw = mengerCurvature(x, y, kLag);   % vectorised, length N-2*kLag
    kappaM     = NaN(N, 1);
    kappaM(kLag + 1 : N - kLag) = kappaM_raw;

    %% ── Instantaneous centripetal demand ─────────────────────────────────
    acRaw = speedNow.^2 .* kappaM;              % NaN at Menger edges; mm/s^2

    %% ── Normalisation scales (trial-level, computed once) ────────────────
    % acScale: median centripetal demand, sets soft-normalisation half-saturation.
    % jScale:  proxy from 1-step Menger slope * median v^2; same units as jPath.
    acValid = acRaw(isfinite(acRaw) & acRaw > 0);
    if isempty(acValid)
        error("PCMv3:degenerateTrajectory", "%s", ...
            "acRaw is all-NaN or all-zero. Check position units and trajectory.");
    end
    acScale = median(acValid);
    if ~(acScale > 0)
        error("PCMv3:degenerateTrajectory", "%s", ...
            "acScale is zero; centripetal demand is degenerate.");
    end

    vMed      = median(speedNow(kLag + 1 : N - kLag));
    dkProxy   = [0; diff(kappaM_raw)] * fs;     % 1-step slope of kappa_M, (1/mm)/s
    jProxy    = abs(dkProxy) * vMed.^2;         % proxy jPath scale, mm/s^3
    jValid    = jProxy(isfinite(jProxy) & jProxy > 0);
    if isempty(jValid)
        jScale = acScale * fs;
        warning("PCMv3:jScaleFallback", "%s", ...
            "jScale proxy is all-zero; using acScale*fs fallback.");
    else
        jScale = median(jValid);
    end

    %% ── Modellable window ────────────────────────────────────────────────
    % tStart: KMemory back from t must not underrun Menger valid range.
    % tEnd:   KPreview ahead from t must not overrun Menger valid range.
    tStart = params.KMemory + kLag + 1;
    tEnd   = N - params.KPreview - kLag;
    if tEnd <= tStart
        error("PCMv3:windowCollapse", "%s", ...
            sprintf("Modellable window is empty (tStart=%d, tEnd=%d). " + ...
            "Reduce KMemory, KPreview, or kMenger, or use a longer trial.", ...
            tStart, tEnd));
    end

    %% ── Pre-allocate outputs as NaN ──────────────────────────────────────
    speedPred    = NaN(N, 1);
    kappaEff     = NaN(N, 1);
    acEff        = NaN(N, 1);
    jPath        = NaN(N, 1);
    acNormTr     = NaN(N, 1);
    jNormTr      = NaN(N, 1);
    reliability  = NaN(N, 1);
    lambdaFutTr  = NaN(N, 1);
    lambdaPastTr = NaN(N, 1);
    alphaFutTr   = NaN(N, 1);

    %% ── Online states ────────────────────────────────────────────────────
    alphaFuture = params.alpha0;
    sigmaPQ_prev = 0.0;    % previous-step total pq noise; drives adaptRate
    floorHits   = 0;

    %% ── Main loop ────────────────────────────────────────────────────────
    for t = tStart:tEnd

        % Update adaptive horizon from previous step's pq noise (no circularity).
        alphaFuture  = alphaFuture + ...
            params.adaptRate * (sigmaPQ_prev - alphaFuture);
        lambdaFuture = params.lambda0 + alphaFuture * sigmaPQ_prev;
        lambdaPast   = params.lambda0 + params.memoryGain * sigmaPQ_prev;

        % --- Integrate kappa_M over past, present, and future ---
        C     = kappaM(t);    % present sample anchors the window
        W     = 1.0;
        jAcc  = 0.0;          % accumulate weighted average slope of kappa_M
        wAcc  = 0.0;

        % Retrospective memory: true past curvature, no forward model.
        for k = 1:params.KMemory
            wk = exp(-lambdaPast * k);
            C  = C + wk * kappaM(t - k);
            W  = W + wk;
        end

        % Prospective preview: kappa_M(t+k) is exact from Menger (no channel
        % fusion needed; all preview positions are available offline).
        % j_path accumulates the k-step average slope of kappa_M, weighted by
        % the same exponential kernel as the horizon integration.
        for k = 1:params.KPreview
            wk   = exp(-lambdaFuture * k);
            kapK = kappaM(t + k);
            C    = C + wk * kapK;
            W    = W + wk;
            % Average rate of curvature change from t to t+k: (Δkappa)/(k/fs).
            % Weighted average over k gives prospective dkappa/dt estimate.
            jAcc = jAcc + wk * (kapK - kappaM(t)) * fs / k;
            wAcc = wAcc + wk;
        end

        kappaEff_t = C / W;

        % j_path: prospective curvature jerk scaled by v(t)^2.
        % Units: (mm/s)^2 * (1/mm)/s = mm/s^3 (centripetal jerk units).
        % Captures kappa_dot contribution to d(a_c)/dt; v_dot contribution
        % is secondary on smooth power-law trajectories near corners.
        jPath_t = speedNow(t).^2 * jAcc / max(wAcc, epsP);

        % Soft-normalised noise signals, both in (0,1).
        acEff_t  = speedNow(t).^2 * kappaEff_t;
        acNorm_t = acEff_t           / (acEff_t           + acScale);
        jNorm_t  = abs(jPath_t)      / (abs(jPath_t)      + jScale);

        % pq reliability: 1 when noise is zero, -> 0 as noise -> inf.
        % At the D16-C optimum (beta=0.5, zero-jerk): jPath -> 0 so jNorm -> 0.
        % The q term then collapses; Type I vestibular channel is silent.
        reliability_t = 1 / (1 + params.noiseGain * ...
            acNorm_t^params.p * abs(jNorm_t)^params.q);

        % Speed prediction.
        kEff = max(kappaEff_t, params.kappaFloor);
        if kappaEff_t < params.kappaFloor
            floorHits = floorHits + 1;
        end
        speedPred(t) = params.vmax * kEff^(-params.beta) * reliability_t;

        % Update pq state for next step's horizon computation.
        sigmaPQ_prev = acNorm_t^params.p * abs(jNorm_t)^params.q;

        % Store traces.
        kappaEff(t)     = kappaEff_t;
        acEff(t)        = acEff_t;
        jPath(t)        = jPath_t;
        acNormTr(t)     = acNorm_t;
        jNormTr(t)      = jNorm_t;
        reliability(t)  = reliability_t;
        lambdaFutTr(t)  = lambdaFuture;
        lambdaPastTr(t) = lambdaPast;
        alphaFutTr(t)   = alphaFuture;
    end

    %% ── Curvature-floor diagnostic ───────────────────────────────────────
    nModelled = tEnd - tStart + 1;
    floorFrac = floorHits / max(nModelled, 1);
    if floorFrac > 0.05
        warning("PCMv3:floorRate", "%s", ...
            sprintf("kappaFloor hit on %.1f%% of modelled samples; " + ...
            "recovered exponent may be biased. Lower kappaFloor or inspect geometry.", ...
            100 * floorFrac));
    end

    %% ── Recover exponent ─────────────────────────────────────────────────
    fit = fitPowerLaw(kappaEff, speedPred, params.kappaFloor, ...
        params.regressType, params.beta);
    if params.adaptRate > 0 && params.adaptRate < 1
        fit.recoveryHalfLife = log(0.5) / log(1 - params.adaptRate) / fs;
    elseif params.adaptRate == 0
        fit.recoveryHalfLife = Inf;
    else
        fit.recoveryHalfLife = 1 / fs;
    end

    %% ── Assemble output ──────────────────────────────────────────────────
    diffInfo.floorFraction = floorFrac;

    out.speedPred    = speedPred;
    out.kappaEff     = kappaEff;
    out.kappaM       = kappaM;
    out.acEff        = acEff;
    out.acRaw        = acRaw;
    out.jPath        = jPath;
    out.acNorm       = acNormTr;
    out.jNorm        = jNormTr;
    out.reliability  = reliability;
    out.lambdaFuture = lambdaFutTr;
    out.lambdaPast   = lambdaPastTr;
    out.alphaFuture  = alphaFutTr;
    out.fit          = fit;
    out.diff         = diffInfo;
    out.params       = params;
end


%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function kM = mengerCurvature(x, y, K)
% MENGERCURVATURE  Vectorised Menger curvature at lag K.
%   Triplets: x[j-K], x[j], x[j+K] for j = K+1 : N-K.
%   Returns vector of length N-2*K. Units: 1/(units of x,y).
%   Zero curvature where triplet area is below numerical threshold (collinear).
%   Verified: returns 1/R for equidistant points on a circle of radius R.
    ax = x(1:end-2*K);   bx = x(K+1:end-K);   cx = x(2*K+1:end);
    ay = y(1:end-2*K);   by = y(K+1:end-K);   cy = y(2*K+1:end);
    area2 = abs((bx-ax).*(cy-ay) - (cx-ax).*(by-ay));  % = 2 * triangle area
    d12   = hypot(bx-ax, by-ay);
    d23   = hypot(cx-bx, cy-by);
    d13   = hypot(cx-ax, cy-ay);
    den   = d12 .* d23 .* d13;
    kM    = zeros(size(ax));
    valid = den > 1e-12;
    kM(valid) = 2 * area2(valid) ./ den(valid);  % = 4*Area/(d12*d23*d13)
end


function [vx, vy, ax, ay, info] = computeKinematics(x, y, fs, ...
        filterType, filterParams, allowFallback)
% COMPUTEKINEMATICS  Velocity via differentiateKinematicsEBR (Fraser et al. 2025).
%   ax, ay returned for signature compatibility but unused in v3_001 (Menger
%   replaces acceleration-based curvature). Hard-fail if EBR absent and
%   allowFallback=false; disclosed central-difference banner if fallback fires.
%
%   differentiateKinematicsEBR signature (verified 2026-05-31):
%     [dx, dy] = differentiateKinematicsEBR(x, y, filterType, filterParams, fs)
%     dx, dy each N-by-4: col1=position, col2=velocity, col3=acceleration, col4=jerk.
    info = struct();
    if exist("differentiateKinematicsEBR", "file") == 2
        try
            [dx, dy]    = differentiateKinematicsEBR(x, y, filterType, filterParams, fs);
            vx          = dx(:, 2);
            vy          = dy(:, 2);
            ax          = dx(:, 3);
            ay          = dy(:, 3);
            info.method = sprintf("differentiateKinematicsEBR(filterType=%d)", filterType);
            fprintf("[PCMv3] Differentiation: %s (Fraser et al. 2025).\n", info.method);
            return
        catch ME
            baseME = MException("PCMv3:ebrFailed", "%s", ...
                sprintf("differentiateKinematicsEBR failed (filterType=%d). " + ...
                "Check filterParams are valid for this filterType.", filterType));
            rethrow(addCause(baseME, ME));
        end
    end
    if ~allowFallback
        error("PCMv3:noDifferentiator", "%s", ...
            "differentiateKinematicsEBR not found. Add src/functions/ to path, " + ...
            "or set allowFallbackDiff=true for the disclosed central-difference fallback.");
    end
    dt = 1 / fs;
    vx = gradient(x, dt);
    vy = gradient(y, dt);
    ax = gradient(vx, dt);
    ay = gradient(vy, dt);
    info.method = "fallback-central-difference";
    fprintf(2, "[PCMv3][WARNING] FALLBACK central-difference differentiation active. " + ...
        "Add src/functions/ to MATLAB path to use differentiateKinematicsEBR.\n");
end


function fit = fitPowerLaw(kappaEff, speedPred, kappaFloor, regressType, beta0)
% FITPOWERLAW  Recover exponent via regressDataEBR (Fraser et al. 2025).
%   Identical failure policy to v2_005:
%     Mode A (EBR absent): disclosed OLS fallback; fit.method = "ols-fallback".
%     Mode B (EBR fails at runtime): fit.betaRecovered = NaN, fit.failed = true.
%   Set noiseGain=0 when calling the parent function to isolate pure beta recovery.
%
%   regressDataEBR signature (verified 2026-05-31):
%     [beta, yGain, stats] = regressDataEBR(speed, kappa, regressType,
%                                            LMseeds, displayGraphs, limitBreak)
    valid = isfinite(kappaEff) & isfinite(speedPred) & ...
            (kappaEff > kappaFloor) & (speedPred > 0);

    fit              = struct();
    fit.nUsed        = nnz(valid);
    fit.failed       = false;
    fit.failureReason = "";

    if fit.nUsed < 10
        warning("PCMv3:tooFewForFit", "%s", ...
            "Fewer than 10 valid samples for power-law fit.");
        fit.betaRecovered = NaN;
        fit.yGain         = NaN;
        fit.r2            = NaN;
        fit.method        = "insufficient-data";
        return
    end

    kappaSub = kappaEff(valid);
    speedSub = speedPred(valid);
    LMseeds  = [1.0, beta0];

    if exist("regressDataEBR", "file") == 2
        try
            [betaOut, yGainOut, ~] = regressDataEBR(speedSub, kappaSub, ...
                regressType, LMseeds, 0, 1);
            fit.betaRecovered = betaOut;
            fit.yGain         = yGainOut;
            fit.method        = sprintf("regressDataEBR(type=%d)", regressType);
            logSpeed    = log(speedSub);
            logKappa    = log(kappaSub);
            logSpeedHat = log(max(yGainOut, eps)) + (-betaOut) * logKappa;
            ssRes = sum((logSpeed - logSpeedHat).^2);
            ssTot = sum((logSpeed - mean(logSpeed)).^2);
            fit.r2 = 1 - ssRes / max(ssTot, eps);
        catch ME
            fit.betaRecovered = NaN;
            fit.yGain         = NaN;
            fit.r2            = NaN;
            fit.failed        = true;
            fit.failureReason = ME.message;
            fit.method        = sprintf("failed-ebr(type=%d)", regressType);
            warning("PCMv3:ebrRegressFail", ...
                "regressDataEBR(type=%d) failed: %s", regressType, ME.message);
        end
        return
    end

    % Mode A: OLS fallback (disclosed; should not appear in sweep runs).
    logSpeed = log(speedSub);
    logKappa = log(kappaSub);
    A        = [ones(numel(logKappa), 1), logKappa];
    c        = A \ logSpeed;
    yHat     = A * c;
    ssRes    = sum((logSpeed - yHat).^2);
    ssTot    = sum((logSpeed - mean(logSpeed)).^2);
    fit.betaRecovered = -c(2);
    fit.yGain         = exp(c(1));
    fit.r2            = 1 - ssRes / max(ssTot, eps);
    fit.method        = "ols-fallback";
    warning("PCMv3:olsFallback", "%s", ...
        "regressDataEBR not on MATLAB path; OLS fallback active. " + ...
        "Add src/functions/ to path before running sweeps.");
end
