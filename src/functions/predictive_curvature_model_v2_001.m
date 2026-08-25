function out = predictive_curvature_model_v2_001(x, y, fs, params)
% PREDICTIVE_CURVATURE_MODEL_V2_001  Anticipatory speed control from integrated curvature.
%
%   out = predictive_curvature_model_v2_001(x, y, fs) runs with defaults.
%   out = predictive_curvature_model_v2_001(x, y, fs, Name, Value, ...) overrides params.
%
%   Successor to predictive_curvature_model_v3 (J. Cook). Three corrections:
%     FIX 1  The integrated curvature estimate now drives speed. In v3 the whole
%            multi-channel estimate was computed (weighted_curvature) then discarded;
%            speed used only 1/sigma and inter-sample distance.
%     FIX 2  The speed law is a POWER law in curvature, v = vmax * kappaEff^(-beta) * R,
%            not an exponential. log v = log vmax - beta*log kappaEff + log R, so the
%            generative law is commensurable with the Lacquaniti exponent beta.
%     FIX 3  Only the SENSORY channel is allowed to read the future (the template is
%            visible). IMAGINATION is a genuine forward model that extrapolates the
%            curvature trajectory from information available at time t and does not peek.
%
%   Jen's central hypothesis is preserved and implemented cleanly: the future time
%   horizon shortens as motor uncertainty rises (lambdaFuture grows with sigmaNorm).
%
%   INPUTS
%     x, y    Position vectors (any orientation; normalised to columns internally).
%     fs      Sampling rate in Hz (> 0). Required for real time units and EBR.
%     params  Name-value pairs (see arguments block) controlling horizons, the
%             exponent, uncertainty dynamics, and fallback behaviour.
%
%   OUTPUT (single struct, to respect the 4-output cap)
%     out.speedPred         N-by-1 predicted speed (NaN at unmodellable edges).
%     out.kappaEff          N-by-1 integrated curvature actually used for speed.
%     out.kappaRaw          N-by-1 instantaneous curvature from EBR kinematics.
%     out.reliability       N-by-1 uncertainty-derived gain in (0,1].
%     out.sigmaNorm         N-by-1 normalised uncertainty in [0,1).
%     out.omega             N-by-1 curvature volatility (EMA of |dkappa/dt|).
%     out.lambdaFuture      N-by-1 future temporal decay (1/sample).
%     out.lambdaPast        N-by-1 past temporal decay (1/sample).
%     out.sensoryWeight     N-by-1 future-fusion weight on sensory preview.
%     out.imaginationWeight N-by-1 future-fusion weight on the forward model.
%     out.fit               Struct: betaRecovered, intercept, r2, nUsed, method.
%     out.diff              Struct: method banner and floor-hit diagnostics.
%     out.params            The resolved parameter struct (provenance).
%
%   FAIL-LOUD CONTRACT (per project policy)
%     - Differentiation uses differentiateKinematicsEBR when wired; otherwise a
%       disclosed central-difference fallback with a printed banner, gated by
%       allowFallbackDiff. With allowFallbackDiff=false and EBR absent, it ERRORS.
%     - Edges that cannot be modelled are left as NaN, never zero-filled (zeros would
%       masquerade as genuine slow movement).
%     - A high curvature-floor hit rate raises a warning rather than silently biasing
%       the recovered exponent.

    arguments
        x  double
        y  double
        fs (1,1) double {mustBePositive}
        params.KPreview          (1,1) double {mustBeInteger, mustBePositive} = 48
        params.KMemory           (1,1) double {mustBeInteger, mustBePositive} = 48
        params.beta              (1,1) double {mustBeNonnegative}              = 1/3
        params.vmax              (1,1) double {mustBePositive}                = 1.0
        params.lambda0           (1,1) double {mustBeNonnegative}             = 0.05
        params.horizonGain       (1,1) double {mustBeNonnegative}             = 2.0
        params.memoryGain        (1,1) double {mustBeNonnegative}             = 0.0
        params.sensoryDecay      (1,1) double {mustBeNonnegative}             = 0.05
        params.forwardModelDecay (1,1) double {mustBeNonnegative}             = 0.20
        params.sigmaEMA          (1,1) double {mustBeInRange(params.sigmaEMA, 0, 1)} = 0.10
        params.omegaEMA          (1,1) double {mustBeInRange(params.omegaEMA, 0, 1)} = 0.10
        params.reliabilityGain   (1,1) double {mustBeNonnegative}             = 1.0
        params.kappaFloor        (1,1) double {mustBePositive}                = 1e-4
        params.allowFallbackDiff (1,1) logical                               = false
        % Pipeline selection for differentiateKinematicsEBR (Fraser et al. 2025).
        % filterType 6 = Savitzky-Golay with fs-scaling (temporal equivalence).
        % filterParams [order, reference_framelen]: reference_framelen at 100 Hz reference.
        % regressType 3 = fitlm OLS with free intercept (log-log); use 5 for IRLS robust.
        params.filterType        (1,1) double {mustBeInteger, mustBePositive} = 6
        params.filterParams      (1,:) double                                 = [4, 17]
        params.regressType       (1,1) double {mustBeInteger, mustBePositive} = 3
    end

    %% Input normalisation and structural validation (fail loud)
    x = x(:);
    y = y(:);
    if numel(x) ~= numel(y)
        error("PCM:lengthMismatch", "%s", "x and y must have the same number of samples.");
    end
    N = numel(x);
    span = params.KMemory + params.KPreview + 1;
    if N <= span
        msg = sprintf("Trajectory too short: N=%d but KMemory+KPreview+1=%d.", N, span);
        error("PCM:tooShort", "%s", msg);
    end
    if any(~isfinite(x)) || any(~isfinite(y))
        error("PCM:nonFinitePosition", "%s", "x and y must be finite (no NaN/Inf).");
    end

    dt = 1 / fs;
    epsP = 1e-9;   % guards a division by accumulated uncertainty, not a fudge factor

    %% Kinematics via EBR (preferred) or disclosed fallback
    [vx, vy, ax, ay, diffInfo] = computeKinematics(x, y, fs, ...
        params.filterType, params.filterParams, params.allowFallbackDiff);

    % Curvature: use curvatureKinematicEBR (Viviani & Stucchi 1992) when available.
    if exist("curvatureKinematicEBR", "file") == 2
        kappaRaw = curvatureKinematicEBR(vx, vy, ax, ay);
        diffInfo.curvatureMethod = "curvatureKinematicEBR";
    else
        speedInst = hypot(vx, vy);
        kappaRaw  = abs(vx.*ay - vy.*ax) ./ (speedInst.^3 + epsP);
        diffInfo.curvatureMethod = "inline-fallback";
        fprintf(2, "[PCM][WARNING] curvatureKinematicEBR not found; using inline formula.\n");
    end

    % Sanitise kappaRaw: curvatureKinematicEBR has no eps protection on the
    % speed^3 denominator, so zero-speed samples (filter transients, stylus
    % lift, BW initialisation artefacts) produce NaN (0/0). These propagate
    % through dkappa into the omega EMA and collapse the precision-weighted
    % channel fusion. Replace with 0: treats the sample as locally straight,
    % conservative for a degenerate velocity estimate, and harmless given
    % temporal integration over the full KMemory/KPreview window.
    nBadKappa = nnz(~isfinite(kappaRaw));
    kappaRaw(~isfinite(kappaRaw)) = 0;
    diffInfo.nKappaSanitised = nBadKappa;

    % Causal curvature rate for the forward model: backward difference uses only the past.
    dkappa = [0; diff(kappaRaw)] / dt;

    % Constant-velocity one-step prediction error (an acceleration-scale quantity).
    xPred = x(1:end-1) + vx(1:end-1)*dt;
    yPred = y(1:end-1) + vy(1:end-1)*dt;
    err   = hypot(x(2:end) - xPred, y(2:end) - yPred);
    err   = [err; err(end)];

    % Robust normalisation scales so horizonGain/reliabilityGain are dimensionless.
    sigmaScale = median(err.^2, "omitnan");
    omegaScale = median(abs(dkappa), "omitnan");
    if ~(sigmaScale > 0)
        error("PCM:degenerateTrajectory", "%s", ...
            "Prediction-error scale is zero; trajectory has no usable dynamics.");
    end
    if ~(omegaScale > 0)
        omegaScale = epsP;   % flat-curvature edge case; volatility term collapses to ~0
    end

    %% Pre-allocate outputs as NaN so unmodellable edges stay visibly empty
    speedPred        = nan(N, 1);
    kappaEff         = nan(N, 1);
    reliability      = nan(N, 1);
    sigmaNormTrace   = nan(N, 1);
    omegaTrace       = nan(N, 1);
    lambdaFutTrace   = nan(N, 1);
    lambdaPastTrace  = nan(N, 1);
    sensoryWeight    = nan(N, 1);
    imaginationWeight = nan(N, 1);

    %% Online states (warmed up from sample 1)
    sigma = err(1)^2;
    omega = abs(dkappa(1));
    floorHits = 0;

    tStart = params.KMemory + 1;
    tEnd   = N - params.KPreview;   % sensory preview needs room ahead; imagination does not

    %% Main loop
    for t = 1:tEnd

        % Update uncertainty and volatility every step so they are warm by tStart.
        sigma = (1 - params.sigmaEMA)*sigma + params.sigmaEMA*err(t)^2;
        omega = (1 - params.omegaEMA)*omega + params.omegaEMA*abs(dkappa(t));

        if t < tStart
            continue
        end

        sigmaNorm = sigma / (sigma + sigmaScale);   % in [0,1)
        omegaNorm = omega / (omega + omegaScale);   % in [0,1)

        % Jen's hypothesis: higher uncertainty -> faster decay -> shorter future horizon.
        lambdaFuture = params.lambda0 + params.horizonGain*sigmaNorm;
        lambdaPast   = params.lambda0 + params.memoryGain*sigmaNorm;

        % --- Integrate curvature over past, present, future ---
        C = kappaRaw(t);   % present anchors the window
        W = 1;
        wsAccum = 0; wiAccum = 0; wFutAccum = 0;

        % Retrospective memory: true past curvature (legitimately available).
        for k = 1:params.KMemory
            wTime = exp(-lambdaPast*k);
            C = C + wTime*kappaRaw(t - k);
            W = W + wTime;
        end

        % Prospective estimate: sensory preview (peeks, visible) fused with a
        % forward model (does not peek), precision-weighted.
        for k = 1:params.KPreview
            wTime = exp(-lambdaFuture*k);

            % Channel 1: sensory preview of the visible template.
            visibility = exp(-params.sensoryDecay*k);
            kSens = kappaRaw(t + k);
            pSens = visibility / (sigma + epsP);

            % Channel 2: forward model. Extrapolate curvature from state at t only.
            kImag = max(kappaRaw(t) + (k*dt)*dkappa(t), 0);
            confFwd = exp(-params.forwardModelDecay*k) * exp(-omegaNorm);
            pImag = confFwd / (sigma + epsP);

            wsum = pSens + pImag;
            if ~(wsum > 0)
                error("PCM:precisionCollapse", "%s", ...
                    "Future-channel precision summed to zero; check decay parameters.");
            end
            ws = pSens / wsum;
            wi = pImag / wsum;

            kFuture = ws*kSens + wi*kImag;
            C = C + wTime*kFuture;
            W = W + wTime;

            wsAccum   = wsAccum + wTime*ws;
            wiAccum   = wiAccum + wTime*wi;
            wFutAccum = wFutAccum + wTime;
        end

        kappaEff(t) = C / W;

        % Disclosed secondary slowing under uncertainty (NOT hidden; see summary doc:
        % it co-varies with curvature, so set reliabilityGain=0 to isolate pure beta).
        reliability(t) = 1 / (1 + params.reliabilityGain*sigmaNorm);

        kEff = max(kappaEff(t), params.kappaFloor);
        if kappaEff(t) < params.kappaFloor
            floorHits = floorHits + 1;
        end

        % FIX 1 + FIX 2: integrated curvature drives a POWER-law speed.
        speedPred(t) = params.vmax * kEff^(-params.beta) * reliability(t);

        % Diagnostics
        sigmaNormTrace(t)  = sigmaNorm;
        omegaTrace(t)      = omega;
        lambdaFutTrace(t)  = lambdaFuture;
        lambdaPastTrace(t) = lambdaPast;
        sensoryWeight(t)     = wsAccum / wFutAccum;
        imaginationWeight(t) = wiAccum / wFutAccum;
    end

    %% Curvature-floor diagnostic (fail loud if it could bias the exponent)
    nModelled = tEnd - tStart + 1;
    floorFrac = floorHits / max(nModelled, 1);
    if floorFrac > 0.05
        wmsg = sprintf("kappaFloor hit on %.1f%% of modelled samples; " + ...
            "the recovered exponent may be biased. Lower kappaFloor or inspect geometry.", ...
            100*floorFrac);
        warning("PCM:floorRate", "%s", wmsg);
    end

    %% Recover the exponent from the model's own output (self-consistency only)
    fit = fitPowerLaw(kappaEff, speedPred, params.kappaFloor, ...
        params.regressType, params.beta);

    %% Assemble output
    out = struct();
    out.speedPred         = speedPred;
    out.kappaEff          = kappaEff;
    out.kappaRaw          = kappaRaw;
    out.reliability       = reliability;
    out.sigmaNorm         = sigmaNormTrace;
    out.omega             = omegaTrace;
    out.lambdaFuture      = lambdaFutTrace;
    out.lambdaPast        = lambdaPastTrace;
    out.sensoryWeight     = sensoryWeight;
    out.imaginationWeight = imaginationWeight;
    out.fit               = fit;
    out.diff              = diffInfo;
    out.diff.floorFraction = floorFrac;
    out.params            = params;
end


function [vx, vy, ax, ay, info] = computeKinematics(x, y, fs, filterType, filterParams, allowFallback)
% COMPUTEKINEMATICS  Velocity and acceleration via differentiateKinematicsEBR (Fraser et al. 2025).
%
%   Real signature (verified against differentiateKinematicsEBR.m 2026-05-31):
%       [dx, dy] = differentiateKinematicsEBR(x, y, filterType, filterParams, fs)
%       x, y         : N-by-1 coordinate vectors (separate, NOT concatenated)
%       filterType   : integer 1-8 selecting the differentiation method
%       filterParams : row vector; content depends on filterType:
%                      case 6 (SG scaled)  -> [order, reference_framelen]  e.g. [4, 17]
%                      case 2 (BWFD)       -> [order, Fc_Hz, zeroLag]      e.g. [4, 10, 1]
%       fs           : sampling rate in Hz
%       Returns dx, dy each N-by-4: col 1=position, col 2=velocity,
%                                   col 3=acceleration, col 4=jerk.
%
%   filterType=6 (SG temporal equivalence, Fraser et al.) is the project default.

    info = struct();

    if exist("differentiateKinematicsEBR", "file") == 2
        try
            [dx, dy] = differentiateKinematicsEBR(x, y, filterType, filterParams, fs);
            vx = dx(:, 2);   % velocity x
            vy = dy(:, 2);   % velocity y
            ax = dx(:, 3);   % acceleration x
            ay = dy(:, 3);   % acceleration y
            info.method = sprintf("differentiateKinematicsEBR(filterType=%d)", filterType);
            fprintf("[PCM] Differentiation: %s (Fraser et al. 2025).\n", info.method);
            return
        catch ME
            msg = sprintf("differentiateKinematicsEBR failed (filterType=%d). " + ...
                "Check filterParams are valid for this filterType.", filterType);
            baseME = MException("PCM:ebrFailed", "%s", msg);
            baseME = addCause(baseME, ME);
            rethrow(baseME);
        end
    end

    if ~allowFallback
        msg = "differentiateKinematicsEBR not found. Add src/functions/ to path, " + ...
              "or set allowFallbackDiff=true for the disclosed central-difference fallback.";
        error("PCM:noDifferentiator", "%s", msg);
    end

    % Disclosed fallback: central differences (gradient). Results are provisional.
    dt = 1 / fs;
    vx = gradient(x, dt);
    vy = gradient(y, dt);
    ax = gradient(vx, dt);
    ay = gradient(vy, dt);
    info.method = "fallback-central-difference";
    fprintf(2, "[PCM][WARNING] FALLBACK central-difference differentiation active. " + ...
        "Add src/functions/ to MATLAB path to use differentiateKinematicsEBR.\n");
end


function fit = fitPowerLaw(kappaEff, speedPred, kappaFloor, regressType, beta0)
% FITPOWERLAW  Recover the exponent via regressDataEBR (Fraser et al. 2025).
%
%   FAILURE POLICY (critical for sweep integrity)
%   Two distinct failure modes are handled differently by design:
%
%   Mode A -- regressDataEBR absent from MATLAB path.
%     Disclosed deployment fallback: OLS in log-log space, fit.method = "ols-fallback".
%     Produces a result; fit.failed = false. Analogous to allowFallbackDiff.
%     Acceptable for single-call diagnostics; should not appear in sweep runs.
%
%   Mode B -- regressDataEBR is on the path but fails at runtime.
%     fit.betaRecovered = NaN, fit.failed = true, fit.failureReason = ME.message.
%     Does NOT fall back to OLS. Rationale: runtime failure (typically optimiser
%     non-convergence for regressType 4/5) is a per-configuration event with
%     methodological content. Non-convergence concentrates in high-noise, high-alpha
%     regions -- precisely where estimator choice matters most. Silent OLS substitution
%     would produce a heterogeneous result set (some cells IRLS, some OLS) with no
%     per-result record of which method applied, contaminating delta_beta estimates
%     and cell-level variance in the pipeline sweep. The sweep script accumulates
%     failed configurations in a failureLog for post-hoc inspection.
%
%   fit.failed is always present (false on success, true on Mode B failure).
%   fit.method records the estimator that produced the result, or "failed-ebr(typeN)".
%
%   Real regressDataEBR signature (verified 2026-05-31):
%       [beta, yGain, stats] = regressDataEBR(speed, kappa, regressType,
%                                              LMseeds, displayGraphs, limitBreak)
%       speed, kappa : raw values (log applied internally)
%       regressType  : 1=OLS-no-intercept, 2=fitlm-no-intercept,
%                      3=fitlm-free-intercept, 4=LM-nonlinear, 5=IRLS-robust
%       LMseeds      : [VGF0, beta0] for types 4 and 5
%       Returns beta positive (~1/3); yGain raw scale factor.
%       R^2 not returned; computed manually here.
%
%   IMPORTANT: betaRecovered is a self-consistency check only. Any empirical or
%   mechanistic claim requires v058 loop-closure / non-identifiability validation.
%   Set reliabilityGain=0 to isolate pure beta recovery.

    valid = isfinite(kappaEff) & isfinite(speedPred) & ...
            (kappaEff > kappaFloor) & (speedPred > 0);

    fit        = struct();
    fit.nUsed  = nnz(valid);
    fit.failed = false;
    fit.failureReason = "";

    if fit.nUsed < 10
        warning("PCM:tooFewForFit", "%s", "Fewer than 10 valid samples for power-law fit.");
        fit.betaRecovered = NaN;  fit.yGain = NaN;  fit.r2 = NaN;
        fit.method = "insufficient-data";
        return
    end

    kappaSub = kappaEff(valid);
    speedSub = speedPred(valid);
    LMseeds  = [1.0, beta0];   % [VGF0, beta0] for nonlinear types 4 and 5

    if exist("regressDataEBR", "file") == 2

        % --- Mode B: EBR present; runtime failure -> NaN, not OLS ---
        try
            [betaOut, yGainOut, ~] = regressDataEBR(speedSub, kappaSub, ...
                regressType, LMseeds, 0, 1);
            fit.betaRecovered = betaOut;   % positive by regressDataEBR convention
            fit.yGain         = yGainOut;
            fit.method        = sprintf("regressDataEBR(type=%d)", regressType);
            % R^2 computed manually; regressDataEBR does not return it.
            logSpeed    = log(speedSub);
            logKappa    = log(kappaSub);
            logSpeedHat = log(max(yGainOut, eps)) + (-betaOut)*logKappa;
            ssRes = sum((logSpeed - logSpeedHat).^2);
            ssTot = sum((logSpeed - mean(logSpeed)).^2);
            fit.r2 = 1 - ssRes / max(ssTot, eps);
        catch ME
            % Record the failure; return NaN. Do NOT substitute OLS.
            % The sweep failureLog captures this configuration for post-hoc analysis.
            fit.betaRecovered = NaN;
            fit.yGain         = NaN;
            fit.r2            = NaN;
            fit.failed        = true;
            fit.failureReason = ME.message;
            fit.method        = sprintf("failed-ebr(type=%d)", regressType);
            warning("PCM:ebrRegressFail", ...
                "regressDataEBR(type=%d) failed: %s. " + ...
                "Result is NaN. Check sweep failureLog for this configuration.", ...
                regressType, ME.message);
        end
        return
    end

    % --- Mode A: regressDataEBR absent from path. Disclosed OLS deployment fallback. ---
    % This path should not appear in sweep runs. If it does, the sweep script will
    % detect fit.method == "ols-fallback" and log a configuration-level warning.
    logSpeed = log(speedSub);
    logKappa = log(kappaSub);
    A      = [ones(numel(logKappa), 1), logKappa];
    coeffs = A \ logSpeed;
    yHat   = A * coeffs;
    ssRes  = sum((logSpeed - yHat).^2);
    ssTot  = sum((logSpeed - mean(logSpeed)).^2);
    fit.betaRecovered = -coeffs(2);       % negate slope: positive beta convention
    fit.yGain         = exp(coeffs(1));
    fit.r2            = 1 - ssRes / max(ssTot, eps);
    fit.method        = "ols-fallback";
    warning("PCM:olsFallback", "%s", ...
        "regressDataEBR not on MATLAB path; using OLS fallback. " + ...
        "Add src/functions/ to path before running sweeps.");
end
