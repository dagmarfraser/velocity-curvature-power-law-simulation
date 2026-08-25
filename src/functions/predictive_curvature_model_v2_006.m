function out = predictive_curvature_model_v2_006(x, y, fs, params)
% PREDICTIVE_CURVATURE_MODEL_V2_006  Adds D17: visual-vestibular dual curvature channel.
%
%   Adds to v2_005:
%     NEW visualChannelWeight  Blends visual (Menger, geometric) and vestibular
%       (kinematic, centripetal force) curvature in the FUTURE SENSORY PREVIEW channel.
%       Past memory and present sample always use κ_kinematic (vestibular).
%       Future sensory preview uses:
%           kSens = (1-w)·κ_kinematic(t+k)  +  w·κ_Menger(t+k)
%       where w = visualChannelWeight ∈ [0,1].
%       visualChannelWeight=0 recovers v2_005 exactly (all kinematic).
%       visualChannelWeight=1 (default): all Menger — models template-tracing where
%       visual system apprehends path geometry ahead without vestibular mediation.
%
%   Physical motivation:
%     κ_kinematic: numerator |vx·ay − vy·ax| ∝ centripetal force. The body FEELS
%       this via Type II vestibular otolith cells sensing v²κ directly.
%     κ_Menger (k=1): computed from three consecutive positions — zero
%       differentiations. The visual system apprehends spatial path geometry
%       without needing to sense velocity or acceleration. The body SEES this.
%
%   DA prediction: PLAC > HALO on fitted visualChannelWeight (reduced DA → impaired
%   precision weighting of visual vs vestibular channels). Separable from adaptRate
%   (horizon width) and volatilityGain (geometric speed modulation).
%
%   Connection to D16: at visualChannelWeight=1.0 and β_gen*=0.5 (D16-C, centripetal
%   noise), jerk j = v³·(κ'/κ)·(1−2β) = 0. Visual preview saturates the future
%   channel; vestibular Type I jerk-sensing is silent. This is the bridge to
%   Huh-Sejnowski (2015): minimum-jerk and minimum-variance-under-centripetal-noise
%   converge at the same trajectory (β=0.5, j=0).
%
%   NEW output fields:
%     out.kappaMenger       N-by-1 Menger curvature trace (k=1 lag, visual channel).
%     out.visualChannelWeight  scalar provenance of the blend weight used.
%
%   For cyclic multi-loop data, pre-segment into single cycles before calling;
%   dGoalNorm is computed over the full trajectory as passed.
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
%     out.alphaFuture       N-by-1 adaptive horizon state (EMA of sigmaNorm).
%     out.lambdaPast        N-by-1 past temporal decay (1/sample).
%     out.dGoalNorm         N-by-1 normalised goal distance (1=start, 0=end).
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
        params.adaptRate         (1,1) double {mustBeInRange(params.adaptRate,0,1)} = 0.05
        params.alpha0            (1,1) double {mustBeNonnegative}             = 0.0
        params.memoryGain        (1,1) double {mustBeNonnegative}             = 0.0
        params.goalDecay         (1,1) double                                 = 0.0
        params.goalDecayMode     (1,:) char                                   = 'kappa_peak'
        params.sensoryDecay      (1,1) double {mustBeNonnegative}             = 0.05
        params.forwardModelDecay (1,1) double {mustBeNonnegative}             = 0.20
        params.sigmaEMA          (1,1) double {mustBeInRange(params.sigmaEMA, 0, 1)} = 0.10
        params.omegaEMA          (1,1) double {mustBeInRange(params.omegaEMA, 0, 1)} = 0.10
        params.volatilityGain    (1,1) double {mustBeNonnegative}             = 0.0
        params.visualChannelWeight (1,1) double {mustBeInRange(params.visualChannelWeight,0,1)} = 1.0
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

    %% Sub-goal distance trace (D2 — v2_004: mode-switchable)
    % goalDecay=0 disables mechanism regardless of mode.
    % 'linear':     monotone trial-arc fraction (was null — kept for reference).
    % 'kappa_peak': Option B — arc to next κ peak, normalised by mean inter-peak arc.
    % 'cycle':      Option A — per-cycle reset; dGoalNorm = 1 − arc_within_cycle/arc_per_cycle.
    switch params.goalDecayMode
        case 'linear'
            arcStep   = hypot(vx, vy) .* dt;
            cumArc    = [0; cumsum(arcStep(1:end-1))];
            totalArc  = cumArc(end) + arcStep(end);
            dGoalNorm = 1 - cumArc ./ max(totalArc, epsP);
        case 'kappa_peak'
            dGoalNorm = pcmKappaPeakGoalNorm(kappaRaw, vx, vy, dt, epsP);
        case 'cycle'
            dGoalNorm = pcmCycleGoalNorm(kappaRaw, vx, vy, dt, epsP);
        otherwise
            error('PCMv4:unknownMode', '%s', ...
                sprintf("Unknown goalDecayMode '%s'. Use 'linear','kappa_peak','cycle'.", ...
                params.goalDecayMode));
    end
    diffInfo.nKappaSanitised = nBadKappa;

    %% Menger curvature (H-S radian K-lag) -- visual apprehension channel (D17)
    % Three anchors: p[t-K] (vestibular prior), p[t] (current), p[t+K] (visual future).
    % K = round(1/lambda0) -- the H-S planning horizon in samples (~29 at 133 Hz, lambda0=0.035).
    % Only p[t+K] requires forward access; p[t-K] is already in vestibular memory.
    % At K=1: noise-dominated (0.15 mm spacing). At K=7: matches kinematic (r=0.963, Gate 2).
    % At K=lambda0 samples: stable (~2-3 cm spacing), forward-looking, distinct from kinematic.
    K_hs = max(1, round(1 / params.lambda0));  % H-S radian window in samples
    i_m = (K_hs + 1 : N - K_hs)';             % valid centre indices
    px  = x(i_m);  py = y(i_m);
    dx_pc = px           - x(i_m - K_hs);  dy_pc = py           - y(i_m - K_hs);  % current - past
    dx_pf = x(i_m+K_hs) - x(i_m - K_hs);  dy_pf = y(i_m+K_hs) - y(i_m - K_hs);  % future  - past
    dx_cf = x(i_m+K_hs) - px;              dy_cf = y(i_m+K_hs) - py;               % future  - current
    area2M = abs(dx_pc.*dy_pf - dx_pf.*dy_pc);  % = 2 x triangle area
    d_pc   = hypot(dx_pc, dy_pc);  % past-to-current distance
    d_cf   = hypot(dx_cf, dy_cf);  % current-to-future distance
    d_pf   = hypot(dx_pf, dy_pf);  % past-to-future distance
    kappaMenger        = zeros(N, 1);
    kappaMenger(i_m)   = 2 * area2M ./ (d_pc .* d_cf .* d_pf + eps);
    diffInfo.K_hs      = K_hs;  % provenance

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
    omegaNormTrace   = nan(N, 1);   % D15: normalised volatility (soft-normalised ω)
    lambdaFutTrace   = nan(N, 1);
    lambdaPastTrace  = nan(N, 1);
    sensoryWeight    = nan(N, 1);
    imaginationWeight = nan(N, 1);
    alphaFutureTrace = nan(N, 1);   % BUG FIX: was growing dynamically in v2_002

    %% Online states (warmed up from sample 1)
    sigma       = err(1)^2;
    omega       = abs(dkappa(1));
    alphaFuture = params.alpha0;   % initial adaptive horizon state
    floorHits   = 0;

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

        % Adaptive horizon (Variant B): alphaFuture tracks sigmaNorm via EMA.
        % recoveryHalfLife = log(0.5)/log(1-adaptRate)/fs seconds.
        % adaptRate=0: static horizon (local controller limit, alphaFuture=alpha0).
        % adaptRate=1: instantaneous tracking (alphaFuture=sigmaNorm each sample).
        alphaFuture  = alphaFuture + params.adaptRate*(sigmaNorm - alphaFuture);
        lambdaFuture = params.lambda0 + alphaFuture*sigmaNorm;
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

            % Channel 1: sensory preview — visual (Menger) blended with vestibular (kinematic).
            % visualChannelWeight=1: all Menger (template tracing, visual apprehension).
            % visualChannelWeight=0: all kinematic (free drawing, recovers v2_005).
            kSensKin    = kappaRaw(t + k);
            kSensMenger = kappaMenger(t + k);
            kSens = (1 - params.visualChannelWeight)*kSensKin + ...
                         params.visualChannelWeight *kSensMenger;
            visibility = exp(-params.sensoryDecay*k);
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

        % D15: curvature-volatility slowing (ω channel, separable from σ channel).
        % volatilityGain=0 recovers v2_004 exactly.
        volatilityFactor = 1 / (1 + params.volatilityGain * omegaNorm);

        % FIX 1 + FIX 2 + D2 + D15: power-law speed with uncertainty, goal, and
        % volatility modulators. Each gain=0 removes its term.
        goalFactor   = exp(-params.goalDecay * dGoalNorm(t));
        speedPred(t) = params.vmax * kEff^(-params.beta) * reliability(t) ...
                       * goalFactor * volatilityFactor;

        % Diagnostics
        sigmaNormTrace(t)  = sigmaNorm;
        omegaTrace(t)      = omega;
        omegaNormTrace(t)  = omegaNorm;   % D15
        lambdaFutTrace(t)  = lambdaFuture;
        lambdaPastTrace(t) = lambdaPast;
        sensoryWeight(t)      = wsAccum / wFutAccum;
        imaginationWeight(t)  = wiAccum / wFutAccum;
        alphaFutureTrace(t)   = alphaFuture;
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

    %% Compute recovery half-life from adaptRate
    if params.adaptRate > 0 && params.adaptRate < 1
        fit.recoveryHalfLife = log(0.5) / log(1 - params.adaptRate) / fs;
    elseif params.adaptRate == 0
        fit.recoveryHalfLife = Inf;
    else
        fit.recoveryHalfLife = 1/fs;
    end

    %% Assemble output
    out = struct();
    out.speedPred         = speedPred;
    out.kappaEff          = kappaEff;
    out.kappaRaw          = kappaRaw;
    out.reliability       = reliability;
    out.sigmaNorm         = sigmaNormTrace;
    out.omega             = omegaTrace;
    out.omegaNorm         = omegaNormTrace;   % D15: normalised volatility trace
    out.kappaMenger       = kappaMenger;      % D17: Menger curvature (visual channel, k=1)
    out.visualChannelWeight = params.visualChannelWeight;  % D17: blend weight (provenance)
    out.lambdaFuture      = lambdaFutTrace;
    out.lambdaPast        = lambdaPastTrace;
    out.dGoalNorm         = dGoalNorm;
    out.goalDecayMode     = params.goalDecayMode;
    out.sensoryWeight     = sensoryWeight;
    out.imaginationWeight = imaginationWeight;
    out.alphaFuture       = alphaFutureTrace;
    out.fit               = fit;
    out.diff              = diffInfo;
    out.diff.floorFraction = floorFrac;
    out.params            = params;
end


%% =========================================================================
%% SUB-FUNCTIONS: Sub-goal geometry (D2, v2_004)
%% =========================================================================

function peakIdx = pcmFindKappaPeaks(kappaRaw, N)
% PCMFINDKAPPAPEAKS  Find prominent curvature peaks for sub-goal computation.
% Operates only on the finite (valid) region of kappaRaw — avoids inflated
% prominence thresholds from NaN-edge extrapolation.
    t_valid = find(isfinite(kappaRaw) & kappaRaw > 0);
    if numel(t_valid) < 10;  peakIdx = [];  return;  end
    % Exclude first/last 2% of samples — SG differentiation edge transients
    % can produce spurious high-κ artifacts that corrupt peak detection.
    nEdge   = max(3, round(0.02 * N));
    t_valid = t_valid(t_valid > nEdge & t_valid < N + 1 - nEdge);
    k_valid   = kappaRaw(t_valid);
    min_sep   = max(5, round(0.05 * numel(t_valid)));
    prom_thr  = 0.25 * (max(k_valid) - median(k_valid));
    try
        [~, pk_local] = findpeaks(k_valid, ...
            'MinPeakProminence', prom_thr, ...
            'MinPeakDistance',   min_sep);
        peakIdx = t_valid(pk_local);   % map back to global sample indices
    catch
        % Fallback: simple local maxima with minimum separation
        lm   = find(k_valid(2:end-1) > k_valid(1:end-2) & ...
                    k_valid(2:end-1) > k_valid(3:end)) + 1;
        peakIdx = [];
        last    = -inf;
        for k = 1:numel(lm)
            if lm(k) - last >= min_sep
                peakIdx(end+1) = t_valid(lm(k)); %#ok<AGROW>
                last = lm(k);
            end
        end
    end
end


function dGN = pcmKappaPeakGoalNorm(kappaRaw, vx, vy, dt, epsP)
% PCMKAPPAPEAKGOALNORM  Option B: normalised arc distance to next kappa peak.
% dGoalNorm(t) = (arc_to_next_peak) / mean_inter_peak_arc, clipped to [0,1].
% Sawtooth: descends toward 0 as each peak approaches, jumps to ~1 just after.
% Reverse-scan: unambiguous and correct for all edge cases.
    N       = numel(kappaRaw);
    arcStep = hypot(vx, vy) .* dt;
    cumArc  = [0; cumsum(arcStep(1:end-1))];
    peakIdx = pcmFindKappaPeaks(kappaRaw, N);
    if numel(peakIdx) < 2;  dGN = zeros(N,1);  return;  end
    peakArcs   = cumArc(peakIdx);
    meanIPArc  = mean(diff(peakArcs));
    if meanIPArc < epsP;  dGN = zeros(N,1);  return;  end
    sentinel       = cumArc(end) + meanIPArc;
    nextPeakArc    = sentinel * ones(N, 1);
    % Reverse scan: for each peak (latest first), assign it as next target
    % to all samples that have not yet reached it.
    for ki = numel(peakArcs):-1:1
        nextPeakArc(cumArc < peakArcs(ki)) = peakArcs(ki);
    end
    dGN = min((nextPeakArc - cumArc) ./ meanIPArc, 1);
    dGN = max(dGN, 0);
end


function dGN = pcmCycleGoalNorm(kappaRaw, vx, vy, dt, epsP)
% PCMCYCLEGOALNORM  Option A: per-cycle reset at kappa peaks.
% dGoalNorm(t) = 1 - arc_within_segment / arc_per_segment. Vectorised.
    N       = numel(kappaRaw);
    arcStep = hypot(vx, vy) .* dt;
    cumArc  = [0; cumsum(arcStep(1:end-1))];
    peakIdx = pcmFindKappaPeaks(kappaRaw, N);
    if numel(peakIdx) < 2;  dGN = zeros(N,1);  return;  end
    peakArcs   = cumArc(peakIdx);
    boundaries = [0; peakArcs(:); cumArc(end) + epsP];
    % For each sample, find which segment it belongs to via histc/discretize
    segIdx     = discretize(cumArc, boundaries);   % in [1, numel(boundaries)-1]
    segIdx     = max(min(segIdx, numel(boundaries)-1), 1);
    segStart   = boundaries(segIdx);
    segLen     = boundaries(segIdx + 1) - segStart;
    segLen     = max(segLen, epsP);
    dGN        = max(min(1 - (cumArc - segStart) ./ segLen, 1), 0);
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
                regressType, LMseeds, 0, 0);   % limitBreak=0: non-convergence → NaN (Fail Loud)
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
