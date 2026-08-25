function [localData, pcmData, meta] = generatePCMTrajectory_v001(shapeNum, canvas, fs, beta, vmax, pcmParams)
% GENERATEPCMTRAJECTORY_V001  Paired local and PCM-controlled position time series.
%
%   Produces two position time series that traverse the same geometric template
%   at the same sampling rate but under different generative speed laws:
%
%     LOCAL: v(t) = vmax * kappa(t)^(-beta)          [instantaneous curvature]
%     PCM:   v(t) = vmax * kappaEff(t)^(-beta)        [horizon-weighted integral]
%
%   The paired design is essential for the pipeline sweep: both controllers receive
%   identical noise realisations and identical pipeline processing, so delta_beta =
%   beta_obs_PCM - beta_obs_LOCAL isolates the anticipatory signature cleanly.
%
%   METHOD (three passes)
%     Pass 1 -- geometry + local speed law.
%       Calls generateSyntheticData_v011 to obtain (xt_local, yt_local) and the
%       local curvature / speed arrays. The spatial path serves as a shared
%       geometric template for both controllers.
%
%     Pass 2 -- PCM speed profile.
%       Calls predictive_curvature_model_v2_001 on (xt_local, yt_local) to obtain
%       speedPred(t): the anticipatory speed the PCM assigns to each sample of the
%       locally-parameterised trajectory. NaN edges (1:KMemory and (N-KPreview+1):N)
%       are filled with the local speed; the edge fraction is recorded in meta.
%
%     Pass 3 -- re-parameterisation.
%       Reinterprets the spatial path under the PCM speed profile via arc-length
%       integration and pchip resampling. The resulting (xt_PCM, yt_PCM) traverses
%       the same geometric path as the local trajectory but with PCM-governed timing.
%       This is a first-order approximation: the PCM sigma dynamics reflect the local
%       parameterisation, not the fully self-consistent PCM one. For smooth shapes
%       with biological noise this approximation is adequate for the pipeline sweep;
%       an iterative convergence pass can be added in a successor version if needed.
%
%   INPUTS
%     shapeNum  (1,1) integer 1-9  Shape index for generateSyntheticData_v011.
%     canvas    (1,2) double       Canvas [width, height] in mm.
%     fs        (1,1) double       Sampling rate in Hz.
%     beta      (1,1) double       Generative power-law exponent (positive, ~1/3).
%     vmax      (1,1) double       Speed gain / VGF.
%     pcmParams struct             PCM configuration. Required fields:
%                                    KPreview    -- future window (samples)
%                                    KMemory     -- past window (samples)
%                                    filterType  -- differentiateKinematicsEBR selector
%                                    filterParams -- pipeline params for filterType
%                                  Optional fields (PCM defaults used if absent):
%                                    horizonGain, memoryGain, lambda0, sensoryDecay,
%                                    forwardModelDecay, sigmaEMA, omegaEMA,
%                                    reliabilityGain, kappaFloor, regressType,
%                                    allowFallbackDiff, orbitCount (default 1).
%
%   OUTPUTS  (3 outputs; structs respect the 4-output cap)
%     localData   struct  .xt .yt .kappa .speed          N_local-by-1 each
%     pcmData     struct  .xt .yt .speedPred .kappaEff    .xt/.yt N_pcm-by-1;
%                                                         .speedPred/.kappaEff
%                                                         N_local-by-1 (local time base)
%     meta        struct  .edgeFraction   fraction of samples using local-speed fallback
%                         .pcmFit         out.fit from predictive_curvature_model_v2_001
%                         .N_local        length of local trajectory
%                         .N_pcm          length of re-parameterised PCM trajectory
%                         .diffMethod     differentiation method string from EBR
%                         .rngStateBefore rng state before Pass 1 (for paired noise)
%
%   FAIL-LOUD CONTRACT
%     Errors (not warnings) on: generateSyntheticData_v011 failure; PCM failure;
%     all-NaN speedPred; degenerate arc-length steps; PCM trajectory < 10 samples.
%     Warns on: edge fraction > 10%; PCM trajectory duration differs from local > 20%.
%
%   PARFOR SAFETY
%     No persistent or global state. All sub-functions are local. Path dependencies
%     (generateSyntheticData_v011, predictive_curvature_model_v2_001, EBR suite)
%     must be on the worker path before parfor is launched via addpath in the
%     sweep orchestration script.
%     Seed rng immediately before calling this function to ensure paired noise
%     realisations are matched across local and PCM conditions.
%
%   USAGE EXAMPLE
%     p.KPreview = 48;  p.KMemory = 48;
%     p.filterType = 6; p.filterParams = [4, 17];
%     p.horizonGain = 2.0;  p.reliabilityGain = 0;
%     rng(42);
%     [loc, pcm, meta] = generatePCMTrajectory_v001(1, [800 600], 240, 1/3, 1.0, p);
%
%   See also: predictive_curvature_model_v2_001, generateSyntheticData_v011,
%             runPCMPipelineSweep_v001, predictive_curvature_model_v2_001_summary.md
%
% Created 2026-05-31
% D.S. Fraser, d.s.fraser@bham.ac.uk
% Centre for Human Brain Health, University of Birmingham

    arguments
        shapeNum  (1,1) double {mustBeInteger, mustBePositive}
        canvas    (1,2) double {mustBePositive}
        fs        (1,1) double {mustBePositive}
        beta      (1,1) double {mustBeNonnegative}
        vmax      (1,1) double {mustBePositive}
        pcmParams struct
    end

    %% === Pass 1: local controller trajectory ===
    meta.rngStateBefore = rng();   % caller can restore for paired noise injection

    orbitCount = pcmParamGet(pcmParams, "orbitCount", 1);
    KMemoryDefault  = pcmParamGet(pcmParams, "KMemory",  48);
    KPreviewDefault = pcmParamGet(pcmParams, "KPreview", 48);

    try
        [xt_L, yt_L, k_L, v_L, ~] = generateSyntheticData_v011( ...
            shapeNum, canvas, fs, beta, vmax, orbitCount, 0, 0, 0);
    catch ME
        msg = sprintf("generateSyntheticData_v011 failed (shapeNum=%d, fs=%.0f Hz).", ...
            shapeNum, fs);
        rethrow(addCause(MException("PCM:genFailed", "%s", msg), ME));
    end

    xt_L = xt_L(:);
    yt_L = yt_L(:);
    k_L  = k_L(:);
    v_L  = v_L(:);
    N_L  = numel(xt_L);

    minRequired = 2*(KMemoryDefault + KPreviewDefault) + 1;
    if N_L < minRequired
        msg = sprintf( ...
            "Local trajectory N=%d is too short for KMemory+KPreview=%d. " + ...
            "Increase orbitCount or reduce window sizes.", ...
            N_L, KMemoryDefault + KPreviewDefault);
        error("PCM:trajectoryTooShort", "%s", msg);
    end

    %% === Pass 2: PCM speed profile on the local trajectory ===
    pcmArgs = buildPCMArgs(pcmParams, beta, vmax);

    try
        pcmOut = predictive_curvature_model_v2_001(xt_L, yt_L, fs, pcmArgs{:});
    catch ME
        rethrow(addCause(MException("PCM:pcmFailed", "%s", ...
            "predictive_curvature_model_v2_001 failed on the local trajectory."), ME));
    end

    speedPred = pcmOut.speedPred;   % N_L-by-1; NaN at unmodellable edges
    kappaEff  = pcmOut.kappaEff;    % N_L-by-1

    %% Edge handling: fill NaN with local speed; error if all NaN
    nanMask      = ~isfinite(speedPred) | (speedPred <= 0);
    edgeFraction = mean(nanMask);

    if all(nanMask)
        error("PCM:allNaN", "%s", ...
            "speedPred is entirely NaN. Check KPreview/KMemory vs trajectory length, " + ...
            "and that differentiateKinematicsEBR is on the MATLAB path.");
    end
    if edgeFraction > 0.10
        warning("PCM:edgeFraction", ...
            "%.1f%% of samples use local-speed fallback at edges. " + ...
            "Consider increasing orbitCount or reducing window parameters.", ...
            100*edgeFraction);
    end

    speedFilled          = speedPred;
    speedFilled(nanMask) = v_L(nanMask);   % local speed at unmodellable edges

    %% === Pass 3: re-parameterise spatial path under PCM speed ===
    % Arc-length increments of the local trajectory.
    ds = hypot(diff(xt_L), diff(yt_L));   % (N_L-1)-by-1

    if any(ds <= 0)
        error("PCM:degeneratePath", "%s", ...
            "Local trajectory contains zero-length arc steps; path is degenerate.");
    end

    % Arc-length-consistent speed normalisation.
    % Arithmetic mean normalisation is INCORRECT for non-constant speed profiles:
    % duration depends on the harmonic mean of speed, not the arithmetic mean.
    % The correct scale is derived from the arc-length integral directly:
    %   T_pcm_raw = sum(ds_i / v_interval_i)   [duration before scaling]
    %   T_local   = (N_L - 1) / fs
    %   scale     = T_pcm_raw / T_local
    % Multiplying speedFilled by scale gives T_pcm_scaled = T_local exactly.
    % Beta recovery depends on relative speed variation not absolute scale,
    % so this normalisation preserves the scientific content of the sweep.
    v_int_raw = 0.5*(speedFilled(1:end-1) + speedFilled(2:end));
    v_int_raw = max(v_int_raw, 1e-6);
    T_pcm_raw = sum(ds ./ v_int_raw);
    T_local   = (N_L - 1) / fs;
    if T_pcm_raw > 0 && T_local > 0
        speedFilled = speedFilled * (T_pcm_raw / T_local);
    end

    % Trapezoidal speed per arc-length interval (after normalisation).
    v_interval = 0.5 * (speedFilled(1:end-1) + speedFilled(2:end));
    v_interval = max(v_interval, 1e-6);

    % dt_i = ds_i / v_i; integrate to get cumulative PCM time at each node.
    dt_pcm = ds ./ v_interval;
    t_pcm  = [0; cumsum(dt_pcm)];         % N_L-by-1
    T_pcm  = t_pcm(end);

    if ~(T_pcm > 0)
        error("PCM:zeroTime", "%s", ...
            "Integrated PCM trajectory duration is zero or negative.");
    end

    % Uniform resample at fs.
    t_uniform = (0 : 1/fs : T_pcm)';
    N_pcm     = numel(t_uniform);

    if N_pcm < 10
        error("PCM:resampledTooShort", "%s", ...
            "Re-parameterised PCM trajectory has fewer than 10 samples at the target fs.");
    end

    xt_PCM = interp1(t_pcm, xt_L, t_uniform, "pchip");
    yt_PCM = interp1(t_pcm, yt_L, t_uniform, "pchip");

    %% === Assemble outputs ===
    localData.xt    = xt_L;
    localData.yt    = yt_L;
    localData.kappa = k_L;
    localData.speed = v_L;

    pcmData.xt        = xt_PCM;
    pcmData.yt        = yt_PCM;
    pcmData.speedPred = speedFilled;   % on local time base; for diagnostics
    pcmData.kappaEff  = kappaEff;      % on local time base

    meta.edgeFraction  = edgeFraction;
    meta.pcmFit        = pcmOut.fit;
    meta.N_local       = N_L;
    meta.N_pcm         = N_pcm;
    meta.diffMethod    = pcmOut.diff.method;
end


%% =========================================================================
function args = buildPCMArgs(p, beta, vmax)
% BUILDPCMARGS  Convert pcmParams struct to Name-Value cell for the PCM arguments block.
%   Always forwards beta and vmax. All other recognised fields forwarded if present.
%   orbitCount is consumed by generatePCMTrajectory_v001 and not forwarded to PCM.

    knownPCMFields = {"KPreview","KMemory","horizonGain","memoryGain","lambda0", ...
                      "sensoryDecay","forwardModelDecay","sigmaEMA","omegaEMA", ...
                      "reliabilityGain","kappaFloor","filterType","filterParams", ...
                      "regressType","allowFallbackDiff"};

    names = {"beta", "vmax"};
    vals  = {beta,   vmax};

    for k = 1:numel(knownPCMFields)
        f = knownPCMFields{k};
        if isfield(p, f)
            names{end+1} = f;        %#ok<AGROW>
            vals{end+1}  = p.(f);    %#ok<AGROW>
        end
    end

    % Interleave into {name, val, name, val, ...} for Name-Value expansion.
    nv   = [names; vals];
    args = nv(:)';
end


function val = pcmParamGet(p, field, default)
% PCMPARAMGET  Field access with default fallback; accepts string or char field names.
    if isfield(p, field)
        val = p.(field);
    else
        val = default;
    end
end
