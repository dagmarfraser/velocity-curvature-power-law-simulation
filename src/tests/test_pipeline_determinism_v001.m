function test_pipeline_determinism_v001()
% test_pipeline_determinism_v001  Identical rng seed -> bit-identical beta_rec.
%
% For each of the 6 pipelines, runs the same trajectory+noise generation twice
% with the same rng seed and verifies that beta_rec is exactly equal (isequal).
% Catches non-determinism from IRLS convergence, filter state initialisation,
% or RNG state leakage.
%
% Run from src/ (tests/ must be on the path).
%
% Dagmar Scott Fraser | PowerLawSimulationPreReg | 2026

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions')));

%% ---- Configuration ----
RNG_SEED    = 42;
FS          = 120;
CANVAS      = [1920; 1080];
SHAPE_NUM   = 6;
ORBIT_COUNT = 10;
RESAMPLE    = 20;
Y_GAIN      = exp(5.2);
BETA_GEN    = 1/3;
SIGMA_PX    = 7.0 * 4.8;   % Zarandi-level (7 mm * 4.8 px/mm)
ALPHA_NOISE = 2.78;
MARGIN_FRAC = 0.10;

pipelines = struct( ...
    'name',        {'BWFD-OLS', 'BWFD-LMLS', 'BWFD-IRLS', 'SG-OLS',  'SG-LMLS',  'SG-IRLS'}, ...
    'filterType',  {2,           2,            2,            6,          6,           6        }, ...
    'filterParams',{[2,10,1],    [2,10,1],     [2,10,1],     [4,17],     [4,17],      [4,17]   }, ...
    'regressType', {3,           4,            5,            3,          4,           5        });

fprintf('=== test_pipeline_determinism_v001 ===\n');
fprintf('rng seed=%d, beta_gen=%.4f, sigma=%.1f px, alpha=%.2f\n\n', ...
    RNG_SEED, BETA_GEN, SIGMA_PX, ALPHA_NOISE);

% Generate trajectory once (shape load is not RNG-dependent)
rng(RNG_SEED, 'twister');
[xt, yt] = generateSyntheticData_v011(SHAPE_NUM, CANVAS, FS, BETA_GEN, ...
    Y_GAIN, ORBIT_COUNT, RESAMPLE, 0, 0);
xt = xt(:); yt = yt(:);
N  = numel(xt);
margin = round(N * MARGIN_FRAC);

allPass = true;
fprintf('  %-14s  %-14s  %-14s  %-12s  %s\n', 'pipeline', 'run1_beta', 'run2_beta', 'delta', 'verdict');
fprintf('  %s\n', repmat('-', 1, 62));

for p = 1:numel(pipelines)
    pip = pipelines(p);

    rng(RNG_SEED + 1, 'twister');
    nx1 = generateCustomNoise_v003(N, ALPHA_NOISE, SIGMA_PX, FS);
    ny1 = generateCustomNoise_v003(N, ALPHA_NOISE, SIGMA_PX, FS);
    beta1 = runOnce(xt + nx1, yt + ny1, pip, FS, margin);

    rng(RNG_SEED + 1, 'twister');     % same seed -> identical noise
    nx2 = generateCustomNoise_v003(N, ALPHA_NOISE, SIGMA_PX, FS);
    ny2 = generateCustomNoise_v003(N, ALPHA_NOISE, SIGMA_PX, FS);
    beta2 = runOnce(xt + nx2, yt + ny2, pip, FS, margin);

    pass  = isequal(beta1, beta2);
    delta = abs(beta1 - beta2);
    if ~pass, allPass = false; end

    fprintf('  %-14s  %-14.8f  %-14.8f  %-12.2e  %s\n', ...
        pip.name, beta1, beta2, delta, verdict(pass));
end

fprintf('\n');
if allPass
    fprintf('PASS: all 6 pipelines are deterministic (bit-identical output)\n');
else
    fprintf('FAIL: non-determinism detected in one or more pipelines\n');
end

end

% =========================================================================
function betaRec = runOnce(xt, yt, pip, fs, margin)
try
    [dx, dy] = differentiateKinematicsEBR(xt, yt, pip.filterType, pip.filterParams, fs);
    vx = dx(:,2); vy = dy(:,2); ax = dx(:,3); ay = dy(:,3);
    v     = sqrt(vx.^2 + vy.^2);
    kappa = curvatureKinematicEBR(vx, vy, ax, ay);
    idx   = (margin + 1):(numel(v) - margin);
    v = v(idx); kappa = kappa(idx);
    ok = (v > 0) & (kappa > 0) & isfinite(v) & isfinite(kappa);
    v = v(ok); kappa = kappa(ok);
    switch pip.regressType
        case 3, [betaRec, ~, ~] = regressDataEBR(v, kappa, 3, [],        0);
        case 4, [betaRec, ~, ~] = regressDataEBR(v, kappa, 4, [1, 1/3], 0, 0);
        case 5, [betaRec, ~, ~] = regressDataEBR(v, kappa, 5, [1, 1/3], 0, 1);
    end
catch
    betaRec = NaN;
end
end

function s = verdict(pass)
if pass, s = 'PASS'; else, s = 'FAIL'; end
end
