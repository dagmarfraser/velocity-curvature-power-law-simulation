function test_pipeline_zero_noise_v001()
% test_pipeline_zero_noise_v001  All 6 pipelines recover beta_gen exactly at sigma=0.
%
% Generates synthetic elliptical trajectories (generateSyntheticData_v011,
% fs=120Hz, cached) with zero additive noise and verifies that each of the
% six analytical pipelines recovers beta_gen to within BETA_TOL.
%
% Six pipelines: {BWFD, SG} x {OLS, LMLS, IRLS}
%   Filter 2 = BWFD  [order=2, cutoff=10, zerolag=1]
%   Filter 6 = SG    [order=4, framelen=17]  (fs-scaled temporal equiv.)
%   Regress 3 = OLS   (fitlm log-log)
%   Regress 4 = LMLS  (fitnlm Levenberg-Marquardt)
%   Regress 5 = IRLS  (fitnlm iteratively reweighted)
%
% Run from src/ (tests/ must be on the path).
%
% Dagmar Scott Fraser | PowerLawSimulationPreReg | 2026

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions')));

%% ---- Configuration ----
FS          = 120;          % Hz -- cached shape file exists for this rate
CANVAS      = [1920; 1080]; % pixels (matches toolchain createConfigSettings)
SHAPE_NUM   = 6;            % ellipse (angular freq '2')
ORBIT_COUNT = 10;
RESAMPLE    = 20;
Y_GAIN      = exp(5.2);     % ~1 Hz orbit frequency
BETA_TOL    = 0.015;        % max |beta_rec - beta_gen| at zero noise
MARGIN_FRAC = 0.10;         % trim this fraction of samples from each edge

betaVals = [0, 1/6, 1/3, 1/2, 2/3];

pipelines = struct( ...
    'name',        {'BWFD-OLS', 'BWFD-LMLS', 'BWFD-IRLS', 'SG-OLS',  'SG-LMLS',  'SG-IRLS'}, ...
    'filterType',  {2,           2,            2,            6,          6,           6        }, ...
    'filterParams',{[2,10,1],    [2,10,1],     [2,10,1],     [4,17],     [4,17],      [4,17]   }, ...
    'regressType', {3,           4,            5,            3,          4,           5        });

fprintf('=== test_pipeline_zero_noise_v001 ===\n');
fprintf('fs=%d Hz, zero noise, tolerance=%.3f\n\n', FS, BETA_TOL);

allPass  = true;
nPipe    = numel(pipelines);
nBetas   = numel(betaVals);

fprintf('  %-14s', 'beta_gen');
for p = 1:nPipe, fprintf('  %-13s', pipelines(p).name); end
fprintf('\n  %s\n', repmat('-', 1, 14 + 15*nPipe));

for bi = 1:nBetas
    betaGen = betaVals(bi);

    try
        [xt, yt] = generateSyntheticData_v011(SHAPE_NUM, CANVAS, FS, betaGen, ...
            Y_GAIN, ORBIT_COUNT, RESAMPLE, 0, 0);
    catch ME
        fprintf('  beta=%.4f: FAIL generating trajectory: %s\n', betaGen, ME.message);
        allPass = false;
        continue;
    end

    xt = xt(:);
    yt = yt(:);
    N  = numel(xt);
    margin = round(N * MARGIN_FRAC);

    fprintf('  %-14.4f', betaGen);
    for p = 1:nPipe
        betaRec = runPipeline(xt, yt, pipelines(p), FS, margin);
        err  = abs(betaRec - betaGen);
        pass = isfinite(betaRec) && (err < BETA_TOL);
        if ~pass, allPass = false; end
        fprintf('  %-13s', sprintf('%.4f %s', betaRec, verdict(pass)));
    end
    fprintf('\n');
end

fprintf('\n');
if allPass
    fprintf('PASS: all 6 pipelines recover beta_gen within %.3f at zero noise\n', BETA_TOL);
else
    fprintf('FAIL: one or more pipelines exceeded tolerance %.3f\n', BETA_TOL);
end

end

% =========================================================================
function betaRec = runPipeline(xt, yt, pip, fs, margin)
try
    [dx, dy] = differentiateKinematicsEBR(xt, yt, pip.filterType, pip.filterParams, fs);
    vx = dx(:,2); vy = dy(:,2); ax = dx(:,3); ay = dy(:,3);
    v     = sqrt(vx.^2 + vy.^2);
    kappa = curvatureKinematicEBR(vx, vy, ax, ay);
    idx   = (margin + 1):(numel(v) - margin);
    v = v(idx); kappa = kappa(idx);
    ok = (v > 0) & (kappa > 0) & isfinite(v) & isfinite(kappa);
    v = v(ok); kappa = kappa(ok);
    if numel(v) < 20, betaRec = NaN; return; end
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
