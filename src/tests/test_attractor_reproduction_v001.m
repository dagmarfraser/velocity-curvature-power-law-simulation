function test_attractor_reproduction_v001()
% test_attractor_reproduction_v001  BWFD-OLS spurious attractor at Zarandi noise.
%
% Tests the central claim of Fraser et al. (2025): at Zarandi empirical noise
% (alpha~2.78, sigma~7mm, fs=100Hz), the BWFD-OLS pipeline does NOT recover
% the true beta_gen. Instead it compresses estimates toward a spurious
% attractor near beta~0.285.
%
% Two claims tested:
%
%   CLAIM A -- Attractor is NOT at 1/3:
%     beta_gen = 1/3 at Zarandi noise -> BWFD-OLS beta_rec << 1/3.
%     Expected: beta_rec in [0.10, 0.25] (claude.md: "~0.15-0.18").
%     The empirical Zarandi result of ~0.33 reflects beta_gen~0.40 in real
%     subjects, not beta_gen = 1/3. A true beta_gen = 1/3 would recover
%     much lower because 1/3 is BELOW the attractor basin.
%
%   CLAIM B -- Compression (not just bias):
%     Across beta_gen in {0, 1/6, 1/3, 1/2, 2/3}, BWFD-OLS beta_rec values
%     span much less than the input range (0.667). The attractor squeezes
%     diverse true values toward the same region.
%     Expected: spread(beta_rec) < COMPRESSION_FACTOR * spread(beta_gen).
%
%   CLAIM C -- SG-IRLS is less compressed:
%     SG-IRLS beta_rec at beta_gen=1/3 stays closer to 1/3 than BWFD-OLS,
%     demonstrating the pipeline-dependence of the bias.
%
% NOTE: Uses fs=120Hz (cached shape file). Zarandi is 100Hz but the BWFD
% filter behaviour at 10Hz cutoff is comparable across these rates.
% Noise sigma in pixels: sigma_mm * pixelScale (4.8 px/mm from toolchain).
%
% Run from src/ (tests/ must be on the path).
%
% Dagmar Scott Fraser | PowerLawSimulationPreReg | 2026

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions')));

%% ---- Configuration ----
FS          = 120;
CANVAS      = [1920; 1080];
SHAPE_NUM   = 6;
ORBIT_COUNT = 10;
RESAMPLE    = 20;
Y_GAIN      = exp(5.2);

% Zarandi empirical noise (claude.md corrected estimates)
SIGMA_MM    = 7.0;
ALPHA_NOISE = 2.78;
PIXEL_SCALE = 4.8;          % px/mm (480px/100mm from createConfigSettings)
SIGMA_PX    = SIGMA_MM * PIXEL_SCALE;   % ~33.6 px

N_REPS      = 30;
MARGIN_FRAC = 0.10;

% Claim A: beta_gen=1/3 must recover well BELOW 1/3 (attractor is not at 1/3)
% claude.md: "True beta_gen = 1/3 would recover as ~0.15-0.18"
CLAIM_A_UPPER = 0.25;   % BWFD-OLS beta_rec at beta_gen=1/3 must be below this

% Claim B: compression across the beta_gen range
% Input span = 2/3 - 0 = 0.667; BWFD-OLS output span should be much smaller
COMPRESSION_FACTOR = 0.6;  % beta_rec span < this * beta_gen span

% Claim C: SG-IRLS less biased than BWFD-OLS at beta_gen=1/3
% SG-IRLS beta_rec should be closer to 1/3 than BWFD-OLS
IRLS_ADVANTAGE_MIN = 0.03; % SG-IRLS must be at least this much closer to 1/3

betaVals  = [0, 1/6, 1/3, 1/2, 2/3];
nBetas    = numel(betaVals);

fprintf('=== test_attractor_reproduction_v001 ===\n');
fprintf('Zarandi noise: sigma=%.1f mm (%.1f px), alpha=%.2f, fs=%d Hz\n', ...
    SIGMA_MM, SIGMA_PX, ALPHA_NOISE, FS);
fprintf('N_REPS=%d per beta_gen\n\n', N_REPS);
fprintf('  %-10s  %-12s  %-12s  %-12s\n', 'beta_gen', 'BWFD-OLS', 'SG-IRLS', 'n');
fprintf('  %s\n', repmat('-', 1, 50));

meanBWFD = nan(nBetas, 1);
meanIRLS = nan(nBetas, 1);

for bi = 1:nBetas
    betaGen = betaVals(bi);

    try
        [xt, yt] = generateSyntheticData_v011(SHAPE_NUM, CANVAS, FS, betaGen, ...
            Y_GAIN, ORBIT_COUNT, RESAMPLE, 0, 0);
    catch ME
        fprintf('  beta=%.4f: FAIL generating trajectory: %s\n', betaGen, ME.message);
        continue;
    end
    xt = xt(:); yt = yt(:);
    N  = numel(xt);
    margin = round(N * MARGIN_FRAC);

    bwfd = nan(N_REPS, 1);
    irls = nan(N_REPS, 1);

    for r = 1:N_REPS
        try
            nx = generateCustomNoise_v003(N, ALPHA_NOISE, SIGMA_PX, FS);
            ny = generateCustomNoise_v003(N, ALPHA_NOISE, SIGMA_PX, FS);
            xn = xt + nx; yn = yt + ny;

            % BWFD-OLS
            [dx, dy] = differentiateKinematicsEBR(xn, yn, 2, [2,10,1], FS);
            [v, k]   = extractVK(dx, dy, margin);
            if ~isempty(v)
                [bwfd(r), ~, ~] = regressDataEBR(v, k, 3, [], 0);
            end

            % SG-IRLS
            [dx, dy] = differentiateKinematicsEBR(xn, yn, 6, [4,17], FS);
            [v, k]   = extractVK(dx, dy, margin);
            if ~isempty(v)
                [irls(r), ~, ~] = regressDataEBR(v, k, 5, [1, 1/3], 0, 1);
            end
        catch
        end
    end

    meanBWFD(bi) = mean(bwfd, 'omitnan');
    meanIRLS(bi) = mean(irls, 'omitnan');
    fprintf('  %-10.4f  %-12.4f  %-12.4f  n=%d\n', ...
        betaGen, meanBWFD(bi), meanIRLS(bi), sum(isfinite(bwfd)));
end

%% ---- Claim A: attractor not at 1/3 ----
idx13    = find(abs(betaVals - 1/3) < 1e-6);
bwfd_at_third = meanBWFD(idx13);
claimA   = isfinite(bwfd_at_third) && (bwfd_at_third < CLAIM_A_UPPER);

fprintf('\n--- Claim A: BWFD-OLS beta_rec at beta_gen=1/3 << 1/3 ---\n');
fprintf('  BWFD-OLS: %.4f  (must be < %.2f)  %s\n', ...
    bwfd_at_third, CLAIM_A_UPPER, verdict(claimA));
fprintf('  [Attractor is not at 1/3; empirical Zarandi ~0.33 reflects true beta_gen~0.40]\n');

%% ---- Claim B: compression across beta_gen range ----
inputSpan  = max(betaVals) - min(betaVals);
outputSpan = max(meanBWFD) - min(meanBWFD);
claimB     = isfinite(outputSpan) && (outputSpan < COMPRESSION_FACTOR * inputSpan);

fprintf('\n--- Claim B: BWFD-OLS compresses beta_rec range ---\n');
fprintf('  beta_gen span = %.4f | beta_rec span = %.4f | ratio = %.2f  (must be < %.2f)  %s\n', ...
    inputSpan, outputSpan, outputSpan/inputSpan, COMPRESSION_FACTOR, verdict(claimB));

%% ---- Claim C: SG-IRLS less biased than BWFD-OLS ----
irls_at_third = meanIRLS(idx13);
bwfd_err      = abs(bwfd_at_third - 1/3);
irls_err      = abs(irls_at_third - 1/3);
claimC        = isfinite(irls_err) && ((bwfd_err - irls_err) > IRLS_ADVANTAGE_MIN);

fprintf('\n--- Claim C: SG-IRLS less biased than BWFD-OLS at beta_gen=1/3 ---\n');
fprintf('  |BWFD-OLS - 1/3| = %.4f | |SG-IRLS - 1/3| = %.4f | advantage = %.4f  (must be > %.2f)  %s\n', ...
    bwfd_err, irls_err, bwfd_err - irls_err, IRLS_ADVANTAGE_MIN, verdict(claimC));

%% ---- Overall ----
allPass = claimA && claimB && claimC;
fprintf('\n');
if allPass
    fprintf('PASS: attractor reproduced -- all three claims confirmed\n');
else
    fprintf('FAIL: one or more attractor claims not met\n');
    fprintf('  (ClaimA=%d, ClaimB=%d, ClaimC=%d)\n', claimA, claimB, claimC);
end

end

% =========================================================================
function [v, k] = extractVK(dx, dy, margin)
vx = dx(:,2); vy = dy(:,2); ax = dx(:,3); ay = dy(:,3);
v  = sqrt(vx.^2 + vy.^2);
k  = curvatureKinematicEBR(vx, vy, ax, ay);
idx = (margin + 1):(numel(v) - margin);
v = v(idx); k = k(idx);
ok = (v > 0) & (k > 0) & isfinite(v) & isfinite(k);
v = v(ok); k = k(ok);
if numel(v) < 20, v = []; k = []; end
end

function s = verdict(pass)
if pass, s = 'PASS'; else, s = 'FAIL'; end
end
