function test_curvature_kinematic_vs_menger_v001()
% test_curvature_kinematic_vs_menger_v001  Kinematic vs Menger curvature on a clean ellipse.
%
% Generates a noise-free analytical ellipse sampled uniformly in arc length,
% then computes curvature by three methods:
%   (1) Analytical formula: kappa = ab / (a^2 sin^2 + b^2 cos^2)^(3/2)
%   (2) Kinematic (cross-product): curvatureKinematicEBR after SG differentiation
%   (3) Menger (three-point circle): curvatureMengerEBR on consecutive triplets
%
% Pass criteria on trimmed interior:
%   RMS(kinematic - analytical) / mean(analytical) < RMS_TOL_REL
%   RMS(menger    - analytical) / mean(analytical) < RMS_TOL_REL
%   RMS(kinematic - menger)     / mean(analytical) < PAIRWISE_TOL_REL
%
% Uses a/b = 2 for wide curvature dynamic range.
% Self-contained -- no dependency on generateSyntheticData_v011.
%
% Run from src/ (tests/ must be on the path).
%
% Dagmar Scott Fraser | PowerLawSimulationPreReg | 2026

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions')));

%% ---- Configuration ----
A                = 80;      % semi-major axis (mm)
B                = 40;      % semi-minor axis (mm, ratio 2:1)
FS               = 120;     % Hz
N_ORBITS         = 5;
V_CONST          = 40;      % mm/s constant tangential speed
SG_ORDER         = 4;
SG_FRAME         = 17;
MARGIN_FRAC      = 0.15;
RMS_TOL_REL      = 0.015;   % vs analytical
PAIRWISE_TOL_REL = 0.020;   % kinematic vs menger

fprintf('=== test_curvature_kinematic_vs_menger_v001 ===\n');
fprintf('Ellipse a=%g mm, b=%g mm (ratio %.1f), fs=%d Hz\n\n', A, B, A/B, FS);

%% ---- Arc-length parametrised trajectory ----
N_HIRES   = 100000;
theta_hr  = linspace(0, 2*pi, N_HIRES + 1)'; theta_hr(end) = [];
x_hr = A * cos(theta_hr);
y_hr = B * sin(theta_hr);
ds   = sqrt(diff([x_hr; x_hr(1)]).^2 + diff([y_hr; y_hr(1)]).^2);
arcLen   = [0; cumsum(ds)];
totalLen = arcLen(end);

totalTime  = N_ORBITS * totalLen / V_CONST;
nSamples   = round(totalTime * FS);
arcSampled = mod(linspace(0, N_ORBITS * totalLen, nSamples)', totalLen);

xt = interp1(arcLen, [x_hr; x_hr(1)], arcSampled, 'pchip');
yt = interp1(arcLen, [y_hr; y_hr(1)], arcSampled, 'pchip');
theta_samp = interp1(arcLen, [theta_hr; theta_hr(1) + 2*pi], arcSampled, 'pchip');

%% ---- (1) Analytical ----
kappa_analytical = (A * B) ./ ...
    ((A^2 * sin(theta_samp).^2 + B^2 * cos(theta_samp).^2).^(3/2));

%% ---- (2) Kinematic (SG + cross-product) ----
[dx, dy] = differentiateKinematicsEBR(xt, yt, 4, [SG_ORDER, SG_FRAME], FS);
kappa_kinematic = curvatureKinematicEBR(dx(:,2), dy(:,2), dx(:,3), dy(:,3));

%% ---- (3) Menger (three-point circles, point-by-point) ----
N = numel(xt);
kappa_menger = nan(N, 1);
for i = 2:(N - 1)
    triplet = [xt(i-1), xt(i), xt(i+1); yt(i-1), yt(i), yt(i+1)];
    kappa_menger(i) = abs(curvatureMengerEBR(triplet));
end

%% ---- Evaluate on trimmed interior ----
margin = round(N * MARGIN_FRAC);
idx = (margin + 1):(N - margin);
ka = kappa_analytical(idx);
kk = kappa_kinematic(idx);
km = kappa_menger(idx);

ok = isfinite(ka) & isfinite(kk) & isfinite(km) & (ka > 0);
ka = ka(ok); kk = kk(ok); km = km(ok);
meanKappa = mean(ka);

rmsKinVsAnal = rms(kk - ka) / meanKappa;
rmsMenVsAnal = rms(km - ka) / meanKappa;
rmsKinVsMen  = rms(kk - km) / meanKappa;

passKin  = rmsKinVsAnal  < RMS_TOL_REL;
passMen  = rmsMenVsAnal  < RMS_TOL_REL;
passPair = rmsKinVsMen   < PAIRWISE_TOL_REL;
allPass  = passKin && passMen && passPair;

fprintf('  %-40s  %-10s  %-10s  %s\n', 'comparison', 'rms_rel', 'threshold', 'verdict');
fprintf('  %s\n', repmat('-', 1, 72));
fprintf('  %-40s  %-10.5f  %-10.4f  %s\n', 'kinematic vs analytical', rmsKinVsAnal, RMS_TOL_REL,      verdict(passKin));
fprintf('  %-40s  %-10.5f  %-10.4f  %s\n', 'menger vs analytical',    rmsMenVsAnal, RMS_TOL_REL,      verdict(passMen));
fprintf('  %-40s  %-10.5f  %-10.4f  %s\n', 'kinematic vs menger',     rmsKinVsMen,  PAIRWISE_TOL_REL, verdict(passPair));
fprintf('\n  mean kappa = %.6f mm^-1, n_interior = %d\n', meanKappa, sum(ok));

fprintf('\n');
if allPass
    fprintf('PASS: kinematic and Menger curvature agree with analytical on clean ellipse\n');
else
    fprintf('FAIL: curvature method discrepancy exceeded tolerance\n');
end

end

function s = verdict(pass)
if pass, s = 'PASS'; else, s = 'FAIL'; end
end
