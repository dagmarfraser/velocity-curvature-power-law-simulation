function test_curvatureKinematicEBR_v001()
% test_curvatureKinematicEBR_v001  Unit test for curvatureKinematicEBR.
%
% For a circle of radius R, curvature = 1/R everywhere.
% Generates a circle uniformly sampled in arc length at constant speed,
% computes analytical velocity and acceleration, then verifies:
%   RMS(kappa_rec - 1/R) / (1/R) < REL_TOL
%
% This is a pure unit test of the cross-product formula
%   kappa = |vx*ay - vy*ax| / (vx^2 + vy^2)^(3/2)
% with exact analytical inputs -- no differentiation filter involved.
%
% A second case uses SG-differentiated inputs to test the combined
% curvatureKinematicEBR + differentiateKinematicsEBR chain on a circle.
%
% Run from src/ (tests/ must be on the path).
%
% Dagmar Scott Fraser | PowerLawSimulationPreReg | 2026

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions')));

RADII       = [10, 50, 200];    % mm -- span the simulation curvature range
FS          = 120;
SPEED       = 30;               % mm/s constant tangential speed
N_ORBITS    = 4;
MARGIN_FRAC = 0.10;

REL_TOL_ANALYTIC = 1e-10;   % analytical inputs -- should be machine precision
REL_TOL_SG       = 0.01;    % SG-differentiated inputs

fprintf('=== test_curvatureKinematicEBR_v001 ===\n\n');
fprintf('  %-12s  %-14s  %-12s  %-14s  %-12s  %s\n', ...
    'radius(mm)', 'analytic_err', 'tol_anal', 'SG_err', 'tol_SG', 'verdict');
fprintf('  %s\n', repmat('-', 1, 72));

allPass = true;

for R = RADII
    trueKappa = 1 / R;

    %% ---- Generate circle (arc-length parametrised, constant speed) ----
    circumference = 2 * pi * R;
    totalTime = N_ORBITS * circumference / SPEED;
    N = round(totalTime * FS);
    theta = linspace(0, N_ORBITS * 2 * pi, N)';

    x = R * cos(theta);
    y = R * sin(theta);

    %% ---- Case 1: analytical velocity and acceleration ----
    omega = SPEED / R;          % angular velocity (rad/s)
    vx_a  = -SPEED * sin(theta);
    vy_a  =  SPEED * cos(theta);
    ax_a  = -omega * SPEED * cos(theta);
    ay_a  = -omega * SPEED * sin(theta);

    kappa_analytic = curvatureKinematicEBR(vx_a, vy_a, ax_a, ay_a);

    margin = round(N * MARGIN_FRAC);
    idx    = (margin + 1):(N - margin);
    err_a  = rms(kappa_analytic(idx) - trueKappa) / trueKappa;
    passA  = err_a < REL_TOL_ANALYTIC;

    %% ---- Case 2: SG-differentiated velocity and acceleration ----
    [dx, dy] = differentiateKinematicsEBR(x, y, 6, [4, 17], FS);
    kappa_sg = curvatureKinematicEBR(dx(idx,2), dy(idx,2), dx(idx,3), dy(idx,3));
    err_sg   = rms(kappa_sg - trueKappa) / trueKappa;
    passSG   = err_sg < REL_TOL_SG;

    pass = passA && passSG;
    if ~pass, allPass = false; end

    fprintf('  %-12.0f  %-14.2e  %-12.2e  %-14.4f  %-12.4f  %s\n', ...
        R, err_a, REL_TOL_ANALYTIC, err_sg, REL_TOL_SG, verdict(pass));
end

fprintf('\n');
if allPass
    fprintf('PASS: curvatureKinematicEBR returns 1/R for circles (analytic and SG inputs)\n');
else
    fprintf('FAIL: curvature formula exceeded tolerance on circle\n');
end

end

function s = verdict(pass)
if pass, s = 'PASS'; else, s = 'FAIL'; end
end
