function test_differentiateKinematicsEBR_v001()
% test_differentiateKinematicsEBR_v001  Unit test for differentiateKinematicsEBR.
%
% Given a signal with analytically known derivatives, both pipeline filter
% types must recover velocity and acceleration within tolerance.
%
% Test signal: quintic polynomial x(t) = t^5 - 2t^3 + t
%   x'(t)  = 5t^4 - 6t^2 + 1    (velocity)
%   x''(t) = 20t^3 - 12t         (acceleration)
%
% A polynomial is recovered exactly by SG (up to filter order) and closely
% by BWFD in the passband. Interior points only (edges have filter artefacts).
%
% Filter types tested:
%   2 = BWFD  [order=2, cutoff=10Hz, zerolag=1]  -- legacy pipeline
%   6 = SG    [order=4, framelen=17]              -- vetted pipeline (fs-scaled)
%
% Pass criteria (interior points):
%   RMS(velocity_rec - velocity_analytic) / RMS(velocity_analytic)   < VEL_TOL_REL
%   RMS(accel_rec    - accel_analytic)    / RMS(accel_analytic)       < ACC_TOL_REL
%
% SG applied to a degree-4 polynomial should be near-exact; BWFD introduces
% some filtering error so its tolerance is looser.
%
% Run from src/ (tests/ must be on the path).
%
% Dagmar Scott Fraser | PowerLawSimulationPreReg | 2026

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions')));

%% ---- Configuration ----
FS          = 120;          % Hz
DURATION    = 3.0;          % seconds -- enough for filter to settle
MARGIN_FRAC = 0.20;         % trim more aggressively for BWFD edge effects

% Tolerances
VEL_TOL_REL_SG   = 0.005;   % SG near-exact on polynomial
VEL_TOL_REL_BWFD = 0.05;    % BWFD attenuates high-freq content
ACC_TOL_REL_SG   = 0.02;
ACC_TOL_REL_BWFD = 0.15;

fprintf('=== test_differentiateKinematicsEBR_v001 ===\n');
fprintf('Signal: quintic polynomial, fs=%d Hz, duration=%.1f s\n\n', FS, DURATION);

%% ---- Generate test signal ----
% Scale time to [-1, 1] so polynomial values are O(1) not O(N^5)
N  = round(DURATION * FS);
t  = linspace(-1, 1, N)';
dt = t(2) - t(1);   % normalised dt (not 1/fs -- derivatives are w.r.t. t, not seconds)

% Signal in normalised time; we'll compare derivatives w.r.t. normalised t
x = t.^5 - 2*t.^3 + t;
y = -t.^4 + 3*t.^2 - t;    % different polynomial for y

% Analytical derivatives w.r.t. normalised time
vx_true  =  5*t.^4 - 6*t.^2 + 1;
vy_true  = -4*t.^3 + 6*t;
ax_true  =  20*t.^3 - 12*t;
ay_true  = -12*t.^2 + 6;

% differentiateKinematicsEBR differentiates w.r.t. real time (1/fs).
% To compare, scale analytic derivatives by chain rule: d/dt_real = d/dt_norm * (dt_norm/dt_real)
% dt_norm = 2/N * (1/fs) ... actually simpler: pass t_real to the function and rescale analytically.
% Easier: just use t_real = (0:N-1)/fs and recompute analytically in real time.
t_real  = (0:N-1)' / FS;
t_norm  = linspace(-1, 1, N)';   % same as t above
scale   = 2 / (DURATION);        % dt_norm/dt_real = 2/(N-1) * fs = 2/DURATION (approx)

vx_analytic = vx_true * scale;
vy_analytic = vy_true * scale;
ax_analytic = ax_true * scale^2;
ay_analytic = ay_true * scale^2;

%% ---- Run both filter types ----
filters = struct( ...
    'name',   {'BWFD (type 2)', 'SG (type 6)'}, ...
    'type',   {2,               6             }, ...
    'params', {[2, 10, 1],      [4, 17]       }, ...
    'velTol', {VEL_TOL_REL_BWFD, VEL_TOL_REL_SG}, ...
    'accTol', {ACC_TOL_REL_BWFD, ACC_TOL_REL_SG});

margin = round(N * MARGIN_FRAC);
idx    = (margin + 1):(N - margin);

allPass = true;
fprintf('  %-18s  %-10s  %-10s  %-10s  %-10s  %s\n', ...
    'filter', 'vel_rms_rel', 'vel_tol', 'acc_rms_rel', 'acc_tol', 'verdict');
fprintf('  %s\n', repmat('-', 1, 72));

for fi = 1:numel(filters)
    f = filters(fi);
    try
        [dx, dy] = differentiateKinematicsEBR(x, y, f.type, f.params, FS);
    catch ME
        fprintf('  %-18s  FAIL (exception: %s)\n', f.name, ME.message);
        allPass = false;
        continue;
    end

    % Interior points only
    vx_rec = dx(idx, 2); vy_rec = dy(idx, 2);
    ax_rec = dx(idx, 3); ay_rec = dy(idx, 3);
    va = vx_analytic(idx); vya = vy_analytic(idx);
    aa = ax_analytic(idx); aya = ay_analytic(idx);

    % RMS error relative to RMS of true signal (combined x+y)
    velErr = rms([vx_rec - va; vy_rec - vya]) / rms([va; vya]);
    accErr = rms([ax_rec - aa; ay_rec - aya]) / rms([aa; aya]);

    velPass = velErr < f.velTol;
    accPass = accErr < f.accTol;
    pass    = velPass && accPass;
    if ~pass, allPass = false; end

    fprintf('  %-18s  %-10.4f  %-10.4f  %-10.4f  %-10.4f  %s\n', ...
        f.name, velErr, f.velTol, accErr, f.accTol, verdict(pass));
end

fprintf('\n');
if allPass
    fprintf('PASS: both filter types recover velocity and acceleration within tolerance\n');
else
    fprintf('FAIL: one or more filter types exceeded derivative tolerance\n');
end

end

function s = verdict(pass)
if pass, s = 'PASS'; else, s = 'FAIL'; end
end
