function out = frenetKinematics_v001(x, y, fs, opts)
% FRENETKINEMATICS_V001  Frenet-Serret kinematics from a planar trajectory.
%
% Implements the Huh & Sejnowski (2015) complex-log method. Two successive
% calls to Smooth_Differentiate yield velocity and angular speed; signed
% curvature follows as kappa = (dtheta/dt) / v.
%
% SYNTAX
%   out = frenetKinematics_v001(x, y, fs)
%   out = frenetKinematics_v001(x, y, fs, N=9, m=1)
%   out = frenetKinematics_v001(x, y, fs, SpeedFloor=50, KRange=[-0.03 0.25])
%
% INPUTS
%   x, y       (L x 1)  position vectors in consistent spatial units (mm or px)
%   fs         (1 x 1)  sampling rate [Hz]
%                       Infer from timestamps: fs = 1 / mean(diff(t))
%
% PARAMETERS (name-value)
%   N          (1 x 1)  Smooth_Differentiate filter width; odd integer (default 9)
%   m          (1 x 1)  accuracy order; 1 <= m < (N-1)/2    (default 1)
%   SpeedFloor (1 x 1)  minimum speed for iValid [units/s]  (default 0)
%   KRange     (1 x 2)  [kMin kMax] curvature window for iValid (default [-Inf Inf])
%
% OUTPUT  struct -- all vectors length M = L - 2*(N-1):
%   pos        (M x 2)  trimmed position [x y]
%   vel        (M x 2)  velocity [xdot ydot]   (units/s)
%   v          (M x 1)  arc-length speed = sqrt(xdot^2 + ydot^2)
%   k          (M x 1)  signed curvature kappa = dtheta/ds  (1/unit)
%   theta      (M x 1)  accumulated Frenet-Serret tangent angle (rad)
%   angSpeed   (M x 1)  angular speed omega = dtheta/dt     (rad/s)
%   logK       (M x 1)  log(|k|)
%   logV       (M x 1)  log(v)
%   iValid     (M x 1)  logical: passes SpeedFloor and KRange thresholds
%   fs         (1 x 1)  sample rate used
%   N          (1 x 1)  filter width used
%   m          (1 x 1)  accuracy order used
%   nEdge      (1 x 1)  samples removed from each end (= N - 1)
%
% ALGORITHM
%   Stage 1  [pos1, vel1] = Smooth_Differentiate([x;y], N, m, fs)
%            logVelCx     = log(vel1(1,:) + 1i*vel1(2,:))
%            speedRaw     = exp(real(logVelCx))       [= |vel|]
%            thetaRaw     = unwrap(imag(logVelCx))    [= atan2(ydot, xdot)]
%   Stage 2  [theta, angSpeed] = Smooth_Differentiate(thetaRaw, N, m, fs)
%            k = angSpeed ./ v_trimmed
%
% TRIMMING
%   Each Smooth_Differentiate call removes M = (N-1)/2 samples per end.
%   Two calls: nEdge = N-1 per end; 2*(N-1) total samples lost.
%   Minimum input length: 4*M + 1  (= 2*(N-1) + 1).
%
% NOTE ON iValid
%   The mask is computed but NOT applied internally. Pass out.iValid to
%   downstream functions (e.g. fresSRpowerLaw_v001) as the quality gate.
%
% DEPENDENCY   Smooth_Differentiate.m  (Huh 2011, UCSD / Salk Institute)
%
% REFERENCE
%   Huh D & Sejnowski TJ (2015) Spectrum of power laws for curved hand
%   movements. PNAS 112(29): E3950-E3958.
%
% SEE ALSO   fresSRpowerLaw_v001, Smooth_Differentiate
%
% PROVENANCE
%   Distilled from Analysis_fnc5.m / Analysis_fnc5_vDag.m (Huh / Cook Lab
%   / Fraser). Part of the 7th toolchain (Frenet-Serret spectral power law)
%   for the PowerLawSimulationPreReg empirical validation pipeline.
%
% AUTHOR  D.S. Fraser, Centre for Human Brain Health, Univ. of Birmingham, 2026.

arguments
    x  (:,1) double
    y  (:,1) double
    fs (1,1) double {mustBePositive}
    opts.N          (1,1) double = 9
    opts.m          (1,1) double = 1
    opts.SpeedFloor (1,1) double = 0
    opts.KRange     (1,2) double = [-inf inf]
end

N = opts.N;
m = opts.m;

%% Guards ------------------------------------------------------------------
if mod(N, 2) == 0
    error('frenetKinematics_v001:badN', '%s', ...
        'N must be an odd integer.');
end

M = (N-1) / 2;

if m < 1 || m >= M
    error('frenetKinematics_v001:badM', '%s', ...
        sprintf('m must satisfy 1 <= m < (N-1)/2 = %d.', M));
end

minLen = 4*M + 1;
if numel(x) < minLen
    error('frenetKinematics_v001:tooShort', '%s', ...
        sprintf('Signal has %d samples; N=%d requires at least %d.', ...
        numel(x), N, minLen));
end

%% Stage 1: position -> velocity -------------------------------------------
% Smooth_Differentiate convention: time along 2nd dimension
posRaw = [x(:)'; y(:)'];                        % 2 x L

[pos1, vel1] = Smooth_Differentiate(posRaw, N, m, fs);
% pos1, vel1 : 2 x (L - 2M)

%% Complex-log decomposition -----------------------------------------------
logVelCx = log(vel1(1,:) + 1i*vel1(2,:));
speedRaw = exp(real(logVelCx));                 % |vel| = sqrt(xdot^2+ydot^2)
thetaRaw = unwrap(imag(logVelCx));              % atan2(ydot, xdot), unwrapped

%% Stage 2: tangent angle -> angular speed ---------------------------------
[thetaTrimmed, angSpeedRow] = Smooth_Differentiate(thetaRaw, N, m, fs);
% both 1 x (L - 4M)

%% Align all signals to output length (L - 4M) ----------------------------
pos      = pos1(:, 1+M : end-M)';              % (L-4M) x 2
vel      = vel1(:, 1+M : end-M)';              % (L-4M) x 2
v        = speedRaw(1+M : end-M)';             % (L-4M) x 1
theta    = thetaTrimmed';                      % (L-4M) x 1
angSpeed = angSpeedRow';                       % (L-4M) x 1

%% Curvature ---------------------------------------------------------------
k = angSpeed ./ v;                             % kappa = dtheta/ds

%% Quality mask (computed, not applied) ------------------------------------
iValid = (v  >  opts.SpeedFloor)  & ...
         (k  >  opts.KRange(1))   & ...
         (k  <  opts.KRange(2));

%% Assemble output ---------------------------------------------------------
out.pos      = pos;
out.vel      = vel;
out.v        = v;
out.k        = k;
out.theta    = theta;
out.angSpeed = angSpeed;
out.logK     = log(abs(k));
out.logV     = log(v);
out.iValid   = iValid;
out.fs       = fs;
out.N        = N;
out.m        = m;
out.nEdge    = 2*M;

end
