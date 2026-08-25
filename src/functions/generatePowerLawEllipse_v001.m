function [x, y, kappa_a, v_a, kappa_m] = generatePowerLawEllipse_v001( ...
        a, b, f0, FS, betaGen, varargin)
% generatePowerLawEllipse_v001  Analytical power-law ellipse generator.
%
% Generates a planar ellipse trajectory in which tangential velocity exactly
% enforces the power law  v = VGF * kappa^(-betaGen)  via the arc-length
% parameterisation used in runLoopClosureFftnoise_v007 processOneTrial_local.
% The generated (x, y) are the ONLY outputs that carry noise after the caller
% adds position-level noise.  The analytical outputs (kappa_a, v_a) are exact
% and noise-free regardless of sigma — they are the noiseless ground truth.
%
% ANALYTICAL OUTPUTS (exact, no differentiation, no edge effects):
%   kappa_a  -- true curvature from ellipse formula at each sample angle phiA
%               kappa = (a*b) / (a^2*sin^2(phiA) + b^2*cos^2(phiA))^(3/2)
%   v_a      -- true speed: VGF * kappa_a^(-betaGen)   (enforced by construction)
%   VGF      -- computed analytically from arc-length conservation:
%               VGF = f0 * ArcLength / <kappa^(-betaGen)>_arc
%
% OPTIONAL OUTPUT (estimate, sampling-density dependent):
%   kappa_m  -- Menger curvature: circumradius of triangle formed by three
%               consecutive position samples.  Uses circular indexing at
%               boundaries (valid for a closed periodic ellipse).
%               Compare kappa_m vs kappa_a to reveal pipeline sampling error
%               at sigma=0, BEFORE any noise is added.
%
% The decomposition this enables (new in session 38, Finding #84):
%   beta_gen - beta_obs_rec  =  (beta_gen - beta_analytic)
%                            +  (beta_analytic - beta_obs_rec)
%   total gap                   geometric/pipeline gap    noise EIV gap
% where beta_analytic is obtained by regressing v_a vs kappa_a with the
% same regression method (OLS/LMLS/IRLS) but NO differentiation noise.
%
% INPUTS (required):
%   a        semi-major axis (mm, or any consistent unit)
%   b        semi-minor axis (mm)
%   f0       fundamental frequency (Hz)
%   FS       sampling rate (Hz)
%   betaGen  power law exponent (0=constant velocity; 1/3=canonical)
%
% NAME-VALUE OPTIONS:
%   'nCycles'  number of complete cycles (default: 10)
%   'M'        override total sample count; if set, nCycles is ignored
%   'theta'    ellipse rotation angle in radians (default: 0)
%   'nPhi'     phi integration grid resolution (default: 10000)
%   'Menger'   logical — compute kappa_m (default: false)
%
% OUTPUTS:
%   x, y     M x 1 position (same units as a, b), ellipse in x-y plane,
%            rotated by theta, centred at origin
%   kappa_a  M x 1 analytical curvature at each sample (exact)
%   v_a      M x 1 analytical speed at each sample (exact)
%   kappa_m  M x 1 Menger curvature estimate (only if 'Menger', true;
%            otherwise empty)
%
% EXAMPLE:
%   % Maoz ellipse, 10 cycles, beta=1/3, compare Menger vs analytical
%   [x,y,ka,va,km] = generatePowerLawEllipse_v001(50, 25, 1, 100, 1/3, ...
%       'nCycles', 10, 'Menger', true);
%   figure; plot(phi, ka, phi, km); legend('Analytical','Menger');
%
% Fraser, D.S. (2026)  v001

    %% --- Parse inputs ---------------------------------------------------
    p = inputParser;
    addRequired(p, 'a',       @(x) isscalar(x) && x > 0);
    addRequired(p, 'b',       @(x) isscalar(x) && x > 0);
    addRequired(p, 'f0',      @(x) isscalar(x) && x > 0);
    addRequired(p, 'FS',      @(x) isscalar(x) && x > 0);
    addRequired(p, 'betaGen', @(x) isscalar(x) && isfinite(x));
    addParameter(p, 'nCycles', 10,    @(x) isscalar(x) && x > 0);
    addParameter(p, 'M',       [],    @(x) isempty(x) || (isscalar(x) && x > 0));
    addParameter(p, 'theta',   0,     @isscalar);
    addParameter(p, 'nPhi',    10000, @(x) isscalar(x) && x >= 100);
    addParameter(p, 'Menger',  false, @islogical);
    parse(p, a, b, f0, FS, betaGen, varargin{:});
    opt = p.Results;

    if ~isempty(opt.M)
        M = round(opt.M);
    else
        M = round(opt.nCycles * FS / f0);
    end
    if M < 4
        error('generatePowerLawEllipse_v001:TooFewSamples', '%s', ...
            sprintf('M=%d is too small (need >= 4).', M));
    end

    %% --- Arc-length parameterisation (same as v007 processOneTrial_local) ---
    nPhi  = opt.nPhi;
    phi   = linspace(0, 2*pi, nPhi)';
    dsdp  = sqrt((a .* sin(phi)).^2 + (b .* cos(phi)).^2);  % |d(ellipse)/dphi|
    den   = (a^2 .* sin(phi).^2 + b^2 .* cos(phi).^2).^1.5;
    kp_phi = (a * b) ./ max(den, eps);                       % kappa(phi), analytic

    % Weight that enforces v = VGF * kappa^(-betaGen)
    wt    = dsdp .* kp_phi .^ betaGen;
    cumT  = cumsum(wt);
    cumT  = cumT / cumT(end);                                 % normalised to [0,1]

    % Sample times mapped to normalised phase within one cycle
    tVec  = (0 : M-1)' / FS;
    tN    = min(mod(tVec * f0, 1), 1 - 1e-9);
    phiA  = interp1(cumT, phi, tN, 'linear', 'extrap');      % M x 1 sample angles

    %% --- Generate (x, y) in local frame, then rotate ----------------------
    xLoc  = a * cos(phiA);
    yLoc  = b * sin(phiA);

    ct = cos(opt.theta);  st = sin(opt.theta);
    x  = ct * xLoc - st * yLoc;
    y  = st * xLoc + ct * yLoc;

    %% --- Analytical curvature at sample angles ----------------------------
    denA    = (a^2 .* sin(phiA).^2 + b^2 .* cos(phiA).^2).^1.5;
    kappa_a = (a * b) ./ max(denA, eps);                     % M x 1, exact

    %% --- Analytical VGF from arc-length conservation ----------------------
    % Total arc length: L = integral(dsdp dphi) over [0, 2pi]
    % Mean curvature term: <kappa^(-beta)>_arc = integral(dsdp * kappa^(-beta) dphi) / L
    % Then VGF = f0 * L / (L / <v>) ... simplifies to:
    % <v> = L * f0  (total arc / period),  VGF = <v> / <kappa^(-beta)>_arc
    dphi     = 2*pi / nPhi;
    arcLen   = sum(dsdp) * dphi;                             % total arc length
    % VGF from arc-length conservation: <v> = f0*L; VGF = <v>/<kappa^(-beta)>_arc
    VGF      = f0 * arcLen / (sum(dsdp .* kp_phi.^(-betaGen)) * dphi);

    %% --- Analytical speed at sample angles --------------------------------
    v_a = VGF .* kappa_a .^ (-betaGen);                     % M x 1, exact

    %% --- Menger curvature (optional) --------------------------------------
    kappa_m = [];
    if opt.Menger
        kappa_m = mengerCurvature_local(x, y);
    end

end

%% =========================================================================
function km = mengerCurvature_local(x, y)
% Menger curvature at each sample using circular indexing for the boundary.
% km(k) = circumradius of triangle (p_{k-1}, p_k, p_{k+1}).
% For a closed curve, indices wrap: p_0 = p_M, p_{M+1} = p_1.
%
% Formula: km = 4*A / (d01 * d12 * d20)
% where A = signed triangle area, dij = distance between points i and j.
% More numerically stable via cross product formulation.

    M  = numel(x);
    km = NaN(M, 1);

    % Circular index helpers
    prev = [M, 1:M-1];
    next = [2:M, 1];

    for k = 1:M
        p0 = [x(prev(k)), y(prev(k))];
        p1 = [x(k),       y(k)      ];
        p2 = [x(next(k)), y(next(k))];

        % Twice the signed triangle area (cross product z-component)
        twoA = (p1(1)-p0(1))*(p2(2)-p0(2)) - (p2(1)-p0(1))*(p1(2)-p0(2));

        d01  = norm(p1 - p0);
        d12  = norm(p2 - p1);
        d20  = norm(p0 - p2);
        denom = d01 * d12 * d20;

        if denom < eps
            km(k) = 0;
        else
            km(k) = 2 * abs(twoA) / denom;
        end
    end
end
