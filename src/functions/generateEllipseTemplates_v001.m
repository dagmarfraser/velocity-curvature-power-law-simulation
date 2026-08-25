function templates = generateEllipseTemplates_v001(a_nat, b_nat, theta, f0, fs, M, betaGenVec)
%GENERATEELLIPSETEMPLATES_V001 Synthetic power-law-compliant ellipse trajectories.
%
% Faithful extraction of runLoopClosureFftnoise_v009.m's inline ellipse-
% template-generation block inside processOneTrial_local (was never even
% a named local function there). Algorithm unchanged, packaged standalone.
%
% Arc-length reparametrisation weighted by curvature^beta_gen: for each
% candidate beta_gen in betaGenVec, generates a trajectory tracing an
% ellipse of the given geometry (a_nat, b_nat semi-axes, theta
% orientation) at tempo f0, such that its velocity-curvature relationship
% follows v = VGF * kappa^(-beta_gen) by construction. This is the
% forward-model generator loop closure inverts: ground-truth trajectories
% with KNOWN beta_gen, built at this trial's own exact geometry and
% tempo -- no grid, no snapping, no generic fixed shape.
%
% SYNTAX:
%   templates = generateEllipseTemplates_v001(a_nat, b_nat, theta, f0, fs, M, betaGenVec)
%
% INPUTS:
%   a_nat, b_nat - Ellipse semi-major/semi-minor axes, in the SAME units
%                  as the trial's own trajectory (native units, i.e.
%                  whatever a_nat/b_nat were fit in -- not necessarily mm;
%                  keep consistent with the noise/residual units used
%                  alongside this template).
%   theta        - Ellipse orientation, radians.
%   f0           - Dominant tempo, Hz (e.g. from estimateF0_v001).
%   fs           - Sampling frequency, Hz.
%   M            - Number of output samples per template (trial length).
%   betaGenVec   - Row vector of candidate beta_gen values to generate
%                  templates for.
%
% OUTPUT:
%   templates - M x 2 x numel(betaGenVec) array. templates(:,1,bi) and
%               templates(:,2,bi) are the x,y trajectory for
%               betaGenVec(bi).
%
% See also: estimateF0_v001, findBothBranches_v001
%
% Extracted from: runLoopClosureFftnoise_v009.m (processOneTrial_local,
% "Synthetic ellipse templates" block)
% Fraser, D.S. (2026)

arguments
    a_nat (1,1) double {mustBePositive}
    b_nat (1,1) double {mustBePositive}
    theta (1,1) double
    f0 (1,1) double {mustBePositive}
    fs (1,1) double {mustBePositive}
    M (1,1) double {mustBeInteger, mustBePositive}
    betaGenVec (1,:) double
end

nPhi   = 10000;
phi    = linspace(0, 2*pi, nPhi)';
dsdp   = sqrt((a_nat.*sin(phi)).^2 + (b_nat.*cos(phi)).^2);
den    = (a_nat^2.*sin(phi).^2 + b_nat^2.*cos(phi).^2).^1.5;
kp_phi = (a_nat*b_nat) ./ max(den, eps);
tVec   = (0:M-1)' / fs;

nBeta     = numel(betaGenVec);
templates = zeros(M, 2, nBeta);
for bi = 1:nBeta
    wt   = dsdp .* kp_phi.^betaGenVec(bi);
    cumT = cumsum(wt); cumT = cumT / cumT(end);
    tN   = min(mod(tVec*f0, 1), 1 - 1e-9);
    phiA = interp1(cumT, phi, tN, 'linear', 'extrap');
    xE   = a_nat*cos(phiA);
    yE   = b_nat*sin(phiA);
    templates(:, 1, bi) = xE*cos(theta) - yE*sin(theta);
    templates(:, 2, bi) = xE*sin(theta) + yE*cos(theta);
end

end
