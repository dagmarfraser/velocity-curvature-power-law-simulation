function [curvature, torsion] = curvatureTorsion3D_v001(velocityX, velocityY, velocityZ, ...
    accelerationX, accelerationY, accelerationZ, jerkX, jerkY, jerkZ)
% CURVATURETORSION3D_V001  3D curvature and torsion from velocity/acceleration/jerk.
%
% 3D generalisation of curvatureKinematicEBR.m (Viviani & Stucchi, 1992;
% adapted from Matic & Gomez-Marin, 2019), extended with torsion via the
% standard Frenet-Serret relation. This is the curvature+torsion pair
% underlying the equi-affine "1/6-power law" of Pollick, Maoz, Handzel,
% Giblin, Sapiro & Flash (2009, Cortex 45(3):325-339), as distinct from
% the naive single-curvature 3D extension tested (and found wanting for
% non-planar movement) by Schaal & Sternad (2001, Exp Brain Res 136:60-72).
%
%   curvature = ||v x a|| / ||v||^3
%   torsion   = (v x a) . j / ||v x a||^2
%
% Verified to reduce exactly (machine precision) to curvatureKinematicEBR.m
% and to torsion == 0 when velocityZ = accelerationZ = jerkZ = 0, i.e. for
% genuinely planar motion embedded in 3D. See
% checkCurvatureTorsion3DSanity_Pilot_v001.m for the live check against a
% real Pilot trial (2026-07-21).
%
%% inputs
% velocity, acceleration, jerk components in x, y, z (each N x 1)
%
%% outputs
% curvature (N x 1)
% torsion   (N x 1) -- NaN/Inf where curvature is zero (torsion is
%                      undefined at zero curvature; this is expected, not
%                      a bug, and callers must not silently zero it out)
%
% Created 2026-07-21
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk

v = [velocityX(:), velocityY(:), velocityZ(:)];
a = [accelerationX(:), accelerationY(:), accelerationZ(:)];
j = [jerkX(:), jerkY(:), jerkZ(:)];

vCrossA = cross(v, a, 2); % N x 3

speedCubed = (sum(v.^2, 2)) .^ 1.5;
curvature = sqrt(sum(vCrossA.^2, 2)) ./ speedCubed;

crossNormSq = sum(vCrossA.^2, 2);
torsion = sum(vCrossA .* j, 2) ./ crossNormSq; % NaN/Inf at zero curvature, by design

end
