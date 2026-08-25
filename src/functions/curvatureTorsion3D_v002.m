function [curvature, torsion] = curvatureTorsion3D_v002(velocity, acceleration, jerk)
% CURVATURETORSION3D_V002  3D curvature and torsion from velocity/acceleration/jerk.
%
% 3D generalisation of curvatureKinematicEBR.m (Viviani & Stucchi, 1992;
% adapted from Matic & Gomez-Marin, 2019), extended with torsion via the
% standard Frenet-Serret relation. This is the curvature+torsion pair
% underlying the equi-affine "1/6-power law" of Pollick, Maoz, Handzel,
% Giblin, Sapiro & Flash (2009, Cortex 45(3):325-339), as distinct from
% the naive single-curvature 3D extension tested (and found wanting for
% non-planar movement) by Schaal & Sternad (2001, Exp Brain Res 136:60-72).
%
%   curvature = norm(v x a) / norm(v)^3
%   torsion   = dot(v x a, j) / norm(v x a)^2
%
% Verified to reduce exactly (machine precision) to curvatureKinematicEBR.m,
% and to torsion == 0, when the z-columns of velocity, acceleration, and
% jerk are all zero, i.e. for genuinely planar motion embedded in 3D. See
% checkCurvatureTorsion3DSanity_Pilot_v002.m for the live check against a
% real Pilot trial (2026-07-21), run through both BWFD-OLS and SG-IRLS.
%
% v002 (2026-07-21): signature changed from 9 separate scalar columns to
% 3 N-by-3 matrices, to respect the project coding guideline limiting
% function inputs to 6. v001 retained on disk, superseded.
%
%% inputs
% velocity, acceleration, jerk - each an N-by-3 matrix, columns [x y z]
%
%% outputs
% curvature - N-by-1
% torsion   - N-by-1 (NaN/Inf where curvature is zero; torsion is
%             undefined at zero curvature by construction, not a bug,
%             and callers must not silently zero it out)
%
% Created 2026-07-21
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk

arguments
    velocity     (:,3) double
    acceleration (:,3) double
    jerk         (:,3) double
end

velocityCrossAcceleration = cross(velocity, acceleration, 2);

speedCubed = (sum(velocity.^2, 2)) .^ 1.5;
curvature = sqrt(sum(velocityCrossAcceleration.^2, 2)) ./ speedCubed;

crossNormSquared = sum(velocityCrossAcceleration.^2, 2);
torsion = sum(velocityCrossAcceleration .* jerk, 2) ./ crossNormSquared;

end
