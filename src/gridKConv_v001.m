function K = gridKConv_v001(beta)
% gridKConv_v001  The v058 grid's tempo conversion constant K(beta), in px^(1-beta).
%   K = gridKConv_v001(beta)   f0 [Hz] = VGF / K(beta),   VGF = f0 * K(beta).
%
%   The generator sets tangential speed v = VGF * kappa^-beta, so one orbit takes
%   T = closed-integral(ds / v) = closed-integral(kappa^beta ds) / VGF, hence
%   f0 = VGF / K(beta) with K(beta) = closed-integral(kappa^beta ds) on the nu = 2
%   shape exactly as the generator builds it (gridShapeGeometry_v002, canvas px).
%
%   BETA IS REQUIRED. There is deliberately no default. K falls steeply with beta
%   (d ln K / d beta about -5.3 at 1/3: K(0.28) is 33% above K(1/3), K(2/3) is 18% of
%   it), so every caller must say which exponent its conversion assumes, and a
%   caller that passes 1/3 out of habit is visible in its own source. Pass a
%   per-trial or per-dataset exponent where one exists; pass 1/3 only as a stated
%   reference convention.
%
%   Replaces the "K_conv" of plotVGFRecovery_v001/v002, checkVGFSEMCentroid_v002/v003,
%   datasetOperatingPointCoverage_v001 and others, computed as perimeter /
%   mean(kappa^-1/3) = 175.2636. That expression is not the period integral: it is
%   5.28% below K(1/3) = 185.02, and it disagrees with the tempo the grid actually
%   realised (checkGridTempoRange_v001; testGridKConv_v001). Do not use it.
%
%   Validity: analytic and exact for any beta >= 0. Agreement with the DB's realised
%   tempo is asserted only for beta <= 0.40 (median 0.04%, max 0.55%, n = 504
%   nodes); above that, nodes with fewer than about 20 samples per orbit carry
%   discretisation error in the realised duration (up to 6.9%), which is not K's.
%
%   INPUT   beta   numeric, real, finite, >= 0; scalar or array.
%   OUTPUT  K      same size as beta.
%   Errors  (Fail Loud): a call without beta, or with a negative, NaN or complex beta.
% Fraser, D.S. (2026)  v001
    arguments
        beta {mustBeNumeric, mustBeReal, mustBeFinite, mustBeNonnegative}
    end
    persistent kap dsV
    if isempty(kap)
        root = fileparts(fileparts(mfilename("fullpath")));
        G    = gridShapeGeometry_v002(root);
        ds   = vecnorm(diff([G.pathPx; G.pathPx(1, :)]), 2, 2);   % segment i -> i+1
        dsV  = (ds + circshift(ds, 1)) / 2;                        % length attributed to vertex i
        kap  = G.kappaPx(:);
        if any(~isfinite(kap)) || any(kap <= 0)
            error("gridKConv:kappa", "%s", "Non-positive or non-finite curvature in the grid shape");
        end
    end
    K = zeros(size(beta));
    for i = 1:numel(beta)
        K(i) = sum(kap .^ beta(i) .* dsV);
    end
end
