function [beta, torsionExponent, yGain, stats] = regressCurvatureTorsion3D_v001(speed, curvature, torsion, regressType, LMseeds)
% REGRESSCURVATURETORSION3D_V001  3D generalisation of regressDataEBR.m
% for the curvature+torsion power law: v = k * curvature^(-beta) * |torsion|^(-torsionExponent).
%
% This is the full equi-affine-style regression (Pollick, Maoz, Handzel,
% Giblin, Sapiro & Flash 2009, Cortex 45(3):325-339), as distinct from the
% naive single-curvature 3D extension (Schaal & Sternad 2001, Exp Brain Res
% 136:60-72) that ignores torsion entirely and is mathematically shown
% there to fail once movement is non-planar. Both exponents are estimated
% freely (not fixed at the classical 1/3, 1/6), consistent with this
% project's own beta_obs convention throughout (regressDataEBR.m).
%
% Population validation on ULTRA-MoCap Crossbody Reach (2026-07-21,
% checkPlanarityTorsionDiagnostics_ULTRA_v002.m, N=217 reaches, 7
% participants): 88% of reaches show a torsion exponent significantly
% different from 0, median ~0.07 (below the classical 1/6, but never
% absent), motivating this as the primary 3D pipeline rather than the
% three-plane stopgap. See checkBeta3DRecovery_ULTRA_v001.m for the
% direct naive-vs-full beta comparison this function is built for.
%
%% inputs
% speed        - tangential speed, N x 1 (matches regressDataEBR's "responses")
% curvature    - N x 1, from curvatureTorsion3D_v002.m
% torsion      - N x 1, from curvatureTorsion3D_v002.m (signed; this
%                function takes abs() internally, matching the |torsion|
%                convention in the equi-affine speed formula)
% regressType  - 3: OLS, linear regression in log-log space (fitlm, both
%                    predictors, matches regressDataEBR case 3 in spirit)
%                4: nonlinear Levenberg-Marquardt (fitnlm, matches
%                    regressDataEBR case 4)
%                5: nonlinear Iteratively Reweighted Least Squares,
%                    bisquare (fitnlm, matches regressDataEBR case 5)
% LMseeds      - [yGain0, beta0, torsionExponent0] initial guess, required
%                for regressType 4 and 5
%
%% outputs
% beta            - curvature exponent (positive value; matches
%                    regressDataEBR's beta sign convention)
% torsionExponent - torsion exponent (positive value)
% yGain           - velocity gain factor k
% stats           - fitted model object (regressType 4/5) or fitlm object
%                    (regressType 3), for downstream CI/significance use
%
% Created 2026-07-21
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk

arguments
    speed       (:,1) double
    curvature   (:,1) double
    torsion     (:,1) double
    regressType (1,1) double
    LMseeds     (1,:) double = []
end

absTorsion = abs(torsion);
stats = [];

switch regressType

    case 3 % OLS: linear regression in log-log space, both predictors
        logResponse = log(speed);
        logPredictors = [log(curvature), log(absTorsion)];
        model = fitlm(logPredictors, logResponse, "linear");
        beta = -model.Coefficients.Estimate(2);
        torsionExponent = -model.Coefficients.Estimate(3);
        yGain = exp(model.Coefficients.Estimate(1));
        stats = model;

    case {4, 5} % nonlinear: Levenberg-Marquardt (4) or IRLS bisquare (5)
        if isempty(LMseeds) || numel(LMseeds) < 3
            error("regressCurvatureTorsion3D:MissingSeeds", "%s", ...
                "regressType 4 and 5 require LMseeds = [yGain0, beta0, torsionExponent0].");
        end

        modelFunction = @(coefficients, predictors) coefficients(1) .* ...
            predictors(:, 1) .^ (-coefficients(2)) .* predictors(:, 2) .^ (-coefficients(3));
        predictors = [curvature, absTorsion];
        initialGuess = LMseeds(1:3);

        lastwarn("");
        if regressType == 4
            model = fitnlm(predictors, speed, modelFunction, initialGuess);
        else
            fitOptions = statset("nlinfit");
            fitOptions.RobustWgtFun = "bisquare";
            model = fitnlm(predictors, speed, modelFunction, initialGuess, Options=fitOptions);
        end

        if isempty(lastwarn)
            beta = model.Coefficients.Estimate(2);
            torsionExponent = model.Coefficients.Estimate(3);
            yGain = model.Coefficients.Estimate(1);
            stats = model;
        else
            beta = NaN;
            torsionExponent = NaN;
            yGain = NaN;
        end

    otherwise
        error("regressCurvatureTorsion3D:UnsupportedRegressType", "%s", ...
            sprintf("regressType %d not supported; use 3 (OLS), 4 (LM), or 5 (IRLS).", regressType));

end

end
