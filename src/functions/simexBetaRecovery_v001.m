function results = simexBetaRecovery_v001(x, y, fs, alpha, sigma, ...
    diffFilterType, diffFilterParams, regressType, LMseeds, options)
%SIMEXBETARECOVERY_V001 SIMEX diagnostic for velocity-curvature beta recovery.
%
% Implements simulation-extrapolation (SIMEX; Cook & Stefanski, 1994) for
% the coupled-measurement-error problem in power law exponent estimation.
% Velocity and curvature are both calculated from the same position
% stream (x, y), so their measurement error is coupled in the sense of
% Archie (1981) and Stratton, Feustel & Newell (1987): the observed
% regression slope is a reliability-weighted mixture of the true slope
% and the slope the coupled measurement error alone would produce.
% Stratton's own closed-form correction requires the calculating function
% to be well approximated by the first two terms of a Taylor expansion
% (true for their multiplicative case, not for two stages of
% differentiation plus a curvature quotient). Holcomb (1999) showed SIMEX
% generalises Stratton's correction to exactly this harder case: add
% synthetic noise of increasing magnitude lambda, re-estimate beta at
% each level, then extrapolate the beta(lambda) curve back to
% lambda = -1, the notional zero-measurement-error point.
%
% THIS IS A DIAGNOSTIC, NOT A CORRECTION. Per this project's own
% noise-matched pointwise-inversion work, the extrapolation is only
% trustworthy where beta(lambda) is smooth and well-behaved (steep, narrow
% bootstrap CI). Where the forward map already folds (non-monotonic
% beta_gen -> beta_obs), beta(lambda) is expected to be flat with an
% exploded bootstrap CI: SIMEX does not rescue non-identifiability, it
% diagnoses it -- the same failure mode as Holcomb's own join-point
% example (his Fig. 2). Do not report betaSimex as a corrected estimate
% without checking betaSimexCI width against the tolerance you are
% willing to accept, and do not treat a narrow CI here as independent
% confirmation of the pointwise-inversion result: both diagnostics are
% sensitive to the same feature of the forward map (local steepness), so
% agreement is convergent validation using an independently-named method,
% not a second, independent line of evidence.
%
% SYNTAX:
%   results = simexBetaRecovery_v001(x, y, fs, alpha, sigma, ...
%                 diffFilterType, diffFilterParams, regressType, LMseeds)
%   results = simexBetaRecovery_v001(..., LambdaGrid=[0 0.5 1 1.5 2], ...
%                 NReplicates=50, NBootstrap=200)
%
% INPUTS:
%   x, y              - Trajectory coordinates, equal-length vectors, mm.
%   fs                - Sampling frequency, Hz.
%   alpha             - Noise colour exponent (1/f^alpha) of the trial's
%                        own residual, e.g. from iraAlphaSigma_v001. Sets
%                        the spectral shape of the synthetic SIMEX noise
%                        (Xu, 2019, via noiseXu).
%   sigma             - Noise magnitude (std, mm) of the trial's own
%                        residual. Synthetic noise added at level lambda
%                        has std sqrt(lambda)*sigma (Cook & Stefanski,
%                        1994).
%   diffFilterType, diffFilterParams - Passed straight to
%                        differentiateKinematicsEBR (Fraser et al., 2025).
%   regressType, LMseeds - Passed straight to regressDataEBR (Fraser et
%                        al., 2025); regressType 4/5 (LMLS/IRLS) require
%                        LMseeds, e.g. [1 1/3].
%   LambdaGrid   - Name-value, noise-inflation levels, default
%                  [0 0.5 1 1.5 2]. Needs >= 3 points for a quadratic fit.
%   NReplicates  - Name-value, independent noise draws per lambda,
%                  default 50.
%   NBootstrap   - Name-value, bootstrap resamples for the extrapolated
%                  CI, default 200. Set to 0 to skip (fast path).
%
% OUTPUTS:
%   results - struct:
%     .lambda      - the lambda grid used
%     .betaMean    - mean beta_obs at each lambda
%     .betaSD      - SD of beta_obs at each lambda
%     .nValid      - valid (non-NaN) replicate count at each lambda.
%                    regressDataEBR can return NaN on non-convergence;
%                    these are surfaced here, never silently dropped or
%                    substituted.
%     .quadCoeffs  - polyfit(lambda, betaMean, 2) coefficients
%     .betaSimex   - extrapolated beta at lambda = -1
%     .betaSimexCI - [lo hi] percentile bootstrap CI (empty if
%                    NBootstrap == 0)
%
% REFERENCES:
%   Archie JP Jr (1981). Mathematic coupling of data: a common source of
%     error. Ann Surg 193:296-303. doi:10.1097/00000658-198103000-00008
%   Stratton HH, Feustel PJ, Newell JC (1987). Regression of calculated
%     variables in the presence of shared measurement error. J Appl
%     Physiol 62:2083-2093. doi:10.1152/jappl.1987.62.5.2083
%   Cook JR, Stefanski LA (1994). Simulation-extrapolation estimation in
%     parametric measurement error models. J Am Stat Assoc 89:1314-1328.
%     doi:10.1080/01621459.1994.10476871
%   Holcomb JP Jr (1999). Regression with covariates and outcome
%     calculated from a common set of variables measured with error:
%     estimation using the SIMEX method. Stat Med 18:2847-2862.
%     doi:10.1002/(SICI)1097-0258(19991115)18:21<2847::AID-SIM240>3.0.CO;2-V
%   Xu C (2019). An easy algorithm to generate colored noise sequences.
%     Astron J 157:127. doi:10.3847/1538-3881/ab037c
%
% EXAMPLE:
%   % Hickman PLAC centroid (EMPIRICAL_DATASETS.md v004: alpha=5.34,
%   % sigma=7.17mm, fs=133Hz), SG-LMLS pipeline:
%   r = simexBetaRecovery_v001(trial.x, trial.y, 133, 5.34, 7.17, ...
%           6, [4 17], 4, [1 1/3]);
%   fprintf('beta_SIMEX = %.3f [%.3f, %.3f], nValid = %s\n', ...
%       r.betaSimex, r.betaSimexCI(1), r.betaSimexCI(2), mat2str(r.nValid));
%
% See also: differentiateKinematicsEBR, curvatureKinematicEBR,
%           regressDataEBR, noiseXu
%
% Created 2026-09-20. Dagmar Scott Fraser, d.s.fraser@bham.ac.uk

arguments
    x (:,1) double
    y (:,1) double
    fs (1,1) double {mustBePositive}
    alpha (1,1) double
    sigma (1,1) double {mustBeNonnegative}
    diffFilterType (1,1) double {mustBeMember(diffFilterType, [1 2 3 4 5 6 7 8])}
    diffFilterParams (1,:) double
    regressType (1,1) double {mustBeMember(regressType, [1 2 3 4 5])}
    LMseeds (1,:) double = [1 1/3]
    options.LambdaGrid (1,:) double {mustBeNonnegative} = [0 0.5 1 1.5 2]
    options.NReplicates (1,1) double {mustBePositive, mustBeInteger} = 50
    options.NBootstrap (1,1) double {mustBeNonnegative, mustBeInteger} = 200
end

if numel(x) ~= numel(y)
    error('simexBetaRecovery_v001:SizeMismatch', ...
        'x (%d) and y (%d) must be the same length.', numel(x), numel(y));
end
if numel(options.LambdaGrid) < 3
    error('simexBetaRecovery_v001:LambdaGridTooShort', '%s', ...
        'LambdaGrid needs at least 3 points to fit a quadratic extrapolant.');
end

N          = numel(x);
lambdaGrid = options.LambdaGrid;
nLambda    = numel(lambdaGrid);
nRep       = options.NReplicates;

betaDraws = NaN(nLambda, nRep);

for iLam = 1:nLambda
    lam = lambdaGrid(iLam);
    for iRep = 1:nRep
        if lam == 0
            xStar = x;
            yStar = y;
        else
            uX = noiseXu(N, alpha, sqrt(lam) * sigma, fs);
            uY = noiseXu(N, alpha, sqrt(lam) * sigma, fs);
            xStar = x + uX;
            yStar = y + uY;
        end

        [dx, dy] = differentiateKinematicsEBR(xStar, yStar, ...
            diffFilterType, diffFilterParams, fs);

        v     = sqrt(dx(:,2).^2 + dy(:,2).^2);
        kappa = curvatureKinematicEBR(dx(:,2), dy(:,2), dx(:,3), dy(:,3));

        valid = isfinite(v) & isfinite(kappa) & (kappa > 0) & (v > 0);
        if nnz(valid) < 10
            continue  % leaves NaN; counted via nValid, never silently substituted
        end

        try
            betaHat = regressDataEBR(v(valid), kappa(valid), regressType, LMseeds, 0, 0);
            betaDraws(iLam, iRep) = betaHat;
        catch
            % Treated the same as non-convergence: leaves NaN, counted
            % via nValid below, never silently substituted. A single bad
            % replicate (e.g. fitnlm returning Inf/NaN internally) no
            % longer discards the other ~249 replicates for this cell.
        end
    end
end

results.lambda   = lambdaGrid;
results.nValid   = sum(isfinite(betaDraws), 2)';
results.betaMean = mean(betaDraws, 2, 'omitnan')';
results.betaSD   = std(betaDraws, 0, 2, 'omitnan')';

if any(results.nValid < nRep)
    warning('simexBetaRecovery_v001:IncompleteReplicates', ...
        ['Some lambda levels had fewer than NReplicates valid fits: ' ...
         'nValid = %s. Extrapolation uses whatever converged; check ' ...
         'before trusting betaSimex.'], mat2str(results.nValid));
end
if any(results.nValid == 0)
    error('simexBetaRecovery_v001:NoValidReplicates', '%s', ...
        ['At least one lambda level produced zero valid beta estimates; ' ...
         'cannot fit an extrapolant. Inspect the pipeline/noise ' ...
         'combination directly rather than proceeding.']);
end

results.quadCoeffs = polyfit(lambdaGrid, results.betaMean, 2);
results.betaSimex  = polyval(results.quadCoeffs, -1);

if options.NBootstrap > 0
    bootBetaSimex = NaN(options.NBootstrap, 1);
    for iBoot = 1:options.NBootstrap
        bootMean = NaN(1, nLambda);
        for iLam = 1:nLambda
            draws = betaDraws(iLam, isfinite(betaDraws(iLam,:)));
            if isempty(draws)
                continue
            end
            resampled     = draws(randi(numel(draws), 1, numel(draws)));
            bootMean(iLam) = mean(resampled);
        end
        if all(isfinite(bootMean))
            bootCoeffs           = polyfit(lambdaGrid, bootMean, 2);
            bootBetaSimex(iBoot) = polyval(bootCoeffs, -1);
        end
    end
    bootBetaSimex = bootBetaSimex(isfinite(bootBetaSimex));
    if numel(bootBetaSimex) < options.NBootstrap / 2
        warning('simexBetaRecovery_v001:BootstrapUnstable', ...
            ['Fewer than half the bootstrap resamples produced a finite ' ...
             'extrapolant (%d of %d); the reported CI is unreliable and ' ...
             'likely reflects exactly the fold/non-identifiability this ' ...
             'diagnostic exists to detect.'], ...
            numel(bootBetaSimex), options.NBootstrap);
    end
    if isempty(bootBetaSimex)
        results.betaSimexCI = [NaN NaN];
    else
        results.betaSimexCI = prctile(bootBetaSimex, [2.5 97.5]);
    end
else
    results.betaSimexCI = [];
end

end
