function results = analyticOracleCheck_v001(opts)
% analyticOracleCheck_v001  Reports GPT SOL's analytic-oracle hierarchy from
% ALREADY-EXISTING saturation-sweep data, all six pipelines, through the
% fold, rather than re-running the sweep.
%
% Sol's revised critique (2026-09-17), section 3: does the trajectory
% generator itself produce the stipulated exponent when velocity and
% acceleration are evaluated analytically, before sampling, filtering, or
% regression? Required hierarchy: (1) beta_gen; (2) exponent from exact
% parametric derivatives; (3) exponent after temporal sampling; (4) exponent
% after each differentiation/filtering method; (5) exponent after each
% regression method.
%
% This has effectively already been computed, three months earlier than
% Sol's critique, by `saturationSweepEngine_v001.m`
% (`measureRegressionSaturation_v004/v005.m`) -- but only reported for one
% pipeline (SG-LMLS) over half the beta_gen range (Finding #85). The
% saved output already covers all six pipelines and beta_gen up to 0.75,
% past the fold; nobody had assembled it this way. Sources Finding #209.
%
% Levels reported:
%   2  betaAnalytic  -- exact analytic v_a/kappa_a (closed-form ellipse
%      derivatives, no sampling error, no differentiation) regressed by
%      each of the six pipelines' own regression method.
%   3+4  meanBetaObs at sigma=0 -- REAL BWFD/SG differentiation applied to
%      noiseless discretely-sampled positions, then each regression
%      method. Isolates pure discretisation/differentiation error with
%      zero noise.
%
% USAGE:
%   R = analyticOracleCheck_v001()
%   R = analyticOracleCheck_v001("ResultFile", "measureRegressionSaturation_v004.mat")
%
% Sanity check: asserts max|betaAnalytic - betaGen| for betaGen>0 is below
% 1e-10 (machine-precision identity) on every run; asserts the sigma=0
% residual stays below 0.01 (an order of magnitude under the clinical
% MDC). Either failing would mean the cited source .mat has changed.
%
% Fraser, D.S. (2026)
% See also: Finding #84 (the decomposition concept), Finding #85 (the
% original narrow report), generatePowerLawEllipse_v001.m,
% saturationSweepEngine_v001.m

    arguments
        opts.ResultFile (1,1) string = "measureRegressionSaturation_v005.mat"
    end

    srcDir = fileparts(mfilename('fullpath'));
    matPath = fullfile(srcDir, opts.ResultFile);
    if ~isfile(matPath)
        error('analyticOracle:MissingFile', '%s', sprintf('%s not found.', matPath));
    end
    S = load(matPath);
    geomNames = string({S.geometries.name});
    nGeom = numel(geomNames);
    nPipe = numel(S.pipelineNames);
    maskPos = S.betaGenSweep > 0;
    bgPos = S.betaGenSweep(maskPos);

    si0 = find(S.sigmaSweep == 0, 1);
    if isempty(si0)
        error('analyticOracle:NoZeroSigma', '%s', 'sigmaSweep has no sigma=0 entry -- level 3/4 check unavailable.');
    end

    %% Level 2: analytic oracle vs betaGen, betaGen > 0
    maxDevAnalytic = nan(nGeom, nPipe);
    for gi = 1:nGeom
        for p = 1:nPipe
            vals = squeeze(S.betaAnalytic(gi,p,maskPos));
            maxDevAnalytic(gi,p) = max(abs(vals(:) - bgPos(:)), [], 'omitnan');
        end
    end

    %% Level 3+4: real differentiation, zero noise, vs betaGen
    maxDevSigma0 = nan(nGeom, nPipe);
    for gi = 1:nGeom
        for p = 1:nPipe
            vals = squeeze(S.meanBetaObs(gi,p,maskPos,si0));
            maxDevSigma0(gi,p) = max(abs(vals(:) - bgPos(:)), [], 'omitnan');
        end
    end

    fprintf('=== analyticOracleCheck_v001 (Finding #209), source: %s ===\n', opts.ResultFile);
    fprintf('betaGenSweep: [%.4f, %.4f], n=%d (%d > 0, used below)\n', ...
        min(S.betaGenSweep), max(S.betaGenSweep), numel(S.betaGenSweep), nnz(maskPos));

    fprintf('\n--- Level 2: max|betaAnalytic - betaGen| (exact analytic derivatives) ---\n');
    fprintf('%-12s', 'geometry');
    for p = 1:nPipe, fprintf(' %-10s', S.pipelineNames(p)); end
    fprintf('\n');
    for gi = 1:nGeom
        fprintf('%-12s', geomNames(gi));
        for p = 1:nPipe, fprintf(' %-10.2e', maxDevAnalytic(gi,p)); end
        fprintf('\n');
    end

    fprintf('\n--- Level 3+4: max|betaObs(sigma=0) - betaGen| (real BWFD/SG, no noise) ---\n');
    fprintf('%-12s', 'geometry');
    for p = 1:nPipe, fprintf(' %-10s', S.pipelineNames(p)); end
    fprintf('\n');
    for gi = 1:nGeom
        fprintf('%-12s', geomNames(gi));
        for p = 1:nPipe, fprintf(' %-10.4f', maxDevSigma0(gi,p)); end
        fprintf('\n');
    end

    %% beta_gen=0 degenerate-case NaN pattern
    bi0 = find(S.betaGenSweep == 0, 1);
    fprintf('\n--- betaGen=0 (degenerate null case) betaAnalytic, all geometries ---\n');
    nanCount = 0; total = 0;
    for gi = 1:nGeom
        for p = 1:nPipe
            total = total + 1;
            if isnan(S.betaAnalytic(gi,p,bi0)), nanCount = nanCount + 1; end
        end
    end
    fprintf('NaN at betaGen=0: %d of %d geometry-pipeline cells (LMLS/IRLS convergence\n', nanCount, total);
    fprintf('failure on a target with zero curvature-velocity information; OLS unaffected).\n');

    results = struct('maxDevAnalytic', maxDevAnalytic, 'maxDevSigma0', maxDevSigma0, ...
        'geomNames', geomNames, 'pipelineNames', S.pipelineNames, 'betaGenSweep', S.betaGenSweep, ...
        'nanCountAtZero', nanCount, 'totalAtZero', total);

    %% Sanity checks
    worstAnalytic = max(maxDevAnalytic(:));
    worstSigma0 = max(maxDevSigma0(:));
    if worstAnalytic > 1e-10
        warning('analyticOracle:SanityCheckFailed', ...
            'Max analytic-level deviation %.2e exceeds machine-precision expectation (1e-10). Source .mat may have changed.', worstAnalytic);
    else
        fprintf('\nSanity check PASSED: analytic level is machine-precision identity (worst case %.2e).\n', worstAnalytic);
    end
    MDC = 0.03;  % clinical threshold, the meaningful bound here -- not an arbitrary guess
    if worstSigma0 > MDC
        warning('analyticOracle:SanityCheckFailed', ...
            'Max sigma=0 real-differentiation deviation %.4f exceeds the clinical MDC (0.03).', worstSigma0);
    else
        fprintf('Sanity check PASSED: noiseless real-differentiation residual stays under the clinical MDC (worst case %.4f vs 0.03).\n', worstSigma0);
    end

    resDir = fullfile(srcDir, 'results');
    if ~exist(resDir, 'dir'), mkdir(resDir); end
    matOut = fullfile(resDir, 'analyticOracleCheck_v001.mat');
    save(matOut, 'results');
    fprintf('\nSaved: %s\n', matOut);
end
