function test_regressDataEBR_v001()
% test_regressDataEBR_v001  Unit test for regressDataEBR.
%
% Given synthetic v = VGF * kappa^(-beta) with zero noise, all three
% regression types must recover beta and VGF to within tight tolerances.
% This is a pure unit test: known input, known exact output, no pipeline,
% no noise, no shape generation.
%
% Regression types tested:
%   3 = OLS  (fitlm log-log with intercept)
%   4 = LMLS (fitnlm Levenberg-Marquardt)
%   5 = IRLS (fitnlm iteratively reweighted)
%
% Pass criteria:
%   |beta_rec  - beta_gen | < BETA_TOL
%   |vgf_rec   - vgf_gen  | / vgf_gen < VGF_REL_TOL
%
% Run from src/ (tests/ must be on the path).
%
% Dagmar Scott Fraser | PowerLawSimulationPreReg | 2026

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions')));

BETA_TOL    = 1e-4;
VGF_REL_TOL = 1e-3;
N_POINTS    = 500;     % enough for stable regression

% Test cases: {beta_gen, vgf_gen, description}
cases = {
    1/3,  1.0,   'canonical 1/3 law, VGF=1';
    1/3,  0.3,   'canonical 1/3 law, VGF=0.3';
    0.5,  2.0,   'beta=0.5, VGF=2';
    0.0,  1.5,   'beta=0 (flat), VGF=1.5';
    2/3,  0.8,   'beta=2/3, VGF=0.8';
};

% Curvature range: span the ellipse curvature range used in the simulation
% Ellipse a=canvas/4, b=canvas/6 roughly gives kappa in ~[0.002, 0.015]
kappa = logspace(-3, -1, N_POINTS)';   % [0.001, 0.1] mm^-1, log-spaced

regressTypes  = [3, 4, 5];
regressNames  = {'OLS', 'LMLS', 'IRLS'};

fprintf('=== test_regressDataEBR_v001 ===\n');
fprintf('N=%d, beta_tol=%.0e, vgf_rel_tol=%.0e\n\n', N_POINTS, BETA_TOL, VGF_REL_TOL);

allPass = true;

for ci = 1:size(cases, 1)
    betaGen = cases{ci, 1};
    vgfGen  = cases{ci, 2};
    descr   = cases{ci, 3};

    % Exact synthetic data: v = VGF * kappa^(-beta)
    v = vgfGen * kappa.^(-betaGen);

    fprintf('  Case %d: %s\n', ci, descr);
    fprintf('  %-6s  %-12s  %-12s  %-12s  %-12s  %s\n', ...
        'type', 'beta_rec', 'beta_err', 'vgf_rec', 'vgf_err_rel', 'verdict');
    fprintf('  %s\n', repmat('-', 1, 70));

    for ri = 1:numel(regressTypes)
        rt   = regressTypes(ri);
        name = regressNames{ri};

        try
            switch rt
                case 3
                    [betaRec, vgfRec, ~] = regressDataEBR(v, kappa, 3, [], 0);
                case 4
                    [betaRec, vgfRec, ~] = regressDataEBR(v, kappa, 4, [vgfGen, betaGen], 0, 0);
                case 5
                    [betaRec, vgfRec, ~] = regressDataEBR(v, kappa, 5, [vgfGen, betaGen], 0, 1);
            end
        catch ME
            fprintf('  %-6s  FAIL (exception: %s)\n', name, ME.message);
            allPass = false;
            continue;
        end

        betaErr  = abs(betaRec - betaGen);
        vgfErr   = abs(vgfRec  - vgfGen) / vgfGen;
        betaPass = betaErr  < BETA_TOL;
        vgfPass  = vgfErr   < VGF_REL_TOL;
        pass     = betaPass && vgfPass;
        if ~pass, allPass = false; end

        fprintf('  %-6s  %-12.6f  %-12.2e  %-12.6f  %-12.2e  %s\n', ...
            name, betaRec, betaErr, vgfRec, vgfErr, verdict(pass));
    end
    fprintf('\n');
end

if allPass
    fprintf('PASS: regressDataEBR recovers beta and VGF exactly on noise-free data\n');
else
    fprintf('FAIL: one or more regression types failed tolerance check\n');
end

end

function s = verdict(pass)
if pass, s = 'PASS'; else, s = 'FAIL'; end
end
