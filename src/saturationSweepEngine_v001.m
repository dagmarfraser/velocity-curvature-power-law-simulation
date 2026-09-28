function saturationSweepEngine_v001(geometries, betaGenSweep, sigmaSweep, N_REPS, outFile)
% saturationSweepEngine_v001  Shared engine for saturation surface sweeps.
%
% Sweeps beta_gen x sigma for each geometry in the input struct array,
% computing both the noisy forward map (beta_obs_rec) and the noiseless
% analytical regression (beta_analytic) that enables the decomposition:
%
%   beta_gen - beta_obs_rec  =  (beta_gen - beta_analytic)   [geometric gap]
%                            +  (beta_analytic - beta_obs_rec) [noise EIV gap]
%
% The geometric gap isolates Wann eccentricity, Schaal-Sternad filter
% artefacts, and LMLS original-space weighting — all at sigma=0.
% The noise EIV gap is the pure sigma-dependent EIV contribution.
%
% INPUTS:
%   geometries   struct array, one entry per geometry. Fields:
%                  .name       char label
%                  .a          semi-major axis (mm)
%                  .b          semi-minor axis (mm)
%                  .f0         fundamental frequency (Hz)
%                  .FS         sampling rate (Hz)
%                  .nCycles    number of cycles (determines M)
%                  .alpha      noise spectral exponent for this geometry
%                  .sigmaEmp   empirical sigma (mm) — stored as reference, not used in sweep
%   betaGenSweep vector of beta_gen values to sweep
%   sigmaSweep   vector of sigma values (mm) to sweep
%   N_REPS       number of replicates (parfor over reps)
%   outFile      full path to output .mat file
%
% OUTPUT FILE variables:
%   betaObs      nGeom x nPP x nBeta x nSig x nReps
%   betaAnalytic nGeom x nPP x nBeta  (noiseless analytical regression)
%   meanBetaObs  nGeom x nPP x nBeta x nSig  (mean over reps)
%   sdBetaObs    nGeom x nPP x nBeta x nSig
%   geometries, betaGenSweep, sigmaSweep, N_REPS, pipelineNames
%
% CALLER SCRIPTS:
%   measureRegressionSaturation_v004.m  (Maoz geometry, LC alpha sweep)
%   measureRegressionSaturation_v005.m  (real dataset geometries, LC alpha)
%
% Fraser, D.S. (2026)  v001

    srcDir = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(srcDir, 'functions')));
    addpath(genpath(fullfile(srcDir, 'req')));

    %% --- Pipeline config (matching runLoopClosureFftnoise_v007) ----------
    filterParamsBW  = [2, 10, 1];
    filterParamsSG  = [4, 17];
    filterTypes     = [2, 4];      % BWFD, SG
    regressTypes    = [3, 4, 5];   % OLS, LMLS, IRLS
    pipelineNames   = ["BWFD-OLS","BWFD-LMLS","BWFD-IRLS","SG-OLS","SG-LMLS","SG-IRLS"];
    nFT  = numel(filterTypes);
    nRT  = numel(regressTypes);
    nPP  = nFT * nRT;
    EDGE_CLIP  = 20;
    ALPHA_MAX  = 5.95;

    nGeom = numel(geometries);
    nBeta = numel(betaGenSweep);
    nSig  = numel(sigmaSweep);

    %% --- Allocate outputs ------------------------------------------------
    betaObs      = NaN(nGeom, nPP, nBeta, nSig, N_REPS);
    betaAnalytic = NaN(nGeom, nPP, nBeta);

    %% --- Main loop -------------------------------------------------------
    for gi = 1:nGeom
        g     = geometries(gi);
        M     = round(g.nCycles * g.FS / g.f0);
        alpha = g.alpha;
        ac    = min(alpha, ALPHA_MAX);

        fprintf('\n=== %s | a=%.1f b=%.1f f0=%.3f FS=%d M=%d alpha=%.2f ===\n', ...
            g.name, g.a, g.b, g.f0, g.FS, M, alpha);

        %% --- Pre-generate ellipse templates for each beta_gen -----------
        templates  = zeros(M, 2, nBeta);
        kappa_a_all = zeros(M, nBeta);
        v_a_all     = zeros(M, nBeta);

        for bi = 1:nBeta
            [xP, yP, ka, va] = generatePowerLawEllipse_v001( ...
                g.a, g.b, g.f0, g.FS, betaGenSweep(bi), 'nCycles', g.nCycles);
            templates(:, 1, bi)  = xP;
            templates(:, 2, bi)  = yP;
            kappa_a_all(:, bi)   = ka;
            v_a_all(:, bi)       = va;
        end

        %% --- Analytical regression (sigma=0 geometric gap) --------------
        for bi = 1:nBeta
            ka  = kappa_a_all(:, bi);
            va  = v_a_all(:, bi);
            % Pre-compute OLS seed once — gives scale-appropriate VGF for LMLS/IRLS
            try
                [bOLS, vOLS] = regressDataEBR(va, ka, 3, [1, -1/3], 0, 0);
                seedGlobal = [vOLS, bOLS];
            catch
                seedGlobal = [1, -1/3];
            end
            pp = 0;
            for fi = 1:nFT
                lmSeed = seedGlobal;   % reset to OLS-derived seed per filter type
                for ri = 1:nRT
                    pp = pp + 1;
                    try
                        [b, vgf] = regressDataEBR(va, ka, regressTypes(ri), lmSeed, 0, 0);
                        betaAnalytic(gi, pp, bi) = b;
                        if regressTypes(ri) == 3 && isfinite(b) && isfinite(vgf)
                            lmSeed = [vgf, b];
                        end
                    catch
                    end
                end
            end
        end
        fprintf('  Analytical regression done.\n');

        %% --- Noisy sweep (parfor over reps) -----------------------------
        for si = 1:nSig
            sigma = sigmaSweep(si);
            fprintf('  sigma=%4gmm ... ', sigma);

            for bi = 1:nBeta
                xPure = templates(:, 1, bi);
                yPure = templates(:, 2, bi);

                parfor rep = 1:N_REPS
                    betaRep = NaN(nPP, 1);

                    % Generate noise
                    if sigma == 0
                        xN = zeros(M, 1);  yN = zeros(M, 1);
                    elseif alpha == 0
                        xN = sigma .* randn(M, 1);
                        yN = sigma .* randn(M, 1);
                    else
                        xN = generateCustomNoise_v003(M, ac, sigma, g.FS);
                        yN = generateCustomNoise_v003(M, ac, sigma, g.FS);
                    end

                    x = xPure + xN;
                    y = yPure + yN;

                    % Differentiate
                    kinX = cell(5,1);  kinY = cell(5,1);
                    for ftD = filterTypes
                        try
                            if ftD == 4
                                [kinX{ftD}, kinY{ftD}] = differentiateKinematicsEBR( ...
                                    x, y, ftD, filterParamsSG, g.FS);
                            else
                                [kinX{ftD}, kinY{ftD}] = differentiateKinematicsEBR( ...
                                    x, y, ftD, filterParamsBW, g.FS);
                            end
                        catch
                            kinX{ftD} = [];  kinY{ftD} = [];
                        end
                    end

                    % Regress
                    pp = 0;
                    for fi = 1:nFT
                        ft = filterTypes(fi);
                        if isempty(kinX{ft}), pp = pp + nRT; continue; end
                        vX = kinX{ft}(EDGE_CLIP:end-EDGE_CLIP, 2);
                        vY = kinY{ft}(EDGE_CLIP:end-EDGE_CLIP, 2);
                        aX = kinX{ft}(EDGE_CLIP:end-EDGE_CLIP, 3);
                        aY = kinY{ft}(EDGE_CLIP:end-EDGE_CLIP, 3);
                        spd = sqrt(vX.^2 + vY.^2);
                        kap = curvatureKinematicEBR(vX, vY, aX, aY);
                        lmSeed = [1, -1/3];
                        for ri = 1:nRT
                            pp = pp + 1;
                            try
                                [b, vgf] = regressDataEBR(spd, kap, regressTypes(ri), lmSeed, 0, 0);
                                betaRep(pp) = b;
                                if regressTypes(ri) == 3 && isfinite(b) && isfinite(vgf)
                                    lmSeed = [vgf, b];
                                end
                            catch
                            end
                        end
                    end
                    betaObs(gi, :, bi, si, rep) = betaRep';
                end
            end
            fprintf('done\n');
        end
    end

    %% --- Summarise -------------------------------------------------------
    meanBetaObs = mean(betaObs, 5, 'omitnan');
    sdBetaObs   = std(betaObs,  0, 5, 'omitnan');

    %% --- Print SG-LMLS (pp=5) table per geometry -------------------------
    for gi = 1:nGeom
        g = geometries(gi);
        fprintf('\n=== %s | SG-LMLS (pp=5) beta_obs ===\n', g.name);
        fprintf('  %-8s', 'bg\\sig');
        for si = 1:nSig, fprintf('  sig=%4g', sigmaSweep(si)); end
        fprintf('  | analytic\n  %s\n', repmat('-',1,12+nSig*10+12));
        for bi = 1:nBeta
            fprintf('  %-8.3f', betaGenSweep(bi));
            for si = 1:nSig
                fprintf('  %7.3f ', meanBetaObs(gi, 5, bi, si));
            end
            fprintf('  | %7.3f\n', betaAnalytic(gi, 5, bi));
        end
    end

    %% --- Save ------------------------------------------------------------
    save(outFile, 'betaObs', 'betaAnalytic', 'meanBetaObs', 'sdBetaObs', ...
        'geometries', 'betaGenSweep', 'sigmaSweep', 'N_REPS', 'pipelineNames', '-v7.3');
    fprintf('\nSaved: %s\n', outFile);

end
