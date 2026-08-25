function results = characteriseNoise_v001(trials, opts)
% characteriseNoise_v001 Estimate noise colour (alpha) and magnitude (sigma)
%
% Characterises the stochastic noise in raw position data from empirical
% movement recordings. Implements the noise characterisation protocol from
% prereg Section 7.3 using multitaper spectral estimation (pmtm) on
% first-differenced position data.
%
% Approach:
%   1. First-difference raw x,y to suppress the movement-dominated low
%      frequencies and expose the noise structure (avoids circular
%      dependence on any smoother choice).
%   2. Multitaper PSD (pmtm) on delta_x, delta_y.
%   3. Fit log10(PSD) vs log10(f) in the noise-dominated frequency band.
%   4. Position noise PSD ~ f^(-alpha), so first-differenced PSD ~
%      f^(2-alpha). Therefore alpha = 2 - fitted_slope.
%   5. Sigma estimated by integrating the back-projected position noise
%      PSD model across the Nyquist bandwidth.
%
% Usage:
%   results = characteriseNoise_v001(trials)
%   results = characteriseNoise_v001(trials, TBP=3, FreqRange=[15 55])
%   results = characteriseNoise_v001(trials, PlotTrials=1:5)
%
% Inputs:
%   trials  - [M x 1] struct array from importDB_*.m with fields:
%             .x, .y, .fs, .subjectID, .trialID, .shape, .database
%
% Name-Value Arguments:
%   TBP        - Time-bandwidth product for pmtm (default: 3; prereg
%                specifies TBP in [2.5, 4], giving 5-7 tapers)
%   FreqRange  - [fLow fHigh] Hz for spectral slope fit (default: [0 0]
%                = auto, which sets fLow = fs/10 and fHigh = fs/2.5)
%   PlotTrials - Vector of trial indices to generate diagnostic plots
%                (default: [] = no plots). Plots PSD with fit overlay.
%   Verbose    - Print progress and summary (default: true)
%
% Output:
%   results - table with one row per trial and columns:
%     trialID, subjectID, shape, database, fs, nSamples,
%     alphaX, alphaY, alphaMean,  (noise colour per axis and mean)
%     sigmaX, sigmaY, sigmaMean,  (noise magnitude in data units)
%     slopeX, slopeY,             (raw fitted slopes on diff spectrum)
%     r2X, r2Y,                   (R-squared of spectral fit)
%     fLow, fHigh                 (frequency range used for fit)
%
% Dependencies:
%   Signal Processing Toolbox (pmtm)
%
% See also: importDB_cook_v002, Stage5_processTrials
%
% Created 2026-03-17  v001
% Dagmar Scott Fraser - d.s.fraser@bham.ac.uk

    arguments
        trials      (1,:) struct
        opts.TBP        (1,1) double {mustBeInRange(opts.TBP, 2, 6)} = 3
        opts.FreqRange  (1,2) double = [0 0]
        opts.PlotTrials (1,:) double = []
        opts.Verbose    (1,1) logical = true
    end

    nTrials = numel(trials);

    % Pre-allocate output columns
    trialIDs   = strings(nTrials, 1);
    subjectIDs = strings(nTrials, 1);
    shapes     = strings(nTrials, 1);
    databases  = strings(nTrials, 1);
    fsVec      = NaN(nTrials, 1);
    nSamps     = NaN(nTrials, 1);
    alphaX     = NaN(nTrials, 1);
    alphaY     = NaN(nTrials, 1);
    alphaMean  = NaN(nTrials, 1);
    sigmaX     = NaN(nTrials, 1);
    sigmaY     = NaN(nTrials, 1);
    sigmaMean  = NaN(nTrials, 1);
    slopeX     = NaN(nTrials, 1);
    slopeY     = NaN(nTrials, 1);
    r2X        = NaN(nTrials, 1);
    r2Y        = NaN(nTrials, 1);
    fLowVec    = NaN(nTrials, 1);
    fHighVec   = NaN(nTrials, 1);

    for ii = 1:nTrials

        t = trials(ii);

        % Metadata
        trialIDs(ii)   = string(t.trialID);
        subjectIDs(ii) = string(t.subjectID);
        shapes(ii)     = string(t.shape);
        databases(ii)  = string(t.database);
        fsVec(ii)      = t.fs;
        nSamps(ii)     = numel(t.x);

        %% Determine frequency range for fit
        fs = t.fs;
        if all(opts.FreqRange == 0)
            fLow  = fs / 10;   % well above typical movement harmonics
            fHigh = fs / 2.5;  % safely below Nyquist, avoids aliasing zone
        else
            fLow  = opts.FreqRange(1);
            fHigh = opts.FreqRange(2);
        end
        fLowVec(ii)  = fLow;
        fHighVec(ii) = fHigh;

        %% First-difference the raw positions
        deltaX = diff(t.x);
        deltaY = diff(t.y);

        % Remove mean (standard practice before spectral estimation)
        deltaX = deltaX - mean(deltaX);
        deltaY = deltaY - mean(deltaY);

        %% Multitaper PSD via pmtm
        % pmtm returns one-sided PSD in (units^2)/Hz
        [psdX, fVec] = pmtm(deltaX, opts.TBP, [], fs);
        [psdY, ~]    = pmtm(deltaY, opts.TBP, [], fs);

        %% Fit log-log slope in the noise-dominated band
        bandIdx = (fVec >= fLow) & (fVec <= fHigh);

        if sum(bandIdx) < 5
            if opts.Verbose
                fprintf("  Trial %d (%s): insufficient points in fit band, skipping\n", ...
                    ii, t.trialID);
            end
            continue
        end

        logF = log10(fVec(bandIdx));
        logPsdX = log10(psdX(bandIdx));
        logPsdY = log10(psdY(bandIdx));

        % OLS fit: log10(PSD) = slope * log10(f) + intercept
        [sX, intX, r2xVal] = fitLogLogSlope(logF, logPsdX);
        [sY, intY, r2yVal] = fitLogLogSlope(logF, logPsdY);

        slopeX(ii) = sX;
        slopeY(ii) = sY;
        r2X(ii)    = r2xVal;
        r2Y(ii)    = r2yVal;

        % Back-calculate alpha from first-differenced slope
        % PSD of diff(noise) ~ f^(2-alpha), so alpha = 2 - slope
        alphaX(ii) = 2 - sX;
        alphaY(ii) = 2 - sY;
        alphaMean(ii) = mean([alphaX(ii), alphaY(ii)]);

        %% Estimate sigma (noise magnitude in position data units)
        % Back-project the noise model to position PSD:
        %   S_pos(f) = C_pos * f^(-alpha)
        % where C_pos = C_diff * (fs / (2*pi))^2  (undoing the diff transfer function)
        %
        % Then sigma^2 = integral of S_pos(f) from fLow to Nyquist.
        % For alpha ~= 1: integral = C/(1-alpha) * [f^(1-alpha)] evaluated at limits
        % For alpha ~= 1: integral = C * ln(fHigh/fLow)
        fNyq = fs / 2;
        sigmaX(ii) = estimateSigma(sX, intX, fs, fLow, fNyq);
        sigmaY(ii) = estimateSigma(sY, intY, fs, fLow, fNyq);
        sigmaMean(ii) = mean([sigmaX(ii), sigmaY(ii)]);

        %% Diagnostic plot (if requested)
        if ismember(ii, opts.PlotTrials)
            plotDiagnostic(fVec, psdX, psdY, bandIdx, ...
                sX, intX, sY, intY, ...
                alphaX(ii), alphaY(ii), sigmaX(ii), sigmaY(ii), ...
                t, ii);
        end

        if opts.Verbose && mod(ii, 20) == 0
            fprintf("  Processed %d/%d trials\n", ii, nTrials);
        end

    end

    %% Assemble results table
    results = table( ...
        trialIDs, subjectIDs, shapes, databases, fsVec, nSamps, ...
        alphaX, alphaY, alphaMean, ...
        sigmaX, sigmaY, sigmaMean, ...
        slopeX, slopeY, r2X, r2Y, ...
        fLowVec, fHighVec, ...
        'VariableNames', { ...
        'trialID', 'subjectID', 'shape', 'database', 'fs', 'nSamples', ...
        'alphaX', 'alphaY', 'alphaMean', ...
        'sigmaX', 'sigmaY', 'sigmaMean', ...
        'slopeX', 'slopeY', 'r2X', 'r2Y', ...
        'fLow', 'fHigh'});

    %% Summary
    if opts.Verbose
        validIdx = ~isnan(alphaMean);
        nValid = sum(validIdx);
        fprintf("\n[characteriseNoise] %d/%d trials characterised\n", nValid, nTrials);

        if nValid > 0
            fprintf("  alpha: %.2f +/- %.2f  (range [%.2f, %.2f])\n", ...
                mean(alphaMean(validIdx)), std(alphaMean(validIdx)), ...
                min(alphaMean(validIdx)), max(alphaMean(validIdx)));
            fprintf("  sigma: %.3f +/- %.3f  (range [%.3f, %.3f])  [data units]\n", ...
                mean(sigmaMean(validIdx)), std(sigmaMean(validIdx)), ...
                min(sigmaMean(validIdx)), max(sigmaMean(validIdx)));
            fprintf("  R2 (fit quality): X=%.3f, Y=%.3f (mean)\n", ...
                mean(r2X(validIdx)), mean(r2Y(validIdx)));
        end
    end

end

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function [slope, intercept, r2] = fitLogLogSlope(logF, logPsd)
% FITLOGLOGSLOPE OLS fit of log10(PSD) vs log10(f)
    X = [ones(numel(logF), 1), logF(:)];
    b = X \ logPsd(:);
    intercept = b(1);
    slope     = b(2);

    fitted = X * b;
    ssRes  = sum((logPsd(:) - fitted).^2);
    ssTot  = sum((logPsd(:) - mean(logPsd)).^2);
    r2     = 1 - ssRes / ssTot;
end

function sigma = estimateSigma(slopeDiff, interceptDiff, fs, fLow, fNyq)
% ESTIMATESIGMA Back-project diff PSD to position PSD, integrate for sigma
%
% The fitted model on first-differenced data is:
%   log10(S_diff) = interceptDiff + slopeDiff * log10(f)
%   => S_diff(f) = 10^interceptDiff * f^slopeDiff
%
% The diff transfer function magnitude squared is:
%   |H(f)|^2 = (2*sin(pi*f/fs))^2  (exact for discrete diff)
%
% For f << fs:  |H|^2 ~ (2*pi*f/fs)^2
%
% So S_pos(f) ~ S_diff(f) / |H(f)|^2
%
% We integrate S_pos(f) from fLow to fNyq to get sigma^2.

    cDiff = 10^interceptDiff;

    % Numerical integration (robust for any alpha)
    nPts = 500;
    fGrid = linspace(fLow, fNyq, nPts);

    % S_diff at each frequency
    sDiff = cDiff * fGrid.^slopeDiff;

    % Exact discrete-diff transfer function
    hSq = (2 * sin(pi * fGrid / fs)).^2;
    hSq(hSq < eps) = eps; % guard against division by zero

    % Position PSD
    sPos = sDiff ./ hSq;

    % Integrate (trapezoidal)
    variance = trapz(fGrid, sPos);

    sigma = sqrt(max(variance, 0));
end

function plotDiagnostic(fVec, psdX, psdY, bandIdx, ...
        sX, intX, sY, intY, aX, aY, sigX, sigY, trial, trialNum)
% PLOTDIAGNOSTIC Generate a 2-panel diagnostic figure for one trial

    figure('Position', [100 100 900 450], 'Name', ...
        sprintf('Noise: %s trial %d', trial.trialID, trialNum));

    fBand = fVec(bandIdx);
    fitX = 10.^(intX + sX * log10(fBand));
    fitY = 10.^(intY + sY * log10(fBand));

    % X axis
    subplot(1,2,1);
    loglog(fVec, psdX, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5); hold on;
    loglog(fBand, psdX(bandIdx), 'b', 'LineWidth', 1);
    loglog(fBand, fitX, 'r--', 'LineWidth', 1.5);
    xline(fBand(1), ':k'); xline(fBand(end), ':k');
    xlabel('Frequency (Hz)');
    ylabel('PSD (units^2/Hz)');
    title(sprintf('\\Deltax:  \\alpha_x = %.2f,  \\sigma_x = %.3f', aX, sigX));
    legend('Full PSD', 'Fit band', sprintf('Slope = %.2f', sX), ...
        'Location', 'southwest');
    grid on;

    % Y axis
    subplot(1,2,2);
    loglog(fVec, psdY, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5); hold on;
    loglog(fBand, psdY(bandIdx), 'b', 'LineWidth', 1);
    loglog(fBand, fitY, 'r--', 'LineWidth', 1.5);
    xline(fBand(1), ':k'); xline(fBand(end), ':k');
    xlabel('Frequency (Hz)');
    ylabel('PSD (units^2/Hz)');
    title(sprintf('\\Deltay:  \\alpha_y = %.2f,  \\sigma_y = %.3f', aY, sigY));
    legend('Full PSD', 'Fit band', sprintf('Slope = %.2f', sY), ...
        'Location', 'southwest');
    grid on;

    sgtitle(sprintf('%s | %s | %s | fs=%d Hz | N=%d', ...
        trial.database, trial.subjectID, trial.shape, ...
        trial.fs, numel(trial.x)), 'Interpreter', 'none');
end
