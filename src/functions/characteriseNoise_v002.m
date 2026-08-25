function results = characteriseNoise_v002(trials, opts)
% characteriseNoise_v002 Estimate noise colour (alpha) and magnitude (sigma)
%
% Characterises the stochastic noise in raw position data from empirical
% movement recordings. Implements the noise characterisation protocol from
% prereg Section 7.3 using multitaper spectral estimation (pmtm) on
% first-differenced position data, with exact discrete-diff transfer
% function correction.
%
% v002 changes from v001:
%   - Undoes the exact discrete-diff transfer function 4*sin^2(pi*f/fs)
%     to recover the corrected position PSD before fitting alpha
%   - Sigma estimated via robust median(PSD)*bandwidth method, not
%     model-dependent integration
%   - R2 of spectral fit reported; alpha flagged as indeterminate when
%     R2 < R2Min (default 0.10)
%   - Fit band defaults to last octave [fs/4, fs/2.2] where the noise
%     floor is cleanest (v001 used [fs/10, fs/2.5] which was still
%     signal-dominated)
%
% Approach:
%   1. First-difference raw x,y to suppress movement-dominated low
%      frequencies (avoids circular dependence on any smoother choice).
%   2. Multitaper PSD (pmtm) on delta_x, delta_y.
%   3. Undo the exact discrete-diff transfer function to recover the
%      corrected position noise PSD.
%   4. Fit log10(PSD_pos) vs log10(f) in the noise-dominated band.
%      alpha = -slope (directly, no 2-minus correction needed).
%   5. Sigma = sqrt(median(PSD_pos) * bandwidth) — robust against
%      spectral outliers and insensitive to alpha uncertainty.
%
% Usage:
%   results = characteriseNoise_v002(trials)
%   results = characteriseNoise_v002(trials, TBP=3, FreqRange=[33 60])
%   results = characteriseNoise_v002(trials, PlotTrials=1:5)
%
% Inputs:
%   trials  - [M x 1] struct array from importDB_*.m with fields:
%             .x, .y, .fs, .subjectID, .trialID, .shape, .database
%
% Name-Value Arguments:
%   TBP        - Time-bandwidth product for pmtm (default: 3; prereg
%                specifies TBP in [2.5, 4], giving 5-7 tapers)
%   FreqRange  - [fLow fHigh] Hz for spectral slope fit (default: [0 0]
%                = auto: [fs/4, fs/2.2], the last octave)
%   R2Min      - Minimum R2 for alpha to be considered determinate
%                (default: 0.10). Below this, alphaReliable = false.
%   SpikeSD    - Speed spike threshold in MAD-scaled units (default: 10).
%                Samples with speed > median + SpikeSD*1.4826*MAD are
%                flagged as spikes. If spikes are found, the longest
%                clean segment is used for PSD estimation.
%   PlotTrials - Vector of trial indices to generate diagnostic plots
%                (default: [] = no plots).
%   Verbose    - Print progress and summary (default: true)
%
% Output:
%   results - table with one row per trial and columns:
%     trialID, subjectID, shape, database, fs, nSamples,
%     alphaX, alphaY, alphaMean,       (noise colour per axis and mean)
%     sigmaX, sigmaY, sigmaMean,       (noise SD in data units)
%     slopeX, slopeY,                  (fitted slopes on corrected PSD)
%     r2X, r2Y,                        (R-squared of spectral fit)
%     alphaReliable,                   (logical: both R2 >= R2Min)
%     spikeCount,                       (number of speed spikes detected)
%     fLow, fHigh                      (frequency range used for fit)
%
% Dependencies:
%   Signal Processing Toolbox (pmtm)
%
% See also: importDB_cook_v002, Stage5_processTrials
%
% Created 2026-03-17  v001
% Updated 2026-03-17  v002 — transfer-function correction, robust sigma,
%                            spike detection with longest-clean-segment
% Dagmar Scott Fraser - d.s.fraser@bham.ac.uk

    arguments
        trials      (1,:) struct
        opts.TBP        (1,1) double {mustBeInRange(opts.TBP, 2, 6)} = 3
        opts.FreqRange  (1,2) double = [0 0]
        opts.R2Min      (1,1) double = 0.10
        opts.PlotTrials (1,:) double = []
        opts.SpikeSD    (1,1) double = 10
        opts.Verbose    (1,1) logical = true
    end

    nTrials = numel(trials);

    % Pre-allocate output columns
    trialIDs      = strings(nTrials, 1);
    subjectIDs    = strings(nTrials, 1);
    shapes        = strings(nTrials, 1);
    databases     = strings(nTrials, 1);
    fsVec         = NaN(nTrials, 1);
    nSamps        = NaN(nTrials, 1);
    alphaX        = NaN(nTrials, 1);
    alphaY        = NaN(nTrials, 1);
    alphaMean     = NaN(nTrials, 1);
    sigmaX        = NaN(nTrials, 1);
    sigmaY        = NaN(nTrials, 1);
    sigmaMean     = NaN(nTrials, 1);
    slopeXout     = NaN(nTrials, 1);
    slopeYout     = NaN(nTrials, 1);
    r2Xout        = NaN(nTrials, 1);
    r2Yout        = NaN(nTrials, 1);
    alphaReliable = false(nTrials, 1);
    spikeCount    = zeros(nTrials, 1);
    fLowVec       = NaN(nTrials, 1);
    fHighVec      = NaN(nTrials, 1);

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
            fLow  = fs / 4;    % last octave before Nyquist
            fHigh = fs / 2.2;  % safely below Nyquist
        else
            fLow  = opts.FreqRange(1);
            fHigh = opts.FreqRange(2);
        end
        fLowVec(ii)  = fLow;
        fHighVec(ii) = fHigh;

        %% First-difference the raw positions
        deltaX = diff(t.x);
        deltaY = diff(t.y);

        %% Spike detection: use longest clean segment if glitches found
        speed = sqrt(deltaX.^2 + deltaY.^2);
        medSpd = median(speed);
        madSpd = median(abs(speed - medSpd));
        spikeIdx = speed > medSpd + opts.SpikeSD * 1.4826 * madSpd;
        spikeCount(ii) = sum(spikeIdx);

        if spikeCount(ii) > 0
            % Find longest run of non-spike samples
            cleanRuns = diff([0; ~spikeIdx(:); 0]);
            runStarts = find(cleanRuns == 1);
            runEnds   = find(cleanRuns == -1) - 1;
            runLens   = runEnds - runStarts + 1;
            [bestLen, bestRun] = max(runLens);

            if bestLen < 100
                if opts.Verbose
                    fprintf("  Trial %d (%s): %d spikes, longest clean segment = %d, skipping\n", ...
                        ii, t.trialID, spikeCount(ii), bestLen);
                end
                continue
            end

            if opts.Verbose
                fprintf("  Trial %d (%s): %d spikes detected, using clean segment (%d/%d samples)\n", ...
                    ii, t.trialID, spikeCount(ii), bestLen, numel(deltaX));
            end

            seg = runStarts(bestRun):runEnds(bestRun);
            deltaX = deltaX(seg);
            deltaY = deltaY(seg);
        end

        % Remove mean (standard practice before spectral estimation)
        deltaX = deltaX - mean(deltaX);
        deltaY = deltaY - mean(deltaY);

        %% Multitaper PSD via pmtm
        [psdDiffX, fVec] = pmtm(deltaX, opts.TBP, [], fs);
        [psdDiffY, ~]    = pmtm(deltaY, opts.TBP, [], fs);

        %% Undo exact discrete-diff transfer function
        % diff transfer: |H(f)|^2 = 4 sin^2(pi f / fs)
        hSq = 4 * sin(pi * fVec / fs).^2;
        hSq(hSq < eps) = eps;
        psdPosX = psdDiffX ./ hSq;
        psdPosY = psdDiffY ./ hSq;

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
        logPosX = log10(psdPosX(bandIdx));
        logPosY = log10(psdPosY(bandIdx));

        % OLS fit: log10(PSD_pos) = slope * log10(f) + intercept
        [sX, ~, r2xVal] = fitLogLogSlope(logF, logPosX);
        [sY, ~, r2yVal] = fitLogLogSlope(logF, logPosY);

        slopeXout(ii) = sX;
        slopeYout(ii) = sY;
        r2Xout(ii)    = r2xVal;
        r2Yout(ii)    = r2yVal;

        % Alpha directly from corrected position PSD: S_pos ~ f^(-alpha)
        alphaX(ii) = -sX;
        alphaY(ii) = -sY;
        alphaMean(ii) = mean([alphaX(ii), alphaY(ii)]);

        % Flag reliability
        alphaReliable(ii) = (r2xVal >= opts.R2Min) && (r2yVal >= opts.R2Min);

        %% Estimate sigma via robust median method
        % sigma^2 = median(PSD_pos) * bandwidth
        % Median is robust to spectral peaks and dips; bandwidth = fHigh - fLow
        bw = fHigh - fLow;
        sigmaX(ii) = sqrt(median(psdPosX(bandIdx)) * bw);
        sigmaY(ii) = sqrt(median(psdPosY(bandIdx)) * bw);
        sigmaMean(ii) = mean([sigmaX(ii), sigmaY(ii)]);

        %% Diagnostic plot (if requested)
        if ismember(ii, opts.PlotTrials)
            plotDiagnostic(fVec, psdPosX, psdPosY, psdDiffX, psdDiffY, ...
                bandIdx, sX, sY, ...
                alphaX(ii), alphaY(ii), sigmaX(ii), sigmaY(ii), ...
                r2xVal, r2yVal, alphaReliable(ii), t, ii);
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
        slopeXout, slopeYout, r2Xout, r2Yout, ...
        alphaReliable, spikeCount, ...
        fLowVec, fHighVec, ...
        'VariableNames', { ...
        'trialID', 'subjectID', 'shape', 'database', 'fs', 'nSamples', ...
        'alphaX', 'alphaY', 'alphaMean', ...
        'sigmaX', 'sigmaY', 'sigmaMean', ...
        'slopeX', 'slopeY', 'r2X', 'r2Y', ...
        'alphaReliable', 'spikeCount', ...
        'fLow', 'fHigh'});

    %% Summary
    if opts.Verbose
        validIdx = ~isnan(alphaMean);
        nValid = sum(validIdx);
        nReliableAlpha = sum(alphaReliable);
        fprintf("\n[characteriseNoise] %d/%d trials characterised\n", nValid, nTrials);
        fprintf("  Alpha reliable (R2 >= %.2f): %d/%d\n", ...
            opts.R2Min, nReliableAlpha, nValid);

        if nValid > 0
            fprintf("  alpha: %.2f +/- %.2f  (range [%.2f, %.2f])", ...
                mean(alphaMean(validIdx)), std(alphaMean(validIdx)), ...
                min(alphaMean(validIdx)), max(alphaMean(validIdx)));
            if nReliableAlpha == 0
                fprintf("  ** all indeterminate **");
            end
            fprintf("\n");

            fprintf("  sigma: %.4f +/- %.4f  (range [%.4f, %.4f])  [data units]\n", ...
                mean(sigmaMean(validIdx)), std(sigmaMean(validIdx)), ...
                min(sigmaMean(validIdx)), max(sigmaMean(validIdx)));
            fprintf("  R2 (fit quality): X=%.3f, Y=%.3f (median)\n", ...
                median(r2Xout(validIdx)), median(r2Yout(validIdx)));
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

function plotDiagnostic(fVec, psdPosX, psdPosY, psdDiffX, psdDiffY, ...
        bandIdx, sX, sY, aX, aY, sigX, sigY, r2x, r2y, reliable, trial, trialNum)
% PLOTDIAGNOSTIC 3-panel diagnostic for one trial

    figure('Position', [100 100 1100 700], 'Name', ...
        sprintf('Noise v002: %s trial %d', trial.trialID, trialNum));

    fBand = fVec(bandIdx);
    fs = trial.fs;

    %% Panel 1: full corrected position PSD
    subplot(2,2,1);
    loglog(fVec, psdPosX, 'b', 'LineWidth', 0.5); hold on;
    loglog(fVec, psdPosY, 'r', 'LineWidth', 0.5);
    xline(fBand(1), '--g', 'LineWidth', 1);
    xline(fBand(end), '--g', 'LineWidth', 1);
    % Reference: white noise level
    medX = median(psdPosX(bandIdx));
    yline(medX, ':k', 'LineWidth', 1);
    xlabel('Frequency (Hz)');
    ylabel('PSD (position, corrected)');
    title('Corrected position PSD');
    legend('x', 'y', 'fit band', '', 'median floor', 'Location', 'southwest');
    grid on;

    %% Panel 2: first-differenced PSD (raw, before correction)
    subplot(2,2,2);
    loglog(fVec, psdDiffX, 'b', 'LineWidth', 0.5); hold on;
    loglog(fVec, psdDiffY, 'r', 'LineWidth', 0.5);
    xline(fBand(1), '--g', 'LineWidth', 1);
    xline(fBand(end), '--g', 'LineWidth', 1);
    xlabel('Frequency (Hz)');
    ylabel('PSD (first-differenced)');
    title('First-differenced PSD (uncorrected)');
    grid on;

    %% Panel 3: zoomed fit band with slope overlay
    subplot(2,2,3);
    loglog(fBand, psdPosX(bandIdx), 'b', 'LineWidth', 1); hold on;
    loglog(fBand, psdPosY(bandIdx), 'r', 'LineWidth', 1);
    % Fit lines
    logF = log10(fBand);
    X = [ones(numel(logF), 1), logF(:)];
    bX = X \ log10(psdPosX(bandIdx));
    bY = X \ log10(psdPosY(bandIdx));
    loglog(fBand, 10.^(X*bX), 'b--', 'LineWidth', 2);
    loglog(fBand, 10.^(X*bY), 'r--', 'LineWidth', 2);
    xlabel('Frequency (Hz)');
    ylabel('PSD (zoomed)');
    reliStr = "RELIABLE";
    if ~reliable, reliStr = "INDETERMINATE"; end
    title(sprintf('\\alpha_x=%.2f  \\alpha_y=%.2f  R^2=%.2f/%.2f  [%s]', ...
        aX, aY, r2x, r2y, reliStr));
    grid on;

    %% Panel 4: text summary
    subplot(2,2,4);
    axis off;
    summaryStr = { ...
        sprintf('Trial: %s', trial.trialID), ...
        sprintf('Database: %s', trial.database), ...
        sprintf('Shape: %s', trial.shape), ...
        sprintf('fs = %d Hz,  N = %d samples', fs, numel(trial.x)), ...
        '', ...
        sprintf('Fit band: [%.1f, %.1f] Hz', fBand(1), fBand(end)), ...
        sprintf('\\alpha_x = %.2f   \\alpha_y = %.2f   mean = %.2f', aX, aY, mean([aX aY])), ...
        sprintf('\\sigma_x = %.4f   \\sigma_y = %.4f   mean = %.4f', sigX, sigY, mean([sigX sigY])), ...
        sprintf('R^2_x = %.3f   R^2_y = %.3f', r2x, r2y), ...
        sprintf('Alpha status: %s (R^2 threshold = 0.10)', reliStr), ...
        '', ...
        'Method: pmtm on diff(pos), exact |H(f)|^2 correction,', ...
        'slope fit on corrected position PSD, sigma via median*BW'};
    text(0.05, 0.95, summaryStr, 'VerticalAlignment', 'top', ...
        'FontName', 'FixedWidth', 'FontSize', 9);

    sgtitle(sprintf('%s | %s | fs=%d Hz', ...
        trial.database, trial.trialID, fs), 'Interpreter', 'none');
end
