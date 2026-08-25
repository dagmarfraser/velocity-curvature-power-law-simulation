function results = characteriseBiologicalNoise_v001(trials, opts)
% characteriseBiologicalNoise_v001 Template-subtraction noise characterisation.
%
% For periodic movement data (ellipse tracing etc.), estimates the biological
% noise component by fitting and subtracting a sliding-window sinusoidal
% template, then characterising the residuals via pmtm.
%
% Generalised from diagnoseSinusoidalResiduals_cook_v3.m. Works with any
% trial struct array conforming to the Stage 5 interface (.x, .y, .fs,
% .subjectID, .trialID, .database).
%
% Returns a table with per-trial noise estimates: alpha (position space),
% sigma, f0, variance explained, and spectral fit R2.
%
% USAGE:
%   R = characteriseBiologicalNoise_v001(trials)
%   R = characteriseBiologicalNoise_v001(trials, NHarmonics=6, FitBand=[1 40])
%
% PREREG REF: Section 7.3 (noise characterisation via pmtm)
% Fraser, D.S. (2025)

    arguments
        trials (:,1) struct
        opts.NHarmonics  (1,1) double = 4        % harmonics in template fit
        opts.NCyclesWin  (1,1) double = 4        % cycles per sliding window
        opts.TBP         (1,1) double = 3        % time-bandwidth product for pmtm
        opts.FitBandLow  (1,1) double = 1        % Hz, low edge of spectral fit band
        opts.FitBandHigh (1,1) double = 0        % Hz, 0 = auto (fs/2.2)
        opts.Verbose     (1,1) logical = true
    end

    nTrials = numel(trials);

    allF0       = NaN(nTrials, 1);
    allVarExpX  = NaN(nTrials, 1);
    allVarExpY  = NaN(nTrials, 1);
    allSigmaX   = NaN(nTrials, 1);
    allSigmaY   = NaN(nTrials, 1);
    allAlphaX   = NaN(nTrials, 1);
    allAlphaY   = NaN(nTrials, 1);
    allR2X      = NaN(nTrials, 1);
    allR2Y      = NaN(nTrials, 1);
    allSubj     = strings(nTrials, 1);
    allTrial    = strings(nTrials, 1);

    for ii = 1:nTrials
        t = trials(ii);
        fs = t.fs;
        N  = numel(t.x);
        allSubj(ii)  = string(t.subjectID);
        allTrial(ii) = string(t.trialID);

        fitHigh = opts.FitBandHigh;
        if fitHigh == 0, fitHigh = fs / 2.2; end

        %% Estimate f0 via FFT peak, cross-validated between x and y
        nfft = 2^nextpow2(4*N);
        xDet = detrend(double(t.x), "linear");
        yDet = detrend(double(t.y), "linear");
        Xf = abs(fft(xDet, nfft));
        Yf = abs(fft(yDet, nfft));
        fAxis = (0:nfft-1)' * fs / nfft;
        searchBand = fAxis > 0.1 & fAxis < 5;

        [~, pkIdxX] = max(Xf(searchBand));
        [~, pkIdxY] = max(Yf(searchBand));
        searchIdx = find(searchBand);
        f0x = fAxis(searchIdx(pkIdxX));
        f0y = fAxis(searchIdx(pkIdxY));

        if abs(f0x - f0y)/f0x < 0.2
            f0 = mean([f0x, f0y]);
        elseif max(Xf(searchBand)) > max(Yf(searchBand))
            f0 = f0x;
        else
            f0 = f0y;
        end
        allF0(ii) = f0;

        %% Sliding window sinusoidal subtraction
        winSamples = round(opts.NCyclesWin / f0 * fs);
        hopSamples = round(winSamples / 2);
        winSamples = min(winSamples, N);
        winSamples = max(winSamples, round(2/f0*fs));

        resX = zeros(N, 1);
        resY = zeros(N, 1);
        winWeight = zeros(N, 1);

        startIdx = 1;
        while startIdx + winSamples - 1 <= N
            endIdx = startIdx + winSamples - 1;
            tWin = (0:winSamples-1)' / fs;
            hannWin = 0.5 * (1 - cos(2*pi*(0:winSamples-1)'/(winSamples-1)));

            % Build design matrix: DC + drift + harmonics
            D = ones(winSamples, 1);
            D(:,2) = tWin;
            for h = 1:opts.NHarmonics
                D(:, end+1) = cos(2*pi*h*f0*tWin); %#ok<AGROW>
                D(:, end+1) = sin(2*pi*h*f0*tWin); %#ok<AGROW>
            end

            xWin = double(t.x(startIdx:endIdx));
            yWin = double(t.y(startIdx:endIdx));
            Dw = D .* hannWin;
            xw  = xWin .* hannWin;
            yw  = yWin .* hannWin;

            bX = Dw \ xw;
            bY = Dw \ yw;

            resX(startIdx:endIdx) = resX(startIdx:endIdx) + (xWin - D*bX) .* hannWin;
            resY(startIdx:endIdx) = resY(startIdx:endIdx) + (yWin - D*bY) .* hannWin;
            winWeight(startIdx:endIdx) = winWeight(startIdx:endIdx) + hannWin;

            startIdx = startIdx + hopSamples;
        end

        % Normalise overlap-add
        validMask = winWeight > 0;
        resX(validMask) = resX(validMask) ./ winWeight(validMask);
        resY(validMask) = resY(validMask) ./ winWeight(validMask);

        % Clip edges
        edgeClip = max(round(winSamples/4), 1);
        cStart = edgeClip;
        cEnd   = max(N - edgeClip, cStart + 100);
        resX = resX(cStart:cEnd);
        resY = resY(cStart:cEnd);
        xClip = double(t.x(cStart:cEnd));
        yClip = double(t.y(cStart:cEnd));

        allVarExpX(ii) = 1 - var(resX)/var(xClip);
        allVarExpY(ii) = 1 - var(resY)/var(yClip);
        allSigmaX(ii)  = std(resX);
        allSigmaY(ii)  = std(resY);

        %% PSD of residuals + spectral slope
        [psdResX, fVec] = pmtm(resX, opts.TBP, [], fs);
        [psdResY, ~]    = pmtm(resY, opts.TBP, [], fs);

        fitBand = fVec >= opts.FitBandLow & fVec <= fitHigh;
        if sum(fitBand) < 3
            continue;  % skip trial if fit band too narrow
        end
        logF  = log10(fVec(fitBand));
        logPX = log10(psdResX(fitBand));
        logPY = log10(psdResY(fitBand));
        Xfit  = [ones(sum(fitBand), 1), logF(:)];
        bFitX = Xfit \ logPX(:);
        bFitY = Xfit \ logPY(:);
        allAlphaX(ii) = -bFitX(2);
        allAlphaY(ii) = -bFitY(2);
        allR2X(ii) = 1 - sum((logPX(:) - Xfit*bFitX).^2) / sum((logPX(:) - mean(logPX)).^2);
        allR2Y(ii) = 1 - sum((logPY(:) - Xfit*bFitY).^2) / sum((logPY(:) - mean(logPY)).^2);

        if opts.Verbose && mod(ii, 20) == 0
            fprintf("  Processed %d/%d\n", ii, nTrials);
        end
    end

    %% Package into table
    results = table(allTrial, allSubj, allF0, ...
        allAlphaX, allAlphaY, mean([allAlphaX, allAlphaY], 2), ...
        allSigmaX, allSigmaY, mean([allSigmaX, allSigmaY], 2), ...
        allR2X, allR2Y, mean([allR2X, allR2Y], 2), ...
        allVarExpX, allVarExpY, ...
        'VariableNames', ["trialID", "subjectID", "f0", ...
        "alphaX", "alphaY", "alphaMean", ...
        "sigmaX", "sigmaY", "sigmaMean", ...
        "r2X", "r2Y", "r2Mean", ...
        "varExpX", "varExpY"]);

    %% Summary
    if opts.Verbose
        db = string(trials(1).database);
        fprintf("\n=== BIOLOGICAL NOISE: %s (%d trials) ===\n", db, nTrials);
        fprintf("f0:     %.3f +/- %.3f Hz\n", mean(allF0, "omitnan"), std(allF0, "omitnan"));
        fprintf("VarExp: X=%.4f Y=%.4f\n", mean(allVarExpX, "omitnan"), mean(allVarExpY, "omitnan"));

        alphaMean = mean([allAlphaX, allAlphaY], 2);
        sigmaMean = mean([allSigmaX, allSigmaY], 2);
        fprintf("\nalpha (position): %.2f +/- %.2f  (after diff: %.2f)\n", ...
            mean(alphaMean, "omitnan"), std(alphaMean, "omitnan"), mean(alphaMean, "omitnan") - 2);
        fprintf("sigma (data units): %.3f +/- %.3f\n", ...
            mean(sigmaMean, "omitnan"), std(sigmaMean, "omitnan"));
        fprintf("R2 spectral fit: %.3f +/- %.3f\n", ...
            mean(results.r2Mean, "omitnan"), std(results.r2Mean, "omitnan"));

        fprintf("\n=== PER-SUBJECT ===\n");
        fprintf("%-8s %6s %6s %8s %8s %6s\n", "Subject", "f0", "alpha", "sigmaX", "sigmaY", "R2");
        uSubj = unique(allSubj);
        for ss = 1:numel(uSubj)
            idx = allSubj == uSubj(ss);
            fprintf("%-8s %6.3f %6.2f %8.3f %8.3f %6.3f\n", ...
                uSubj(ss), mean(allF0(idx), "omitnan"), ...
                mean(alphaMean(idx), "omitnan"), ...
                mean(allSigmaX(idx), "omitnan"), ...
                mean(allSigmaY(idx), "omitnan"), ...
                mean(results.r2Mean(idx), "omitnan"));
        end
    end

end
