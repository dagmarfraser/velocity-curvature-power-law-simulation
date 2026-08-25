function results = characteriseBiologicalNoise_v004(trials, opts)
% characteriseBiologicalNoise_v004  Biological noise characterisation.
%
% Iterates from v003. Two targeted changes to the IRASA method:
%
%   (1) IRASA now operates on TEMPLATE-SUBTRACTED RESIDUALS, not raw data.
%       v003 ran IRASA on raw x/y, which underestimated alpha by 1-2 units
%       because dense ellipse harmonics (f0=0.836 Hz, 23 harmonics in
%       [1,20 Hz]) defeat IRASA resampling cancellation: shifted harmonics
%       land on adjacent originals and do not average out in the geometric
%       mean. Template subtraction removes the harmonics before IRASA runs.
%       (Finding #97)
%
%   (2) IRASA fitting band upper bound: min(20, fs/2-1) replacing fs/2.2.
%       This matches the in-runner IRASA in runLoopClosureFftnoise_v007
%       exactly, making stored alpha and in-runner alpha directly comparable.
%       For Pilot (240 Hz) the old bound was 109 Hz, spanning the electronics
%       floor and dragging alpha down further. The fixed 20 Hz cap targets
%       the biological motor noise regime below the floor transition.
%
%   (3) TBP default: 4 (was 3), matching iraAlphaSigma_v001 in v007.
%
% These three changes together explain the 1-2 unit gap between v003 stored
% alpha values and the in-runner alpha used for noise generation. After v004
% the stored alpha and in-runner alpha are computed identically.
%
% IMPACT ON RESULTS:
%   Core results (beta_gen*, HKSJ, CCC, R7): UNAFFECTED.
%   runLoopClosureFftnoise_v007 independently recomputes per-trial alpha
%   from residuals inside processOneTrial_local and uses that for all noise
%   generation. Stored alpha in noiseCharacterisation mats is a reference
%   label only (used for per_subject trial selection; not used in all mode).
%
% METHODS (unchanged from v003):
%   'template' : Sliding-window sinusoidal template subtraction.
%                Best sigma estimate. Alpha underestimates above ~3.
%   'detrend'  : Staged polynomial detrending (prereg Section 7.3).
%                Better alpha than template up to ~4; collapses above 4.5.
%   'irasa'    : IRASA (Wen & Liu 2016) via iraAlphaSigma_v001.
%                Operates on RESIDUALS (v004 change). Best alpha estimate.
%                Sigma output is UNRELIABLE : do not use ira_sigma.
%                Use template sigma instead.
%
% METHOD OPTIONS:
%   'template'  : template only
%   'detrend'   : detrend only
%   'both'      : template + detrend  [backward-compatible default]
%   'irasa'     : template + IRASA    (template for sigma, IRASA for alpha)
%   'all'       : all three methods
%
% USAGE:
%   R = characteriseBiologicalNoise_v004(trials)
%   R = characteriseBiologicalNoise_v004(trials, Method="irasa")
%   R = characteriseBiologicalNoise_v004(trials, Method="all", NHarmonics=6)
%
% OUTPUT TABLE COLUMNS (identical to v003):
%   Common:      trialID, subjectID, f0
%   Template:    alphaX, alphaY, alphaMean, sigmaX, sigmaY, sigmaMean,
%                r2X, r2Y, r2Mean, varExpX, varExpY
%   Detrend:     det_alphaX, det_alphaY, det_alphaMean,
%                det_sigmaX, det_sigmaY, det_sigmaMean,
%                det_r2X, det_r2Y, det_r2Mean, det_polyOrd
%   IRASA:       ira_alphaX, ira_alphaY, ira_alphaMean
%                (no ira_sigma -- unreliable, see note above)
%
% FINDINGS: #97 (raw-data artefact in v003), #98 (corrected calibration slope)
% Fraser, D.S. (2026)  v004
% See also: iraAlphaSigma_v001, characteriseBiologicalNoise_v003

    arguments
        trials           (:,1) struct
        opts.Method      (1,1) string {mustBeMember(opts.Method, ...
            ["template","detrend","both","irasa","all"])} = "both"
        opts.NHarmonics  (1,1) double = 4
        opts.NCyclesWin  (1,1) double = 4
        opts.TBP         (1,1) double = 4        % v004: 4 (was 3), matches v007
        opts.FitBandLow  (1,1) double = 1.0      % Hz
        opts.FitBandHigh (1,1) double = 0        % 0 = auto: min(20, fs/2-1)
        opts.MaxPolyOrder (1,1) double = 3
        opts.IrasaHset   (1,:) double = 1.1:0.05:1.9
        opts.Verbose     (1,1) logical = true
    end

    doTemplate = ismember(opts.Method, ["template","both","irasa","all"]);
    doDetrend  = ismember(opts.Method, ["detrend", "both","all"]);
    doIRASA    = ismember(opts.Method, ["irasa",   "all"]);

    if doIRASA && ~doTemplate
        error('characteriseBiologicalNoise_v004:IRASARequiresTemplate', ...
            '%s', ['IRASA requires template to run concurrently (residuals needed). ' ...
            'Use Method="irasa" or Method="all".']);
    end

    nTrials = numel(trials);

    allF0    = NaN(nTrials, 1);
    allSubj  = strings(nTrials, 1);
    allTrial = strings(nTrials, 1);

    tpl = struct( ...
        'alphaX', NaN(nTrials,1), 'alphaY', NaN(nTrials,1), ...
        'sigmaX', NaN(nTrials,1), 'sigmaY', NaN(nTrials,1), ...
        'r2X',    NaN(nTrials,1), 'r2Y',    NaN(nTrials,1), ...
        'varExpX',NaN(nTrials,1), 'varExpY',NaN(nTrials,1));

    det = struct( ...
        'alphaX', NaN(nTrials,1), 'alphaY', NaN(nTrials,1), ...
        'sigmaX', NaN(nTrials,1), 'sigmaY', NaN(nTrials,1), ...
        'r2X',    NaN(nTrials,1), 'r2Y',    NaN(nTrials,1), ...
        'polyOrd',NaN(nTrials,1));

    ira = struct( ...
        'alphaX', NaN(nTrials,1), ...
        'alphaY', NaN(nTrials,1));

    for ii = 1:nTrials
        t  = trials(ii);
        fs = t.fs;
        allSubj(ii)  = string(t.subjectID);
        allTrial(ii) = string(t.trialID);

        % v004: fitting band upper bound matches runLoopClosureFftnoise_v007
        fitHigh = opts.FitBandHigh;
        if fitHigh == 0
            fitHigh = min(20.0, fs/2 - 1);   % v004: was fs/2.2 in v003
        end

        %% ---- Estimate f0 ----
        f0 = estimateF0_local(t, fs);
        allF0(ii) = f0;

        %% ---- Template subtraction ----
        % Pre-initialise so IRASA block can always reference resX/resY.
        resX = [];  resY = [];

        if doTemplate
            minWinSamples = round(2 / max(f0, 0.1) * fs);
            if numel(t.x) < minWinSamples
                if opts.Verbose
                    fprintf('  Trial %d (%s): too short for template (%d < %d samples)\n', ...
                        ii, allTrial(ii), numel(t.x), minWinSamples);
                end
            else
                [resX, resY, vexpX, vexpY] = templateSubtract_local( ...
                    t, fs, f0, opts.NHarmonics, opts.NCyclesWin);
                tpl.varExpX(ii) = vexpX;
                tpl.varExpY(ii) = vexpY;
                tpl.sigmaX(ii)  = std(resX);
                tpl.sigmaY(ii)  = std(resY);
                [tpl.alphaX(ii), tpl.r2X(ii)] = fitSpectralSlope_local( ...
                    resX, fs, opts.TBP, opts.FitBandLow, fitHigh);
                [tpl.alphaY(ii), tpl.r2Y(ii)] = fitSpectralSlope_local( ...
                    resY, fs, opts.TBP, opts.FitBandLow, fitHigh);
            end
        end

        %% ---- Staged detrending ----
        if doDetrend
            [dResX, dResY, polyOrd] = stagedDetrend_local( ...
                t, fs, f0, opts.MaxPolyOrder);
            det.polyOrd(ii) = polyOrd;
            det.sigmaX(ii)  = std(dResX);
            det.sigmaY(ii)  = std(dResY);
            [det.alphaX(ii), det.r2X(ii)] = fitSpectralSlope_local( ...
                dResX, fs, opts.TBP, opts.FitBandLow, fitHigh);
            [det.alphaY(ii), det.r2Y(ii)] = fitSpectralSlope_local( ...
                dResY, fs, opts.TBP, opts.FitBandLow, fitHigh);
        end

        %% ---- IRASA ----
        % v004: operates on template-subtracted RESIDUALS (not raw signal).
        % This eliminates harmonic contamination that caused v003 to
        % underestimate alpha by 1-2 units. Fitting band now min(20,fs/2-1)
        % matching runLoopClosureFftnoise_v007. (Findings #97, #98)
        if doIRASA
            if ~isempty(resX) && numel(resX) >= 64
                [ira.alphaX(ii), ~] = iraAlphaSigma_v001( ...
                    resX(:), fs, opts.FitBandLow, fitHigh, opts.IrasaHset);
                [ira.alphaY(ii), ~] = iraAlphaSigma_v001( ...
                    resY(:), fs, opts.FitBandLow, fitHigh, opts.IrasaHset);
            end
        end

        if opts.Verbose && mod(ii, 20) == 0
            fprintf('  Processed %d/%d\n', ii, nTrials);
        end
    end

    %% ---- Build output table ----
    results = buildResultsTable_local(allTrial, allSubj, allF0, ...
        tpl, det, ira, doTemplate, doDetrend, doIRASA);

    %% ---- Summary ----
    if opts.Verbose
        printSummary_local(results, trials, doTemplate, doDetrend, doIRASA);
    end

end

%% ========== LOCAL FUNCTIONS ==========

function f0 = estimateF0_local(t, fs)
    N    = numel(t.x);
    nfft = 2^nextpow2(4*N);
    xDet = detrend(double(t.x), 'linear');
    yDet = detrend(double(t.y), 'linear');
    Xf   = abs(fft(xDet, nfft));
    Yf   = abs(fft(yDet, nfft));
    fAxis = (0:nfft-1)' * fs / nfft;
    band  = fAxis > 0.1 & fAxis < 5;
    [~, pkX] = max(Xf(band));
    [~, pkY] = max(Yf(band));
    idx = find(band);
    f0x = fAxis(idx(pkX));
    f0y = fAxis(idx(pkY));
    if abs(f0x - f0y) / max(f0x, 0.01) < 0.2
        f0 = mean([f0x, f0y]);
    elseif max(Xf(band)) > max(Yf(band))
        f0 = f0x;
    else
        f0 = f0y;
    end
end

function [resX, resY, varExpX, varExpY] = templateSubtract_local( ...
        t, fs, f0, nHarmonics, nCyclesWin)
    N = numel(t.x);
    winSamples = round(nCyclesWin / f0 * fs);
    hopSamples = round(winSamples / 2);
    winSamples = min(winSamples, N);
    winSamples = max(winSamples, round(2/f0*fs));

    resX      = zeros(N, 1);
    resY      = zeros(N, 1);
    winWeight = zeros(N, 1);

    startIdx = 1;
    while startIdx + winSamples - 1 <= N
        endIdx  = startIdx + winSamples - 1;
        tWin    = (0:winSamples-1)' / fs;
        hannWin = 0.5 * (1 - cos(2*pi*(0:winSamples-1)'/(winSamples-1)));

        D = ones(winSamples, 1);
        D(:,2) = tWin;
        for h = 1:nHarmonics
            D(:, end+1) = cos(2*pi*h*f0*tWin); %#ok<AGROW>
            D(:, end+1) = sin(2*pi*h*f0*tWin); %#ok<AGROW>
        end
        xWin = double(t.x(startIdx:endIdx));
        yWin = double(t.y(startIdx:endIdx));
        Dw   = D .* hannWin;
        bX   = Dw \ (xWin .* hannWin);
        bY   = Dw \ (yWin .* hannWin);

        resX(startIdx:endIdx) = resX(startIdx:endIdx) + (xWin - D*bX) .* hannWin;
        resY(startIdx:endIdx) = resY(startIdx:endIdx) + (yWin - D*bY) .* hannWin;
        winWeight(startIdx:endIdx) = winWeight(startIdx:endIdx) + hannWin;
        startIdx = startIdx + hopSamples;
    end

    valid = winWeight > 0;
    resX(valid) = resX(valid) ./ winWeight(valid);
    resY(valid) = resY(valid) ./ winWeight(valid);

    edgeClip = max(round(winSamples/4), 1);
    cStart   = edgeClip;
    cEnd     = min(max(N - edgeClip, cStart + 100), N);
    xClip    = double(t.x(cStart:cEnd));
    yClip    = double(t.y(cStart:cEnd));
    resX     = resX(cStart:cEnd);
    resY     = resY(cStart:cEnd);

    varExpX = 1 - var(resX) / var(xClip);
    varExpY = 1 - var(resY) / var(yClip);
end

function [resX, resY, polyOrd] = stagedDetrend_local(t, fs, f0, maxPolyOrd)
    xRaw = double(t.x(:));
    yRaw = double(t.y(:));
    N    = numel(xRaw);

    resX    = detrend(xRaw, 'linear');
    resY    = detrend(yRaw, 'linear');
    polyOrd = 1;

    periodSamples = round(fs / max(f0, 0.01));
    if periodSamples > 1 && periodSamples < N/2
        acfX = xcorr(resX, periodSamples, 'coeff');
        acfY = xcorr(resY, periodSamples, 'coeff');
        periodicityDetected = max(abs(acfX(end)), abs(acfY(end))) > 0.3;
    else
        periodicityDetected = false;
    end

    if periodicityDetected && maxPolyOrd > 1
        polyOrd = maxPolyOrd;
        tVec    = (0:N-1)' / fs;
        V       = zeros(N, polyOrd + 1);
        for k = 0:polyOrd
            V(:, k+1) = tVec.^k;
        end
        resX = xRaw - V * (V \ xRaw);
        resY = yRaw - V * (V \ yRaw);
    end

    edgeClip = min(20, floor(N/10));
    if edgeClip > 0
        resX = resX(edgeClip+1:end-edgeClip);
        resY = resY(edgeClip+1:end-edgeClip);
    end
end

function [alpha, r2] = fitSpectralSlope_local(residuals, fs, tbp, fLow, fHigh)
    alpha = NaN;  r2 = NaN;
    if numel(residuals) < 64, return; end
    [psd, fVec] = pmtm(residuals, tbp, [], fs);
    band = fVec >= fLow & fVec <= fHigh;
    if sum(band) < 3, return; end
    logF   = log10(fVec(band));
    logPSD = log10(psd(band));
    X      = [ones(sum(band), 1), logF(:)];
    b      = X \ logPSD(:);
    alpha  = -b(2);
    ssRes  = sum((logPSD(:) - X*b).^2);
    ssTot  = sum((logPSD(:) - mean(logPSD)).^2);
    r2     = 1 - ssRes / max(ssTot, eps);
end

function results = buildResultsTable_local( ...
        trialIDs, subjIDs, f0, tpl, det, ira, doTpl, doDet, doIRA)
    results = table(trialIDs, subjIDs, f0, ...
        'VariableNames', ["trialID","subjectID","f0"]);

    if doTpl
        alphaMean = mean([tpl.alphaX, tpl.alphaY], 2);
        sigmaMean = mean([tpl.sigmaX, tpl.sigmaY], 2);
        r2Mean    = mean([tpl.r2X,    tpl.r2Y],    2);
        tplT = table( ...
            tpl.alphaX, tpl.alphaY, alphaMean, ...
            tpl.sigmaX, tpl.sigmaY, sigmaMean, ...
            tpl.r2X,    tpl.r2Y,    r2Mean, ...
            tpl.varExpX, tpl.varExpY, ...
            'VariableNames', [ ...
            "alphaX","alphaY","alphaMean", ...
            "sigmaX","sigmaY","sigmaMean", ...
            "r2X","r2Y","r2Mean","varExpX","varExpY"]);
        results = [results, tplT];
    end

    if doDet
        dAlphaMean = mean([det.alphaX, det.alphaY], 2);
        dSigmaMean = mean([det.sigmaX, det.sigmaY], 2);
        dR2Mean    = mean([det.r2X,    det.r2Y],    2);
        detT = table( ...
            det.alphaX, det.alphaY, dAlphaMean, ...
            det.sigmaX, det.sigmaY, dSigmaMean, ...
            det.r2X,    det.r2Y,    dR2Mean, det.polyOrd, ...
            'VariableNames', [ ...
            "det_alphaX","det_alphaY","det_alphaMean", ...
            "det_sigmaX","det_sigmaY","det_sigmaMean", ...
            "det_r2X","det_r2Y","det_r2Mean","det_polyOrd"]);
        results = [results, detT];
    end

    if doIRA
        iraAlphaMean = mean([ira.alphaX, ira.alphaY], 2);
        iraT = table( ...
            ira.alphaX, ira.alphaY, iraAlphaMean, ...
            'VariableNames', ["ira_alphaX","ira_alphaY","ira_alphaMean"]);
        results = [results, iraT];
    end
end

function printSummary_local(results, trials, doTpl, doDet, doIRA)
    db = string(trials(1).database);
    nT = height(results);
    fprintf('\n=== BIOLOGICAL NOISE v004: %s (%d trials) ===\n', db, nT);
    fprintf('f0: %.3f +/- %.3f Hz\n', ...
        mean(results.f0,'omitnan'), std(results.f0,'omitnan'));

    if doTpl
        fprintf('\n--- Template subtraction ---\n');
        fprintf('VarExp: X=%.4f Y=%.4f\n', ...
            mean(results.varExpX,'omitnan'), mean(results.varExpY,'omitnan'));
        fprintf('alpha:  %.3f +/- %.3f\n', ...
            mean(results.alphaMean,'omitnan'), std(results.alphaMean,'omitnan'));
        fprintf('sigma:  %.3f +/- %.3f mm\n', ...
            mean(results.sigmaMean,'omitnan'), std(results.sigmaMean,'omitnan'));
        fprintf('R2:     %.3f +/- %.3f\n', ...
            mean(results.r2Mean,'omitnan'), std(results.r2Mean,'omitnan'));
    end

    if doDet
        fprintf('\n--- Staged detrending ---\n');
        fprintf('alpha:  %.3f +/- %.3f\n', ...
            mean(results.det_alphaMean,'omitnan'), std(results.det_alphaMean,'omitnan'));
        fprintf('sigma:  %.3f +/- %.3f mm  (inflated by signal residuals)\n', ...
            mean(results.det_sigmaMean,'omitnan'), std(results.det_sigmaMean,'omitnan'));
    end

    if doIRA
        fprintf('\n--- IRASA on residuals (v004) ---\n');
        fprintf('alpha:  %.3f +/- %.3f  [residual IRASA, fHigh=min(20,fs/2-1)]\n', ...
            mean(results.ira_alphaMean,'omitnan'), std(results.ira_alphaMean,'omitnan'));
        fprintf('sigma: not reported (use template sigmaMean)\n');
    end

    if doTpl && doIRA
        fprintf('\n--- Template vs IRASA alpha ---\n');
        dAlpha = results.ira_alphaMean - results.alphaMean;
        fprintf('IRASA - Template: %.3f +/- %.3f  (positive = IRASA higher)\n', ...
            mean(dAlpha,'omitnan'), std(dAlpha,'omitnan'));
    end
end
