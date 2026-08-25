function results = characteriseBiologicalNoise_v002(trials, opts)
% characteriseBiologicalNoise_v002 Biological noise characterisation.
%
% Estimates the biological noise component in periodic movement data
% (ellipse tracing etc.) using one or both of two methods:
%
%   'template' — Sliding-window sinusoidal template subtraction (task-
%                geometry-informed, avoids smoother circularity).
%   'detrend'  — Staged polynomial detrending per prereg Section 7.3
%                (linear first, polynomial if systematic structure remains).
%   'both'     — Both methods, side by side. Default for ellipse databases.
%
% For 'both', the output table contains template columns (alphaX etc.) and
% parallel detrend columns (det_alphaX etc.) for direct comparison.
%
% Generalised from characteriseBiologicalNoise_v001.m.
%
% USAGE:
%   R = characteriseBiologicalNoise_v002(trials)
%   R = characteriseBiologicalNoise_v002(trials, Method="detrend")
%   R = characteriseBiologicalNoise_v002(trials, Method="both", NHarmonics=6)
%
% PREREG REF: Section 7.3 (noise characterisation via pmtm)
% Fraser, D.S. (2025)

    arguments
        trials (:,1) struct
        opts.Method      (1,1) string {mustBeMember(opts.Method, ...
            ["template", "detrend", "both"])} = "both"
        opts.NHarmonics  (1,1) double = 4
        opts.NCyclesWin  (1,1) double = 4
        opts.TBP         (1,1) double = 3
        opts.FitBandLow  (1,1) double = 1        % Hz
        opts.FitBandHigh (1,1) double = 0        % 0 = auto (fs/2.2)
        opts.MaxPolyOrder (1,1) double = 3       % ceiling for staged detrend
        opts.Verbose     (1,1) logical = true
    end

    doTemplate = ismember(opts.Method, ["template", "both"]);
    doDetrend  = ismember(opts.Method, ["detrend",  "both"]);
    nTrials    = numel(trials);

    % --- Pre-allocate: common ---
    allF0      = NaN(nTrials, 1);
    allSubj    = strings(nTrials, 1);
    allTrial   = strings(nTrials, 1);

    % --- Pre-allocate: template ---
    tpl = struct( ...
        'alphaX', NaN(nTrials,1), 'alphaY', NaN(nTrials,1), ...
        'sigmaX', NaN(nTrials,1), 'sigmaY', NaN(nTrials,1), ...
        'r2X',    NaN(nTrials,1), 'r2Y',    NaN(nTrials,1), ...
        'varExpX',NaN(nTrials,1), 'varExpY',NaN(nTrials,1));

    % --- Pre-allocate: detrend ---
    det = struct( ...
        'alphaX', NaN(nTrials,1), 'alphaY', NaN(nTrials,1), ...
        'sigmaX', NaN(nTrials,1), 'sigmaY', NaN(nTrials,1), ...
        'r2X',    NaN(nTrials,1), 'r2Y',    NaN(nTrials,1), ...
        'polyOrd',NaN(nTrials,1));

    for ii = 1:nTrials
        t  = trials(ii);
        fs = t.fs;
        allSubj(ii)  = string(t.subjectID);
        allTrial(ii) = string(t.trialID);

        fitHigh = opts.FitBandHigh;
        if fitHigh == 0, fitHigh = fs / 2.2; end

        %% ---- Estimate f0 (needed by both methods) ----
        f0 = estimateF0(t, fs);
        allF0(ii) = f0;

        %% ---- Template subtraction path ----
        if doTemplate
            minWinSamples = round(2 / max(f0, 0.1) * fs);
            if numel(t.x) < minWinSamples
                % Trial too short for even one template window; leave NaN
                if opts.Verbose
                    fprintf("  Trial %d (%s): too short for template (%d < %d samples)\n", ...
                        ii, allTrial(ii), numel(t.x), minWinSamples);
                end
            else
            [resX, resY, vexpX, vexpY] = templateSubtract( ...
                t, fs, f0, opts.NHarmonics, opts.NCyclesWin);
            tpl.varExpX(ii) = vexpX;
            tpl.varExpY(ii) = vexpY;
            tpl.sigmaX(ii)  = std(resX);
            tpl.sigmaY(ii)  = std(resY);
            [tpl.alphaX(ii), tpl.r2X(ii)] = fitSpectralSlope( ...
                resX, fs, opts.TBP, opts.FitBandLow, fitHigh);
            [tpl.alphaY(ii), tpl.r2Y(ii)] = fitSpectralSlope( ...
                resY, fs, opts.TBP, opts.FitBandLow, fitHigh);
            end  % minWinSamples guard
        end

        %% ---- Staged detrending path (prereg Section 7.3) ----
        if doDetrend
            [dResX, dResY, polyOrd] = stagedDetrend( ...
                t, fs, f0, opts.MaxPolyOrder);
            det.polyOrd(ii) = polyOrd;
            det.sigmaX(ii)  = std(dResX);
            det.sigmaY(ii)  = std(dResY);
            [det.alphaX(ii), det.r2X(ii)] = fitSpectralSlope( ...
                dResX, fs, opts.TBP, opts.FitBandLow, fitHigh);
            [det.alphaY(ii), det.r2Y(ii)] = fitSpectralSlope( ...
                dResY, fs, opts.TBP, opts.FitBandLow, fitHigh);
        end

        if opts.Verbose && mod(ii, 20) == 0
            fprintf("  Processed %d/%d\n", ii, nTrials);
        end
    end

    %% ---- Build output table ----
    results = buildResultsTable(allTrial, allSubj, allF0, ...
        tpl, det, doTemplate, doDetrend);

    %% ---- Summary ----
    if opts.Verbose
        printSummary(results, trials, doTemplate, doDetrend);
    end

end

%% ========== LOCAL FUNCTIONS ==========

function f0 = estimateF0(t, fs)
% Estimate fundamental tracing frequency via FFT peak, cross-validated
% between x and y channels.
    N    = numel(t.x);
    nfft = 2^nextpow2(4*N);
    xDet = detrend(double(t.x), "linear");
    yDet = detrend(double(t.y), "linear");
    Xf   = abs(fft(xDet, nfft));
    Yf   = abs(fft(yDet, nfft));
    fAxis = (0:nfft-1)' * fs / nfft;
    band  = fAxis > 0.1 & fAxis < 5;

    [~, pkX] = max(Xf(band));
    [~, pkY] = max(Yf(band));
    idx = find(band);
    f0x = fAxis(idx(pkX));
    f0y = fAxis(idx(pkY));

    if abs(f0x - f0y)/max(f0x, 0.01) < 0.2
        f0 = mean([f0x, f0y]);
    elseif max(Xf(band)) > max(Yf(band))
        f0 = f0x;
    else
        f0 = f0y;
    end
end

function [resX, resY, varExpX, varExpY] = templateSubtract( ...
        t, fs, f0, nHarmonics, nCyclesWin)
% Sliding-window sinusoidal template subtraction (H+N model).
    N = numel(t.x);
    winSamples = round(nCyclesWin / f0 * fs);
    hopSamples = round(winSamples / 2);
    winSamples = min(winSamples, N);
    winSamples = max(winSamples, round(2/f0*fs));

    resX = zeros(N, 1);
    resY = zeros(N, 1);
    winWeight = zeros(N, 1);

    startIdx = 1;
    while startIdx + winSamples - 1 <= N
        endIdx = startIdx + winSamples - 1;
        tWin   = (0:winSamples-1)' / fs;
        hannWin = 0.5 * (1 - cos(2*pi*(0:winSamples-1)'/(winSamples-1)));

        % Design matrix: DC + linear drift + harmonics
        D = ones(winSamples, 1);
        D(:,2) = tWin;
        for h = 1:nHarmonics
            D(:, end+1) = cos(2*pi*h*f0*tWin); %#ok<AGROW>
            D(:, end+1) = sin(2*pi*h*f0*tWin); %#ok<AGROW>
        end

        xWin = double(t.x(startIdx:endIdx));
        yWin = double(t.y(startIdx:endIdx));
        Dw = D .* hannWin;

        bX = Dw \ (xWin .* hannWin);
        bY = Dw \ (yWin .* hannWin);

        resX(startIdx:endIdx) = resX(startIdx:endIdx) + (xWin - D*bX) .* hannWin;
        resY(startIdx:endIdx) = resY(startIdx:endIdx) + (yWin - D*bY) .* hannWin;
        winWeight(startIdx:endIdx) = winWeight(startIdx:endIdx) + hannWin;

        startIdx = startIdx + hopSamples;
    end

    % Normalise overlap-add
    valid = winWeight > 0;
    resX(valid) = resX(valid) ./ winWeight(valid);
    resY(valid) = resY(valid) ./ winWeight(valid);

    % Clip edges
    edgeClip = max(round(winSamples/4), 1);
    cStart = edgeClip;
    cEnd   = min(max(N - edgeClip, cStart + 100), N);
    xClip  = double(t.x(cStart:cEnd));
    yClip  = double(t.y(cStart:cEnd));
    resX   = resX(cStart:cEnd);
    resY   = resY(cStart:cEnd);

    varExpX = 1 - var(resX)/var(xClip);
    varExpY = 1 - var(resY)/var(yClip);
end

function [resX, resY, polyOrd] = stagedDetrend(t, fs, f0, maxPolyOrd)
% Staged polynomial detrending per prereg Section 7.3.
%   1. Linear detrend
%   2. Check residuals for systematic periodic structure
%   3. Escalate to polynomial if structure remains
%
% For periodic movement (ellipses), linear detrend removes DC+slope but
% leaves sinusoidal structure intact. Polynomial detrend (order 3) removes
% cubic trends but still cannot remove periodicity. The pmtm spectrum of
% these residuals will therefore be dominated by movement harmonics rather
% than noise — this is the expected and informative outcome.

    xRaw = double(t.x(:));
    yRaw = double(t.y(:));
    N    = numel(xRaw);

    % Stage 1: linear detrend
    resX = detrend(xRaw, "linear");
    resY = detrend(yRaw, "linear");
    polyOrd = 1;

    % Stage 2: check for systematic periodic structure
    % Compute autocorrelation at the expected period (1/f0 samples)
    periodSamples = round(fs / max(f0, 0.01));
    if periodSamples > 1 && periodSamples < N/2
        acfX = xcorr(resX, periodSamples, "coeff");
        acfY = xcorr(resY, periodSamples, "coeff");
        % Peak at lag = periodSamples (index periodSamples+1 in the xcorr output)
        peakLagX = acfX(end);
        peakLagY = acfY(end);
        periodicityDetected = max(abs(peakLagX), abs(peakLagY)) > 0.3;
    else
        periodicityDetected = false;
    end

    % Stage 3: escalate to polynomial if periodic structure detected
    if periodicityDetected && maxPolyOrd > 1
        polyOrd = maxPolyOrd;
        tVec = (0:N-1)' / fs;
        % Fit polynomial via Vandermonde matrix
        V = zeros(N, polyOrd + 1);
        for k = 0:polyOrd
            V(:, k+1) = tVec.^k;
        end
        resX = xRaw - V * (V \ xRaw);
        resY = yRaw - V * (V \ yRaw);
    end

    % Clip 20 samples from each edge (same as PowerLawEmpirical.m convention)
    edgeClip = min(20, floor(N/10));
    if edgeClip > 0
        resX = resX(edgeClip+1:end-edgeClip);
        resY = resY(edgeClip+1:end-edgeClip);
    end
end

function [alpha, r2] = fitSpectralSlope(residuals, fs, tbp, fLow, fHigh)
% Fit log-log spectral slope via pmtm. Returns alpha = -slope and R2.
    alpha = NaN;
    r2    = NaN;
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

function results = buildResultsTable(trialIDs, subjIDs, f0, tpl, det, doTpl, doDet)
% Assemble output table depending on which methods were run.
    results = table(trialIDs, subjIDs, f0, ...
        'VariableNames', ["trialID", "subjectID", "f0"]);

    if doTpl
        alphaMean = mean([tpl.alphaX, tpl.alphaY], 2);
        sigmaMean = mean([tpl.sigmaX, tpl.sigmaY], 2);
        r2Mean    = mean([tpl.r2X,    tpl.r2Y],    2);
        tplT = table(tpl.alphaX, tpl.alphaY, alphaMean, ...
            tpl.sigmaX, tpl.sigmaY, sigmaMean, ...
            tpl.r2X, tpl.r2Y, r2Mean, ...
            tpl.varExpX, tpl.varExpY, ...
            'VariableNames', ["alphaX","alphaY","alphaMean", ...
            "sigmaX","sigmaY","sigmaMean", ...
            "r2X","r2Y","r2Mean","varExpX","varExpY"]);
        results = [results, tplT];
    end

    if doDet
        dAlphaMean = mean([det.alphaX, det.alphaY], 2);
        dSigmaMean = mean([det.sigmaX, det.sigmaY], 2);
        dR2Mean    = mean([det.r2X,    det.r2Y],    2);
        detT = table(det.alphaX, det.alphaY, dAlphaMean, ...
            det.sigmaX, det.sigmaY, dSigmaMean, ...
            det.r2X, det.r2Y, dR2Mean, det.polyOrd, ...
            'VariableNames', ["det_alphaX","det_alphaY","det_alphaMean", ...
            "det_sigmaX","det_sigmaY","det_sigmaMean", ...
            "det_r2X","det_r2Y","det_r2Mean","det_polyOrd"]);
        results = [results, detT];
    end
end

function printSummary(results, trials, doTpl, doDet)
% Print human-readable comparison summary.
    db = string(trials(1).database);
    nT = height(results);
    fprintf("\n=== BIOLOGICAL NOISE: %s (%d trials) ===\n", db, nT);
    fprintf("f0: %.3f +/- %.3f Hz\n", ...
        mean(results.f0, "omitnan"), std(results.f0, "omitnan"));

    if doTpl
        fprintf("\n--- METHOD: Template subtraction ---\n");
        fprintf("VarExp: X=%.4f Y=%.4f\n", ...
            mean(results.varExpX, "omitnan"), mean(results.varExpY, "omitnan"));
        fprintf("alpha (position): %.2f +/- %.2f  (after diff: %.2f)\n", ...
            mean(results.alphaMean, "omitnan"), ...
            std(results.alphaMean, "omitnan"), ...
            mean(results.alphaMean, "omitnan") - 2);
        fprintf("sigma: %.3f +/- %.3f\n", ...
            mean(results.sigmaMean, "omitnan"), std(results.sigmaMean, "omitnan"));
        fprintf("R2 spectral fit: %.3f +/- %.3f\n", ...
            mean(results.r2Mean, "omitnan"), std(results.r2Mean, "omitnan"));
    end

    if doDet
        fprintf("\n--- METHOD: Staged detrending (prereg Section 7.3) ---\n");
        fprintf("Polynomial order used: %.1f +/- %.1f\n", ...
            mean(results.det_polyOrd, "omitnan"), std(results.det_polyOrd, "omitnan"));
        fprintf("alpha (position): %.2f +/- %.2f  (after diff: %.2f)\n", ...
            mean(results.det_alphaMean, "omitnan"), ...
            std(results.det_alphaMean, "omitnan"), ...
            mean(results.det_alphaMean, "omitnan") - 2);
        fprintf("sigma: %.3f +/- %.3f\n", ...
            mean(results.det_sigmaMean, "omitnan"), std(results.det_sigmaMean, "omitnan"));
        fprintf("R2 spectral fit: %.3f +/- %.3f\n", ...
            mean(results.det_r2Mean, "omitnan"), std(results.det_r2Mean, "omitnan"));
    end

    if doTpl && doDet
        fprintf("\n--- COMPARISON: Template vs Detrend ---\n");
        dAlpha = results.alphaMean - results.det_alphaMean;
        dSigma = results.sigmaMean - results.det_sigmaMean;
        fprintf("Delta alpha (template - detrend): %.2f +/- %.2f\n", ...
            mean(dAlpha, "omitnan"), std(dAlpha, "omitnan"));
        fprintf("Delta sigma (template - detrend): %.3f +/- %.3f\n", ...
            mean(dSigma, "omitnan"), std(dSigma, "omitnan"));
        fprintf("Detrend sigma / Template sigma: %.1f +/- %.1f\n", ...
            mean(results.det_sigmaMean ./ results.sigmaMean, "omitnan"), ...
            std(results.det_sigmaMean ./ results.sigmaMean, "omitnan"));
    end
end
