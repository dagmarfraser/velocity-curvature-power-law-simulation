function [alphaEst, sigmaEst, pFractal, fVec, ftInfo, pmInfo] = estimateIRASA(x, fs, fLow, fHigh, options)
%ESTIMATEIRASA IRASA spectral-exponent and noise-magnitude estimator.
%
% Toolbox wrap of PowerLawSimulationPreReg's iraAlphaSigma_v003.m (not
% v001 -- v003 adds the pFractal/fVec outputs needed for Lorentzian-knee
% work at no downside; same core estimator otherwise). Core algorithm
% (resampling-median IRASA on pmtm PSDs, Wen & Liu 2016 mechanism) is
% unchanged from source.
%
% SYNTAX:
%   [alphaEst, sigmaEst] = estimateIRASA(x, fs, fLow, fHigh)
%   [alphaEst, sigmaEst, pFractal, fVec] = estimateIRASA(x, fs, fLow, fHigh)
%   [___] = estimateIRASA(___, Hset=customHset)
%   [___, ftInfo] = estimateIRASA(___, CompareFieldTrip=true)
%   [___, ftInfo, pmInfo] = estimateIRASA(___, ComparePMTM=true)
%
% INPUTS:
%   x     - Time series, real column vector
%   fs    - Sampling rate in Hz
%   fLow  - Lower bound of the log-log fitting band, Hz
%   fHigh - Upper bound of the log-log fitting band, Hz
%   Hset  - Name-value. Resampling ratio set for the IRASA median.
%           Default 1.1:0.05:1.9 (matches source default).
%   CompareFieldTrip - Name-value logical, default false. If true, also
%           runs FieldTrip's ft_freqanalysis(method='irasa') on the same
%           signal and band for comparison, returned in ftInfo.
%   ComparePMTM - Name-value logical, default false. If true, also fits
%           a direct log-log slope to the un-resampled pmtm PSD (no
%           IRASA resampling-median step at all) over the same [fLow,
%           fHigh] band, returned in pmInfo. This is the comparison that
%           actually tests the resampling-median step in isolation:
%           CompareFieldTrip varies the PSD estimator (pmtm vs Hanning
%           FFT) while keeping resampling-median fixed, so HB-vs-FT
%           agreement cannot rule out both sharing a resampling-median
%           blind spot. ComparePMTM holds the PSD estimator fixed
%           (pmtm, same as HB) and removes resampling-median entirely,
%           isolating that step as the variable under test -- see
%           README_IRASAThreeWayAgreement_v001.md Section 2.
%
% OUTPUTS:
%   alphaEst - Spectral exponent (homebrew IRASA estimate)
%   sigmaEst - Noise magnitude, sqrt(integrated fractal PSD * 2), where
%              the integral is taken over the full available frequency
%              range of the fractal PSD estimate (all valid bins), NOT
%              restricted to [fLow, fHigh] -- that band is used only for
%              the alpha log-log fit. Do not expect sigmaEst to equal
%              std(x); it is a spectral-magnitude estimate, not a
%              reconstruction of the time-domain standard deviation.
%   pFractal - Aperiodic PSD, linear units, geometric median across
%              h-pairs (same length as fVec)
%   fVec     - Frequency axis, Hz
%   ftInfo   - Struct, populated only when CompareFieldTrip=true:
%       .available - logical, whether FieldTrip was found on the path
%       .alphaFT   - FieldTrip's recovered alpha (NaN if unavailable or
%                    the fit failed)
%       .gapHBFT   - abs(alphaEst - alphaFT) (NaN if either is NaN)
%       .message   - string, empty on success; explains why
%                    available=false or alphaFT=NaN otherwise
%       When CompareFieldTrip=false (default), ftInfo is returned as a
%       struct with available=false and message="CompareFieldTrip not requested".
%   pmInfo   - Struct, populated only when ComparePMTM=true:
%       .available - logical, always true when requested (no external
%                    dependency; pmtm ships with Signal Processing
%                    Toolbox, already required by the core estimator)
%       .alphaPM   - direct pmtm log-log slope estimate, no resampling
%                    (NaN if fewer than 3 valid bins in [fLow, fHigh])
%       .gapHBPM   - abs(alphaEst - alphaPM) (NaN if either is NaN)
%       .message   - string, empty on success; explains why alphaPM=NaN
%                    otherwise
%       When ComparePMTM=false (default), pmInfo is returned as a
%       struct with available=false and message="ComparePMTM not requested".
%
% VERIFIED RANGE (Finding #100, PowerLawSimulationPreReg):
%   Independent estimator recoverability -- i.e. this function's alpha
%   estimate agreeing with a known ground truth -- is verified for
%   alpha in [-2, 6]. This is a NARROWER claim than the generator's
%   numerical stability range (alpha in [-2, 7], see noiseXu); at
%   alpha=7 three independent estimators (this homebrew IRASA, FieldTrip
%   IRASA, and a resampling-free direct pmtm slope fit) diverge from one
%   another, consistent with an information limit at typical trial
%   lengths rather than a specific defect in any one estimator. Do not
%   treat this function's output as independently verified ground truth
%   outside alpha in [-2, 6]. ComparePMTM is the estimator most directly
%   responsible for establishing this range in the original finding --
%   prefer enabling it over CompareFieldTrip alone when checking whether
%   a given alphaEst is trustworthy, since HB-vs-FT agreement alone
%   cannot rule out a shared resampling-median artefact.
%
% FIELDTRIP DEPENDENCY POLICY:
%   CompareFieldTrip defaults to false, so estimateIRASA has no hard
%   FieldTrip dependency for its core (homebrew) estimate. Requesting
%   CompareFieldTrip=true without FieldTrip on the path does NOT error --
%   it returns ftInfo.available=false with an explanatory message and a
%   visible warning banner (Fail Loud, Never Fake: the degraded mode is
%   disclosed, not silently skipped). The `which('ft_defaults')` gate
%   (not `which('ft_freqanalysis')`) is used to detect FieldTrip
%   availability, following the corrected pattern in
%   benchmarkIRASA_v004.m: ft_freqanalysis.m sits in the FieldTrip root
%   and `which` can find it before ft_defaults has run, but FieldTrip's
%   internal private/ utility calls will still fail until ft_defaults has
%   added them to the path. Checking for ft_defaults itself is the
%   correct availability test.
%
% See also: noiseXu, checkSpectralKnee
%
% Toolbox: fractalnoise
% Extracted from: src/functions/iraAlphaSigma_v003.m (algorithm inlined
% below as a local function -- the toolbox does not depend on the
% original file existing outside its own tree)
% (PowerLawSimulationPreReg, session 52)
% Citation: Wen, H., & Liu, Z. (2016). Separating Fractal and
%   Oscillatory Components in the Power Spectrum of Neurophysiological
%   Signal. Brain Topography, 29(1), 13-26.
%   https://doi.org/10.1007/s10548-015-0448-0

arguments
    x (:, 1) double
    fs (1, 1) double {mustBePositive}
    fLow (1, 1) double {mustBePositive}
    fHigh (1, 1) double {mustBePositive}
    options.Hset (1, :) double = 1.1:0.05:1.9
    options.CompareFieldTrip (1, 1) logical = false
    options.ComparePMTM (1, 1) logical = false
end

[alphaEst, sigmaEst, pFractal, fVec, pRaw] = estimateIRASACore(x, fs, fLow, fHigh, options.Hset);

if options.CompareFieldTrip
    ftInfo = estimateIRASAFieldTripCompare(x, fs, fLow, fHigh, options.Hset, alphaEst);
else
    ftInfo = struct("available", false, "alphaFT", NaN, "gapHBFT", NaN, ...
        "message", "CompareFieldTrip not requested");
end

if options.ComparePMTM
    pmInfo = estimateIRASAPmtmCompare(fVec, pRaw, fLow, fHigh, alphaEst);
else
    pmInfo = struct("available", false, "alphaPM", NaN, "gapHBPM", NaN, ...
        "message", "ComparePMTM not requested");
end

end

function pmInfo = estimateIRASAPmtmCompare(fVec, pRaw, fLow, fHigh, alphaEst)
%ESTIMATEIRASAPMTMCOMPARE Direct pmtm log-log slope, no resampling-median.
%   Uses the same un-resampled pmtm PSD (pRaw, fVec) that
%   estimateIRASACore already computes once at h=1 before entering the
%   resampling-median loop -- no extra pmtm call. Isolates the
%   resampling-median step as the variable under test against alphaEst
%   (see README_IRASAThreeWayAgreement_v001.md Section 2, "PM").

fitMask = fVec >= fLow & fVec <= fHigh & ~isnan(pRaw) & pRaw > 0;
if sum(fitMask) < 3
    pmInfo = struct("available", true, "alphaPM", NaN, "gapHBPM", NaN, ...
        "message", "PMTM fit failed: fewer than 3 valid bins in fitting band");
    return
end

pFit = polyfit(log10(fVec(fitMask)), log10(pRaw(fitMask)), 1);
alphaPM = -pFit(1);

if isfinite(alphaEst) && isfinite(alphaPM)
    gapHBPM = abs(alphaEst - alphaPM);
else
    gapHBPM = NaN;
end

pmInfo = struct("available", true, "alphaPM", alphaPM, "gapHBPM", gapHBPM, ...
    "message", "");

end

function ftInfo = estimateIRASAFieldTripCompare(x, fs, fLow, fHigh, hset, alphaEst)
%ESTIMATEIRASAFIELDTRIPCOMPARE Optional FieldTrip IRASA comparison panel.
%   Availability gate follows benchmarkIRASA_v004.m: check ft_defaults,
%   not ft_freqanalysis (see FIELDTRIP DEPENDENCY POLICY above).

if isempty(which("ft_defaults"))
    warning("estimateIRASA:NoFieldTrip", "%s", ...
        "CompareFieldTrip=true requested but FieldTrip not found on path " + ...
        "(checked for ft_defaults). Returning ftInfo.available=false. " + ...
        "Install from https://www.fieldtriptoolbox.org and addpath the " + ...
        "FieldTrip root folder to enable this comparison.");
    ftInfo = struct("available", false, "alphaFT", NaN, "gapHBFT", NaN, ...
        "message", "FieldTrip not found on path (ft_defaults absent)");
    return
end

ft_defaults;
if isempty(which("ft_freqanalysis"))
    warning("estimateIRASA:FieldTripIncomplete", "%s", ...
        "ft_defaults ran but ft_freqanalysis still not found -- check your " + ...
        "FieldTrip installation. Returning ftInfo.available=false.");
    ftInfo = struct("available", false, "alphaFT", NaN, "gapHBFT", NaN, ...
        "message", "ft_defaults ran but ft_freqanalysis not found");
    return
end

N = numel(x);
ftData.trial = {x(:)'};
ftData.time = {(0:N - 1) / fs};
ftData.fsample = fs;
ftData.label = {'ch1'};
ftData.sampleinfo = [1, N]; %#ok<STRNU> ftData used inside evalc() string below; analyzer cannot trace it. Also avoids FieldTrip's fixsampleinfo reconstruction warning.
% FieldTrip's internal validation requires char vectors here, not string
% objects -- {"ch1"} throws "Cell array input must be a cell array of
% character vectors" inside ft_freqanalysis. Single-quoted char literal
% is required, not a style choice.

ftCfg = [];
ftCfg.method = 'irasa';
ftCfg.output = 'fractal';
ftCfg.foilim = [fLow fHigh];
ftCfg.taper = 'hanning';
ftCfg.pad = 'nextpow2';
ftCfg.hset = hset;
ftCfg.keeptrials = 'no';

alphaFT = NaN;
message = "";
try
    [~, ftFreq] = evalc("ft_freqanalysis(ftCfg, ftData)");
    fFT = ftFreq.freq;
    pFT = squeeze(ftFreq.powspctrm);
    fitMask = fFT >= fLow & fFT <= fHigh & ~isnan(pFT) & pFT > 0;
    if sum(fitMask) >= 3
        pFit = polyfit(log10(fFT(fitMask)), log10(pFT(fitMask)), 1);
        alphaFT = -pFit(1);
    else
        message = "FieldTrip fit failed: fewer than 3 valid bins in fitting band";
    end
catch ME
    message = "FieldTrip call failed: " + string(ME.message);
end

if isfinite(alphaEst) && isfinite(alphaFT)
    gapHBFT = abs(alphaEst - alphaFT);
else
    gapHBFT = NaN;
end

ftInfo = struct("available", true, "alphaFT", alphaFT, "gapHBFT", gapHBFT, ...
    "message", message);

end

function [alphaEst, sigmaEst, pFractal, fVec, pRaw] = estimateIRASACore(x, fs, fLow, fHigh, hset)
%ESTIMATEIRASACORE Inlined copy of iraAlphaSigma_v003's IRASA algorithm.
%   Bundled here, not called externally, so the toolbox has no dependency
%   on PowerLawSimulationPreReg's src/functions/ directory existing on
%   the path. Algorithm unchanged from source: resampling-median IRASA
%   on pmtm PSDs (Wen & Liu 2016 mechanism).
%
%   pRaw (fifth output) is the un-resampled h=1 pmtm PSD, computed here
%   regardless (needed for fVec) but previously discarded. Kept and
%   returned so ComparePMTM can reuse it without a second pmtm call --
%   this is exactly the "PM" estimator from
%   README_IRASAThreeWayAgreement_v001.md: same PSD estimator as HB,
%   no resampling-median step.

[pRaw, fVec] = pmtm(x, 4, [], fs);
nF = numel(fVec);
geoMeans = nan(nF, numel(hset));

for hi = 1:numel(hset)
    h = hset(hi);

    xUp = resample(x, round(h * 1000), 1000);
    [pU, fU] = pmtm(xUp(:), 4, [], fs * h);
    pUI = interp1(fU, pU, fVec, 'linear', NaN);

    xDn = resample(x, 1000, round(h * 1000));
    [pD, fD] = pmtm(xDn(:), 4, [], fs / h);
    pDI = interp1(fD, pD, fVec, 'linear', NaN);

    geoMeans(:, hi) = sqrt(pUI .* pDI);
end

pFractal = median(geoMeans, 2, 'omitnan');

fitMask = fVec >= fLow & fVec <= fHigh & ~isnan(pFractal) & pFractal > 0;
if sum(fitMask) < 3
    alphaEst = NaN;
    sigmaEst = NaN;
    return
end
pFit = polyfit(log10(fVec(fitMask)), log10(pFractal(fitMask)), 1);
alphaEst = -pFit(1);

validMask = ~isnan(pFractal) & pFractal > 0;
noiseVar = trapz(fVec(validMask), pFractal(validMask)) * 2;
sigmaEst = sqrt(max(noiseVar, 0));

end
