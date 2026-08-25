function result = checkSpectralKnee(resid, fs, options)
%CHECKSPECTRALKNEE Lorentzian-vs-fixed spectral knee diagnostic (single trial).
%
% Toolbox generalisation of PowerLawSimulationPreReg's
% checkLorentzianKnee_v001.m per-trial fitting core (~lines 90-140:
% lorentz_log10 model + per-trial fminsearch fit). The source script's
% dataset-loop, import, and figure-generation scaffolding is NOT carried
% over -- this function takes one residual series in and returns one
% diagnostic result out; call it per trial from your own loop.
%
% BACKGROUND (Preston, Smith & Voytek, 2026, Nature Human Behaviour):
%   Neural field potential spectra often follow a Lorentzian (kneed) form
%   rather than a pure power law: power is approximately flat below a
%   "knee frequency" f_knee, then decays as 1/f^alpha above it. IRASA
%   assumes scale-free (fractal) input; if a knee falls within the
%   fitting band, the returned alpha is biased by the transition region.
%   This function tests whether a given residual's spectrum is better
%   explained by a Lorentzian (kneed) model than a fixed power law, via
%   AIC comparison, and whether any fitted knee falls inside the
%   requested fitting band (bias risk).
%
% MODELS (fitted in log10-space to the IRASA fractal PSD):
%   Fixed:  log10(P) = b - alpha * log10(f)                [2 parameters]
%   Knee:   log10(P) = b - (alpha/2) * log10(k^2 + f^2)    [3 parameters]
%           where k = knee_hz, fitted unconstrained via th(2) = log10(k)
%
% AIC comparison: AIC = N*ln(RSS/N) + 2*nParams
%   deltaAIC = AIC_knee - AIC_fixed
%   deltaAIC < -10: strong support for knee model
%   deltaAIC < -2:  marginal support
%   deltaAIC >= -2: no meaningful support for knee (fixed model adequate)
%
% SYNTAX:
%   result = checkSpectralKnee(resid, fs)
%   result = checkSpectralKnee(resid, fs, FitBandLow=1.0, WideBandLow=0.15)
%
% INPUTS:
%   resid - Residual time series, real column vector. Internally
%           demeaned before analysis (matches source: sig = x - mean(x)).
%   fs    - Sampling rate in Hz.
%   FitBandLow  - Name-value. Lower bound of the standard (narrow)
%                 fitting band, Hz. Default 1.0.
%   WideBandLow - Name-value. Lower bound of the wide band used for the
%                 knee search itself. Default 0.15.
%   FitBandHigh - Name-value. Upper bound of both bands, Hz. Default 0,
%                 meaning auto = fs/2.2.
%   Hset        - Name-value. IRASA resampling ratio set, passed through
%                 to estimateIRASA. Default 1.1:0.05:1.9.
%   ExtraKneeStarts - Name-value. Additional knee-frequency (Hz) starting
%                 guesses for the fminsearch fit, e.g. a generator's
%                 known theoretical knee if you are validating a
%                 surrogate. Default [] (empty). Always combined with two
%                 built-in starts (WideBandLow and the geometric mean of
%                 the wide band) regardless of this argument.
%
% OUTPUT (struct):
%   .alphaStd    - Power-law exponent fit on [FitBandLow, fHigh]
%   .alphaWide   - Power-law exponent fit on [WideBandLow, fHigh]
%                  (the fixed-model comparator for the AIC test)
%   .kneeHz      - Best-fit knee frequency, Hz
%   .deltaAIC    - AIC_knee - AIC_fixed (more negative = stronger knee support)
%   .verdict     - "strong" | "marginal" | "none"
%   .kneeInBand  - logical; true if kneeHz falls within [FitBandLow, fHigh]
%                  (bias risk: an IRASA alpha fit over that band would be
%                  biased by the knee transition)
%   .fHigh       - Resolved upper band bound, Hz
%   .fBandLow    - FitBandLow, echoed for convenience
%   .fWideLow    - WideBandLow, echoed for convenience
%   .fRes        - Approximate frequency resolution, fs/numel(resid), Hz
%   .nSamples    - numel(resid)
%
% NOTES:
%   - Fail Loud, Never Fake: if the residual is too short or the fractal
%     PSD yields fewer than 5 valid bins in the wide band, this function
%     errors rather than returning a degraded/placeholder result. The
%     source script's per-trial loop silently skipped such trials
%     (acceptable in a many-trials batch context); a single-call function
%     has no such context to skip within, so it surfaces the failure.
%
% See also: estimateIRASA
%
% Toolbox: fractalnoise
% Extracted from: src/checkLorentzianKnee_v001.m (per-trial fitting core)
% (PowerLawSimulationPreReg, session 52)
% Citation: Preston, M., Smith, S. & Voytek, B. (2026). Potential
%   mechanisms and functional significance of aperiodic neural activity.
%   Nature Human Behaviour.

arguments
    resid (:, 1) double
    fs (1, 1) double {mustBePositive}
    options.FitBandLow (1, 1) double {mustBePositive} = 1.0
    options.WideBandLow (1, 1) double {mustBePositive} = 0.15
    options.FitBandHigh (1, 1) double {mustBeNonnegative} = 0
    options.Hset (1, :) double = 1.1:0.05:1.9
    options.ExtraKneeStarts (1, :) double = []
end

sig = resid - mean(resid);
nSamples = numel(sig);
if nSamples < 80
    error("checkSpectralKnee:TooShort", ...
        "Residual has %d samples; at least 80 are required for a " + ...
        "meaningful spectral knee fit.", nSamples);
end

if options.FitBandHigh > 0
    fHigh = options.FitBandHigh;
else
    fHigh = fs / 2.2;
end
fBandLow = options.FitBandLow;
fWideLow = options.WideBandLow;

[~, ~, pFrac, fVec] = estimateIRASA(sig, fs, fWideLow, fHigh, Hset=options.Hset);
fVec = fVec(:);
pFrac = pFrac(:);

% Standard-band alpha (reference power-law fit; not itself part of the
% AIC comparison, which uses the wide-band fixed model as its comparator)
maskStd = fVec >= fBandLow & fVec <= fHigh & ~isnan(pFrac) & pFrac > 0;
if sum(maskStd) < 3
    error("checkSpectralKnee:StdBandInsufficient", "%s", ...
        "Fewer than 3 valid PSD bins in the standard fitting band " + ...
        "[" + fBandLow + ", " + fHigh + "] Hz -- cannot fit alphaStd.");
end
pStd = polyfit(log10(fVec(maskStd)), log10(pFrac(maskStd)), 1);
alphaStd = -pStd(1);

% Wide-band fixed model (2 parameters): the knee model's comparator
maskWide = fVec >= fWideLow & fVec <= fHigh & ~isnan(pFrac) & pFrac > 0;
if sum(maskWide) < 5
    error("checkSpectralKnee:WideBandInsufficient", "%s", ...
        "Fewer than 5 valid PSD bins in the wide fitting band " + ...
        "[" + fWideLow + ", " + fHigh + "] Hz -- cannot fit the knee model.");
end
fWide = fVec(maskWide);
logPWide = log10(pFrac(maskWide));
nWide = numel(fWide);
pWide = polyfit(log10(fWide), logPWide, 1);
bWide = pWide(2);
alphaWide = -pWide(1);
rssFixed = sum((logPWide - (bWide - alphaWide * log10(fWide))) .^ 2);
aicFixed = nWide * log(rssFixed / nWide) + 2 * 2;

% Knee model (3 parameters), multi-start fminsearch
lorentzLog10 = @(th, f) th(1) - (th(3) / 2) .* log10(10 ^ (2 * th(2)) + f .^ 2);
sseFcn = @(th) sum((logPWide - lorentzLog10(th, fWide)) .^ 2);

kneeStartsHz = unique([fWideLow, sqrt(fWideLow * fHigh), fBandLow, options.ExtraKneeStarts]);
starts = [repmat(bWide, numel(kneeStartsHz), 1), log10(kneeStartsHz(:)), repmat(alphaWide, numel(kneeStartsHz), 1)];

fminOpts = optimset("Display", "off", "TolFun", 1e-10, "TolX", 1e-10, "MaxFunEvals", 10000);
bestRss = Inf;
thBest = starts(1, :);
for si = 1:size(starts, 1)
    [thTry, rssTry] = fminsearch(sseFcn, starts(si, :), fminOpts);
    if rssTry < bestRss
        bestRss = rssTry;
        thBest = thTry;
    end
end

kneeHz = 10 ^ thBest(2);
aicKnee = nWide * log(bestRss / nWide) + 2 * 3;
deltaAIC = aicKnee - aicFixed;

if deltaAIC < -10
    verdict = "strong";
elseif deltaAIC < -2
    verdict = "marginal";
else
    verdict = "none";
end

result.alphaStd = alphaStd;
result.alphaWide = alphaWide;
result.kneeHz = kneeHz;
result.deltaAIC = deltaAIC;
result.verdict = verdict;
result.kneeInBand = kneeHz >= fBandLow && kneeHz <= fHigh;
result.fHigh = fHigh;
result.fBandLow = fBandLow;
result.fWideLow = fWideLow;
result.fRes = fs / nSamples;
result.nSamples = nSamples;

end
