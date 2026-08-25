function [alphaEst, sigmaEst] = iraAlphaSigma_v001(x, fs, fLow, fHigh, hset)
% IRAALPHA_V001  Lean IRASA spectral exponent and sigma estimator.
%
% Implements Wen & Liu (2016) Irregular-Resampling Auto-Spectral Analysis.
% For each resampling factor h and its reciprocal 1/h:
%   - Resample signal by h, compute PSD via pmtm, re-interpolate onto
%     the original frequency axis
%   - Geometric mean of up- and down-sampled PSDs
% Median across all h pairs isolates the fractal (aperiodic) component.
% Periodic peaks at f0 and harmonics shift in frequency under resampling
% and are averaged away without requiring any knowledge of f0.
%
% INPUTS:
%   x      - signal (column vector, position in mm)
%   fs     - sampling rate (Hz)
%   fLow   - lower bound of spectral fitting band (Hz)
%   fHigh  - upper bound of spectral fitting band (Hz)
%   hset   - resampling factors, e.g. 1.1:0.05:1.9 (17 pairs)
%
% OUTPUTS:
%   alphaEst - spectral exponent (slope of log PSD vs log f, negated)
%   sigmaEst - noise sigma (mm): sqrt of integrated fractal PSD x 2
%
% HURST / HAUSDORFF EQUIVALENCES (for reference):
%   H = (alpha + 1) / 2   continuous extension valid for alpha > -1
%   D = 2 - H             Hausdorff / box-counting dimension of trajectory
%
% Author:  d.s.fraser@bham.ac.uk
% Created: 2026-04-18

arguments
    x      (:,1) double
    fs     (1,1) double {mustBePositive}
    fLow   (1,1) double {mustBePositive}
    fHigh  (1,1) double {mustBePositive}
    hset   (1,:) double = 1.1:0.05:1.9
end

[pOrig, fOrig] = pmtm(x, 4, [], fs);
nF = numel(fOrig);

geoMeans = nan(nF, numel(hset));

for hi = 1:numel(hset)
    h = hset(hi);

    xUp = resample(x, round(h * 1000), 1000);
    [pU, fU] = pmtm(xUp(:), 4, [], fs * h);
    pUI = interp1(fU, pU, fOrig, 'linear', NaN);

    xDn = resample(x, 1000, round(h * 1000));
    [pD, fD] = pmtm(xDn(:), 4, [], fs / h);
    pDI = interp1(fD, pD, fOrig, 'linear', NaN);

    geoMeans(:, hi) = sqrt(pUI .* pDI);
end

pFractal = median(geoMeans, 2, 'omitnan');

fitMask = fOrig >= fLow & fOrig <= fHigh & ~isnan(pFractal) & pFractal > 0;
if sum(fitMask) < 3
    alphaEst = NaN;
    sigmaEst = NaN;
    return
end
pFit     = polyfit(log10(fOrig(fitMask)), log10(pFractal(fitMask)), 1);
alphaEst = -pFit(1);

validMask = ~isnan(pFractal) & pFractal > 0;
noiseVar  = trapz(fOrig(validMask), pFractal(validMask)) * 2;
sigmaEst  = sqrt(max(noiseVar, 0));

end
