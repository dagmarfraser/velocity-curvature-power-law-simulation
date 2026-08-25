function [alphaEst, sigmaEst, pFractal, fVec] = iraAlphaSigma_v003(x, fs, fLow, fHigh, hset)
% IRAALPHA_V003  IRASA estimator — exposes fractal PSD for downstream use.
%
% Identical to iraAlphaSigma_v001 (pmtm-based) but adds optional 3rd and
% 4th outputs: pFractal (geometric-median aperiodic PSD) and fVec (Hz).
% These are required by checkLorentzianKnee_v001 for Lorentzian model
% fitting (Preston et al. 2026 knee observation diagnostic).
%
% INPUTS / OUTPUTS: identical to iraAlphaSigma_v001 except:
%   [alphaEst, sigmaEst]             — as v001
%   [alphaEst, sigmaEst, pFractal]   — also returns aperiodic PSD (linear)
%   [alphaEst, sigmaEst, pFractal, fVec] — also returns frequency axis (Hz)
%
% pFractal: column vector, same length as fVec, geometric median across h-
%           pairs, in linear PSD units (same as pmtm output). The alpha
%           slope is estimated from the log10-log10 fit in [fLow, fHigh].
%
% Author:  d.s.fraser@bham.ac.uk
% Created: 2026-06-24
% See also: iraAlphaSigma_v001, checkLorentzianKnee_v001

arguments
    x      (:,1) double
    fs     (1,1) double {mustBePositive}
    fLow   (1,1) double {mustBePositive}
    fHigh  (1,1) double {mustBePositive}
    hset   (1,:) double = 1.1:0.05:1.9
end

[~, fVec]  = pmtm(x, 4, [], fs);
nF         = numel(fVec);
geoMeans   = nan(nF, numel(hset));

for hi = 1:numel(hset)
    h   = hset(hi);

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
pFit     = polyfit(log10(fVec(fitMask)), log10(pFractal(fitMask)), 1);
alphaEst = -pFit(1);

validMask = ~isnan(pFractal) & pFractal > 0;
noiseVar  = trapz(fVec(validMask), pFractal(validMask)) * 2;
sigmaEst  = sqrt(max(noiseVar, 0));

end  % iraAlphaSigma_v003
