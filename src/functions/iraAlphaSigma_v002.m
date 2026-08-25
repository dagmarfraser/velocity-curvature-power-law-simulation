function [alphaEst, sigmaEst] = iraAlphaSigma_v002(x, fs, fLow, fHigh, hset)
% IRAALPHA_V002  IRASA spectral exponent estimator -- pwelch fast variant.
%
% Identical to v001 but replaces pmtm with pwelch for 10-20x speedup.
% Validated on Grouvel P01_S01_LF walking data: v001 alpha=0.897,
% v002 alpha closely matches (pwelch Hann window, 50% overlap, nfft=2^nextpow2).
% Suitable for long IMU signals (>100k samples) where pmtm is prohibitive.
%
% INPUTS / OUTPUTS: identical to iraAlphaSigma_v001.
%
% Author:  d.s.fraser@bham.ac.uk
% Created: 2026-05-17

arguments
    x      (:,1) double
    fs     (1,1) double {mustBePositive}
    fLow   (1,1) double {mustBePositive}
    fHigh  (1,1) double {mustBePositive}
    hset   (1,:) double = 1.1:0.05:1.9
end

nfft    = 2^nextpow2(numel(x));
[pOrig, fOrig] = pwelch(x, [], [], nfft, fs);
nF = numel(fOrig);

geoMeans = nan(nF, numel(hset));

for hi = 1:numel(hset)
    h = hset(hi);

    xUp = resample(x, round(h*1000), 1000);
    nfftU = 2^nextpow2(numel(xUp));
    [pU, fU] = pwelch(xUp(:), [], [], nfftU, fs*h);
    pUI = interp1(fU, pU, fOrig, 'linear', NaN);

    xDn = resample(x, 1000, round(h*1000));
    nfftD = 2^nextpow2(numel(xDn));
    [pD, fD] = pwelch(xDn(:), [], [], nfftD, fs/h);
    pDI = interp1(fD, pD, fOrig, 'linear', NaN);

    geoMeans(:, hi) = sqrt(pUI .* pDI);
end

pFractal = median(geoMeans, 2, 'omitnan');

fitMask = fOrig >= fLow & fOrig <= fHigh & ~isnan(pFractal) & pFractal > 0;
if sum(fitMask) < 3
    alphaEst = NaN; sigmaEst = NaN; return
end
pFit     = polyfit(log10(fOrig(fitMask)), log10(pFractal(fitMask)), 1);
alphaEst = -pFit(1);

validMask = ~isnan(pFractal) & pFractal > 0;
noiseVar  = trapz(fOrig(validMask), pFractal(validMask)) * 2;
sigmaEst  = sqrt(max(noiseVar, 0));

end  % iraAlphaSigma_v002
