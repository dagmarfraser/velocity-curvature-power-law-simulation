function [nMaj, nMin, info] = endMatchedNoise_v001(segMaj, segMin, FS, alphaMaj, alphaMin, opts)
%ENDMATCHEDNOISE_V001 One realisation of the end-matched random-phase noise for both axes (the primary noise, 3b, 2026-11-14).
%   [nMaj, nMin, info] = endMatchedNoise_v001(segMaj, segMin, FS, alphaMaj, alphaMin) returns two noise series of the
%   segments' length L. segMaj/segMin are the residual axes cut to ONE end-matched segment (endToEndSegment_v003 on
%   [resid_maj, resid_min]); alphaMaj/alphaMin are alphaPWLW_v001 of those segments over [1, min(20, FS/2 - 1)] Hz.
%   Per axis: |FFT(segment)| on uniformly random phases (Schreiber & Schmitz 2000; the endToEndWhite / e2eFair construction
%   of endMismatchSensitivity_v002/v006), mean removed, scaled to std(segment); redrawn until |alpha(noise) - alpha| <= Tol
%   or MaxDraw draws (alpha selection, as in v006). An axis with no accepted draw is returned as NaN (Fail Loud: the caller
%   skips that replicate and counts it), never as the last unaccepted draw.
%   Draws come from the global stream: the caller sets the cell's substream first (docs/README_SeedPractice_v001.md).
%   The synthetic trajectory this noise is added to must have length L too (Dagmar, 2026-11-14: shorten the trajectory).
%   info: nDrawMaj, nDrawMin, acceptedMaj, acceptedMin, dAlphaMaj, dAlphaMin (of the returned draw; NaN if none accepted).
%   Validated on real residuals in docs/EndMismatch_FairShareRule_FullLength_20261114_v001.md (median filtered velocity and
%   acceleration noise 1.00 in Hickman and Cook, 1.01-1.02 Fraser).
arguments
    segMaj (:,1) double {mustBeFinite}
    segMin (:,1) double {mustBeFinite}
    FS (1,1) double {mustBePositive}
    alphaMaj (1,1) double
    alphaMin (1,1) double
    opts.Tol (1,1) double {mustBePositive} = 0.3
    opts.MaxDraw (1,1) double {mustBeInteger, mustBePositive} = 20
end
if numel(segMaj) ~= numel(segMin), error('endMatchedNoise:length', 'the two axes must share one segment (%d vs %d samples).', numel(segMaj), numel(segMin)); end
if ~isfinite(alphaMaj) || ~isfinite(alphaMin), error('endMatchedNoise:alpha', 'target alpha is not finite (%g, %g).', alphaMaj, alphaMin); end
band = [1 min(20, FS/2 - 1)];
[nMaj, kMaj, okMaj, dMaj] = drawAxis(segMaj, alphaMaj, FS, band, opts.Tol, opts.MaxDraw);
[nMin, kMin, okMin, dMin] = drawAxis(segMin, alphaMin, FS, band, opts.Tol, opts.MaxDraw);
info = struct('nDrawMaj', kMaj, 'nDrawMin', kMin, 'acceptedMaj', okMaj, 'acceptedMin', okMin, 'dAlphaMaj', dMaj, 'dAlphaMin', dMin);
end

function [u, k, ok, dA] = drawAxis(ref, aRef, FS, band, tol, maxDraw)
L = numel(ref);  A = abs(fft(ref - mean(ref)));  h = floor((L - 1) / 2);  sR = std(ref);  ok = false;  dA = NaN;
for k = 1:maxDraw
    Z = zeros(L, 1);  Z(2:h + 1) = A(2:h + 1) .* exp(2i * pi * rand(h, 1));  Z(L:-1:L - h + 1) = conj(Z(2:h + 1));
    if mod(L, 2) == 0, Z(L/2 + 1) = A(L/2 + 1) * sign(randn); end
    u = real(ifft(Z));  u = u - mean(u);  u = u / std(u) * sR;
    d = NaN;  try, d = alphaPWLW_v001(u, FS, band) - aRef; catch, end      % a failed estimate is a rejected draw
    if abs(d) <= tol, ok = true;  dA = d;  return, end
end
u = NaN(L, 1);
end
