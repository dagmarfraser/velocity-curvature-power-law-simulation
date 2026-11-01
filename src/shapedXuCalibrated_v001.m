function [n, info] = shapedXuCalibrated_v001(resid, FS, alphaIn, fLow, fHigh, hset, tol, maxIter)
%SHAPEDXUCALIBRATED_V001  shaped_xu with a closed-loop correction of the carrier alpha (spec 10.60).
%   Plain shaped_xu (generateLoopClosureNoise_v003) delivers a lower IRASA alpha than the residual it is built from (about 0.13 to 0.18 lower in Fraser, spec 10.57 C).
%   Here the carrier alpha is adjusted until the IRASA alpha of the OUTPUT equals the IRASA alpha of the input residual under the SAME estimator and band
%   (iraAlphaSigma_v001). The Xu carrier draw is held fixed across iterations so only alpha moves. Output keeps the residual's SD.
%   [n, info] = shapedXuCalibrated_v001(resid, FS, alphaIn, fLow, fHigh, hset, tol, maxIter)  info: target, achieved, carrierAlpha, nIter, converged.
%   Fails loud: no convergence within maxIter returns info.converged = false (the caller counts and reports it); non-finite alpha is an error.
if nargin < 7, tol = 0.02; end
if nargin < 8, maxIter = 8; end
resid = resid(:); M = numel(resid); sig = std(resid);
target = iraAlphaSigma_v001(resid, FS, fLow, fHigh, hset);
if ~isfinite(target), error('shapedXuCalibrated:target','IRASA alpha of the input residual is not finite.'); end
Eemp = abs(fft(resid)); aC = alphaIn;
seed = randi(2^31-1);          % fixed carrier seed per call (drawn from the caller's stream so calls differ)
best = struct('err',Inf,'n',[],'aC',aC,'meas',NaN); info = struct('target',target,'converged',false);
for it = 1:maxIter
    sc = rng; rng(seed);
    base = generateCustomNoise_v004(M, aC, sig, FS);
    rng(sc);
    n = real(ifft(Eemp .* exp(1i*angle(fft(base)))));
    n = n * (sig / std(n));
    meas = iraAlphaSigma_v001(n, FS, fLow, fHigh, hset);
    if ~isfinite(meas), error('shapedXuCalibrated:meas','IRASA alpha of the surrogate is not finite at carrier alpha %.3f.', aC); end
    err = target - meas;
    if abs(err) < best.err, best = struct('err',abs(err),'n',n,'aC',aC,'meas',meas); end
    if abs(err) <= tol, info.converged = true; break, end
    aC = aC + err;
end
n = best.n; info.achieved = best.meas; info.carrierAlpha = best.aC; info.nIter = it; info.err = best.err;
end
