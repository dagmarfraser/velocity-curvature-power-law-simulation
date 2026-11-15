function [alpha, k, info] = alphaPWLW_v001(x, fs, band, k)
%ALPHAPWLW_V001 Prewhitened, band-limited local Whittle estimate of a spectral exponent (knife item K-IRASA, route B).
%   [alpha, k] = alphaPWLW_v001(x, fs, band)      chooses the differencing order k from a DFA-3 pre-estimate
%   [alpha, k] = alphaPWLW_v001(x, fs, band, k)   uses the given k (oracle arm)
%   Steps: difference x k times (removes the low-frequency power that leaks into the band), Hann-tapered
%   periodogram, divide by |1 - e^{-i lam}|^{2k} to recolour to the original series, then the Robinson (1995)
%   local Whittle objective for S ~ G lam^-alpha over the fitting band, minimised by fminbnd on [-2, 10].
%   Targets the in-band (band-effective) slope, like IRASA over the same band. k = floor((a~ + 1)/2) clamped
%   to 0..3, a~ = dfaAlpha_v001(x, 3) (the rule validated in fractal-noise stage 1). info.alphaDFA3 = a~.
arguments
    x double {mustBeVector}
    fs (1,1) double {mustBePositive}
    band (1,2) double {mustBePositive}
    k double = []
end
x = x(:) - mean(x);
info.alphaDFA3 = dfaAlpha_v001(x, 3);
if isempty(k), k = min(max(floor((info.alphaDFA3 + 1) / 2), 0), 3); end
xd = x;  for j = 1:k, xd = diff(xd); end
N = numel(xd);  w = 0.5 * (1 - cos(2*pi*(0:N-1)' / N));
I = abs(fft((xd - mean(xd)) .* w)).^2 / sum(w.^2);
f = (0:N-1)' * fs / N;  m = f >= band(1) & f <= band(2) & f < fs/2;
if nnz(m) < 10, error('alphaPWLW:band', 'Only %d periodogram bins in [%g %g] Hz (N = %d, fs = %g).', nnz(m), band, N, fs); end
lam = 2*pi*f(m)/fs;  Ib = I(m) ./ (2 - 2*cos(lam)).^k;
obj = @(a) log(mean(lam.^a .* Ib)) - a * mean(log(lam));
[alpha, ~, flag] = fminbnd(obj, -2, 10, optimset('TolX', 1e-6));
if flag ~= 1 || abs(alpha + 2) < 1e-3 || abs(alpha - 10) < 1e-3
    error('alphaPWLW:bound', 'Local Whittle optimiser on its bound or not converged (alpha = %.4f, flag %d).', alpha, flag);
end
info.nBins = nnz(m);
end
