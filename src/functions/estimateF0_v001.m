function f0 = estimateF0_v001(x, y, fs)
%ESTIMATEF0_V001 Estimate a trial's dominant movement tempo from raw x,y.
%
% Faithful extraction of runLoopClosureFftnoise_v009.m's local function
% estimateF0_local -- algorithm unchanged, packaged standalone so it can
% be reused outside that dataset-specific runner (e.g. by a standalone
% per-trial worker taking raw trajectory input directly). v009's own
% internal copy deliberately left untouched, matching the same
% duplication-over-refactor rationale as findBothBranches_v001.m.
%
% Peak-picking on a 4x zero-padded FFT of the linearly-detrended x and y
% channels separately, restricted to a plausible drawing-tempo band
% (0.1-5 Hz). If the two channels' peaks agree within 20%, returns their
% mean; otherwise returns whichever channel has the larger peak magnitude
% (the more confidently periodic axis).
%
% SYNTAX:
%   f0 = estimateF0_v001(x, y, fs)
%
% INPUTS:
%   x, y - Raw trajectory coordinates, real column vectors, same length.
%          Units do not matter (detrended before any spectral estimate).
%   fs   - Sampling frequency in Hz.
%
% OUTPUT:
%   f0   - Estimated dominant tempo in Hz.
%
% See also: templateSubtract_v001
%
% Extracted from: runLoopClosureFftnoise_v009.m (estimateF0_local)
% Fraser, D.S. (2026)

arguments
    x (:,1) double
    y (:,1) double
    fs (1,1) double {mustBePositive}
end

N    = numel(x); nfft = 2^nextpow2(4*N);
Xf   = abs(fft(detrend(x(:), 'linear'), nfft));
Yf   = abs(fft(detrend(y(:), 'linear'), nfft));
fAx  = (0:nfft-1)' * fs / nfft;
band = fAx > 0.1 & fAx < 5; idx = find(band);
[~, pkX] = max(Xf(band)); [~, pkY] = max(Yf(band));
f0x = fAx(idx(pkX)); f0y = fAx(idx(pkY));
if abs(f0x - f0y) / max(f0x, 0.01) < 0.2
    f0 = mean([f0x, f0y]);
elseif max(Xf(band)) > max(Yf(band))
    f0 = f0x;
else
    f0 = f0y;
end

end
