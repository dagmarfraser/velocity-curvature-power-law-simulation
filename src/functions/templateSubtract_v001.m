function [resX, resY, fitX, fitY] = templateSubtract_v001(x, y, fs, f0, nH, nCW)
%TEMPLATESUBTRACT_V001 Windowed harmonic-regression template subtraction.
%
% Faithful extraction of runLoopClosureFftnoise_v009.m's local function
% templateSubtract_local -- algorithm unchanged, packaged standalone.
% v009's own internal copy deliberately left untouched, matching the same
% duplication-over-refactor rationale as findBothBranches_v001.m.
%
% Removes the periodic drawing template from a trajectory by fitting, in
% overlapping Hann-windowed segments, a harmonic regression (linear trend
% + nH harmonics of f0) and subtracting it -- NOT a PCA-ellipse residual,
% a Fourier/harmonic one. This distinction matters: it is what makes
% IRASA on the residual meaningful (removes the dense drawing-frequency
% harmonics that otherwise defeat IRASA's resampling-cancellation
% mechanism -- Finding #97) and is a different, complementary step from
% the PCA geometry fit (which supplies ellipse axes a,b,theta only).
%
% SYNTAX:
%   [resX, resY, fitX, fitY] = templateSubtract_v001(x, y, fs, f0, nH, nCW)
%
% INPUTS:
%   x, y - Raw trajectory coordinates, real column vectors, same length.
%   fs   - Sampling frequency in Hz.
%   f0   - Dominant tempo in Hz (e.g. from estimateF0_v001).
%   nH   - Number of harmonics of f0 included in the regression (v009
%          default: 4).
%   nCW  - Window width in cycles of f0 (v009 default: 4).
%
% OUTPUTS:
%   resX, resY - Residual (template-subtracted) trajectory, edge-clipped
%                by roughly win/4 samples at each end (win = round(nCW/f0*fs),
%                clamped to [round(2/f0*fs), N]). Shorter than the input
%                by construction -- this is the residual IRASA/noise
%                characterisation should run on, not a same-length signal.
%   fitX, fitY - The fitted (template) trajectory over the same clipped
%                range, i.e. x(clipped)-resX = fitX exactly. Used for
%                deterministic smoke-testing (add the original residual's
%                own spectral content back onto the fit and confirm the
%                pipeline recovers what it started with).
%
% See also: estimateF0_v001
%
% Extracted from: runLoopClosureFftnoise_v009.m (templateSubtract_local)
% Fraser, D.S. (2026)

arguments
    x (:,1) double
    y (:,1) double
    fs (1,1) double {mustBePositive}
    f0 (1,1) double {mustBePositive}
    nH (1,1) double {mustBeInteger, mustBePositive} = 4
    nCW (1,1) double {mustBePositive} = 4
end

N   = numel(x);
win = round(nCW/f0*fs);
hop = round(win/2);
win = min(win, N);
win = max(win, round(2/f0*fs));

if win < 4
    error('templateSubtract_v001:WindowTooSmall', '%s', ...
        sprintf(['Computed window (%d samples) is too small to fit a harmonic ' ...
        'regression -- trial too short relative to f0=%.3fHz at fs=%.1fHz. ' ...
        'Discard this trial rather than proceed.'], win, f0, fs));
end

resX = zeros(N,1); resY = zeros(N,1); ww = zeros(N,1);
s = 1;
while s + win - 1 <= N
    e   = s + win - 1;
    tw  = (0:win-1)'/fs;
    han = 0.5*(1 - cos(2*pi*(0:win-1)'/(win-1)));
    D   = [ones(win,1), tw];
    for h = 1:nH
        D(:,end+1) = cos(2*pi*h*f0*tw); %#ok<AGROW>
        D(:,end+1) = sin(2*pi*h*f0*tw); %#ok<AGROW>
    end
    Dw        = D .* han;
    bX        = Dw \ (x(s:e) .* han);
    bY        = Dw \ (y(s:e) .* han);
    resX(s:e) = resX(s:e) + (x(s:e) - D*bX) .* han;
    resY(s:e) = resY(s:e) + (y(s:e) - D*bY) .* han;
    ww(s:e)   = ww(s:e) + han;
    s         = s + hop;
end

if ~any(ww > 0)
    error('templateSubtract_v001:NoWindowsFit', '%s', ...
        'No windows fit inside the trial -- trial too short for this f0/nCW combination.');
end

ok = ww > 0;
resX(ok) = resX(ok) ./ ww(ok);
resY(ok) = resY(ok) ./ ww(ok);
cl   = max(round(win/4), 1);
cs   = cl;
ce   = min(max(N - cl, cs + 100), N);
fitX = x(cs:ce) - resX(cs:ce);
fitY = y(cs:ce) - resY(cs:ce);
resX = resX(cs:ce);
resY = resY(cs:ce);

end
