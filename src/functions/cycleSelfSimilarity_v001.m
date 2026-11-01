function out = cycleSelfSimilarity_v001(x, y, fs, opts)
%CYCLESELFSIMILARITY_V001 Cycle-to-cycle self-similarity of one repeated-ellipse trial.
%
% Tangential speed is cut into windows one speed cycle long and every window
% is compared with its nearest non-overlapping neighbour using MATLAB's
% matrixProfile (Predictive Maintenance Toolbox, z-normalised Euclidean
% distance). The distance is converted to 1 - rho (rho = Pearson correlation
% between the two windows), which does not depend on window length or
% sampling rate and so can be compared across trials and datasets:
%
%   1 - rho = d^2 / (2 L)         (exact for z-normalised windows of length L)
%
% 0 means the cycle recurs exactly; larger values mean a cycle has no close
% relative elsewhere in the trial. Because windows have a fixed length, the
% measure reflects shape AND tempo drift: a wandering tempo moves windows out
% of alignment.
%
% SYNTAX:
%   out = cycleSelfSimilarity_v001(x, y, fs)
%   out = cycleSelfSimilarity_v001(x, y, fs, WindowsPerRev=1, MinWindows=4)
%
% INPUTS:
%   x, y  - Raw trajectory, finite column vectors of equal length.
%   fs    - Sampling rate (Hz).
% NAME-VALUE:
%   WindowsPerRev - Speed windows per revolution (default 2: an ellipse's
%                   tangential speed repeats twice per revolution; 1 = one
%                   full revolution).
%   WindowSamples - Window length in samples. Default NaN: derive from f0
%                   (estimateF0_v001). Set it to compare a manipulated trial
%                   on the clean trial's window.
%   MinWindows    - Minimum windows needed after trimming (default 3).
%   FilterParams  - [order, reference frame length] for SG case 6 of
%                   differentiateKinematicsEBR (default [4 17]).
%
% OUTPUT (struct):
%   f0, L, trim - tempo (Hz), window length and trimmed samples at each end
%   speed       - trimmed tangential speed that was profiled
%   dissim      - 1 - rho profile, one value per window start (index i is the
%                 window starting at original sample i + trim)
%   typical     - median(dissim): the headline self-similarity number
%   worst       - max(dissim): the least repeatable cycle (discord)
%   nWindows    - trimmed length expressed in windows
%
% Errors (identifier cycleSelfSimilarity:tooShort) if fewer than MinWindows
% windows remain after trimming one window at each end; no fallback is used.
%
% See also: estimateF0_v001, differentiateKinematicsEBR, matrixProfile
% Fraser, D.S. (2026)

arguments
    x (:,1) double {mustBeFinite}
    y (:,1) double {mustBeFinite}
    fs (1,1) double {mustBePositive}
    opts.WindowsPerRev (1,1) double {mustBePositive} = 2
    opts.WindowSamples (1,1) double = NaN
    opts.MinWindows (1,1) double {mustBePositive} = 3
    opts.FilterParams (1,2) double = [4 17]
end

if ~any(exist('matrixProfile') == [2 6]) %#ok<EXIST>
    error('cycleSelfSimilarity:noToolbox', '%s', ...
        'matrixProfile not found: Predictive Maintenance Toolbox (R2024b or later) is required.');
end
if numel(y) ~= numel(x)
    error('cycleSelfSimilarity:sizeMismatch', '%s', 'x and y must have the same length.');
end

f0 = estimateF0_v001(x, y, fs);
L  = opts.WindowSamples;
if isnan(L), L = round(fs / (opts.WindowsPerRev * f0)); end
if L < 2 || L ~= round(L)
    error('cycleSelfSimilarity:badWindow', 'Window length %g samples is not a usable integer.', L);
end

nTrim = numel(x) - 2*L;
if nTrim < opts.MinWindows * L
    error('cycleSelfSimilarity:tooShort', ...
        'Only %.2f windows (L = %d, f0 = %.2f Hz) remain after trimming; need %g.', ...
        max(nTrim, 0) / L, L, f0, opts.MinWindows);
end

[dx, dy] = differentiateKinematicsEBR(x, y, 6, opts.FilterParams, fs);
speed = hypot(dx(:,2), dy(:,2));
s  = speed(L+1 : end-L);
mp = matrixProfile(s, L, ExclusionZoneLength=L-1);   % L-1 excludes windows overlapping by any sample
d1 = mp(:).^2 / (2*L);

out = struct('f0', f0, 'L', L, 'trim', L, 'speed', s, 'dissim', d1, ...
    'typical', median(d1), 'worst', max(d1), 'nWindows', numel(s) / L);
end
