function cmap = errorColormap_v001(nColors, threshold, maxErr)
% errorColormap_v001  Colormap for |beta_rec - beta_gen| with a solid
% "faithful" band for errors below a specified threshold.
%
% Maps [0, threshold] to a single saturated green ("within clinical
% detection noise") and (threshold, maxErr] to a warm gradient from
% yellow through red to dark maroon ("how wrong"). Values above maxErr
% saturate to the darkest shade.
%
% Default threshold 0.03 is the minimum detectable change (MDC) for
% ASD vs CTRL beta differences (Cook et al. 2023; cf. Fraser et al. 2025).
%
% USAGE:
%   cmap = errorColormap_v001();                   % 256x3, threshold=0.03, maxErr=0.5
%   cmap = errorColormap_v001(256, 0.03, 0.5);     % explicit
%
% Fraser, D.S. (2025)

arguments
    nColors   (1,1) double = 256
    threshold (1,1) double = 0.03
    maxErr    (1,1) double = 0.5
end

if threshold >= maxErr
    error("errorColormap:badRange", "threshold (%g) must be < maxErr (%g)", ...
        threshold, maxErr);
end

% Proportional split: faithful band occupies (threshold / maxErr) of the map
nFaithful = max(1, round(nColors * threshold / maxErr));
nGradient = nColors - nFaithful;

% Faithful band: solid saturated green
faithfulCol  = [0.18 0.72 0.33];
cmapFaithful = repmat(faithfulCol, nFaithful, 1);

% Error gradient: yellow -> orange -> red -> dark maroon
anchorErr = [threshold; 0.10; 0.25; maxErr];
anchorRGB = [1.00 0.88 0.15;   % yellow
             1.00 0.55 0.10;   % orange
             0.85 0.15 0.10;   % red
             0.30 0.05 0.05];  % dark maroon

gradErrPts   = linspace(threshold, maxErr, nGradient);
cmapGradient = [interp1(anchorErr, anchorRGB(:,1), gradErrPts(:), "linear"), ...
                interp1(anchorErr, anchorRGB(:,2), gradErrPts(:), "linear"), ...
                interp1(anchorErr, anchorRGB(:,3), gradErrPts(:), "linear")];

cmap = [cmapFaithful; cmapGradient];

end
