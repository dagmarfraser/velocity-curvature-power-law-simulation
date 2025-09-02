function paramSpace = defineParameterSpace(debug)
% DEFINEPARAMETERSPACE Creates parameter space configuration for power law analysis
%
% This function defines the parameter space for multiverse analysis of
% velocity-curvature power law calculations based on the debug level.
%
% Input:
%   debug - Debug level (0: maximal run, 1: minimal debug, 2: rebuild shapes, 3: parallel debug)

% delete all the shapes x sampling freq files before running debug 2
%
% Output:
%   paramSpace - Structure containing all parameter space definitions
%
% Created May 2025
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk

% Initialize the parameter space structure
paramSpace = struct();

% ===== SAMPLING RATES =====
% Sampling rates based on common digitizers (Fraser et al., 2025)
if debug == 0
    % Full run
    paramSpace.samplingRates = [60, 120, 240]; % Hz
elseif debug == 1
    % Minimal debug run
    paramSpace.samplingRates = [120]; % Hz
elseif debug == 3
    % Parallel debug run (medium size for testing parallel performance)
    paramSpace.samplingRates = [120, 240]; % Hz
else % redo the shapes at at frequencies
    % Comprehensive debug run
    paramSpace.samplingRates = [60, 120, 240]; % Hz
end

% ===== SHAPES =====
% Shape choices based on angular frequency φ (Huh & Sejnowski, 2015)
% Shape indices map to {'0' '2/33' '2/5' '4/5' '4/3' '2' '3' '4' '6'}
% where index 6 = angular frequency '2' (ellipse)
if debug == 0
    % Full run - Range around the ellipse (φ = 2)
    % paramSpace.shapes = 3:7; 
    paramSpace.shapes = [6]; % Just ellipse (φ = 2)
elseif debug == 1
    % Minimal debug run
    paramSpace.shapes = [6]; % Just ellipse (φ = 2)
elseif debug == 3
    % Parallel debug run
    paramSpace.shapes = [6]; % Ellipse (φ = 2) and φ = 3
else
    % shapes rebuild
    paramSpace.shapes = [6]; % Ellipse (φ = 2) and φ = 3
end

% ===== POWER LAW EXPONENTS =====
% Enforced beta values in trajectory creation
if debug == 0
    % Full run - range of values centered on 1/3
    paramSpace.generatedBetas = 0:(2/3)/20:2/3;
elseif debug == 1
    % Minimal debug run
    paramSpace.generatedBetas = [1/3]; % Classic one-third only
elseif debug == 3
    % Parallel debug run
    paramSpace.generatedBetas = [0, 1/3, 2/3]; % Full range but fewer points
else
    % Comprehensive debug run
    paramSpace.generatedBetas = [1/3]; % Classic and two deviations
end

% ===== VELOCITY GAIN FACTORS =====
% VGF values to give appropriate tempos
if debug == 0
    % Full run - ~10s (0.5Hz) to ~2.5s (2Hz) per orbit for a standard
    % ellipse 
    paramSpace.vgfValues = exp(4.5:0.1:5.8);
elseif debug == 1
    % Minimal debug run
    paramSpace.vgfValues = [exp(5.2)]; % ~5s tempo per orbit (1Hz)
elseif debug == 3
    % Parallel debug run
    paramSpace.vgfValues = [exp(5.0), exp(5.2), exp(5.5)]; % Three different tempos
else
    % Comprehensive debug run
    paramSpace.vgfValues = [exp(5.2)]; % ~1Hz and ~0.8Hz tempos
end

% ===== NOISE TYPES =====
% Noise spectral types (noise color) as in v023
% -1: blue, 0: white, 1: pink, 2: red/brown, 3: black
if debug == 0
    % Full run
    paramSpace.noiseTypes = 0:0.1:3; % White, pink, brown/red, black
elseif debug == 1
    % Minimal debug run
    paramSpace.noiseTypes = [0]; % White only
elseif debug == 3
    % Parallel debug run
    paramSpace.noiseTypes = [0, 2]; % White and brown/red
else
    % Comprehensive debug run
    paramSpace.noiseTypes = [1]; % White, pink, brown/red
end

% ===== NOISE MAGNITUDE =====
% Noise magnitude (standard deviation in mm)
if debug == 0
    % Full run
    paramSpace.noiseMagnitudes = [0:0.025:0.1 0.25:0.25:2.25 4  6  8  10]; % The regressions plateau at < 2mm
elseif debug == 1
    % Minimal debug run
    paramSpace.noiseMagnitudes = [0, 1]; % No noise and 1mm
elseif debug == 3
    % Parallel debug run
    paramSpace.noiseMagnitudes = [0, 1, 2]; % No, medium, high
else
    % Comprehensive debug run
    paramSpace.noiseMagnitudes = [0]; % No, low, medium, high
end

% ===== FILTER TYPES =====
% Filter types for differentiation
if debug == 0
    % Full run
    %paramSpace.filterTypes = [1, 2, 3, 4];
    paramSpace.filterTypes = [2, 6]; % Butterworth and SG as per the paper.
    % 6 is the revised Sg which scales the width from the 17 for 100Hz
elseif debug == 1
    % Minimal debug run
    paramSpace.filterTypes = [4]; % SG only
elseif debug == 3
    % Parallel debug run
    paramSpace.filterTypes = [2, 6]; % Butterworth and SG-FS
else
    % Comprehensive debug run
    paramSpace.filterTypes = [4]; %  SG
end

% Filter parameters
% Store as cell array of structs for flexibility
paramSpace.filterParams = {
    struct(), % No params for simple diff
    struct('order', 2, 'cutoff', 10, 'zerolag', 1), % Standard BW params
    struct('order', 2, 'cutoff', 10, 'zerolag', 1), % Same for DFBW
    struct('order', 4, 'width', 17) % Standard SG params
};

% ===== REGRESSION TYPES =====
% Regression types for beta calculation
if debug == 0
    % Full run
    %paramSpace.regressTypes = 1:5; % All regression types
    paramSpace.regressTypes = 3:5; % fitlm, LMLS, IRLS
elseif debug == 1
    % Minimal debug run
    paramSpace.regressTypes = [5]; % IRLS only
elseif debug == 3
    % Parallel debug run
    paramSpace.regressTypes = [3, 5]; % fitlm and IRLS
else
    % regen shapes x fs
    paramSpace.regressTypes = 5; % fitlm, LMLS, IRLS
end

% ===== TRIALS =====
% Number of trials for repeatability
if debug == 0
    paramSpace.repeatTrial = 5;
elseif debug == 1
    paramSpace.repeatTrial = 2;
elseif debug == 3
    paramSpace.repeatTrial = 2;
else % regen shapes * trials
    paramSpace.repeatTrial = 1;
end

% Generate trials array
paramSpace.trialNums = 1:paramSpace.repeatTrial;

% Calculate and attach total configuration count for info purposes
paramSpace.totalConfigCount = length(paramSpace.samplingRates) * ...
                             length(paramSpace.shapes) * ...
                             length(paramSpace.generatedBetas) * ...
                             length(paramSpace.vgfValues) * ...
                             length(paramSpace.noiseTypes) * ...
                             length(paramSpace.noiseMagnitudes) * ...
                             length(paramSpace.filterTypes) * ...
                             length(paramSpace.regressTypes) * ...
                             paramSpace.repeatTrial

% Attach debug level to the parameter space for reference
paramSpace.debugLevel = debug;

end
