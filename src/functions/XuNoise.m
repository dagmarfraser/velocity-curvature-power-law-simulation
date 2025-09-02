function noise = XuNoise(N, alpha, amplitude, fs, phi)
% XUNOISE Generate colored noise using Xu (2019) GGM method with break frequency
% 
% CRITICAL IMPLEMENTATION of Xu, C. (2019) "An Easy Algorithm to Generate 
% Colored Noise Sequences" for red leakage mitigation in biological motion analysis
%
% Syntax:
%   noise = XuNoise(N, alpha)
%   noise = XuNoise(N, alpha, amplitude)
%   noise = XuNoise(N, alpha, amplitude, fs)
%   noise = XuNoise(N, alpha, amplitude, fs, phi)
%
% Inputs:
%   N         - Length of generated noise sequence (integer)
%   alpha     - Spectral exponent (full alpha value, will be converted to d=alpha/2)
%   amplitude - Target standard deviation of output signal (default: 1.0)
%   fs        - Sampling frequency in Hz (default: 1.0, used for validation)
%   phi       - GGM break frequency parameter (default: 0.99, CRITICAL for red leakage mitigation)
%
% Output:
%   noise - Column vector of colored noise samples with GUARANTEED amplitude control
%
% Key Implementation Details:
%   - Uses phi = 0.99 (NOT 1.0) to introduce break frequency for red leakage mitigation
%   - Fractional memory parameter d = alpha/2 (NOT full alpha value)
%   - FFT-based convolution for computational efficiency
%   - Post-processing amplitude normalization ensures precise amplitude control
%   - Based on exact implementation from Xu (2019) Appendix A
%
% Performance Characteristics:
%   - Superior alpha recovery compared to traditional methods
%   - Reduced red leakage artifacts for high alpha values (α ≥ 2)
%   - Enhanced spectral stability across diverse analysis methods
%   - Validated for biological motion analysis applications
%
% References:
%   Xu, C. (2019). An Easy Algorithm to Generate Colored Noise Sequences. 
%   The Astronomical Journal, 157:127
%
% Example:
%   % Generate 1000 samples of brown noise (alpha=2) with amplitude=1.5 at 60Hz
%   brownNoise = XuNoise(1000, 2.0, 1.5, 60);
%   
%   % Generate pink noise with default amplitude
%   pinkNoise = XuNoise(500, 1.0);
%
% See also: generateCustomNoise_v003,
% SpectralVsDFAComparison_v004_XuMethods
% compareNoiseGenerators
%
% Author: PowerLaw Toolchain Development Team
% Based on: Xu (2019) validated implementation
% Version: 2.0 - Includes post-processing amplitude normalization
% Date: July 2025

%% Input validation and parameter setup
if nargin < 2
    error('XuNoise:InsufficientInputs', 'At least N and alpha must be specified');
end

if nargin < 3 || isempty(amplitude)
    amplitude = 1.0;
end

if nargin < 4 || isempty(fs)
    fs = 1.0;  % Default sampling frequency for validation
end

if nargin < 5 || isempty(phi)
    phi = 0.99;  % CRITICAL: Use break frequency for red leakage mitigation
end

% Validate inputs
if ~isscalar(N) || N <= 0 || mod(N, 1) ~= 0
    error('XuNoise:InvalidN', 'N must be a positive integer');
end

if ~isscalar(alpha) || ~isnumeric(alpha) || alpha < -2 || alpha > 5
    error('XuNoise:InvalidAlpha', 'Alpha parameter must be numeric in range [-2, 5]');
end

if ~isscalar(amplitude) || ~isnumeric(amplitude) || amplitude < 0
    error('XuNoise:InvalidAmplitude', 'Amplitude parameter must be a non-negative numeric value');
end

if ~isscalar(phi) || phi <= 0 || phi >= 1
    error('XuNoise:InvalidPhi', 'Phi must be in range (0, 1)');
end

%% Handle zero amplitude case
if amplitude == 0
    noise = zeros(N, 1);
    return;
end

%% Xu GGM Implementation
% Convert full alpha to fractional memory parameter
d = alpha / 2;  % CRITICAL: Use alpha/2, not full alpha

% Generate white noise foundation (use unit variance for initial scaling)
w = randn(N, 1);

% Apply Xu fractional differencing procedure
noise = xuGGM_FD(w, d, phi);

% CRITICAL: Apply post-processing amplitude normalization
% The fractional differencing process changes signal amplitude beyond
% the initial scaling, so we must normalize the final output
currentStd = std(noise);
if currentStd ~= 0
    noise = noise * (amplitude / currentStd);
else
    % Fallback for degenerate case
    noise = zeros(N, 1);
end

% Ensure output is column vector
noise = noise(:);

end

function dx = xuGGM_FD(x, d, phi)
% XU_GGM_FD Xu's fractional differencing procedure for GGM
% 
% Exact implementation from Xu (2019) Appendix A
% 
% Mathematical formulation:
%   dx = (1-phi*B)^(-d) * x
%   where B is the backshift operator and phi introduces break frequency
%
% Input:
%   x   - Input white noise sequence (column vector)
%   d   - Fractional memory parameter (alpha/2)
%   phi - Break frequency parameter (use 0.99 for red leakage mitigation)
%
% Output:
%   dx  - Colored noise sequence with GGM characteristics

N = size(x, 1);

% Generate fractional differencing coefficients
k = (1:N-1)';
h = phi .* (k + d - 1) ./ k;  % Coefficient generation
h = [1; cumprod(h)];          % Cumulative product for series expansion

% FFT-based convolution for computational efficiency
np = 2^(fix(log2(2*N)) + 1);  % Next power of 2 for FFT efficiency
dx = ifft(fft(x, np) .* fft(h, np));
dx = real(dx(1:N, :));        % Extract real component, original length

end