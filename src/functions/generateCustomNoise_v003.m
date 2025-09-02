function noise = generateCustomNoise_v003(N, alpha, amplitude, fs)
% GENERATECUSTOMNOISE_V003 Clean interface wrapper for XuNoise generation
%
% SYNTAX:
%   noise = generateCustomNoise_v003(N, alpha, amplitude, fs)
%
% INPUTS:
%   N         - Signal length (number of samples)
%   alpha     - Spectral exponent (-2 to 5, where 0 = white noise)
%   amplitude - Target standard deviation of output signal (>=0, zero allowed)
%   fs        - Sampling frequency in Hz
%
% OUTPUT:
%   noise     - Generated colored noise signal (column vector)
%                GUARANTEED to have std(noise) ≈ amplitude
%
% IMPLEMENTATION:
%   This function serves as a clean interface wrapper around XuNoise,
%   which handles all noise generation and amplitude control internally.
%   The amplitude normalization logic has been moved into XuNoise for
%   improved architecture and code reusability.
%
% ARCHITECTURE NOTES:
%   - XuNoise: Core noise generation with post-processing amplitude control
%   - generateCustomNoise_v003: Clean interface wrapper for compatibility
%   - All validation and amplitude normalization handled by XuNoise
%
% See also: XuNoise, SpectralVsDFAComparison_v004_XuMethods
%
% Author: PowerLaw Toolchain Development Team
% Based on: XuNoise wrapper architecture
% Version: 3.1 - Clean wrapper implementation
% Date: July 2025

% Input validation
if nargin < 4
    error('generateCustomNoise_v003:InsufficientArguments', ...
          'Function requires 4 input arguments: N, alpha, amplitude, fs');
end

% Generate noise using XuNoise (handles all validation and amplitude control)
try
    noise = XuNoise(N, alpha, amplitude, fs, 0.99);
catch ME
    % Re-throw with function-specific error context
    error('generateCustomNoise_v003:XuNoiseFailure', ...
          'Noise generation failed: %s', ME.message);
end

end