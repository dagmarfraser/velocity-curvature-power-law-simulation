function noise = noiseXu(N, alpha, amplitude, fs, options)
%NOISEXU Generate 1/f^alpha colored noise via Xu (2019) GGM fractional differencing.
%
% Toolbox extraction of PowerLawSimulationPreReg's XuNoise_v002 (core
% generator) with generateCustomNoise_v004's wrapper folded into this
% signature. There is no separate wrapper function in the toolbox; Phi
% is exposed as an optional name-value argument instead of a positional
% one, so the common 4-argument call is unchanged from the old wrapper.
%
% SYNTAX:
%   noise = noiseXu(N, alpha)
%   noise = noiseXu(N, alpha, amplitude)
%   noise = noiseXu(N, alpha, amplitude, fs)
%   noise = noiseXu(N, alpha, amplitude, fs, Phi=phiValue)
%
% INPUTS:
%   N         - Length of generated noise sequence (positive integer)
%   alpha     - Spectral exponent (full alpha; converted internally to
%               fractional memory parameter d = alpha/2)
%   amplitude - Target standard deviation of output signal (>=0; default 1.0)
%   fs        - Sampling frequency in Hz, retained for interface
%               compatibility with generateCustomNoise_v004; not used in
%               the generation maths itself (default 1.0)
%   Phi       - Name-value. GGM break frequency parameter, range (0,1)
%               exclusive. Default 0.99 (red-leakage mitigation per Xu
%               2019 Appendix A). Only override this if you know why.
%
% OUTPUT:
%   noise     - Column vector, length N, with std(noise) numerically
%               equal to amplitude (post-hoc normalised; see Notes)
%
% VERIFIED RANGE (Finding #100, PowerLawSimulationPreReg):
%   This generator is numerically stable (no NaN/Inf output) for alpha
%   in the range [-2, 7]. That is NOT the same claim as "alpha=7 is
%   independently measurable" -- three-way estimator agreement
%   (homebrew IRASA, FieldTrip IRASA, direct pmtm slope fit) verified
%   INDEPENDENT RECOVERABILITY only up to alpha=6. At alpha=7 the
%   three estimators diverge from one another, consistent with an
%   information-limit at typical trial lengths rather than a generator
%   defect. Treat alpha=7 output as usable for stress-testing pipeline
%   robustness, not as ground truth for alpha-recovery validation.
%
% NOTES:
%   - No alpha-range validation is enforced here by design: this is a
%     stress-test-capable generator (matching XuNoise_v002 behaviour).
%     Numerical breakdown (NaN/Inf) is the empirical failure signal and
%     is the caller's responsibility to catch.
%   - amplitude=0 returns an exact zero vector, not a degenerate
%     near-zero-variance draw.
%   - fs is accepted but unused in the generation maths; it exists
%     purely so callers migrating from generateCustomNoise_v004 do not
%     need to change their call sites.
%
% REFERENCES:
%   Xu, C. (2019). An Easy Algorithm to Generate Colored Noise Sequences.
%   The Astronomical Journal, 157:127.
%
% EXAMPLE:
%   brownNoise = noiseXu(1000, 2.0, 1.5, 60);
%   pinkNoise  = noiseXu(500, 1.0);
%   customPhi  = noiseXu(1000, 3.0, 1.0, 100, Phi=0.995);
%
% See also: shapeNoise, estimateIRASA
%
% Toolbox: fractalnoise
% Extracted from: src/functions/XuNoise_v002.m + generateCustomNoise_v004.m
% (PowerLawSimulationPreReg, session 51/52)

arguments
    N (1, 1) double {mustBeInteger, mustBePositive}
    alpha (1, 1) double
    amplitude (1, 1) double {mustBeNonnegative} = 1.0
    fs (1, 1) double {mustBePositive} = 1.0 %#ok<INUSA> intentionally unused; kept for call-site compatibility with generateCustomNoise_v004
    options.Phi (1, 1) double {mustBeGreaterThan(options.Phi, 0), mustBeLessThan(options.Phi, 1)} = 0.99
end

phi = options.Phi;

%% Zero-amplitude short circuit
if amplitude == 0
    noise = zeros(N, 1);
    return
end

%% Xu GGM generation
% Full alpha to fractional memory parameter (Xu 2019 convention: d = alpha/2)
d = alpha / 2;

whiteSeed = randn(N, 1);
noise = xuGgmFracDiff(whiteSeed, d, phi);

% Post-hoc amplitude normalisation: fractional differencing changes
% signal amplitude beyond the initial unit-variance seed, so the final
% std must be corrected explicitly rather than trusted from the seed.
currentStd = std(noise);
if currentStd ~= 0
    noise = noise * (amplitude / currentStd);
else
    noise = zeros(N, 1);
end

noise = noise(:);

end

function dx = xuGgmFracDiff(x, d, phi)
%XUGGMFRACDIFF Xu's fractional differencing procedure for GGM noise.
%
% Exact implementation from Xu (2019) Appendix A:
%   dx = (1 - phi*B)^(-d) * x
% where B is the backshift operator and phi introduces a break
% frequency relative to the unmodified (phi=1) fractional differencing
% series.
%
% INPUTS:
%   x   - Input white noise sequence, column vector
%   d   - Fractional memory parameter (alpha/2)
%   phi - Break frequency parameter, (0,1) exclusive
%
% OUTPUT:
%   dx  - Colored noise sequence with GGM spectral characteristics

N = size(x, 1);

k = (1:N - 1)';
h = phi .* (k + d - 1) ./ k;
h = [1; cumprod(h)];

np = 2 ^ (fix(log2(2 * N)) + 1);
dx = ifft(fft(x, np) .* fft(h, np));
dx = real(dx(1:N, :));

end
