%[text] # Getting Started with `fractalnoise`
%[text] A four-function tour: generate, shape, estimate, and diagnose 1/f^alpha noise. Each section is runnable on its own.
%[text] ## Installation
%[text] Add the toolbox folder to your path. That's it -- no compiled components, no install script.
%[text] ```matlabCodeExample
%[text] addpath('/path/to/fractalnoise/toolbox');
%[text] ```
%[text] Everything below runs with nothing else on path. Two of the functions (`shapeNoise`, `estimateIRASA`) have **optional** extra features that need extra MathWorks toolboxes or FieldTrip -- covered at the end of this document -- but none of that is required for normal use.
%[text] ## 1. Generate: `noiseXu`
%[text] Generates 1/f^alpha noise via Xu's (2019) fractional-differencing method. Four arguments: length, spectral exponent alpha, target amplitude (standard deviation), sampling rate.
rng(1);
n = noiseXu(2000, 2.0, 1.5, 60);
fprintf("Generated %d samples, std=%.3f (target 1.5)\n", numel(n), std(n));
%%
%[text] alpha=0 is white noise, alpha=1 is pink, alpha=2 is Brownian/red, alpha=3 is black. Non-integer and negative values are also valid -- see `noiseXu`'s docstring for the full verified range.
%[text] ## 2. Estimate: `estimateIRASA`
%[text] Recovers alpha (and a noise-magnitude estimate) from a time series via IRASA. Takes the signal, sampling rate, and a frequency band to fit over.
[alphaEst, sigmaEst] = estimateIRASA(n, 60, 1.0, 15.0);
fprintf("Recovered alpha=%.3f (target 2.0)\n", alphaEst);
fprintf("Noise magnitude sigma=%.3f\n", sigmaEst);
%%
%[text] Note `sigmaEst` integrates the fractal PSD over the full available frequency range (not restricted to the `[fLow, fHigh]` band used for the alpha fit), so it won't generally equal the generator's amplitude parameter exactly -- it's a spectral-magnitude estimate, not a reconstruction of the original time-domain standard deviation. Use `std()` directly on your signal if you want the latter.
%[text] Ask for a second opinion from FieldTrip (if installed) with `CompareFieldTrip=true`:
[~, ~, ~, ~, ftInfo] = estimateIRASA(n, 60, 1.0, 15.0, CompareFieldTrip=true);
if ftInfo.available
    fprintf("FieldTrip agrees: alpha=%.3f (gap from homebrew: %.3f)\n", ftInfo.alphaFT, ftInfo.gapHBFT);
else
    fprintf("FieldTrip comparison unavailable: %s\n", ftInfo.message);
end
%%
%[text] See `doc/WhyHomebrewIRASA.mlx` for why this function has its own IRASA implementation rather than only calling FieldTrip's.
%[text] ## 3. Shape: `shapeNoise`
%[text] Imposes an **empirical** residual's amplitude spectrum onto a fresh phase carrier -- useful for generating surrogate data that matches a real signal's spectral shape without reusing its exact phase (and therefore its exact realisation).
myResidual = cumsum(randn(500, 1));
myResidual = myResidual - mean(myResidual);

surrogate = shapeNoise(myResidual, 100, 3.0);
fprintf("Surrogate std=%.3f, original std=%.3f (should match)\n", std(surrogate), std(myResidual));
%%
%[text] Five phase-carrier choices are available via `PhaseSource`: `"xu"` (default), `"white"`, `"pinknoise"`, `"dsp"`, or a custom function handle. See `shapeNoise`'s docstring for when each applies.
%[text] ## 4. Diagnose: `checkSpectralKnee`
%[text] Tests whether a residual's spectrum is better explained by a Lorentzian (kneed) model than a plain power law -- relevant because IRASA-based alpha estimates (like `estimateIRASA` above) assume a scale-free, non-kneed spectrum.
result = checkSpectralKnee(myResidual, 100);
fprintf("Verdict: %s (deltaAIC=%.2f, knee at %.3f Hz, in fit band: %d)\n", ...
    result.verdict, result.deltaAIC, result.kneeHz, result.kneeInBand);
%%
%[text] See `examples/KneeCheckWorkedExample.mlx` for this diagnostic applied to real kinematic data, where three of five datasets came back clean and two showed a genuine but harmless below-the-fit-band knee.
%[text] ## Optional toolbox dependencies
%[text] Two functions have optional features that degrade **visibly**, not silently, when the underlying toolbox isn't installed: requesting the feature without the dependency returns a clearly-flagged unavailable result and a warning, never a silent fallback or an outright crash.
%[text] - `shapeNoise(..., PhaseSource="pinknoise")` needs MATLAB's built-in  `pinknoise` function (Signal Processing Toolbox).
%[text] - `shapeNoise(..., PhaseSource="dsp")` needs `dsp.ColoredNoise` (DSP System Toolbox).
%[text] - `estimateIRASA(..., CompareFieldTrip=true)` needs FieldTrip ([https://www.fieldtriptoolbox.org](https://www.fieldtriptoolbox.org)) on the path, specifically  `ft_defaults` having been run. \
%[text] None of the four core functions' default behaviour needs any of these.
%[text] ## Where to go next
%[text] - `doc/WhyHomebrewIRASA.mlx` -- the validation story behind  `estimateIRASA`'s verified alpha range
%[text] - `examples/KneeCheckWorkedExample.mlx` -- `checkSpectralKnee` against real data
%[text] - `examples/CompareNoiseGenerators.mlx` -- a side-by-side look at the different noise-generation and phase-shaping options
%[text] - `tests/` -- the test suite doubles as executable documentation of expected behaviour for every function \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
