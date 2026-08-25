# fractalnoise

A small, standalone MATLAB toolbox for generating, shaping, estimating,
and diagnosing 1/f^alpha ("fractal" or "colored") noise. Extracted and
generalised from internal functions of
[velocity-curvature-power-law-simulation](https://github.com/dagmarfraser/velocity-curvature-power-law-simulation)
(University of Birmingham, School of Psychology), where these tools were
built to characterise biological and instrumental noise in kinematic
drawing-task data. `fractalnoise` itself is not about power laws,
kinematics, or drawing tasks -- it is the noise-handling machinery that
project needed, offered here as a general-purpose, dependency-light
utility.

**Status:** v1.1.0 (2026-07-15). Packaged as `.mltbx` -- see
`release/Fractal Noise Toolbox.mltbx`. See [CHANGELOG.md](CHANGELOG.md).
Source: [github.com/dagmarfraser/fractal-noise](https://github.com/dagmarfraser/fractal-noise).

## Requirements

**MATLAB R2025a or later.** This is driven entirely by the documentation
format: `doc/` and `examples/` use the plain-text Live Code (`.m`)
format, which did not exist before R2025a. The four core functions
(`noiseXu`, `shapeNoise`, `estimateIRASA`, `checkSpectralKnee`) have no
R2025a-specific syntax themselves -- if you only need the functions and
not the rendered documentation, copying the `.m` files directly (rather
than installing the packaged `.mltbx`, which enforces the declared
minimum release) may work on considerably older MATLAB releases, but
this is untested and not officially supported.

Two optional features need additional MathWorks toolboxes: `shapeNoise`'s
`"dsp"` phase carrier needs DSP System Toolbox; `estimateIRASA`'s
`CompareFieldTrip` option needs [FieldTrip](https://www.fieldtriptoolbox.org)
on the path separately (not a MathWorks toolbox). Neither is required for
normal use -- see the degradation note below.

---

## What's in it

| Function | Purpose | Toolbox dependencies |
|---|---|---|
| `noiseXu.m` | Generate 1/f^alpha noise via Xu (2019) GGM fractional differencing | None |
| `shapeNoise.m` | Impose an empirical amplitude spectrum onto a chosen phase carrier | None for default carrier; optional Signal Processing / DSP System Toolbox for two of five carrier choices |
| `estimateIRASA.m` | IRASA spectral-exponent and noise-magnitude estimator | **Signal Processing Toolbox required** (`pmtm` is called directly in the core estimate, not gated behind an optional-feature check); optional FieldTrip for a `CompareFieldTrip` cross-validation panel; `ComparePMTM` needs no toolbox beyond Signal Processing Toolbox (reuses the same PSD the core estimate already computes) |
| `checkSpectralKnee.m` | Lorentzian-vs-fixed spectral knee diagnostic (Preston, Smith & Voytek, 2026) | **Signal Processing Toolbox required** (built on `estimateIRASA`, inherits its `pmtm` dependency) |

**Signal Processing Toolbox is a genuine hard dependency for `estimateIRASA.m`
and `checkSpectralKnee.m`** (both call `pmtm` directly and unconditionally) --
this is a real requirement, not something that degrades gracefully if absent.
`noiseXu.m` and `shapeNoise.m`'s default path need no toolbox beyond base MATLAB.
Beyond that hard dependency, every OPTIONAL feature (FieldTrip, `dsp.ColoredNoise`,
built-in `pinknoise`) degrades visibly, not silently: requesting an optional
feature without the underlying toolbox installed returns a clearly-flagged
unavailable result with a warning, rather than erroring outright or pretending
to succeed.

## Where `noiseXu` sits relative to MathWorks' own options

MATLAB ships several ways to generate noise, and nearly all of them are
fixed at a single alpha or capped at a narrow range:

| | Method | Alpha range | Requires |
|---|---|---|---|
| `randn` | White noise -- no alpha parameter at all | `0`, fixed | Base MATLAB |
| `wgn` / `awgn` | White Gaussian noise, power specified in dBW/dBm/linear -- `randn` under the hood with a power convention layered on top | `0`, fixed | Communications Toolbox |
| `pinknoise` | Built-in pink-noise generator | `1`, fixed | Signal Processing Toolbox |
| `dsp.ColoredNoise` | Kasdin (1995) IIR/FIR filter-coefficient recursion | `[-2, 2]`, closed | DSP System Toolbox |
| `fractalcoef` | Same Kasdin (1995) method -- hands you the raw filter coefficients instead of a System object | `(0, 2)`, open (no negative/anti-persistent alpha at all) | Sensor Fusion and Tracking Toolbox |
| `noiseXu` (this toolbox) | Xu (2019) GGM fractional differencing -- a different method, not a repackaging of Kasdin | `[-2, 7]`, numerically stable (see the verified-range caveat below for what "stable" does and doesn't mean) | None |

Two independent MathWorks toolboxes (`dsp.ColoredNoise`, `fractalcoef`)
cap out at the same boundary because they implement the same 1995
algorithm, not because alpha=2 is some fundamental limit on 1/f^alpha
noise generation -- and everything else in the table above is fixed at a
single alpha with no range at all. `noiseXu` exists specifically because
this project's own empirical noise (Cook/Hickman kinematic residuals,
alpha=3.2-5.4) sits well past where any built-in MATLAB option can go.

## Installation

**Option A -- MATLAB Toolbox file (recommended):** double-click
`release/Fractal Noise Toolbox.mltbx`, or install it via
`matlab.addons.toolbox.installToolbox("release/Fractal Noise Toolbox.mltbx")`.
This also handles `external/fftnoise/` automatically (added to path along
with everything else), which plain `addpath` below does not.

**Option B -- source folder:**
```matlab
addpath('/path/to/fractalnoise/toolbox');
```
If you need `external/fftnoise/` too (only `examples/KneeCheckWorkedExample.m`
calls it directly), add that explicitly:
```matlab
addpath('/path/to/fractalnoise/toolbox/external/fftnoise');
```

## Quick example

```matlab
% Generate 2000 samples of alpha=2 (Brownian) noise at 60 Hz, amplitude 1.5mm
n = noiseXu(2000, 2.0, 1.5, 60);

% Estimate it back out
[alphaEst, sigmaEst] = estimateIRASA(n, 60, 1.0, 15.0);

% Impose an empirical residual's amplitude spectrum onto a fresh Xu carrier
shaped = shapeNoise(myResidual, 100, 3.5);

% Check whether a residual's spectrum has a Lorentzian knee (Preston et al. 2026)
result = checkSpectralKnee(myResidual, 100);
fprintf('%s knee support (deltaAIC=%.2f)\n', result.verdict, result.deltaAIC);
```

See `doc/WhyHomebrewIRASA.m` for why `estimateIRASA` exists as a
homebrew implementation rather than relying solely on FieldTrip, and
`examples/` for worked examples against real data, including
`examples/CompareFftnoiseVsShapedXu.m` for a direct demonstration of
why `fftnoise` is the wrong surrogate above alpha~3 and `shapeNoise`
is the fix.

## A verified-range caveat, stated once here and repeated in every relevant docstring

`noiseXu`'s generator is numerically stable (no NaN/Inf output) for alpha
in **[-2, 7]**. That is a *narrower* claim than whether `estimateIRASA`
can independently *recover* a given alpha back out, which is verified
only for alpha in **[-2, 6]** (three mechanistically distinct estimators
-- homebrew IRASA, FieldTrip IRASA (`CompareFieldTrip`), and a
resampling-free direct pmtm slope fit (`ComparePMTM`) -- diverge from one
another at alpha=7, consistent with an information limit at typical
trial lengths rather than a defect in any one estimator). Do not conflate
generator stability with estimator reliability; see `estimateIRASA.m`'s
docstring for the full account, and `examples/CompareFftnoiseVsShapedXu.m`
for a worked demonstration of all three estimators disagreeing with a
broken generator (`fftnoise`) in the same way, and agreeing with a
correct one (`shapeNoise`) throughout.

## Testing

```matlab
results = runtests('/path/to/fractalnoise/tests');
```

32 tests across five files, covering source-equivalence (each function's
original extraction origin checked against the toolbox version, same RNG
seed, where a genuine legacy source exists to compare against),
independent-estimator cross-checks (`ComparePMTM`, `CompareFieldTrip`),
edge cases, error-path guards, and self-containment (the toolbox works
with only its own folder on path -- see `tests/testToolboxSelfContained.m`,
added after this exact class of bug was caught during development; see
`CHANGELOG.md`). Verified on a genuinely fresh path
(`restoredefaultpath` + only the toolbox directory added, no legacy
project code, no FieldTrip): 23 passed, 0 failed, 9 correctly filtered
(legacy source-equivalence tests self-filter with no legacy code
present, by design; exactly one side of each optional-dependency pair
filters depending on what ships as part of the base MATLAB install on
the test machine).

## Credits

- `noiseXu.m` implements the fractional-differencing method of Xu, C.
  (2019), *An Easy Algorithm to Generate Colored Noise Sequences*, The
  Astronomical Journal, 157(3), 127 -- see "Where `noiseXu` sits relative
  to MathWorks' own options" above for how this differs from, and
  extends past, `dsp.ColoredNoise` and `fractalcoef`, both built on
  Kasdin, N. J. (1995), *Discrete Simulation of Colored Noise and
  Stochastic Processes and 1/f^alpha Power Law Noise Generation*,
  Proceedings of the IEEE.
- `estimateIRASA.m` implements the resampling-median IRASA algorithm of
  Wen, H., & Liu, Z. (2016), *Separating Fractal and Oscillatory
  Components in the Power Spectrum of Neurophysiological Signal*, Brain
  Topography, 29(1), 13-26. https://doi.org/10.1007/s10548-015-0448-0
- `checkSpectralKnee.m` implements a diagnostic motivated by Preston, M.,
  Smith, S. & Voytek, B. (2026), *Potential mechanisms and functional
  significance of aperiodic neural activity*, Nature Human Behaviour.
- `external/fftnoise/fftnoise.m` is bundled unmodified, (c) Aslak
  Grinsted 2011, BSD 2-clause license (see
  `external/fftnoise/license.txt`). It is a phase-surrogate generator
  (randomises an existing FFT's phase, preserving magnitude) -- a
  different tool from `noiseXu`'s parametric generator. Offered as an
  alternative phase source for `shapeNoise`'s `PhaseSource` argument
  (via a custom function handle), and used directly in
  `examples/KneeCheckWorkedExample.m` to build synthetic Lorentzian-
  knee test signals with an exact known ground truth, by handing it an
  analytically constructed target magnitude spectrum instead of an
  empirical one.

## Origin

Every function here was extracted from a specific source file or script
in [velocity-curvature-power-law-simulation](https://github.com/dagmarfraser/velocity-curvature-power-law-simulation);
each docstring states its extraction origin and what, if anything,
changed in the generalisation. `fractalnoise` has no dependency back on
that project -- power-law-specific logic (regression, curvature, VGF)
was deliberately left behind. (Developed under the working name
`noisetools`; see `CHANGELOG.md` for the rename rationale.)

## License

MIT. See `license.txt`. Except `external/fftnoise/`, which retains its
own BSD 2-clause license and copyright notice (Aslak Grinsted, 2011).
