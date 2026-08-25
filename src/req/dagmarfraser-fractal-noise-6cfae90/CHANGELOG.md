# Changelog

`fractalnoise` follows semantic versioning. The internal `PowerLawSimulationPreReg`
`_vNNN` suffix convention (e.g. `XuNoise_v002`, `iraAlphaSigma_v003`) does NOT
carry into this toolbox's releases -- that convention exists to track
in-place iteration of a single file across a long-running research project's
session history, which is not this toolbox's situation. Each function's
docstring states which specific internal vintage it was extracted from, so
that history is not lost, just relocated to prose rather than the filename.

## v1.1.0 -- 2026-07-15

### Naming rationale confirmed against a second MathWorks option: `fractalcoef`

The original v1.0.0 naming rationale (below) positioned `fractalnoise`
against `dsp.ColoredNoise`. Checked directly this release whether
`fractalcoef` (Sensor Fusion and Tracking Toolbox) offered a genuinely
different alpha range, since MathWorks documentation itself uses the
term "fractal noise" for this function. It does not: `fractalcoef`'s own
citation is Kasdin, N. J. (1995), *Discrete Simulation of Colored Noise
and Stochastic Processes and 1/f^alpha Power Law Noise Generation*,
Proceedings of the IEEE -- the identical algorithm `dsp.ColoredNoise` is
built on, just returned as raw filter coefficients (`Numerator`/
`Denominator` for `filter()`) rather than wrapped as a System object.
Its alpha range is `(0, 2)`, open -- narrower than `dsp.ColoredNoise`'s
`[-2, 2]`, not wider or different. This strengthens rather than
complicates the original naming case: two independent MathWorks
toolboxes cap out at the same boundary because they implement the same
1995 method, confirming `noiseXu`'s Xu (2019) GGM fractional-
differencing approach (range `[-2, 7]`, numerically stable) fills a
genuine gap rather than duplicating an existing MATLAB option under a
different name. See `README.md`, "Where `noiseXu` sits relative to
MathWorks' own options."

### `ComparePMTM` added to `estimateIRASA.m`

New name-value option, sixth output `pmInfo` matching `ftInfo`'s shape
(`available`, `alphaPM`, `gapHBPM`, `message`). Motivated by a direct
challenge: validating an unsteady noise generator against an estimator's
unverified range is dangerous unless the check is genuinely independent,
not just agreement between two variants of the same mechanism.
`CompareFieldTrip` shares the resampling-median step with the homebrew
estimator, so HB-vs-FT agreement alone cannot rule out both sharing a
blind spot; `ComparePMTM` is a direct log-log pmtm slope fit with no
resampling step at all -- the genuinely independent check. Reuses the
PSD `estimateIRASACore` already computes internally (previously
discarded), so it costs no extra `pmtm` call. Three new tests added,
including a source-equivalence check confirming the underlying PSD is
bit-identical (AbsTol 1e-12) to a direct `pmtm(x,4,[],fs)` call.

### Test suite decoupled from a hardcoded legacy machine path

Every test method except the one genuine source-equivalence test
(`testMatchesSourceForSameSeed`) was calling `generateCustomNoise_v004`
purely to generate a synthetic signal with known alpha -- reaching
across the toolbox's own boundary to a `PowerLawSimulationPreReg`
function hardcoded to one specific machine's path, when the toolbox's
own `noiseXu` (Item 1's extraction of that exact function) already
provides the same thing. Swapped every such call to `noiseXu`. Verified
via a full fresh-path run (`restoredefaultpath`, only the toolbox
directory added): 23 passed, 0 failed, 9 correctly filtered (7 legacy
source-equivalence tests self-filtering as designed, one FieldTrip and
one DSP-System-Toolbox optional-dependency test filtering on the side
that doesn't apply to this particular machine's base install).

### `.mlx` documentation converted to the plain-text Live Code (`.m`) format

Since R2025a, the Live Editor supports a genuinely plain-text alternative
to `.mlx` (custom `%[text]`/`%[output]`/`%[appendix]` comment markup),
solving the git-diffability problem `.mlx`'s binary zip format always
had. All four existing `.mlx` docs/examples (`GettingStarted`,
`WhyHomebrewIRASA`, `CompareNoiseGenerators`, `KneeCheckWorkedExample`)
converted via the Live Editor's own Save As (`export()` does NOT produce
this format -- confirmed empirically, it silently falls back to the old
publish-style `%%` markup instead). Old `.mlx` files removed.

**This raises the toolbox's minimum MATLAB release from `R2021a` to
`R2025a`** -- the plain-text Live Code format did not exist before
R2025a. The four core functions themselves have no new syntax
requirements; the bump is driven entirely by the documentation format.

### New example: `examples/CompareFftnoiseVsShapedXu.m`

Demonstrates the fftnoise-vs-shaped_xu divergence directly: above target
alpha~3, `fftnoise` (phase-randomised surrogate) plateaus near alpha~2.7
regardless of target, while `shapeNoise(...,PhaseSource="xu")` continues
tracking the target. Checked with all three of `estimateIRASA`'s methods
simultaneously (HB/FT/PM), citing Brookshire (2022, *Nature Human
Behaviour*) as the rationale for not trusting a single method's
agreement with itself, especially in a noise regime (alpha 3-5+) well
outside `dsp.ColoredNoise`'s supported range. Authored directly in the
plain-text Live Code format (no existing `.mlx` to convert).

### `functionSignatures.json` updated

`estimateIRASA`'s entry was missing `ComparePMTM` from inputs and
`pmInfo` from outputs after the above change -- caught by inspection,
not automatically. Fixed; validated as syntactically correct JSON.

### Repackaged

`release/Fractal Noise Toolbox.mltbx` rebuilt: version 1.0.0 -> 1.1.0,
`MinimumMatlabRelease` R2021a -> R2025a, `ToolboxGettingStartedGuide`
repointed from `doc/GettingStarted.mlx` to `doc/GettingStarted.m`, same
GUID preserved (`f070d0f9-6828-49e2-9178-1861e71f6fe2`) so this is not
treated as a different toolbox on reinstall. Verified genuinely
functional post-repackage, same discipline as v1.0.0: dev folder removed
from path, installed into a clean environment, all four core functions
smoke-tested from the installed Add-On location (not the dev copy,
confirmed via `which`), `ComparePMTM` specifically exercised and
returning sane output, `installedToolboxes` confirmed correct
name/version/GUID, then uninstalled and dev path explicitly restored.

## v1.0.0 -- 2026-07-02

First release. Four core functions, delivered from `PowerLawSimulationPreReg`
session 52 (named `noisetools` during development; renamed to `fractalnoise`
before release -- see "Naming" below):

- `noiseXu.m` -- extracted from `XuNoise_v002.m`, with
  `generateCustomNoise_v004.m`'s wrapper interface folded in
- `shapeNoise.m` -- generalised from `shapeXu_local` (a private local
  function inside `generateLoopClosureNoise_v003.m`), with the phase
  carrier made pluggable (`PhaseSource`: `"xu"`, `"white"`, `"pinknoise"`,
  `"dsp"`, or a custom function handle)
- `estimateIRASA.m` -- wraps `iraAlphaSigma_v003.m`, with an optional
  FieldTrip cross-validation panel gated behind `CompareFieldTrip`
- `checkSpectralKnee.m` -- generalised from `checkLorentzianKnee_v001.m`'s
  per-trial fitting core into a clean single-trial function

Plus `external/fftnoise/` (Aslak Grinsted, 2011, BSD-2, unmodified), a
full `doc/`/`examples/` set (four `.mlx` files), `functionSignatures.json`
for editor tab-completion, and a 29-test suite across five files
(`tests/test*.m`).

### Naming

Developed under the working name `noisetools`. Renamed to `fractalnoise`
before release: the original name implied a narrower, `dsp.ColoredNoise`-
adjacent scope than the toolbox actually has (generation AND estimation
AND diagnostics of fractal/1-over-f noise, not just an extension of one
MathWorks class that only one of four functions optionally touches).
`fractalnoise` also ties directly to `pFractal`, the term already used
throughout `estimateIRASA`/`checkSpectralKnee`, and matches
`external/fftnoise/`'s lowercase-concatenated naming convention.

### Three real bugs found and fixed

All three were only caught once real dependencies were actually tested
against, rather than just their absence -- every test had passed against
the "tool not installed" degraded-mode paths before these were found:

1. `shapeNoise.m`'s `"dsp"` carrier referenced the class as
   `dsp.colorednoise`; the correct name is `dsp.ColoredNoise` (case-
   sensitive) -- the lowercase name silently doesn't exist, so this
   branch always fell through to the "unavailable" error even with the
   toolbox installed. Also redesigned to use `dsp.ColoredNoise`'s
   continuous `InverseFrequencyPower` property (range [-2,2], clamped
   with a warning outside it) rather than snapping to three fixed noise
   colours.
2. `estimateIRASA.m`'s FieldTrip comparison panel built its
   `ftData`/`ftCfg` structs with double-quoted MATLAB `string` values;
   FieldTrip's internal validation requires char vectors specifically
   and errors on `string` input. Also added `ftData.sampleinfo` to avoid
   two benign-but-noisy reconstruction warnings FieldTrip otherwise
   emits every call.
3. **Most serious:** `estimateIRASA.m` was never actually self-contained.
   It called `iraAlphaSigma_v003` directly -- a `PowerLawSimulationPreReg`
   function documented as the extraction origin but never actually
   copied into the toolbox, directly contradicting the toolbox's stated
   design goal of not depending on that project. Every test had passed
   because the test file's setup always added `PowerLawSimulationPreReg`'s
   legacy functions directory to path alongside the toolbox directory,
   for its own source-equivalence comparisons -- masking the hidden
   dependency. Fixed by inlining the algorithm as a local function.
   `tests/testToolboxSelfContained.m` was added specifically to catch
   this class of bug going forward: it adds ONLY the toolbox directory
   to path and confirms all four public functions still work.

### Citations added

- Xu, C. (2019). An Easy Algorithm to Generate Colored Noise Sequences.
  The Astronomical Journal, 157(3), 127.
- Wen, H., & Liu, Z. (2016). Separating Fractal and Oscillatory
  Components in the Power Spectrum of Neurophysiological Signal. Brain
  Topography, 29(1), 13-26. https://doi.org/10.1007/s10548-015-0448-0
- Preston, M., Smith, S. & Voytek, B. (2026). Potential mechanisms and
  functional significance of aperiodic neural activity. Nature Human
  Behaviour.
- Grinsted, A. (2011). `fftnoise.m` (bundled unmodified, BSD-2).

### `examples/KneeCheckWorkedExample.mlx` strengthened

Initially demonstrated the diagnostic against `noiseXu(..., Phi=0.99)` as
an approximation of the real Hickman-like sub-floor knee -- reasonable,
but only qualitative, and covered just one of the three scenarios
`checkSpectralKnee` distinguishes. Rebuilt to construct synthetic
residuals with an *exact, known* Lorentzian spectrum: the target
magnitude `(kneeHz^2+f^2)^(-alpha/4)` is computed analytically at every
FFT bin, then handed to the bundled `fftnoise.m` to randomise phase
while preserving that exact magnitude. This lets recovered `kneeHz` and
`alpha` be checked against genuine ground truth (not possible with real
data, where the true values are unknown), and covers all three
scenarios: sub-floor knee (matches the real Hickman result), in-band
knee (the actual bias-risk case, not observed in any of the five real
datasets, now quantified: a naive fit straddling a known knee at
alpha=4.5 recovers ~4.05-4.08, a real ~0.4-unit bias), and no knee at
all.

### Packaged as `.mltbx`

Built via `matlab.addons.toolbox.ToolboxOptions` + `packageToolbox`
(fully scripted, no GUI needed) to `release/Fractal Noise Toolbox.mltbx`.
Caught `.DS_Store` sneaking into the auto-detected file list before
packaging -- filtered out, and added to `.gitignore` (built from
MathWorks' own `mathworks/gitignore` MATLAB template plus
`fractalnoise`-specific entries) so it can't happen again silently.
Verified genuinely functional, not just "a file exists": installed into
a clean environment (dev folder removed from path first, so functions
resolved from the installed Add-On location, not the dev copy), ran all
four core functions successfully, confirmed correct name/version/GUID
via `installedToolboxes`, then uninstalled to leave the environment as
found. Installing via `.mltbx` also resolves the `external/fftnoise/`
path gap automatically (the packager adds all subfolders to path on
install) -- unlike the raw dev `addpath(toolbox/)` workflow, which needs
the manual extra `addpath` `KneeCheckWorkedExample.mlx` adds for itself.
