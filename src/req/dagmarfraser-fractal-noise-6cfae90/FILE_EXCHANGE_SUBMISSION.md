# File Exchange Submission Content (draft)

Status: draft, not yet submitted. GitHub repo now exists:
https://github.com/dagmarfraser/fractal-noise -- see TODO Item 9. Next
decision: connect live vs cut a Release first and link that (below),
then complete the submission form.

## Route decision: GitHub-connected, not direct upload

**File Exchange only supports BSD license for direct "Upload Files"
submissions** -- other licenses (including this toolbox's MIT) are only
permitted via "Connect to GitHub" or "Link to an External Website".
Direct upload would silently force a BSD relicense, or get flagged/
rejected against our stated `license.txt`. This makes the GitHub-first
route (already `mathworks/toolboxdesign`'s top recommendation for other
reasons -- version control over what ships) effectively required here,
not just preferred.

**Practical implication, now actionable:** the repo exists
(`dagmarfraser/fractal-noise`) -- either connect it directly (latest
code always used) or create a GitHub Release first and link that (lets
us pin exactly what ships -- recommended, matches
`mathworks/toolboxdesign`'s guidance). Decide this before filling in the
submission form below, since it's a field in that form, not a separate
step afterward.

## Auto-populated from the `.mltbx` (already done via `ToolboxOptions`)

File Exchange extracts Name, Version, Summary, Description, and Author
info directly from the packaged toolbox's project metadata. These were
already set when building `release/Fractal Noise Toolbox.mltbx`:

- **Name:** Fractal Noise Toolbox
- **Version:** 1.1.0
- **Summary:** Generate, shape, estimate, and diagnose 1/f^alpha
  (fractal) noise in MATLAB.
- **Author:** Dagmar Scott Fraser, University of Birmingham,
  d.s.fraser@bham.ac.uk

**Minimum MATLAB release: R2025a**, declared in the packaged `.mltbx`.
This gates File Exchange install eligibility on older releases, and is
worth stating explicitly in the submission text (not just relying on the
auto-populated metadata) so users on pre-R2025a MATLAB understand why
before attempting to install and hitting a version-check failure. Driven
by the `doc/`/`examples/` plain-text Live Code (`.m`) format, not by the
four core functions themselves.

If any of these need editing for a public File Exchange audience
specifically (as opposed to the internal packaging metadata), edit them
in the File Exchange submission form directly -- doing so triggers a
re-package using the updated info (confirmed MathWorks behaviour), it
does not require touching the `.mltbx` again.

## Description (expand/replace the auto-populated summary text)

A small, standalone MATLAB toolbox for generating, shaping, estimating,
and diagnosing 1/f^alpha ("fractal" or "colored") noise -- the kind of
noise that follows a power-law relationship between spectral power and
frequency, common in biological, neurophysiological, and instrumental
signals.

**Four functions:**

- `noiseXu` -- generate 1/f^alpha noise via Xu's (2019) fractional-
  differencing method (an exact spectral-synthesis method, not a filter
  approximation)
- `shapeNoise` -- impose an empirical residual's amplitude spectrum onto
  a chosen phase carrier (5 carrier options, including MATLAB's built-in
  `dsp.ColoredNoise` and a custom-function-handle option)
- `estimateIRASA` -- recover the spectral exponent (alpha) and noise
  magnitude from a signal via IRASA (Wen & Liu, 2016), with two optional
  independent cross-validation checks: a FieldTrip comparison panel, and
  a resampling-free direct pmtm slope fit (no extra toolbox needed for
  the latter -- the genuinely independent check between the two, since
  FieldTrip's IRASA shares a resampling-median mechanism with the
  homebrew estimator)
- `checkSpectralKnee` -- test whether a signal's spectrum is better
  explained by a Lorentzian (kneed) model than a plain power law
  (Preston, Smith & Voytek, 2026) -- relevant because IRASA-based alpha
  estimates assume a scale-free spectrum, and a knee inside your fitting
  band biases the result

**No hard dependencies.** Every function that has an optional feature
requiring another toolbox (MATLAB's Signal Processing Toolbox, DSP
System Toolbox, or FieldTrip) degrades visibly, not silently: requesting
the feature without the dependency returns a clearly-flagged unavailable
result with a warning, never a silent fallback.

**Verified, not just documented.** 32-test suite covering source-
equivalence, independent-estimator cross-checks, edge cases, error
paths, and self-containment (confirmed on a genuinely fresh MATLAB path
via `restoredefaultpath`, not just a dedicated single test -- no hidden
dependencies anywhere in the suite). See `doc/GettingStarted.m` to get
started, `doc/WhyHomebrewIRASA.m` for the validation story behind the
alpha-recovery range this toolbox states explicitly, and
`examples/CompareFftnoiseVsShapedXu.m` for a direct demonstration of why
naive phase-randomised surrogates (`fftnoise`) fail above alpha~3 and
what to use instead: alpha in [-2, 6] verified for estimation, [-2, 7]
for generation stability -- these are different claims and the toolbox
does not conflate them.

## Tags

colored noise, fractal noise, "1/f noise", "power law noise", IRASA,
"pink noise", "spectral analysis", "spectral exponent", "aperiodic
activity", "signal processing", "noise generation", neurophysiology

(File Exchange convention: comma-separated, quotation marks around
multiword tags.)

## Category

Not directly selectable -- File Exchange auto-assigns via text analytics
on tags/title/word-density. Likely landing zone: Signal Processing >
Spectral Analysis, or similar. If it lands wrong, correct via File
Exchange support (confirmed possible, per MathWorks Answers threads) --
not worth trying to game via tag ordering.

## Toolbox image

Not yet created. Optional but recommended (shows in FEX Overview page
and the Add-On Explorer). MathWorks' stated required aspect ratio:
3/4 (height/width). Would live in `images/toolboxIcon.jpg` per
`mathworks/toolboxdesign`'s convention, named to match the root folder.
Not blocking submission -- can be added later without a version bump if
edited directly on the File Exchange Overview page (confirmed possible
for non-GitHub-derived fields), though for a GitHub-connected submission
specifically, verify this still holds (the edit-without-version-bump
behaviour was confirmed for direct/other submissions, not explicitly
tested here for the GitHub-Release-linked case).

## Acknowledgements

- Xu, C. (2019). An Easy Algorithm to Generate Colored Noise Sequences.
  The Astronomical Journal, 157(3), 127.
- Wen, H., & Liu, Z. (2016). Separating Fractal and Oscillatory
  Components in the Power Spectrum of Neurophysiological Signal. Brain
  Topography, 29(1), 13-26. https://doi.org/10.1007/s10548-015-0448-0
- Preston, M., Smith, S. & Voytek, B. (2026). Potential mechanisms and
  functional significance of aperiodic neural activity. Nature Human
  Behaviour.
- Grinsted, A. (2011). fftnoise.m (bundled unmodified, BSD-2), used
  directly in `examples/KneeCheckWorkedExample.m` and
  `examples/CompareFftnoiseVsShapedXu.m`.
- Brookshire, G. (2022). Putative rhythms in attentional switching can
  be explained by aperiodic temporal structure. Nature Human Behaviour.
  https://doi.org/10.1038/s41562-022-01364-0 (cited in
  `examples/CompareFftnoiseVsShapedXu.m` as the rationale for
  cross-checking with multiple independent estimators rather than
  trusting one method's self-agreement).

## Default citation instruction

File Exchange auto-adds a default citation instruction to new entries.
Confirmed deletable/editable if the default text doesn't fit (e.g. if
Dagmar wants a specific citation format pointing to this toolbox
distinct from [velocity-curvature-power-law-simulation](https://github.com/dagmarfraser/velocity-curvature-power-law-simulation),
the research project it was extracted from).

## Still open before this can actually be submitted

1. ~~GitHub repo~~ -- done: `https://github.com/dagmarfraser/fractal-noise`.
   Naming mismatch vs the toolbox's own `fractalnoise` name flagged in
   TODO Item 9, not yet resolved -- confirm intended before this URL
   gets baked into more places.
2. Verify this new repo's remote auth isn't a plaintext-embedded PAT
   (the concern flagged for `PowerLawSimulationPreReg`'s own repo) --
   not yet re-checked here.
3. Decide: connect repo directly, or cut a GitHub Release first and link
   that (recommended, per `mathworks/toolboxdesign`)
4. Toolbox image (optional, not blocking)
5. Dagmar's MathWorks/File Exchange account confirmed and logged in
