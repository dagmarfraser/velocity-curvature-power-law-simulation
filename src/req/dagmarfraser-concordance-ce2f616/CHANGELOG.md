# Changelog

`concordance` follows semantic versioning.

## Unreleased

- **Fixed:** `tests/testLinCCC.m` had a hardcoded personal absolute path
  (`LegacyFunctionsDir`, pointing at a specific machine's Dropbox
  folder) for its optional source-equivalence check against the
  original `linCCC_v001.m`. This should never have shipped publicly --
  it exposed a personal folder structure/username for no functional
  benefit, and was already stale relative to the machine it was written
  on (the test had accordingly been silently filtering on every machine
  in every session since v1.0.0, never actually verified). Replaced with
  an opt-in `CONCORDANCE_LEGACY_DIR` environment variable: unset (the
  default for everyone, including now-current project machines) filters
  cleanly with the same clear message as before; set to a real path, the
  test now genuinely runs and has been confirmed passing (10/10, bit-
  identical match to `linCCC_v001.m` across all 11 output fields) for
  the first time since extraction.

## v1.0.0 -- 2026-07-02

First release. One function, delivered from `PowerLawSimulationPreReg`
session 52:

- `linCCC.m` -- extracted from `linCCC_v001.m`, `_vNNN` suffix dropped.
  Core computation (CCC, Pearson r, Cb bias correction, Z-transform
  confidence interval) unchanged from source and verified bit-identical,
  including against the canonical Bland-Altman (1986) PEFR dataset.

New in this extraction, not present in the source:

- **`CompareR` opt-in feature.** Cross-validates against R's
  `DescTools::CCC()` live, on request, gated behind an availability check
  (same design pattern as `fractalnoise`'s `estimateIRASA` /
  `CompareFieldTrip`). Requesting it without R + DescTools installed does
  not error -- warns visibly and returns a clearly-flagged unavailable
  result. Both branches (available and degraded-mode) were genuinely
  verified during development: the success path by confirming exact
  agreement with R (gap=0.000000 across all 7 validation cases), and the
  degraded-mode path by actually stripping the `PATH` environment
  variable to simulate R's absence and confirming the real function
  still behaves correctly, not just inferred from reading the code.
  The R-calling logic (adapted from `computeCCC_R.m`) is inlined as a
  local function within `linCCC.m` itself, not calling an external file
  outside the toolbox's own tree -- applying `fractalnoise`'s
  self-containment lesson proactively rather than discovering the same
  bug class reactively.
- **A caught-and-fixed regression.** MATLAB's `arguments` block cannot
  express a cross-parameter validation (y1 and y2 must be the same
  length), so initially removing the source's explicit manual check for
  this left it falling through to a confusing generic MATLAB
  broadcasting error instead of a clear message. Restored the explicit
  check.
- **A disclosed validation improvement.** `alpha` (significance level)
  now requires `mustBeInRange(alpha, 0, 1, "exclusive")` -- the source
  let `alpha=0` or `alpha=1` through silently, producing `Inf`
  confidence bounds rather than a clear error.
- **Naming decision:** kept `alpha` as the parameter name (matches Lin
  1989's own notation and standard statistical terminology), with an
  explicit docstring note disambiguating it from the unrelated spectral-
  exponent "alpha" used throughout the separate `fractalnoise` toolbox.

Plus `doc/GettingStarted.mlx`, `doc/WhyValidatedAgainstR.mlx` (the R
cross-validation story, reproducible live via `CompareR`),
`examples/AgreementAnalysisWorkedExample.mlx` (a realistic method-
comparison scenario distinguishing good agreement, systematic bias, and
proportional/scale bias -- deliberately new content, not reusing the
PEFR dataset already covered in the validation doc), and a 10-test suite.

See `docs/README_MakingToolboxes_v001.md` (in the parent
`PowerLawSimulationPreReg` project) for the general lessons this
extraction drew on, from building `fractalnoise` first.
