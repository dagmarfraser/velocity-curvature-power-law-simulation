# concordance

A small, standalone MATLAB toolbox for Lin's Concordance Correlation
Coefficient (CCC): a measure of agreement between two continuous
variables, distinct from (and stricter than) Pearson correlation -- two
variables can correlate perfectly while disagreeing completely, and CCC
catches that where Pearson's r doesn't.

Originally written for a kinematics research project
([velocity-curvature-power-law-simulation](https://github.com/dagmarfraser/velocity-curvature-power-law-simulation)),
then generalised into this standalone, dependency-light package.

**Status:** v1.0.0 (2026-07-02).

---

## What's in it

One function: `linCCC.m` -- computes CCC, its confidence interval, and a
precision/accuracy decomposition (Pearson r and bias correction Cb).
Optionally cross-validates against R's `DescTools::CCC()` for an
independent second opinion, via `CompareR=true`.

No optional-toolbox dependencies for the core computation. `CompareR`
needs R (with the DescTools package) installed and degrades visibly, not
silently, when it isn't: a warning and a clearly-flagged unavailable
result, never an error and never a silent skip.

## Installation

```matlab
addpath('/path/to/concordance/toolbox');
```

## Quick example

```matlab
methodA = randn(30, 1) * 10 + 100;
methodB = methodA + randn(30, 1) * 3;

result = linCCC(methodA, methodB);
fprintf('CCC = %.4f\n', result.ccc);

% Or with no output argument, for a printed interpretation:
linCCC(methodA, methodB);
```

See `doc/GettingStarted.m` for a full tour, and
`examples/AgreementAnalysisWorkedExample.m` for a realistic
method-comparison scenario distinguishing good agreement from systematic
bias and proportional (scale) bias.

## Why CCC and not just Pearson r

Pearson r measures whether two variables move together linearly -- it's
blind to systematic offsets or scale differences. Two measurement methods
can correlate at r=0.99 while one consistently reads higher than the
other; Pearson won't flag this, but CCC will, because CCC specifically
penalises deviation from the 45-degree identity line, not just linear
association. See `doc/GettingStarted.m` for a worked demonstration of
exactly this.

## Validation

Cross-validated against R's `DescTools::CCC()` across 7 test cases
(|MATLAB - R| < 0.001 for every case), including the canonical
Bland-Altman (1986) Peak Expiratory Flow Rate dataset -- the same worked
example used in `DescTools`' and `epiR`'s own documentation. See
`doc/WhyValidatedAgainstR.m` for the full comparison, reproducible live
if you have R installed.

## Testing

```matlab
results = runtests('/path/to/concordance/tests');
```

10 tests, covering source-equivalence against the original internal
function, the canonical PEFR dataset, `arguments`-block validation, error
paths, NaN-pair handling, and both branches of the `CompareR` feature
(genuinely verified -- both the R-available success path and the
R-unavailable degraded-mode path were tested directly against real
environment states, not just inferred from reading the code).

## Credits

- Lin LI-K (1989). A concordance correlation coefficient to evaluate
  reproducibility. Biometrics 45(1):255-268.
- Lin LI-K (2000). A note on the concordance correlation coefficient.
  Biometrics 56(1):324-325. [Correction to the variance formula]
- McBride GB (2005). A proposal for strength-of-agreement criteria for
  Lin's concordance correlation coefficient. NIWA Client Report.
- Bland JM, Altman DG (1986). Statistical methods for assessing agreement
  between two methods of clinical measurement. Lancet 327:307-310.
  (Source of the canonical PEFR validation dataset.)
- Signorell A (2025). DescTools: Tools for Descriptive Statistics. R
  package version 0.99.60. (Reference implementation for `CompareR`;
  version as verified installed and tested against on this machine,
  2026-07-20.)

## Origin

Extracted and generalised from `linCCC_v001.m`, an internal function of
[velocity-curvature-power-law-simulation](https://github.com/dagmarfraser/velocity-curvature-power-law-simulation)
(University of Birmingham, School of Psychology). Each
extraction/generalisation change is documented in `linCCC.m`'s own
docstring and in `CHANGELOG.md`. `concordance` has no dependency back on
that project.

## License

MIT. See `license.txt`.
