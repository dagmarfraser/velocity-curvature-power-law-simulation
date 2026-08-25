# mantel

A small, standalone MATLAB toolbox for the Mantel test -- a permutation-
based test for correlation between two distance/dissimilarity matrices
on the same set of items. Extracted and generalised from an internal
function of the `PowerLawSimulationPreReg` project (University of
Birmingham, School of Psychology). `mantel` itself has no connection to
that project's power-law or kinematics work -- it's a general-purpose
statistical test, offered here as a standalone, dependency-light utility.

**Status:** v1.0.1 (2026-07-27).

---

## Why this toolbox

MathWorks' Statistics and Machine Learning Toolbox has no built-in
Mantel test. File Exchange search results for "Mantel test" are mostly
the **Mantel-Haenszel test** -- a different statistic by the same first
author (Nathan Mantel), for stratified contingency tables, not distance-
matrix correlation. The actual matrix-correlation Mantel test (Mantel,
1967) exists as `bramila_mantel.m` (Glerean, 2013), but is hosted on
figshare rather than File Exchange, is single-author with no visible
independent validation, and is not widely known in the MATLAB community.
`mantel` aims to fill that gap properly: File-Exchange-hosted, clearly
disambiguated from Mantel-Haenszel, general-purpose, and cross-validated
against R's authoritative `vegan::mantel()` implementation.

## What's in it

One function: `mantelTest.m` -- computes the Mantel correlation
statistic (Pearson or Spearman) and its permutation-based significance,
either from two pre-computed distance matrices or from raw observations
(converted via Euclidean distance). Automatically uses exact permutation
enumeration when the sample size makes it tractable (n<=8 by default),
giving an exact p-value rather than a random-sampling approximation.
Optionally cross-validates against R's `vegan::mantel()`, via
`CompareR=true`.

No hard dependency on Statistics and Machine Learning Toolbox: Pearson
correlation uses `corrcoef` (base MATLAB), Spearman is implemented via a
hand-rolled average-rank tie handler, and the raw-observations input path
computes Euclidean distance directly -- none of `corr`, `pdist`,
`squareform`, or `tiedrank` are used.

`CompareR` needs R (with the `vegan` package) installed and degrades
visibly, not silently, when it isn't: a warning and a clearly-flagged
unavailable result, never an error and never a silent skip.

## Installation

```matlab
addpath('/path/to/mantel/toolbox');
```

## Quick example

```matlab
D1 = [0 12 8; 12 0 6; 8 6 0];   % e.g. geographic distance
D2 = [0 0.4 0.3; 0.4 0 0.2; 0.3 0.2 0];   % e.g. a dissimilarity measure

[mantelR, pValue, info] = mantelTest(D1, D2);
fprintf('Mantel r = %.3f, p = %.3f (n=%d, exact=%d)\n', ...
    mantelR, pValue, info.n, info.exact);
```

See `doc/GettingStarted.m` for both input forms (distance matrix, raw
observations) and `examples/IsolationByDistanceWorkedExample.m` for a
classic population-genetics application at a larger sample size.

## Validation

Cross-validated against R's `vegan::mantel()`: the correlation statistic
matched to ~1e-16 (machine precision), and the exact-enumeration p-value
matched precisely (not just approximately) once `mantelTest`'s
`Exact="auto"` behaviour was built to replicate `vegan`'s own automatic
exact-enumeration threshold. See `doc/WhyValidatedAgainstR.m` for the
full comparison, including an explanation of why a p-value computed with
random-sampling permutation (`Exact="off"`) will legitimately differ
from an exact one at small n -- expected behaviour, not a bug.

## Testing

```matlab
results = runtests('/path/to/mantel/tests');
```

15 tests across two test classes. `testMantelTest.m` covers source-
equivalence, the known exact-mode p-value against a validated `vegan`
result, both input forms, both correlation methods, `Exact` auto/on/off
behaviour including the hard safety limit, validation error paths,
degenerate-input handling, both branches of `CompareR` (genuinely
verified via real environment manipulation, not just inferred from
reading the code), and -- new in v1.0.1 -- `mantelRank`'s average-rank
tie handling cross-validated against R's own `rank(ties.method="average")`
at both a small (n=30) and moderate (n=300) tie-rich scale.
`testToolboxSelfContained.m` verifies the toolbox resolves and runs
correctly with nothing but its own folder on path -- no hidden
dependency on convenience paths added by other test setup.

## Credits

- Mantel, N. (1967). The detection of disease clustering and a
  generalized regression approach. Cancer Research, 27(2), 209-220.
- Legendre, P., Fortin, M.-J., & Borcard, D. (2015). Should the Mantel
  test be used in spatial analysis? Methods in Ecology and Evolution,
  6(11), 1239-1247.
- Somers, K. M., & Jackson, D. A. (2022). Putting the Mantel test back
  together again. Ecology, 103(10), e3780.
- Oksanen, J. et al. (2024). vegan: Community Ecology Package. R package
  version 2.6-8. (Reference implementation for `CompareR`.)

## Origin

Generalised from a private local function (`mantel_local`) in
`src/constellationMetrics_v002.m` of `PowerLawSimulationPreReg` -- a
narrow, 6-scalar-pipeline-value construction, extended here to accept
arbitrary distance matrices or raw observations, matching the standard
textbook Mantel test interface. See `CHANGELOG.md` for the full list of
extraction changes. `mantel` has no dependency back on
`PowerLawSimulationPreReg`.

## License

MIT. See `license.txt`.
