# Changelog

`mantel` follows semantic versioning.

## v1.0.1 -- 2026-07-27

- **Fixed: O(n^2) tie-handling in `mantelRank` at scale.** The average-
  rank tie correction (looping once per unique value, each iteration
  doing a full-length vectorized comparison against the whole array)
  was correct but not scale-aware -- fine at the toolbox's own
  validation sizes (n=6, n=12) and the narrow originating 6-scalar use
  case, but caught during first real production use:
  `PowerLawSimulationPreReg`'s trial-by-trial RSA calls `mantelTest` on
  the lower triangle of a 4237x4237 trial dissimilarity matrix (~9
  million elements), where O(n^2) tie-handling doesn't just run
  slowly, it never completes in practical time. Manifested as a silent
  multi-day `parfor` hang with no error output, since the failure was a
  single non-terminating call rather than an exception. Diagnosed by a
  minimal probe using a simpler no-tie-correction rank (fast, 3.8s at
  NPerm=10) to isolate that ranking itself, not permutation count or
  worker allocation, was the bottleneck, then confirmed by direct
  inspection of the loop. Fixed via `sort` + group-boundary detection
  (`diff`) + `accumarray(..., @mean)`: mathematically identical
  average-rank output, O(n log n) instead of O(n^2). No interface
  change. Existing n=6/n=12 validation results were re-verified
  unaffected.

## v1.0.0 -- 2026-07-02

First release. One function, delivered from `PowerLawSimulationPreReg`
session 52:

- `mantelTest.m` -- generalised from `mantel_local` (a private local
  function in `src/constellationMetrics_v002.m`), which computed a
  narrow special case: correlating `abs(scalar_p - scalar_q)`
  "distances" derived from 6 raw pipeline values. Generalised to accept
  either arbitrary distance/dissimilarity matrices directly, or raw
  observations converted via Euclidean distance -- matching the
  standard textbook Mantel test interface (Mantel, 1967) rather than
  the narrower originating use case.

New in this extraction, not present in the source:

- **`Exact` auto-enumeration.** The source always used random-sampling
  permutation regardless of sample size. R's `vegan::mantel()` was found
  during validation to automatically switch to exact permutation
  enumeration when the sample size makes it tractable (n! small enough),
  giving an exact p-value instead of an approximation. `mantelTest`'s
  `Exact="auto"` (default) replicates this: exact enumeration triggers
  automatically at n<=8, is forceable via `Exact="on"` up to a hard
  safety limit of n<=10 (errors clearly above that rather than
  attempting an intractable allocation), and can be disabled entirely
  via `Exact="off"` to match the source's original always-random
  behaviour. Verified to reproduce `vegan`'s exact p-value precisely
  (0.001388889 on the validation dataset, not just approximately).
- **No hard Statistics and Machine Learning Toolbox dependency.**
  Checked whether `corr`/`pdist`/`squareform`/`tiedrank` (all Statistics
  Toolbox) were needed and rewrote the design to avoid them: Pearson
  correlation via `corrcoef` (base MATLAB), Spearman via a hand-rolled
  average-rank tie handler followed by Pearson-on-ranks, and the raw-
  observations input path computes Euclidean distance directly.
  `perms` (used for exact enumeration) is also base MATLAB. Matches the
  "dependency-light" ethos already established by `fractalnoise` and
  `concordance`.
- **`CompareR` opt-in feature.** Cross-validates against R's
  `vegan::mantel()` live, on request, gated behind an availability check
  (`Rscript` on path, then the `vegan` package specifically) -- same
  design pattern as `concordance`'s `linCCC`/`CompareR` and
  `fractalnoise`'s `estimateIRASA`/`CompareFieldTrip`. Requesting it
  without R + vegan installed does not error -- warns visibly and
  returns a clearly-flagged unavailable result. Both branches (available
  and degraded-mode) were genuinely verified during development: the
  success path by confirming exact agreement with R (statistic gap
  1.11e-16, p-value gap 0 once exact enumeration matched on both sides),
  and the degraded-mode path by actually stripping the `PATH`
  environment variable to simulate R's absence. The R-calling logic is
  inlined as a local function within `mantelTest.m` itself, applying
  `fractalnoise`'s self-containment lesson proactively.
- **A caught-and-fixed bug.** `isempty("")` is `false` for MATLAB
  `string` type (a 1x1 array with zero-length content, not an empty
  array) -- unlike `isempty('')` for `char`, which is `true`. This
  silently prevented `rInfo.message`'s explanatory disclaimer text from
  ever populating on the `CompareR` success path. Fixed via
  `strlength(message) == 0`. Checked `concordance`'s `linCCC.m` for the
  same bug pattern -- confirmed absent there.
- **Fail-Loud input validation added beyond source.** The source assumed
  well-formed 6-element input implicitly; `mantelTest` explicitly
  validates that distance matrices are square, symmetric (within 1e-9),
  and have a zero diagonal (within 1e-9), erroring clearly rather than
  silently proceeding with malformed input.

Plus `doc/GettingStarted.mlx`, `doc/WhyValidatedAgainstR.mlx` (the R
cross-validation story, including the exact-vs-random p-value
distinction explained, reproducible live via `CompareR`),
`examples/IsolationByDistanceWorkedExample.mlx` (a classic population-
genetics application at n=12, deliberately larger than the validation
example's n=6 to also demonstrate the random-sampling fallback when
exact enumeration is intractable), and a 12-test suite.

See `docs/README_MakingToolboxes_v001.md` (in the parent
`PowerLawSimulationPreReg` project) for the general lessons this
extraction drew on, from building `fractalnoise` and `concordance` first.
