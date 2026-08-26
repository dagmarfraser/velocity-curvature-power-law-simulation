# Velocity-curvature power law: protocol characterisation and identifiability framework

How reliably can we measure the exponent β in the velocity-curvature power law, and does the answer depend on which analytical protocol we use — and can a recovered β ever be trusted to reveal the true underlying value?

This repository contains the pre-registration, simulation code, and analysis pipeline for a systematic evaluation of six protocols for estimating β. The work follows from Fraser et al. (2025), who found that conventional protocols compress β estimates toward 1/3 in ways that could mask genuine biological variation. That matters because reported β divergences of about |0.03| between autistic and neurotypical populations (Cook et al., 2026; Fourie et al., 2024) sit right at the edge of what these methods can resolve.

## Headline finding

Characterising *how much* a protocol biases β turned out not to be the whole question. The real question is whether the mapping from the true generative exponent (β_gen) to the value a protocol recovers (β_rec) is even invertible — whether, given a recovered value, the generator can be recovered at all.

Two methods answer that question differently:

- **Grid inversion** — the pre-registered design. Build one simulation grid spanning the full parameter space, then look up any recovered β against it directly. This assumes the β_gen → β_rec mapping is monotonic everywhere the grid covers.
- **Pointwise inversion** — a correction, not a parallel design. Checking that assumption directly (a local-monotonicity gate applied per trial) shows it fails for most real data: across seven datasets and six protocols, only 4 of 42 dataset-protocol combinations pass cleanly, 6 are conditional, and 32 fail outright. Where grid inversion isn't safe, pointwise inversion instead checks, trial by trial, whether a local, noise-calibrated inversion is possible — and reports "not recoverable" rather than a false answer when it isn't.

Grid inversion still produces the headline pooled estimate (β_gen* ≈ 0.30, k=3 core datasets). Pointwise inversion is what makes that number honest: it's the reason the paper reports invertibility status alongside any recovered value, rather than treating every β_rec as automatically meaningful.

## What's in this release

Everything through simulation, model adequacy, empirical noise characterisation, per-trial inversion (both methods), constellation validation, and structural (RSA/Mantel) validation is complete and has been run against seven real movement datasets. What follows is a guide to reproducing it, in the order it was actually built.

## Reproducing this pipeline

The whole pipeline is one lineage, not two parallel tracks: **grid**, then a **deviation** into **pointwise** correction once the grid's own core assumption (global monotonicity) was tested and found not to hold everywhere. Some stages regenerate cleanly from what's checked into this repository; others depend on artefacts too large for GitHub (100MB hard limit) that will be hosted on UBIRA eData, with a placement script provided so a fresh clone can pull them into the right place.

### 1. Grid — build the simulation

`Toolchain_caller_v057.m` + `src/functions/defineParameterSpace.m` (the `debug==5` "expanded noise grid" branch) generate the production parameter space and populate `powerlaw_debug_v058.db`.

| | |
|---|---|
| Configurations | 18,045,720 (3 sampling rates × 1 shape × 22 generative β × 14 VGF × 31 noise colours × 21 noise magnitudes × 2 kinematic-derivation methods × 3 regression methods × 5 repetitions) |
| Noise colour (α) | 0 to 6.0 (extended past the pre-registered 0–3.0 ceiling; empirical Cook/Hickman α sits at 4.1–4.6, above the original ceiling) |
| Noise magnitude (σ) | 0 to 20mm |
| Compute | ~200+ hours, single BlueBEAR HPC node (72 cores) |
| Output | `powerlaw_debug_v058.db` (~5GB) |
| **Status** | **Not regenerable from a clone alone** — needs eData download (see below) or the full HPC run |

### 2. Grid — characterise it (Model Adequacy Framework)

`ModelAdequacy_Master_v2_001.m` orchestrates four stages against `powerlaw_debug_v058.db`:

1. **Stage 1** (`ModelAdequacy_Stage1_KitchenSink_v2_001.m`) — global 7-way factorial linear mixed-effects model, all parameter interactions. Produces `stage1_results_*.mat` (~15GB) and the smaller `L9_coefficients_v004.csv`.
2. **Stage 2** (`ModelAdequacy_Stage2_Assessment_v003.m`) — adequacy check (SEM < 0.011, derived from MDC/2.77 where MDC = 0.03).
3. **Stage 3** (`ModelAdequacy_Stage3_Conditional_v005.m`) — conditional models for parameter regions the global model doesn't cover.
4. **Stage 4** (`ModelAdequacy_Stage4_Integration_v003.m`) — decision framework: given a noise profile and sampling rate, which protocol to trust.

Also produced: `perCoordinateSEM_*.mat` (measurement-error surface across the grid).

**Status:** the small outputs (`L9_coefficients_v004.csv`, `perCoordinateSEM_*.mat`) are in this repo. The stage outputs themselves (`stage1_results_*.mat` and larger) are not — they require stage 1's own multi-GB inputs, in turn requiring `powerlaw_debug_v058.db`.

Figures 1–3 (`plotBetaRecovery_v005.m`, `attractorLocation_v001.m`, `plotLMMCoefficients_v002.m`, `extractSimpleEffects_v003.m` → `plotSimpleEffects_v001.m`) plot the small Stage 1 outputs above and run cleanly from what's in this repo.

### 3. Grid inversion — apply the grid to real data

The pre-registered method: for each real trial, look up its recovered β against the simulation grid built in step 1, to recover an estimate of the true generator. This is `runLoopClosureFftnoise_v007`/`_v008` (Fraser uses `_v008`, a newer runner version, not a different tier) — sort-and-`interp1` lookup, no monotonicity check. It produces the headline pooled β_gen* (≈0.30, k=3 core datasets: Fraser, Cook CTRL, Hickman PLAC).

**Status:** needs `powerlaw_debug_v058.db` (step 1) plus each dataset's own empirical noise characterisation (step 4 below) and raw per-trial recovered β values. Not regenerable from this repo alone.

### 4. The deviation point

Grid inversion assumes the β_gen → β_rec map is monotonic everywhere. Checking that assumption directly — a local-monotonicity gate, per trial, per protocol — found it isn't: **32 of 42 dataset-protocol combinations fail invertibility outright** (4 pass, 6 conditional). A single global lookup table cannot be trusted as a general-purpose inversion tool. That finding is the reason a second method exists.

### 5. Pointwise inversion — the correction

`runLoopClosureFftnoise_v012.m` is the corrected production runner: same forward-map generation as steps 1–3, but adds a genuine local-monotonicity gate (`findBothBranches_v008.m`, part of the standalone [`monotonicity`](https://github.com/dagmarfraser/monotonicity) toolbox) before inverting. Where the gate passes, it inverts; where it doesn't, it reports that honestly rather than returning a number.

Upstream of this step:

- **Noise characterisation** (`characteriseBiologicalNoise_v004.m`) — pmtm/IRASA-based α and σ extraction from each dataset's own raw recordings. Produces `noiseCharacterisation_<dataset>.mat` (in this repo). **Needs each dataset's own raw trajectory data**, which is not redistributed here (see `data/README.md`).
- **shaped_xu noise synthesis** (`generateLoopClosureNoise_v002.m`) — matched surrogate noise from a dataset's own (α, σ), used to build the local forward map around each trial. Needs `powerlaw_debug_v058.db`.

Output: `loopClosureResults_<dataset>_all_shaped_xu_v012.mat`, one per dataset (all in this repo — each under GitHub's 100MB limit).

### 6. Constellation validation

`constellationMetrics_v004.m` — consumes step 5's outputs only. **The first point in the whole pipeline that regenerates cleanly from this repository alone**, no eData download required. Produces per-trial CCC/Pearson/MAE/Mantel contrasts between predicted and observed pipeline behaviour, across all seven datasets.

Downstream of this (all regenerate cleanly from what's checked in):

- Bland-Altman limits of agreement (`analyzeBlandAltman_v001.m`)
- Pattern preservation (`analyzePatternPreservation_v002.m`)
- SEM agreement (`analyzeSEMAgreement_v002.m`)
- Coverage probability / TDI (`analyzeTDICoverage_v001.m`)
- Pooling (`loopClosureVarDecomp_v010.m`/`_v011.m`) — the headline β_gen* pool
- Structural (RSA/Mantel) validation (`buildConstellationRDM_v002.m`, `runConstellationRSA_HPC_v007.m`) — whether the *pattern* of pipeline behaviour across noise regimes matches between simulation and real data, not just individual values

### Large files not yet in this repository

`powerlaw_debug_v058.db`, `stage1_results_*.mat`, and related multi-GB Stage 1–4 intermediates are hosted on [UBIRA eData](https://edata.bham.ac.uk) rather than GitHub (over the 100MB limit by two to three orders of magnitude). A placement script and `data/README.md`-style index for the eData-hosted files is planned but not yet built — check back here, or the eData record directly, for the current state.

## Why this exists

The velocity-curvature power law, *v(t) = VGF × κ(t)^(-β)*, is widely cited as a kinematic invariant, with β close to 1/3 for typical movement. Fraser et al. (2025) showed that conventional analytical protocols (Butterworth filtering followed by log-transformed OLS regression) introduce systematic biases that push β estimates toward 1/3 regardless of the true value. If the measurement tool itself is biased toward the expected answer, you can't tell when the law genuinely holds and when it doesn't — and, as this project's own headline finding shows, you can't always tell whether the true value is recoverable from the measurement at all.

## The six protocols

Two kinematic derivation methods crossed with three regression approaches.

**Kinematic derivation:**

BWFD (Butterworth + finite differences) applies a second-order low-pass filter at 10 Hz with zero-phase correction, then differentiates numerically. Each differentiation step amplifies high-frequency noise, and the filter alters the noise structure.

SG (Savitzky-Golay) fits local polynomials with a sampling-rate-scaled window, extracting smoothed derivatives in one operation. This sidesteps the repeated-differentiation issue noted above.

**Regression:**

OLS: log-transform both sides, fit a line. Simple and ubiquitous, but the log transformation can distort the underlying relationship.

LMLS (Levenberg-Marquardt): nonlinear fit on untransformed data. Avoids the log distortion but remains sensitive to outliers.

IRLS (iteratively reweighted least squares): robust nonlinear fit with bisquare weighting. Resists outliers at the cost of longer computation.

**What we know:**

BWFD-OLS is the dominant protocol in the literature and shows the most severe compression toward 1/3 under noise. SG-IRLS is the best-performing single protocol across the constellation validation; SG-LMLS is a close second. No single protocol is adequate everywhere in the noise/sampling-rate parameter space — the four-stage model adequacy framework (step 2 above) exists to map exactly where each one is and isn't trustworthy.

## Empirical validation datasets

Stage 5 (post-simulation empirical validation) tests the framework against seven trial sets from five studies, spanning clinical populations and recording technologies. Data are linked rather than redistributed; see `data/README.md` for full references and access details.

| Dataset | Description | Sampling rate |
|---|---|---|
| Fraser | iPad Pro 13" M4 + Apple Pencil, StimuliApp, shape tracing (N=45 usable of 48 recruited) | 240 Hz |
| Zarandi et al. 2023 | WACOM Bamboo Slate, practised ellipses | 100 Hz |
| Dhieb et al. 2022 | Graphical tablet, ellipse drawing, ages 19–85 | 100 Hz |
| Cook et al. 2026 | WACOM, autistic vs neurotypical tracing (CTRL / ASD) | 133 Hz |
| Hickman et al. 2024 | WACOM, haloperidol (D2 antagonist) vs placebo, neurotypical adults (PLAC / HALO) | 133 Hz |

## Toolboxes developed alongside this project

Three general-purpose MATLAB toolboxes were extracted from this project's own machinery and shipped independently (MIT licensed, MATLAB File Exchange + GitHub):

- [`monotonicity`](https://github.com/dagmarfraser/monotonicity) — the local-monotonicity segmenter behind pointwise inversion's own gate.
- [`concordance`](https://github.com/dagmarfraser/concordance) — Lin's Concordance Correlation Coefficient, cross-validated against R's `DescTools::CCC`.
- [`mantel`](https://github.com/dagmarfraser/mantel) — Mantel test for structural (RDM-based) validation.
- [`fractal-noise`](https://github.com/dagmarfraser/fractal-noise) — 1/f^α coloured noise generation (Xu, 2019 fractional-differencing method), used throughout the simulation.

## Notation

| Symbol | Meaning |
|---|---|
| β_gen | True (generative) β |
| β_rec | Recovered β from a pipeline |
| β_gen* | Pooled estimate of the true generator, recovered via inversion |
| bias | β_gen − β_rec |
| SEM | Standard error of measurement |
| MDC | Minimal detectable change (0.03, the clinically-relevant autism/neurotypical divergence) |
| α | Noise colour exponent (0 = white, 1 = pink, 2 = red, 3 = black) |
| σ | Noise magnitude in mm |
| VGF | Velocity gain factor |
| CCC | Lin's Concordance Correlation Coefficient |
| Grid inversion | Pre-registered whole-grid lookup method (`v007`/`v008`) |
| Pointwise inversion | Per-trial, noise-calibrated, monotonicity-gated correction method (`v009`–`v012`) |

## Lin's CCC

MATLAB has no built-in CCC. The `concordance` toolbox's implementation agrees with R's `DescTools::CCC` to six decimal places across all test cases.

```matlab
result = linCCC_v001(y1, y2);
```

McBride (2005) thresholds: > 0.99 almost perfect, 0.95–0.99 substantial, 0.90–0.95 moderate, < 0.90 poor.

## Usage

Requires: Database Toolbox, Parallel Computing Toolbox, Curve Fitting Toolbox, Statistics and Machine Learning Toolbox.

```matlab
% Check your installation has everything it needs
cd src
validate_prereg_readiness

% Quick validation run with ersatz data (no full simulation)
ModelAdequacy_Master_v2_001(2, 5, false)

% Production (requires powerlaw_debug_v058.db — see "Large files" above)
Toolchain_caller_v057;                       % Grid: build the simulation
ModelAdequacy_Master_v2_001(9, 5, true)      % Grid: characterise it (full parameter space)

% Empirical validation, once loopClosureResults_*.mat exist for a dataset
constellationMetrics_v004()                  % Pointwise: constellation validation
```

`ModelAdequacy_Master_v2_001(level, nObs, useDatabase, resumeStage)` takes tractability levels 1–9. Levels 1–3 are for quick prototyping; level 9 (Full-Original, 14.8M of the 18.0M total configurations) needs HPC. The optional `resumeStage` argument picks up after interruption.

## Project links

- OSF: https://osf.io/dwxa2/
- GitHub: https://github.com/dagmarfraser/velocity-curvature-power-law-simulation
- Large datasets: https://edata.bham.ac.uk (UBIRA eData)
- Target journal: Behavior Research Methods

## Xiao et al. reimplementation

The `power_analysis.m` function reimplements the general guidelines from Xiao et al. (2011) for choosing between log-transformed regression (LR), nonlinear regression (NLR), and AICc-weighted model averaging when fitting power law data. The `useRevised` parameter (default: `true`) controls whether bootstrap confidence intervals recompute AICc weights per resample (fixing a bug in the original R code where global weights were reused).

`test_xiao_crossvalidate_v002.m` validates all six conditions (three method-selection paths times two switch positions) against the original R source, confirming agreement to machine precision for LR, ~10⁻⁸ for NLR (optimiser tolerance), and comparable bootstrap CIs for model averaging.

```matlab
run('src/test_xiao_crossvalidate_v002.m')  % requires R on PATH
```

## File structure

```
├── README.md
├── prereg_v101.docx
├── CITATION.cff
├── LICENSE
├── src/
│   ├── Toolchain_caller_v057.m                Grid: simulation orchestration
│   ├── Toolchain_func_v032.m                  Grid: simulation processing engine
│   ├── ModelAdequacy_Master_v2_001.m          Grid: model adequacy orchestration
│   ├── ModelAdequacy_Stage[1-4]_*             Grid: model adequacy stages
│   ├── runLoopClosureFftnoise_v012.m          Pointwise: production inversion runner
│   ├── constellationMetrics_v004.m            Pointwise: constellation validation
│   ├── analyze{BlandAltman,PatternPreservation,SEMAgreement,TDICoverage}_*.m
│   ├── buildConstellationRDM_v002.m           Pointwise: structural (RSA) validation
│   ├── runConstellationRSA_HPC_v007.m
│   ├── loopClosureVarDecomp_v0{10,11}.m       Pointwise: headline pooling
│   ├── functions/
│   │   ├── differentiateKinematicsEBR.m
│   │   ├── regressDataEBR.m
│   │   ├── curvatureKinematicEBR.m
│   │   ├── generateSyntheticData_v011.m
│   │   ├── linCCC_v001.m
│   │   ├── iraAlphaSigma_v001.m               IRASA spectral exponent (pmtm-based)
│   │   ├── generateCustomNoise_v004.m         Coloured noise via XuNoise (α 0-6)
│   │   ├── findBothBranches_v008.m            Local-monotonicity gate (monotonicity toolbox)
│   │   ├── defineParameterSpace.m
│   │   └── ...
│   ├── req/                                   Vendored third-party dependencies
│   ├── test_xiao_crossvalidate_v002.m
│   ├── test_R_CCC.m
│   └── validate_prereg_readiness.m
├── data/                                      See data/README.md for sources
├── results/                                    Populated after simulation
├── figures/
└── reports/                                    Model adequacy HTML reports (generated)
```
