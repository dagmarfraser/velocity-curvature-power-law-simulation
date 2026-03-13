# Velocity-curvature power law: protocol characterisation framework

How reliably can we measure the exponent β in the velocity-curvature power law, and does the answer depend on which analytical protocol we use?

This repository contains the pre-registration, simulation code, and analysis pipeline for a systematic evaluation of six protocols for estimating β. The work follows from Fraser et al. (2025), who found that conventional protocols compress β estimates toward 1/3 in ways that could mask genuine biological variation. That matters because reported β divergences of about |0.03| between autistic and neurotypical populations (Cook et al., 2026; Fourie et al., 2024) sit right at the edge of what these methods can resolve.

## What's in this release (and what isn't)

This is the pre-registration release (tagged `v1.0.2-prereg`). The simulation and model adequacy code are complete and ready to run. Empirical validation tooling (Stage 5) will follow after the simulation itself is done.

**Present:**

- Pre-registration document (`prereg_v101.docx` - and registered here https://doi.org/10.17605/OSF.IO/CKR8H)
- Phase 1 simulation pipeline (`Toolchain_caller_v057.m`, `Toolchain_func_v032.m`)
- Phase 2 model adequacy framework, all four stages
- Core analytical functions: kinematic derivation, regression, curvature, synthetic data generation
- Lin's CCC implementation, validated against R's DescTools to six decimal places
- Xiao et al. power law analysis: MATLAB reimplementation (`power_analysis.m`) cross-validated against the original R code across all three method-selection paths (LR, NLR, model averaging) and both bootstrap variants
- Parameter space definition with tractability levels 1 through 9
- Database infrastructure (SQLite backend with checkpointing)

**Not present yet:**

- Noise characterisation wrapper (pmtm-based α and σ extraction from empirical recordings)
- Empirical validation scripts (Stage 5, Section 7 of the pre-registration)
- RSA/Mantel test implementation (exploratory supplementary analysis)
- Bias correction and invertibility assessment (Section 7.7)
- Simulation results themselves (these need roughly 200 hours of HPC time and will likely be too large for github)
- Empirical datasets (see `data/README.md` for sources and access; not redistributed here)

## Why this exists

The velocity-curvature power law, *v(t) = VGF × κ(t)^(-β)*, is widely cited as a kinematic invariant, with β close to 1/3 for typical movement. Fraser et al. (2025) showed that conventional analytical protocols (Butterworth filtering followed by log-transformed OLS regression) introduce systematic biases that push β estimates toward 1/3 regardless of the true value. If the measurement tool itself is biased toward the expected answer, you can't tell when the law genuinely holds and when it doesn't.

This framework maps measurement precision across 14.7 million parameter configurations, so that researchers can make informed choices about which protocol to trust given their data's characteristics.

## The six protocols

Two kinematic derivation methods crossed with three regression approaches.

**Kinematic derivation:**

BWFD (Butterworth + finite differences) applies a second-order low-pass filter at 10 Hz with zero-phase correction, then differentiates numerically. Each differentiation step amplifies high-frequency noise, and the filter alters the noise structure.

SG (Savitzky-Golay) fits local polynomials with a sampling-rate-scaled window, extracting smoothed derivatives in one operation. This sidesteps the repeated-differentiation issue noted above.

**Regression:**

OLS: log-transform both sides, fit a line. Simple and ubiquitous, but the log transformation can distort the underlying relationship.

LMLS (Levenberg-Marquardt): nonlinear fit on untransformed data. Avoids the log distortion but remains sensitive to outliers.

IRLS (iteratively reweighted least squares): robust nonlinear fit with bisquare weighting. Resists outliers at the cost of longer computation.

**What we know and what we don't:**

BWFD-OLS is the dominant protocol in the literature. Fraser et al. (2025) showed it introduces systematic compression. SG-LMLS and SG-IRLS are the "vetted" protocols from that paper, which don't compress to 1/3.  However this was demonstrated only for a narrow simulation.

The four hybrid combinations (BWFD-LMLS, BWFD-IRLS, SG-OLS, SG-LMLS) haven't been simulated, indeed none of the 6 protocols has been systematically characterised across a broad spectrum of biologically plausible paramters. The interesting question is whether robust regression can compensate for the artefacts that Butterworth filtering introduces. If BWFD-IRLS performs comparably to SG-OLS, it suggests regression choice dominates overall performance. If SG-based protocols consistently outperform BWFD-based ones regardless of regression method, the derivation step is what matters and can't be rescued downstream.

## Implementation

### Phase 1: simulation (`Toolchain_caller_v057.m`)

Synthetic elliptical trajectories with known β and VGF are contaminated with coloured noise and processed through all six protocols. The parameter space covers:

- Sampling rates: 60, 120, 240 Hz
- Generative β: 0 to 2/3 (21 levels)
- Velocity gain factor: 14 levels
- Noise colour α: 0 to 3.0 (31 levels, white through black)
- Noise magnitude σ: 0 to 10 mm (18 levels)
- Repetitions: 5 per condition

That gives 14,764,680 total configurations, 5 GB+ of results, and roughly 200 hours on a single BlueBEAR HPC node (72 cores).

### Phase 2: model adequacy (`ModelAdequacy_Master_v002.m`)

Four stages, run sequentially:

Stage 1 fits global linear mixed-effects models (via `fitlme`) with all parameter interactions, following an interaction-first analysis strategy with simple effects.

Stage 2 checks whether the global model is adequate everywhere, using bootstrap stability, residual diagnostics, and cross-validation. The adequacy criterion is SEM < 0.011, which comes from MDC/2.77 where MDC = 0.03 is the clinically relevant difference for autism research.

Stage 3 builds conditional models for parameter regions where the global model falls short.

Stage 4 pulls it all together into a decision framework: given your data's noise profile and sampling rate, which protocol should you use?

## Usage

Requires: Database Toolbox, Parallel Computing Toolbox, Curve Fitting Toolbox, Statistics and Machine Learning Toolbox, Bioinformatics Toolbox.

```matlab
% Check your installation has everything it needs
cd src
validate_prereg_readiness

% Quick validation run with ersatz data (random response values, no simulation)
ModelAdequacy_Master_v002(2, 5, false)

% Production
Toolchain_caller_v057;                  % Phase 1
ModelAdequacy_Master_v002(9, 5, true)   % Phase 2, full parameter space
```

`Toolchain_caller_v057` accepts debug levels 0-3. Level 0 is full production; higher levels reduce the parameter space for development and testing.

`ModelAdequacy_Master_v002(level, nObs, useGroundTruth, resumeStage)` takes tractability levels 1-9. Levels 1-3 are for quick prototyping; 7-9 need HPC. The optional `resumeStage` argument picks up after interruption.

## Notation

| Symbol | Meaning |
|--------|---------|
| β_gen | True (generative) β |
| β_rec | Recovered β from a pipeline |
| bias | β_gen - β_rec |
| SEM | Standard error of measurement |
| MDC | Minimal detectable change |
| α | Noise colour exponent (0 = white, 1 = pink, 2 = red) |
| σ | Noise magnitude in mm |
| CCC | Lin's Concordance Correlation Coefficient |

## Lin's CCC

MATLAB has no built-in CCC. The custom implementation (`src/functions/linCCC_v001.m`) agrees with R's DescTools::CCC to six decimal places across all test cases.

```matlab
result = linCCC_v001(y1, y2);
```

McBride (2005) thresholds: > 0.99 almost perfect, 0.95-0.99 substantial, 0.90-0.95 moderate, < 0.90 poor.

## Empirical validation databases

Stage 5 (post-simulation) will test the framework against seven movement databases spanning species, clinical populations, and recording technologies. Data are linked rather than redistributed; see `data/README.md` for full references and access details.

| Dataset | Description |
|---------|-------------|
| Fraser (pilot) | iPad Pro 240 Hz shape tracing |
| Zarandi et al. 2023 | WACOM 100 Hz, practised ellipses |
| Dhieb et al. 2022 | Tablet drawing, ages 19-85 |
| Cook et al. 2026 | WACOM 133 Hz, autistic vs neurotypical |
| Dagenais et al. 2021 | Elephant trunk 3D trajectories |
| James et al. 2020 | Bumblebee locomotion |
| Hickman et al. 2024 | Parkinson's ON/OFF medication |

## Project links

- OSF: https://osf.io/dwxa2/
- GitHub: https://github.com/dagmarfraser/velocity-curvature-power-law-simulation
- Large datasets: https://edata.bham.ac.uk (UBIRA eData)

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
│   ├── Toolchain_caller_v057.m       Phase 1 orchestration
│   ├── Toolchain_func_v032.m         Phase 1 processing engine
│   ├── ModelAdequacy_Master_v002.m   Phase 2 orchestration
│   ├── ModelAdequacy_Stage[1-4]_*    Phase 2 stages
│   ├── functions/
│   │   ├── differentiateKinematicsEBR.m
│   │   ├── regressDataEBR.m
│   │   ├── curvatureKinematicEBR.m
│   │   ├── generateSyntheticData_v011.m
│   │   ├── linCCC_v001.m
│   │   ├── defineParameterSpace.m
│   │   └── ...
│   ├── req/
│   │   └── xiaoxiao/3551973/
│   │       ├── power_analysis.m           MATLAB reimplementation of Xiao et al.'s guidelines
│   │       ├── Sup_2_Guidelines.r         Original R source (Xiao et al.)
│   │       └── Sup_2_Guidelines_revised.r Revised R source (bootstrap bug fix)
│   ├── test_xiao_crossvalidate_v002.m    Cross-validates MATLAB vs R across all 3 paths × 2 switch positions
│   ├── test_R_CCC.m                      Validates linCCC_v001 against R's DescTools::CCC
│   └── validate_prereg_readiness.m       Pre-submission verification checks
├── data/                             See data/README.md for sources
├── docs/                             CCC validation demos
├── results/                          Populated after simulation
├── figures/
└── reports/                          Model adequacy HTML reports (generated)
```

