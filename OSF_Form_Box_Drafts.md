# OSF Registration Form: Trimmed Box Content

**Strategy:** Each box gives enough for a reviewer to understand the study without opening the PDF, then redirects with a consistent closing line. The attached pre-registration document (prereg_v101.docx) is the authoritative specification.

---

## Metadata

### Title
Systematic Characterisation of Velocity-Curvature Power Law Analysis Protocols Across Biologically Informed Parameter Space

### Description
Systematic Characterisation of Velocity-Curvature Power Law Analysis Protocols Across Biologically Informed Parameter Space

### Contributors
Dagmar Scott Fraser, Jennifer Cook, Massimiliano Di Luca

---

## Overview

### Research questions or hypotheses

Fraser et al. (2025) demonstrated that conventional analytical protocols for the velocity-curvature power law (v ∝ κ^(−β)) systematically bias β toward ~1/3 in specific simulations, potentially masking genuine biological variation. However, neither conventional nor alternative protocol elements have been characterised across a comprehensive range of biologically plausible conditions.

This study asks: across a full factorial space of biologically plausible β, noise profiles, and sampling rates, which analytical protocols can reliably recover the true β exponent with the measurement precision (SEM < 0.011) needed to detect clinically meaningful divergences (|Δβ| = 0.03) between populations?

We systematically evaluate six analytical protocols in 14.7 million parameter combinations, then validate the resulting characterisation framework against seven empirical movement databases spanning clinical populations, developmental trajectories, and species.

Full rationale and theoretical framework: attached pre-registration document, Sections 1-2.

### Foreknowledge of data or evidence

Data does not yet exist. No part of the data that will be used for this analysis plan exists, and no part will be generated until after this plan is registered.

---

## Research Design

### Study type
Simulation study: Using synthetic data to assess performance of analytical protocols, characterise measurement precision boundaries, and validate predictions against empirical databases. This includes demonstration of methods, modelling, and prediction of protocol performance across parameter space.

### Intention for causal interpretation
No causal relationship inferred: This study is not intended to inform a causal relationship.

### Blinding
No blinding is involved.

---

## Sampling

### Data collection procedures

This study uses fully synthetic data generated in MATLAB. Elliptical trajectories with known ground truth β and velocity gain factor (VGF) are admixed with 1/f^α noise spanning instrumental white noise through biological pink noise, through red and black noise, using the fractional differencing approach of Xu (2019). The parameter space crosses five simulation factors (sampling rate, generative β, generated VGF, noise colour, noise magnitude) with six analytical pipelines, yielding 14,764,680 total configurations. A hierarchical tractability framework (Levels 1-9) manages computational feasibility.

Full parameter specifications and generation procedures: attached pre-registration document, Section 4.

### Sample size

The primary unit of analysis is the unique parameter combination (fs × β_gen × VGF × α × σ), of which there are 3 × 21 × 14 × 31 × 18 = 494,604. Each combination is processed by all 6 pipelines with 5 repetitions, yielding 14,764,680 individual trajectory analyses. Sample size is fully determined by the factorial design.

### Sample size rationale

Five repetitions per condition balance statistical reliability for SEM estimation against computational tractability (~200 CPU-hours at Level 9). Scaling to ten or twenty repetitions would exhaust available HPC allocation without increasing parameter space coverage, where the primary scientific contribution lies. The adequacy criterion SEM < 0.011 derives from the MDC framework: given |Δβ| = 0.03 as the smallest clinically meaningful divergence (Cook et al., 2026; Fourie et al., 2024), SEM < MDC/2.77 at 95% confidence yields this threshold. Where within-condition SEM estimates prove unstable across five repetitions, that instability is treated as a substantive result indicating sensitivity to individual noise realisations.

Full justification: attached pre-registration document, Section 4.1.

### Starting and stopping rules

This study generates simulated data programmatically; all 14,764,680 configurations will be generated in full. Prior to this, the LMM framework has been explored with ersatz test data (randomly generated response values populating the full factorial table structure, without running the simulation pipeline) to demonstrate feasibility of convergence at scale. If the LMM does not converge on the simulated data, we propose employing smaller subsets of the data. These tractability Levels, 1 to 9, at Level 9 (the complete parameter space), successively smaller subsets (Levels 8 through 1) provide fallback positions until convergence is achieved. The highest completed tractability level will be reported. No data will be excluded post-hoc on the basis of results.

---

## Variables

### Manipulated variables

Seven variables are manipulated in a fully crossed factorial design. Five are simulation parameters: generative β (21 levels, 0 to 2/3), velocity gain factor (14 levels), noise colour α (31 levels, 0 to 3.0), noise magnitude σ (18 levels, 0 to 10 mm), and sampling rate (60, 120, 240 Hz). Two are analytical pipeline factors: kinematic derivation method (Low Pass Butterworth Filter + finite differences differentiation vs Savitzky-Golay smoothing differential filter) crossed with regression method (LMLS, IRLS and log linear OLS), yielding six pipelines applied to every combination with five repetitions.

Full variable specifications: attached pre-registration document, Sections 3 and 4.

### Measured variables

Two primary outcomes per pipeline per parameter combination: recovered β (β_rec), entering models as β_bias = β_gen − β_rec; and recovered VGF (VGF_rec), entering as VGF_bias = VGF_gen − VGF_rec conditional on adequate β recovery. In the empirical validation stage (Stage 5), where ground truth is unavailable, the primary outcomes shift to the 15 pairwise pipeline differences per trial (Δβ_observed), compared against simulation-derived predictions.

Full variable definitions: attached pre-registration document, Sections 3 and 7.

### Indices

Three derived indices are central: β_bias (β_gen − β_rec) as the LMM response variable, with within-cell SD giving SEM; the SEM adequacy classification (adequate < 0.011, marginal 0.011-0.020, inadequate > 0.020); and Δβ_PQ, the pairwise pipeline difference used in Stage 5 constellation comparisons, quantified via Lin's Concordance Correlation Coefficient (CCC), MAE, and supplementary RSA.

Full formulae and implementation details: attached pre-registration document, Sections 5 and 7.

---

## Analysis Plan

### Statistical models

The study employs a four-stage progressive modelling framework using Linear Mixed-effects Models (LMM) fitted via MATLAB's fitlme.

Stage 1 fits a seven-way fully factorial LMM to the complete simulation output, with β_bias as the response and (1|combination) random intercepts. The analysis strategy is interaction-first, with simple effects as the primary reporting unit. Stage 2 applies diagnostics (bootstrap CIs, residual pattern analysis, cross-validation) to identify parameter regions where the global model fails. Stage 3 develops region-specific conditional LMMs for those regions. Stage 4 synthesises global and conditional models into a unified decision framework mapping each (α, σ, fs) coordinate to recommended pipelines with SEM bounds.

A subsequent empirical validation phase then tests the framework: all six pipelines are applied to trials from seven movement databases, trial-level noise profiles are estimated via multitaper spectral analysis (pmtm), and the 15-element predicted and observed pipeline constellation vectors are compared using CCC (McBride, 2005 thresholds), MAE, and supplementary Mantel tests.

Full model specifications, interaction decomposition, and validation logic: attached pre-registration document, Sections 5-7.

### Transformations

The primary response β_bias requires no transformation; it is continuous and directly interpretable in the metric of the MDC framework. For the OLS pipeline, log-transformation of velocity and curvature is part of that pipeline's definition, not a pre-processing step applied to LMM inputs. Noise colour (α) and noise magnitude (σ) enter models on their natural scales; if Stage 1 residual diagnostics reveal non-linearity, predictor transformations will be considered and documented. Categorical factors (derivation method, regression method, sampling rate) enter as nominal predictors.

### Inference criteria

The primary criterion is SEM < 0.011, derived from MDC/2.77 where MDC = 0.03. This is a pre-specified adequacy criterion, not a significance threshold. For LMM fixed effects, the analysis is interaction-first: omnibus interaction terms are tested before main effects. Stage 5 validation uses CCC thresholds from McBride (2005): excellent > 0.90, good 0.70-0.90, moderate 0.50-0.70, failure < 0.50, with MAE < 0.03 as a co-criterion. Bootstrap CIs in Stage 2 flag coefficients as unstable when spanning ±0.03.

### Data inclusion and exclusion

All 14,764,680 parameter combinations are included by design. Computational failures are logged and reported as a substantive result. In Stage 5, trials whose estimated noise profile falls outside the simulation parameter space (α > 3 or σ > 10 mm) are reported separately rather than silently excluded. No subject-level exclusions apply; for empirical databases, original published exclusion criteria govern data quality.

### Missing data

Synthetic data is complete by construction. Computational failures producing missing β_rec values are treated as convergence failures and excluded from LMM fitting only for the affected pipeline. For empirical databases, missing trials reflect original dataset structure and are not imputed; systematic missingness patterns are noted.

### Other planned analysis

Bias correction framework: for pipelines with systematic, invertible bias, trial-level corrected estimates are computed conditional on Stage 5 validation success (Section 7.7). VGF recovery analysis parallels the primary β model but is fitted only where β SEM < 0.011. RSA structural validation via Mantel tests is reported as supplementary material for all databases.

Full details: attached pre-registration document, Sections 5-7.

---

## Other

### Context and additional information

The attached pre-registration document (prereg_v101.docx) is the authoritative and complete specification for this study. It contains full theoretical framework, glossary, parameter space definitions, model specifications, empirical validation logic, expected deliverables, timeline, and open science commitments. The OSF form fields above provide an overview; for implementation-level detail, please consult the attached document.
