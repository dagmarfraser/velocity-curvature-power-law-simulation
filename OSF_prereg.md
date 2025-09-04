# OSF Pre-Registration: Systematic Characterisation of Velocity-Curvature Power Law Analysis Protocols Across Parameter Space

**Authors:** Fraser, D. S., Di Luca, M., & Cook, J. L.  
**Institution:** University of Birmingham, School of Psychology

**OSF:** <https://osf.io/dwxa2/>  
**GitHub:** <https://github.com/dagmarfraser/velocity-curvature-power-law-simulation>  
**Data Repository:** <https://edata.bham.ac.uk>

## Research Questions or Hypotheses

The power law established by Lacquaniti et al (1983), relating tangential velocity to curvature (v ∝ κ^(-β), where β ≈ 1/3), potentially encodes a fundamental principle governing biological movement. Fraser et al. (2025) demonstrated that legacy analytical protocols systematically bias towards β ≈ 1/3, potentially distorting measurements of the underlying generative mechanism. Fraser and colleagues subsequently introduced a vetted protocol designed to eliminate this bias. However, neither legacy nor vetted protocols have been characterised across the comprehensive range of biologically plausible parameters. This methodological lacuna undermines potential diagnostic applications of the power law, particularly where subtle divergences of β (|Δβ| = 0.03–0.04) distinguish between autistic and neurotypical populations (Cook et al., 2023; Fourie et al., 2024).

**Primary Research Question:** Which velocity-curvature power law analytical protocols - legacy or vetted of Fraser et al. (2025) - enable reliable detection of β, and |Δβ| = 0.03–0.04 divergences documented in preliminary autism research. How do protocol-parameter interactions determine β calculation precision across experimental methodology choices and biologically informed parameter combinations?

## Factors in the Data-Generating Mechanism

This systematic parameter space exploration employs seven key factors that span experimental methodology with analytical variations, and biologically informed parameter combinations:

**Sampling Rates**: Consumer tablet through professional motion capture (60, 120, 240 Hz)

**Generative β Values**: Empirical range observed in Huh and Sejnowski (2015);from 0 to 6 (21 levels: 0 to 2/3 in increments of 2/3 / 20)

**Velocity Gain Factors**: Corresponding to value swhich give rise ellipse tracing frequencies, of approximately 0.5–2 Hz (14 levels: exp(4.5:0.1:5.8)) 

**Noise Colour Characteristics**: White through black noise; accounting for white instrumental noise, pink 'biological' noise, and differentiation transformations whereby red and black noise become white and pink respectively (31 levels: 1/f^α where α ranges 0:0.1:3.0)

**Noise Magnitudes**: Instrument precision (sub-millimetre) through average human error to challenging measurement conditions (18 levels: 0:0.025:0.1, 0.25:0.25:2.25, 4, 6, 8, 10 mm)

**Filter Types**: Legacy versus vetted approaches (2 levels: Butterworth + finite differences vs. Savitzky-Golay smoothing differentiators)

**Regression Methods**: Traditional versus robust estimation approaches (3 levels: Linear, Levenberg-Marquardt Keast Squares (LMLS), Iteratively Weighted Least Squares (IRLS) regression)

## Factor Level Combinations and Number of Conditions

**Total parameter combinations**: 14,764,680 configurations including 5 repetitions per condition. This comprehensive factorial design systematically examines all possible combinations across the seven-factor space to ensure complete characterisation of protocol performance boundaries.

## Estimands or Targets of Interest

**Primary Estimand**: Δβ = β_generated - β_recovered, representing the deviation between ground truth and recovered power law exponents.

**Precision Boundaries**: ±0.03 precision threshold based on autism research findings, with additional empirically derived thresholds through bootstrap confidence interval analysis.

**Secondary Estimand**: ΔVGF = VGF_generated - VGF_recovered, representing the deviation between ground truth and recovered velocity gain factors. Whilst the literature establishes that VGF should remain constant within movement segments (Lacquaniti et al., 1983), a specifc range is not established.

**Protocol Adequacy Metrics**: Coefficient stability, residual pattern analysis, and prediction accuracy assessment across parameter regions.

## Included Methods and Extracted Quantities

**Four-Stage Progressive Analysis Framework** (necessitated by anticipated heterogeneous protocol performance across parameter space):

**Stage 1 - Global Interaction Modelling**: Comprehensive Linear Mixed Effects Models capturing all parameter interactions: 


- Primary model: δβ ~ β_generated × VGF_generated × sampling_rate × filter_type × regression_type × noise_magnitude × noise_colour + (1|parameter combination)


- Secondary model (conditional on adequate β performance): δVGF ~ β_generated × VGF_generated × sampling_rate × filter_type × regression_type × noise_magnitude × noise_colour + (1|parameter combination)

**Stage 2 - Systematic Failure Detection**: Adequacy assessment designed to identify parameter regions where the global model break down due to heterogeneous noise effects, protocol-specific breakdowns, or nonlinear interactions. Assessment employs: (a) coefficient stability via bootstrap confidence intervals, (b) residual pattern analysis (Cohen's d > 0.5 for systematic deviations), and (c) changepoint analysis of cross-validation R² distributions to detect performance degradation.

**Stage 3 - Regional Specialisation**: Conditional analysis triggered by Stage 2 observations, developing specialised models for parameter regions where the global approach demonstrates inadequate performance. This stage addresses the anticipated reality that legacy and vetted protocols will exhibit different results across the noise spectrum and sampling rate combinations.

**Stage 4 - Unified Decision Framework**: Integration of global and conditional models into a coherent decision support system that optimises protocol and methodolgy selection based on experimental conditions, accounting for the heterogeneous performance landscape identified in Stages 1-3.

**Rationale**: The 4-stage approach is necessitated by our anticipation that no analytical model will adequately describe protocol performance across the entire 14.7M parameter space. High noise conditions, specific filter-regression combinations, and variant sampling rates are expected to create distinct performance regimes requiring specialised treatment.

## Performance Measures

**Primary Performance Criteria**:

- β recovery accuracy within ±0.03 for optimal conditions (white noise 0–0.1mm standard deviation), indicating that β can be recovered with commercially available instrumentation in the absence of additional errors such as physiological noise
- Noise resilience across physiological ranges (pink through black noise spectra), indicating that β can be recovered from standard experimental tasks and reflect the underlying generator  
- Computational stability across all parameter combinations, indicating the regions within parameter space where β can be recovered faithfully, thereby guiding future research in experimental and analytical methodology

**Secondary Performance Criteria** (conditional on β performance within ±0.03):

- VGF recovery consistency across parameter combinations, assessing whether the velocity gain factor remains stable across repetitions as theoretically predicted, or exhibits systematic variation across samping rates and noise conditions.

**Model Adequacy Expectations**:

- **Expected Global Model Success**: Coefficient stability (bootstrap CI ≤ ±0.03) and residual patterns (Cohen's d ≤ 0.5) in optimal conditions (low noise, adequate sampling rates)
- **Anticipated Regional Failures**: Global model inadequacy expected in high noise conditions (>2mm), extreme sampling rates (<100Hz), and specific legacy protocol combinations, necessitating Stage 3 conditional modelling
- **Protocol Differentiation**: Legacy and vetted protocols expected to exhibit different failure boundaries, with vetted protocols maintaining adequacy across broader parameter ranges

## Number of Simulation Repetitions per Condition

**Five repetitions per parameter combination** to ensure statistical reliability for precision assessment whislt balancing computational tractablity. Totalling 14,764,680 individual simulations across the complete parameter space. 

## Handling Stochastic Noise and Simulation Variability

**Replication Strategy**: Statistical reliability ensured through five repetitions per parameter combination with random intercepts for parameter combinations in mixed-effects models to account for between-replication variability introduced by stochastic noise injection.

**Bootstrap Procedures**: Bootstrap confidence interval analysis (minimum 1000 bootstraps) employed for threshold determination and uncertainty quantification of model parameters.

**Noise Control**: Simulated trajectories generated and contamination with noise (Gaussian white noise and 1/f^α coloured noise) with known characteristics, enabling separation of method-induced bias from noise-related uncertainty.

**Reproducibility**: Fixed random seeds, documented in crash-safe toolchain implementation, ensure exact replication of noise and downstream results.


## Statistical Software and Packages

**MATLAB 2023B** with Statistics and Machine Learning Toolbox for mixed-effects modelling, bootstrap analysis, and comprehensive simulation framework implementation.

## Planned Computational Environment

BlueBEAR High-performance computing infrastructure, Code capable of handling 200+ compute hours with crash-safe simulated groundtruth generation implementation, and synthetic data validated 4-stage Progressive Analysis Framework ensuring computational feasibility and convergence.  Tractability levels to permit subset evalaution on desktiop hardware, less powerful HPC clusters, or in the event of non convergence of Linear Mixed Effect Models.

## Reproducibility Measures

The complete MATLAB simulation toolchain, with characterisation metadata will be publicly available via OSF under MIT licence. All simulation parameters, ground truth generation procedures, and analysis protocols are documented with accompanying precision stratification modules and methodological decision support algorithms. Fixed random seed management ensures exact reproducibility of all stochastic elements.

---

## Extended Background and Rationale

### Theoretical Foundation

The velocity-curvature power law, codified by Lacquaniti et al. (1983), describes a remarkably consistent relationship governing biological movement across scales and species. This kinematic invariant relates tangential velocity to trajectory curvature through a power law:

*v(t) = VGF × κ(t)^(-β)*

where *v* represents tangential velocity, *κ* denotes curvature, *VGF* is the velocity gain factor (theoretically constant within movement segments), and β typically approximates 1/3. The VGF represents the baseline movement speed and should remain stable across different curvature values within a single movement trajectory, making it a critical parameter for validating power law calculations. This "one-third power law" emerges across diverse motor behaviours, from eye saccades to drawing movements to bodily trajectories. The law emerges in non-human species, suggesting a fundamental computational principle underlying motor control (Flash, 2021).

### Methodological Stakes: Small Effects with Large Implications

Recent autism research has elevated the urgency of characterising the methodological precision of power law analyses. Cook et al. (2023) and Fourie et al. (2024) independently documented that autistic individuals demonstrate subtle but consistent deviations from canonical power law values:

- **Fourie et al. (2024):** Autistic children demonstrated β = 0.36 versus neurotypical β = 0.33 (|Δβ| = 0.03, Cohen's d = 0.46)
- **Cook et al. (2023):** Adult autistic populations exhibited similar effect sizes (|Δβ| = 0.03–0.04)

This β exponent can purportedly be extracted from any curved movement, regardless of scale or duration. Therefore, these numerically small divergences offer potential non-invasive kinematic markers for neurodevelopmental conditions, though their broader significance requires systematic validation. However, detecting such subtle effects demands analytical precision. Characterising precision within legacy or vetted protocols across the full range of biologically plausible parameters constitutes the aim of this study.

### Ground Truth Parameter Space Justification

The comprehensive characterisation of analytical protocols requires simulated elliptical trajectories with known β and VGF ground truth values, contaminated with noise. This ground truth approach enables definitive assessment of protocol performance.

| **Parameter** | **Values** | **Count** | **Methodological Justification** |
|---------------|------------|-----------|-----------------------------------|
| Sampling Rates (fs) | [60, 120, 240] Hz | 3 | Consumer tablet through professional motion capture spectrum |
| Generative β | 0:(2/3)/20:(2/3) | 21 | Empirical range observed by Huh and Sejnowski (2015) for shapes with angular frequency φ [0...6] |
| VGF Values | exp(4.5:0.1:5.8) | 14 | Velocity gain factors corresponding to ellipse tracing frequencies of ~0.5–2 Hz |
| Noise Colour (1/f^α) | 0:0.1:3.0 | 31 | White (α=0) replicating Maoz et al. (2005), Pink (α≈1) physiological noise, through to Black (α=3) accounting for α transformations after differentiation |
| Noise Magnitudes | [0:0.025:0.1, 0.25:0.25:2.25, 4, 6, 8, 10] mm | 18 | Instrument precision <0.1 mm, regression-safe <2.25 mm, to challenging measurement conditions up to 10 mm |
| Filter Types | Legacy vs. Vetted | 2 | Butterworth + finite differences vs. Savitzky-Golay comparison |
| Regression Types | Linear, LMLS, IRLS | 3 | Legacy versus nonlinear regression methods |
| Parameter Combination | 5 repetitions | 5 | Statistical reliability for precision assessment |

### Expected Outcomes and Research Translation

Based on Fraser et al. (2025), we anticipate systematic protocol-parameter interactions, with vetted protocols maintaining precision across broader parameter ranges whilst legacy approaches demonstrate progressive degradation with increasing noise magnitude and decreasing noise colour (i.e. as α approaches zero). Global modelling will characterise these performance boundaries, with conditional analysis providing targeted optimisation for challenging measurement regimes.

**Deliverables include:**

- **Protocol selection flowchart** based on experimental conditions
- **Performance boundary maps** delineating protocol operational ranges
- **Precision lookup tables** for common paradigms
- **Uncertainty quantification tools** for individual assessments

This investigation provides the first systematic characterisation of velocity-curvature power law analysis protocols across comprehensive parameter space. By mapping where methodological approaches succeed and fail, we enable evidence-based protocol selection for detecting small-effect kinematic markers essential to research translation.
