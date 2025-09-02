# Pre-Registration: Systematic Characterisation of Velocity-Curvature Power Law Analysis Protocols Across Parameter Space

**Authors:** Fraser, D. S., Di Luca, M., & Cook, J. L.  
**Institution:** University of Birmingham, School of Psychology

NB OSF <https://osf.io/dwxa2/> github
<https://github.com/dagmarfraser/velocity-curvature-power-law-simulation>
UBIRA/edata <https://edata.bham.ac.uk> tested but not uploaded yet

## Abstract

The power law of Lacquaniti et al (1983), relating the tangential
velocity to curvature (*v* ∝ *κ*^(-β), where β ≈ 1/3) may encode a
fundamental principle governing biological movement. Fraser et al.
(2025) demonstrated that legacy protocols systematically bias towards β
≈ 1/3, potentially skewing measurement of the underlying generator.
Fraser and colleagues introduced a vetted protocol that avoids this
bias. However, neither legacy nor vetted protocol is characterised
across a wide range of biologically plausible parameters. This
methodological gap undermines potential diagnostic applications of the
power law where divergences of β, |Δβ| = 0.03-0.04, are calculated
between autistic and neurotypical populations (Cook et al., 2023; Fourie
et al., 2024).

We present a comprehensive characterisation framework, systematically
evaluating performance across 14.7 million parameter combinations; across
analysis protocols and kinematic parameters. Simulated data with known ground
truth, allows
Linear Mixed Effects Models to map the effective limits of
different analytical approaches, identifying optimal strategies for
specific measurement contexts.

This investigation establishes the first evidence-based protocol
selection framework for velocity-curvature power law analysis, providing
essential guidance for translating kinematic principles into reliable
kinematic assessments.

## 1. The Velocity-Curvature Power Law: From Fundamental Principle to Kinematic Assessment

### 1.1 Theoretical Foundation

The velocity-curvature power law, codified by Lacquaniti et al. (1983),
describes a remarkably consistent relationship governing biological
movement across scales and species. This kinematic invariant relates
tangential velocity to trajectory curvature through a power law:

*v(t) = VGF × κ(t)^(-β)*

where *v* represents tangential velocity, *κ* denotes curvature, *VGF*
is the constant velocity gain factor, and β typically approximates 1/3.
This "one-third power law" emerges across diverse motor behaviours;
from eye saccades, to drawing movements, to bodily trajectories. The law
emerges in non-human species suggesting a fundamental computational
principle underlying motor control (Flash, 2021).

### 1.2 Methodological Concerns in Power Law Analysis

Fraser et al. (2025) reported that conventional analytical protocols
contain systematic biases that artificially reduced deviations from β ≈
1/3 regardless of underlying signal characteristics. Their investigation
identified three critical failure modes:

1. **Differentiation artefacts:** Finite differences differentiation
   amplifies noise, creating artificial power law confirmations at
   higher magnitudes

2. **Filtering distortions:** Butterworth filters, ostensibly employed
   to address noise alter the underlying signal, creating spurious
   power law results of one-third

3. **Log-transform / linear regression complications:** Logarithmic
   projection, required for linear regression, fundamentally alters
   noise distributions, minimising deviations from one-third

The methodological implications of this uncertainty are profound: the β
range compression around one-third, introduced by legacy protocols may
mask genuine biological variation, at best rendering conservative the
observed divergences in existing literature, at worst explaining the
prior ubiquity of the one-third power law.

Fraser and colleagues demonstrated these biases for β = 0 (i.e. constant
velocity regardless of curvature), yielding spurious β ≈ 1/3 through the
legacy protocol. The authors in turn proposed a vetted protocol
employing Savitzky-Golay smoothing differential filters and non-linear
regression to minimise these distortions. This vetted protocol was demonstrated to avoid
spurious one-third confirmation. However yet a full characterisation of protocol
behaviour across a spectrum of β is required to assess if and where any
protocol perform adequately, not merely avoiding spurious 1/3 compliance.

### 1.3 Methodological Stakes: Small Effects with Large Implications

Recent autism research has elevated the urgency of characterising the
methodological precision of power law analyses. Cook et al. (2023) and
Fourie et al. (2024) independently documented that autistic individuals
demonstrate subtle but consistent deviations from canonical power law
values:

- **Fourie et al. (2024):** Autistic children showed β = 0.36 versus
  neurotypical β = 0.33 (|Δβ| = 0.03, Cohen's d = 0.46)

- **Cook et al. (2023):** Adult autistics populations exhibited similar effect
  sizes (|Δβ| = 0.03-0.04) 

This β exponent can purportedly be extracted from any
curved movement, regardless of scale or duration. Therefore, these
numerically small divergences offer potential non-invasive kinematic
markers for neurodevelopmental conditions, though their broader
significance requires systematic validation. However, detecting such
subtle effects demands analytical precision. Characterizing
precision within legacy or vetted protocols, across the
full range of biologically plausible parameters, is the aim of this
study.

## 2. Research Questions and Innovation

### 2.1 Primary Investigation

**Which analytical protocols enable reliable detection of the |Δβ| =
0.03-0.04 deviations documented in preliminary autism research, and how
do protocol-parameter interactions determine measurement precision?**

### 2.2 Methodological Framework

We introduce a comprehensive characterisation framework that:

- Maps protocol performance boundaries across parameter space through
  systematic evaluation

- Identifies optimal analytical strategies for specific measurement
  contexts

- Quantifies uncertainty boundaries for each protocol-parameter
  combination

- Provides evidence-based guidance for future protocol selection

## 3. Ground Truth Parameter Space

### 3.1 Comprehensive Parameter Space

The comprehensive characterisation of the analytical protocols requires
simulated elliptical trajectories with known β and VGF ground truth
values, contaminated with noise. This
ground truth approach enables a definitive assessment of protocol
performance.

Synthetic trajectories are generated, and in turn examined, with the
following parameters;

- Generative β values covering the empirical spectrum observed in biological movements by Huh and Sejnowski (2015)

- Velocity gain factors represent different movement speeds, from
  slower deliberate actions to rapid fluid movements; 0.5Hz to 2Hz ellipse trajectories, at β = 1/3. These values straddle a value under which
  power law may be more consistently recovered empirically i.e. when
  ellipse tracing approximates or exceeds ~1Hz (Matic and
  Gomez-Marin, 2023)

- Noise characteristics span from instrument precision (submillimetre)
  through average human motor variability for such tasks (<4mm;
  Madirolas et al., 2022) to high measurement degradation (10mm). For
  computational tractability this wide range is sparse above 2mm.
  Noise is added as injected as per Maoz et al., (2005)

- Noise colours range from the white noise of instrumentation, the
  presumed pink of biological motion.  Red and black noise are also considered, as after differentiation these noises are transformed into white and pink
  (e.g. coordinates to velocities). Noise generation employs the
  fractional differencing approach of Xu et al. (2022) for accurate
  1/f^α noise across the complete spectrum from white (α=0) through
  black (α=3)

- Sampling rates span widely available instrumented tablets, e.g. 60Hz
  Android tablets, 120Hz for WACOM devices, to iPad Pro with
  proprietary stylus at 240Hz

- Protocol variations compare established legacy approaches against
  the improved vetted methods introduced by Fraser et al. (2025). See
  attached table

| **Parameter** | **Values** | **Count** | **Methodological Justification** |
|---------------|------------|-----------|-----------------------------------|
| Sampling Rates (fs) | [60, 120, 240] Hz | 3 | Consumer tablet through professional motion capture spectrum |
| Generative β | 0:(2/3)/20:(2/3) | 21 | Empirical range observed by Huh and Sejnowski (2015) for shapes with angular frequency φ [0...6] |
| VGF Values | exp(4.5:0.1:5.8) | 14 | Velocity gain factors corresponding to ellipse tracing frequencies of ~0.5-2 Hz |
| Noise Colour (1/f^α) | 0:0.1:3.0 | 31 | White (α=0) replicating Maoz et al. (2005), Pink (α≈1) physiological noise, through to Black (α=3) accounting for α transformations after differentiation |
| Noise Magnitudes | [0:0.025:0.1, 0.25:0.25:2.25, 4, 6, 8, 10] mm | 18 | Instrument precision < 0.1mm, regression-safe < 2.25mm, to challenging measurement conditions up to 10mm |
| Filter Types | Legacy vs. Vetted | 2 | Butterworth + finite differences vs. Savitzky-Golay comparison |
| Regression Types | Linear, LMLS, IRLS | 3 | Legacy versus non-linear regression methods |
| Parameter Combination | 5 repetitions | 5 | Statistical reliability for precision assessment |

**Total configurations:** 14,764,680 (including 5 repetitions per condition)

### 3.2 Precision Stratification

Systematic assessment of the adequacy of protocols will employ the ±0.03
precision boundary based on preliminary findings from autism research
(Cook et al., 2023; Fourie et al., 2024). However it should be noted this precision
boundary derives from legacy protocol calculations which compress the
expressed range of β. Therefore, these thresholds potentially represent conservative
estimates i.e. true biological differences may be larger. Additional
precision boundaries will be empirically derived through bootstrap
confidence interval analysis and cross-validation performance assessment
across the 14.7M parameter combinations.

## 4. Validation Criteria

### 4.1 Parameter Space Validation

Prior to simulating ground truth trajectories, the computational tractability of such a large parameter space has been comprehensively validated through our 4-stage Model Adequacy Framework (ModelAdequacy_Master_v002). The seven-way factorial mixed-effects model with full interaction structure and random intercepts for parameter combinations has been demonstrated to converge with test data. However convergence of is not guaranteed with the ground truth simulations. In anticipation this framework also demonstrates computational feasibility across tractability levels (1-9) which employ progressively smaller subsets of the parameter space. The true subset framework ensures that all tractability levels 1-8 use identical parameter values as exact subsets of Level 9, maintaining methodological consistency. Additionally this enables scalable analysis across different computational budgets. 

### 4.2 Ground Truth Validation Criteria

This controlled ground truth generation approach provides the empirical foundation
necessary to translate the observations of Δβ in clinical populations
into a diagnostic tool.

- Verification of β recovery accuracy within ±0.03 for optimal
  conditions (white noise 0-0.1mm standard deviation), in which only
  instrument error is observed. Failure in this zone will suggest the
  method can never accurately recover β.

- Demonstration of noise resilience across physiological ranges (pink
  through black noise spectra). Failures in this zone will indicate
  that only biased estimate of β can be recovered for biological
  motion.

- Computational stability assessment across all parameter
  combinations. This will delineate the boundaries within which β can
  be recovered faithfully.


## 5. Four-Stage Progressive Analysis Framework

### Stage 1: Global Interaction Modelling

Fit comprehensive Linear Mixed Effects Model capturing all parameter
interactions across the complete parameter space, relating them to
divergences from the groundtruth generated beta:

δβ ~ β_generated × VGF × sampling_rate × filter_type × regression_type
× noise_magnitude × noise_colour + (1|parameter combination)

Where δβ is β_generated - β_recovered.

This global model quantifies how each parameter and their interactions
systematically affect β recovery accuracy. The random effect accounts
for the five simulation repetitions per parameter combination, whilst
the fixed effects reveal which measurement conditions introduce
systematic bias versus random error. The comprehensive 7-way interaction
structure identifies the specific parameter combinations where protocols
fail and maintain precision.

### Stage 2: Systematic Adequacy Assessment

Stage 2 determines whether the global model adequately captures
parameter recovery patterns or requires conditional analysis for
specific parameter regions. Three empirical criteria identify model
inadequacy:

**Coefficient Stability:** Bootstrap confidence intervals exceeding
±0.03 indicate unstable parameter estimates that cannot reliably detect
autism-relevant effect sizes (Cook et al., 2023; Fourie et al., 2024).

**Residual Pattern Analysis:** Cohen's d > 0.5 for systematic
deviations reveals parameter regions where the global model
systematically under- or over-predicts β recovery accuracy.

**Prediction Accuracy Assessment:** Changepoint analysis of
cross-validation R² distributions will allow identification of
inflection points where model performance degrades significantly, thus
indicating inadequate representation of underlying parameter
relationships.

Parameter regions failing any criterion proceed to Stage 3 conditional
analysis. Regions meeting all criteria confirm global model adequacy for
those measurement conditions.

### Stage 3: Conditional Analysis 

Stage 3 develops specialized models for parameter regions where the
global model demonstrates inadequate performance. This targeted approach
addresses the specific sources of model failure identified in Stage 2
rather than imposing a single analytical framework across the entire
parameter space.

Conditional analysis isolates problematic parameter combinations and
develops region-specific models that account for the unique interaction
patterns within those measurement contexts. For example, high-noise
conditions may require different analytical structures than low-noise
scenarios, or specific combinations of sampling rate and filter type may
exhibit non-linear relationships not captured by the global model.

Each conditional model undergoes rigorous validation through regional
cross-validation to ensure improved performance over the global
approach. Only conditional models demonstrating statistically
significant improvement proceed to Stage 4 integration.

### Stage 4: Integrated Assessment Framework

Stage 4 synthesizes the global and conditional models into a unified
decision framework that optimizes analytical approach selection based on
specific experimental conditions. This integration provides researchers
with evidence-based guidance for choosing between global model
predictions and conditional model recommendations depending on their
measurement parameters.

The framework generates precision-stratified protocol recommendations
that account for the uncertainty inherent in different measurement
contexts. These recommendations include quantified confidence bounds
that enable researchers to assess whether their experimental conditions
support reliable detection of the small effect sizes relevant to
neurodevelopmental research applications.

## 6. Expected Outcomes and Research Translation

### 6.1 Anticipated Results

Based on Fraser et al. (2025), we anticipate systematic
protocol-parameter interactions: vetted protocols maintaining precision
across broader parameter ranges, whilst legacy approaches demonstrate
progressive degradation with increasing noise magnitude and decreasing
noise colour (i.e. as α trends to 0). Global modelling will characterise
these performance boundaries, with conditional analysis providing
targeted optimisation for challenging measurement regimes.

Specifically, we expect:

- For white noise (α=0) at magnitudes exceeding instrument precision
  (>0.1mm), progressive β estimation degradation

- For pink noise (α≈1) beyond average human error (>4mm),
  reliability boundaries for both protocols

- For red noise (α=2) in challenging ranges (>10mm), limited protocol
  utility

- For black noise (α=3), potential adequacy across all magnitudes due
  to spectral characteristics

These patterns will manifest as improved noise magnitude tolerance and
lower noise colour resilience in vetted protocols compared to legacy
methods.

### 6.2 Methodological Decision Support

Deliverables include:

- **Protocol selection flowchart** based on experimental conditions

- **Performance boundary maps** delineating protocol operational
  ranges

- **Precision lookup tables** for common paradigms

- **Uncertainty quantification tools** for individual assessments

## 7. Significance and Future Directions

This investigation provides the first systematic characterisation of
velocity-curvature power law analysis protocols across comprehensive
parameter space. By mapping where methodological approaches succeed and
fail, we enable evidence-based protocol selection for detecting
small-effect kinematic markers essential to research translation.

Our framework extends beyond autism applications to any condition
affecting movement kinematics. The open-source implementation ensures
reproducibility whilst the parametric characterisation approach
establishes a methodological paradigm for validating complex biomarker
protocols under realistic measurement conditions.

## 8. Timeline

**Weeks 1-2:** Finalize framework validation ✓ COMPLETE - 4-stage Model Adequacy Framework validated (computational feasibility confirmed across tractability levels 1-9)  
**Weeks 3-5:** Execute 14.7M parameter simulations using validated crash-safe toolchain (210 compute hours)  
**Weeks 6-8:** Progressive precision analysis and threshold investigation  
**Week 9:** Integration and documentation of methodological translation framework

## 9. Open Science Commitment

All materials will be publicly available via OSF (CC BY 4.0):

- Complete MATLAB simulation toolchain

- Simulation database with characterisation metadata

- Precision stratification modules

- Methodological decision support algorithms

## References
