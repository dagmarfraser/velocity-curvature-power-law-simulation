# Velocity-Curvature Power Law Simulation and Analysis

Comprehensive framework for testing the velocity-curvature power law in biological motion analysis, implementing the vetted and legacy protocols from Fraser et al. (2025) "Biological kinematics: a detailed review of the velocity-curvature power law calculation" (Experimental Brain Research).

This repository implements a sequential two-phase analytical pipeline that generates comprehensive ground truth data and subsequently performs systematic model adequacy assessment across 14,764,680 parameter configurations, providing the first evidence-based protocol selection framework for velocity-curvature power law analysis.

---

## Documentation Hierarchy and Current Status

This repository contains multiple documentation files serving different purposes. Understanding their status and relevance helps navigate the project effectively.

### Primary Documentation
- **README.md** (this file): Central documentation hub providing comprehensive project overview, installation instructions, and current methodologies
- **revised\_prereg\_v071.md**: **CURRENT** - Pre-registration with complete 4-stage progressive analysis methodology

Project Hosted https://osf.io/dwxa2/?view_only=193d498afc6441e990a7ebec77e16f1f
and on https://github.com/dagmarfraser/velocity-curvature-power-law-simulation
and eventually UBIRA eData (Large Datasets)

Multiverse Database: [DOI when available] - 5.2GB SQLite database with complete parameter space results - Generated using Toolchaincallerv057.m - 14M+ configurations across biological parameter ranges---

## Core Workflow Architecture

### Analytical Pipeline

The research framework follows a sequential two-phase approach where ground truth data generation enables comprehensive model adequacy assessment.

**Phase 1: Ground Truth Database Generation**
**Primary Tool:** `Toolchain_caller_v057.m`
**Purpose:** Generate comprehensive ground truth database of synthetic trajectories with known parameters
**Output:** Complete parameter space database (14,764,680 configurations) with ground truth β values for comparison against recovered β values, enabling definitive assessment of protocol performance across biologically plausible measurement conditions
**Computational Requirements:** 200+ hours on BlueBEAR HPC (72-core node)

**Phase 2: Model Adequacy Assessment** 
**Primary Tool:** `ModelAdequacy_Master_v002.m`
**Purpose:** Systematic evaluation of analytical protocol performance to answer the primary research question: Which analytical protocols enable reliable detection of the |Δβ| = 0.03-0.04 deviations documented in autism research (Cook et al., 2023; Fourie et al., 2024), and how do protocol-parameter interactions determine measurement precision?
**Validation Status:** Framework extensively tested with synthetic data prior to ground truth database analysis (Weeks 1-2)
**Scope:** Four-Stage Progressive Analysis Framework with evidence-based conditional modeling
**Clinical Translation:** Quantified adequacy boundaries for robust method selection, essential for detecting small-effect kinematic markers in neurodevelopmental research

---

## Phase 1: Ground Truth Database Generation

### Ground Truth Parameter Space

Toolchain_caller_v057 generates comprehensive ground truth ellipse trajectories across the complete empirical parameter space documented by Huh and Sejnowski (2015). The system creates systematic coverage of biologically plausible measurement conditions with known ground truth values:

**Parameter Space Specification (14,764,680 total configurations):**

**Sampling Rates:** 60Hz (Android tablets), 120Hz (WACOM devices), 240Hz (iPad Pro) spanning consumer through professional motion capture spectrum

**Generative β Values:** 21 values from 0 to 2/3, covering the empirical range observed by Huh and Sejnowski (2015) for shapes with angular frequency φ ranging from 0 to 6

**Velocity Gain Factors (VGF):** 14 exponentially spaced values (exp(4.5:0.1:5.8)) corresponding to ellipse tracing frequencies of approximately 0.5Hz to 2Hz, straddling the empirical threshold where power law recovery becomes more consistent

**Noise Characteristics:** 31 noise colors (1/f^α where α ranges from 0 to 3.0) spanning white noise of instrumentation (α=0) through pink noise of biological motion (α≈1) to black noise (α=3), accounting for spectral transformations after differentiation

**Noise Magnitudes:** 18 levels from instrument precision (<0.1mm) through average human motor variability (<4mm based on Madirolas et al., 2022) to challenging measurement degradation (10mm)

**Protocol Variations:** Legacy (Butterworth + finite differences) versus Vetted (Savitzky-Golay) filtering approaches with Linear, Levenberg-Marquardt Least Squares (LMLS), and Iteratively Reweighted Least Squares (IRLS) regression methods

**Statistical Reliability:** 5 repetitions per parameter combination ensuring robust precision assessment across the complete parameter space

**Key Configuration Switches** (edit at top of Toolchain_caller_v057.m):

**Debug Level Control:**
```matlab
debug = 0;  % Full production run (default for complete analysis)
debug = 1;  % Minimal debug run (small subset, no parallel computing)
debug = 2;  % Comprehensive debug run (rebuild shapes, no parallel computing) 
debug = 3;  % Parallel debug run (tests parallel functionality)
```

**Parameter Generation Method:**
```matlab
useFastBatch = true;   % Fast streaming batch method (recommended, ~100x faster)
useFastBatch = false;  % Original individual INSERT method (slow but verified)
```

**Version 057 Features:**
- **Optimized Checkpoint Interval**: 200 configurations per checkpoint for high-core efficiency, reducing database write contention
- **Just-in-time Parallel Pool Creation**: Pool created only when needed after config generation, preventing timeout errors
- **HPC Performance**: Eliminates SQLite serialization bottlenecks on high-core-count systems
- **MATLAB 2022a HPC Compatibility**: All internal functions positioned for optimal HPC performance
- **Fault Tolerance**: Checkpoint-based recovery with improved reliability

**Requirements:**
- Database Toolbox
- Parallel Computing Toolbox  
- Curve Fitting Toolbox
- Statistics and Machine Learning Toolbox

**Automatic Features:**
- **Job Resumption**: Automatically detects and offers to resume incomplete jobs
- **Environment Detection**: Automatically configures for HPC (ProcessPool) vs Local (ThreadPool) environments
- **Resource Monitoring**: Built-in system resource monitoring with warnings
- **Method Fallback**: Automatically falls back to verified method if fast batch fails

### Phase 1 Implementation Architecture

The ground truth database generation implements a hierarchical structure optimized for massive parameter space exploration.

**Toolchain_caller_v057.m** serves as the master orchestrator managing parallel computation across the complete parameter space. This system initializes SQLite database infrastructure for result storage, defines the comprehensive parameter space encompassing shapes, noise types, magnitudes, filters, and regression methods, creates optimized parallel pools for high-performance computing environments, generates parameter configurations using fast batch methods, and orchestrates parallel execution with checkpoint-based fault tolerance.

**Toolchain_func_v032.m** functions as the core processing engine called for each parameter configuration. This component generates synthetic trajectories using pureCurveGenerator, applies controlled noise via generateCustomNoise_v003, implements filtering methods including Butterworth and Savitzky-Golay variants, calculates kinematics using differentiateKinematicsEBR, performs regression analysis via regressDataEBR, and returns beta values, velocity gain factors, and comprehensive error metrics.

---

## Phase 2: Four-Stage Progressive Analysis Framework

### Comprehensive Model Adequacy Assessment

The framework systematically evaluates model adequacy through progressive analysis of the ground truth database generated in Phase 1, implementing the methodology validated in Weeks 1-2 of the research timeline.

**Stage 1: Global Interaction Modelling**
Fits comprehensive Linear Mixed Effects Model capturing all parameter interactions across the complete parameter space: `δβ ~ βgenerated × VGF × samplingRate × filterType × regressionType × noiseMagnitude × noiseColour + (1|parameterCombination)` where δβ represents the difference between generated and recovered beta values. This global model quantifies how each parameter and their interactions systematically affect beta recovery accuracy across approximately 192 coefficients.

**Stage 2: Systematic Adequacy Assessment**
Determines whether the global model adequately captures parameter recovery patterns through three empirical criteria. Coefficient Stability assessment identifies bootstrap confidence intervals exceeding ±0.03 (the clinical significance threshold from autism research). Residual Pattern Analysis employs Cohen’s d > 0.5 threshold for systematic deviations revealing parameter regions where the global model systematically under- or over-predicts beta recovery accuracy. Prediction Accuracy Assessment uses changepoint analysis of cross-validation R² distributions to identify inflection points where model performance degrades significantly.

**Stage 3: Conditional Parameter Analysis**
Develops specialized models for parameter regions where the global model demonstrates inadequate performance, addressing specific sources of model failure identified in Stage 2. This targeted approach isolates problematic parameter combinations and develops region-specific models that account for unique interaction patterns within those measurement contexts. Each conditional model undergoes rigorous validation through regional cross-validation to ensure statistically significant improvement over the global approach.

**Stage 4: Integrated Assessment Framework**
Synthesizes the global and conditional models into a unified decision framework that optimizes analytical approach selection based on specific experimental conditions. This integration provides researchers with evidence-based guidance for choosing between global model predictions and conditional model recommendations, generating precision-stratified protocol recommendations with quantified confidence bounds.

### ModelAdequacy_Master_v002 Usage Framework
```matlab
% Complete 4-stage progressive analysis framework
ModelAdequacy_Master_v002()                    % Default: Level 2, small subset analysis
ModelAdequacy_Master_v002(2, 5, false)        % Level 2, 5 obs/combo, synthetic test data
ModelAdequacy_Master_v002(2, 5, true)         % Level 2, ground truth database analysis
ModelAdequacy_Master_v002(2, 5, false, 2)     % Resume from Stage 2 (crash recovery)

% Tractability levels (computational feasibility versus comprehensiveness):
% Level 1: Conservative (28.8K obs, 99.0% reduction) - Initial testing
% Level 2: Focused (472K obs, 96.8% reduction) - Framework validation
% Level 7-9: Comprehensive (4.4M-14.8M obs) - Complete analysis
```

### Expected Outcomes and Research Translation

Based on Fraser et al. (2025), the investigation anticipates systematic protocol-parameter interactions with vetted protocols maintaining precision across broader parameter ranges while legacy approaches demonstrate progressive degradation with increasing noise magnitude and decreasing noise colour. Global modelling will characterise these performance boundaries, with conditional analysis providing targeted optimisation for challenging measurement regimes.

**Anticipated Protocol Performance Patterns:**

For white noise at magnitudes exceeding instrument precision (>0.1mm), progressive beta estimation degradation is expected. Pink noise beyond average human error (>3.3mm) will establish reliability boundaries for both protocols. Red noise in challenging ranges (>10mm) will demonstrate limited protocol utility, while black noise may show adequacy across all magnitudes due to spectral characteristics. These patterns will manifest as improved noise magnitude tolerance and enhanced noise colour resilience in vetted protocols compared to legacy methods.

**Methodological Decision Support Deliverables:**

Protocol selection flowchart based on experimental conditions will enable evidence-based analytical approach selection. Performance boundary maps will delineate protocol operational ranges across the complete parameter space. Precision lookup tables for common paradigms will provide immediate guidance for standard measurement contexts. Uncertainty quantification tools will support individual assessments with confidence bounds appropriate for detecting small-effect kinematic markers essential to neurodevelopmental research applications.

---

## Model Adequacy Assessment Results Guide

### Primary Adequacy Criteria (prereg_v071)
**Residual Pattern Threshold**: Cohen's d > 0.5 indicates systematic bias requiring conditional analysis
**Coefficient Instability**: Bootstrap CI > ±0.03 (clinical significance from Cook et al. 2023)
**Cross-Validation Degradation**: >15% performance loss indicates poor regional fit
**Minimum Region Size**: n ≥ 200 for reliable conditional modeling

### Framework Execution Outcomes
**Stage 1 → Stage 2 Assessment**: Complete 192-coefficient structure successfully fitted, demonstrating computational tractability of massive linear mixed effects models with robust optimization achieved across full parameter space.

**Stage 2 → Stage 3 Decision**: Framework successfully detects both adequate and inadequate scenarios through systematic adequacy assessment, with synthetic testing validating methodology using artificially generated poor fit regions.

**Stage 3 → Stage 4 Integration**: Region-specific modeling with enhanced performance metrics, systematic evaluation of hierarchical integration benefits, and synthetic enhancement methodology testing with artificial integration-worthy regions.

**Framework Status Interpretation**: Complete 4-stage framework execution providing evidence-based method selection with quantified adequacy boundaries.

---

## Data Sources and Open Science Commitment

**Ground Truth Database**: Generated through Toolchain_caller_v057 requiring 200+ computational hours on 72-core HPC systems with checkpoint interval optimization, implementing the complete 14,764,680 parameter configuration space.

**Framework Validation Status**: Weeks 1-2 - Four-Stage Progressive Analysis Framework validated with synthetic data, demonstrating computational feasibility across tractability levels 1-9 prior to ground truth database analysis.

**Open Science Materials**: All materials are publicly available via OSF (CC BY 4.0) including complete MATLAB simulation toolchain, simulation database with characterisation metadata, precision stratification modules, and methodological decision support algorithms. The repository implements reproducible science principles with comprehensive documentation and version control.

### System Requirements
- MATLAB (R2020b or later recommended)
- Statistics and Machine Learning Toolbox
- Curve Fitting Toolbox
- Database Toolbox
- Parallel Computing Toolbox
- PsychoPhysics Toolbox (PTB) - required only for generating baseline shape files

### File Structure
- **/src**: Complete Model Adequacy Framework (all 4 stages), and Toolchain_caller implementations
- **/functions**: Core analysis functions and standalone tools
- **/data**: Raw empirical datasets and processed results
- **/results**: Analysis outputs and databases
- Main scripts and enhanced diagnostics in root directory

---

## Getting Started: Complete Sequential Pipeline

### Prerequisites: Framework Testing and Validation
```matlab
% Test the model adequacy framework with synthetic data before ground truth analysis
ModelAdequacy_Master_v002(2, 5, false)  % Level 2, 5 obs/combo, synthetic test data
```

### Phase 1: Ground Truth Database Generation
```matlab
% Generate comprehensive ground truth database (requires substantial HPC resources)
Toolchain_caller_v057;     % Generates complete parameter space database
```

### Phase 2: Model Adequacy Assessment of Ground Truth Database
```matlab
% Execute complete 4-stage model adequacy assessment on ground truth data
ModelAdequacy_Master_v002(9, 5, true)   % Level 9 comprehensive analysis, ground truth database
```

### Framework Validation Notes
The model adequacy framework has been extensively tested using synthetic data to validate methodology before application to the computationally expensive ground truth database. This testing approach ensures robust analytical procedures prior to committing substantial computational resources to the complete pipeline execution.

## License

MIT License - see LICENSE file for details.

## Contact

Dagmar Scott Fraser - d.s.fraser@bham.ac.uk

---
