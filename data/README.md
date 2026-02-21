# Data sources

This folder documents the seven movement databases used for Stage 5 empirical validation. **No data are redistributed through this repository**, except potentially the Fraser pilot dataset collected by the first author. Everything else must be obtained from the original authors or publication supplements at the links below.

## Databases

**Fraser (pilot)** -- iPad Pro 240 Hz shape tracing via StimuliApp. Low-noise laboratory conditions where all six pipelines should converge. May be included directly in a future release.
Contact: D. S. Fraser, University of Birmingham.

**Zarandi et al. 2023** -- Practised ellipse drawing on WACOM at 100 Hz. Clean, well-controlled data from experienced participants.
Zarandi, Z., Stucchi, N. A., Fadiga, L., & Pozzo, T. (2023). The effect of the preferred hand on drawing movement. *Scientific Reports*, 13, 8264.
https://doi.org/10.1038/s41598-023-34861-x

**Dhieb et al. 2022** -- Tablet ellipse drawing across the lifespan (ages 19 to 85). Gives us biological variability without clinical group differences.
Dhieb, T., et al. (2022). Do individual characteristics influence the beta-elliptic modelling errors during ellipse drawing movements? *CMBBE*, 25(7), 783-793.
https://doi.org/10.1080/10255842.2021.1978434

**Cook et al. 2026** -- WACOM 133 Hz shape tracing, autistic vs neurotypical comparison. This is the dataset behind the 0.03 MDC threshold that motivates the whole project.
Cook, J. L., Fraser, D. S., Hickman, L. J., Brewer, R., & Huh, D. (2026). Autistic kinematics diverge from the power laws that typically govern movement.
https://doi.org/10.1101/2023.03.23.532745
Contact: J. Cook, University of Birmingham.

**Dagenais et al. 2021** -- 3D elephant trunk tip trajectories (Chishuru, African elephant) during reaching tasks. About as far from human tablet drawing as you can get while still obeying a power law, which is exactly the point.
Dagenais, P., Hensman, S., Haechler, V., & Milinkovitch, M. C. (2021). Elephants evolved strategies reducing the biomechanical complexity of their trunk. *Current Biology*, 31(21), 4727-4737.
https://doi.org/10.1016/j.cub.2021.08.029
Contact original authors for redistribution rights.

**James et al. 2020** -- Bumblebee locomotion trajectories. High-noise conditions that should push all six pipelines hard.
James, L., et al. (2020). Do bumblebees have signatures? *PLoS ONE*, 15(1), e0226393.
https://doi.org/10.1371/journal.pone.0226393

**Hickman et al. 2024** -- Parkinson's ON/OFF medication and haloperidol manipulation. Pharmacologically induced kinematic changes rather than developmental ones.
Hickman, L. J., et al. (2024). Dopaminergic manipulations affect the modulation and meta-modulation of movement speed. *BBR*, 474, 115213.
https://doi.org/10.1016/j.bbr.2024.115213
Contact: J. Cook, University of Birmingham.

## Processing (planned, not yet implemented)

Each dataset will go through coordinate extraction, timestamp verification, unit standardisation (mm, seconds), and quality checks. Noise colour (α) and magnitude (σ) will be estimated via multitaper spectral analysis (pmtm) to locate each trial within the simulation parameter space. The noise characterisation wrapper will be added alongside the Stage 5 validation scripts after the simulation run.

## Contact

For data access queries on Birmingham datasets:
- Jennifer Cook: j.l.cook@bham.ac.uk
- Dagmar Fraser: d.s.fraser@bham.ac.uk
