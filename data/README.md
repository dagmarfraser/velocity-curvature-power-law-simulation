# Data

This project analyses fourteen movement datasets. **No raw data is redistributed through this repository** — every subfolder below has its own `README.md` with the citation, DOI, and how to obtain the dataset yourself. Raw data files, extracted archives, and internal working notes are excluded from version control (see `.gitignore`); only the citation READMEs are tracked.

## Core empirical validation datasets (Stage 5 of the pre-registration)

| Dataset | Folder |
|---|---|
| Zarandi et al. (2023) | `zarandi/` |
| Cook et al. (2026) | `cook/` |
| Dhieb (RETED) | `dhieb_RETED/` |
| Hickman et al. (2024, 2026) | `hickman/` |
| Fraser | `fraser/` |
| Dagenais et al. (2021) | `dagenais/` |
| James et al. (2020) | `james/` |

`fraser/` supersedes an earlier internal pilot cohort, retained for provenance at `pilot/` but no longer part of the active analysis.

## Side-project datasets (Coherence Gap sub-project — not part of the primary pre-registered analysis)

| Dataset | Folder |
|---|---|
| ULTRA-MoCap | `ULTRA-MoCap/` |
| Santos et al. (2022) | `santos/` |
| Grouvel et al. (2023) | `grouvel/` |

## Exploratory, closed

| Dataset | Folder | Status |
|---|---|---|
| Peters, Wang & Limanowski | `peters/` | Closed pending a retest that never became available; not part of the primary analysis. Has its own participant-privacy restriction — see its README. |

## Downloaded, not currently adopted

| Dataset | Folder | Status |
|---|---|---|
| Niño-Tejada et al. | `nino-tejada/` | Candidate alongside ULTRA-MoCap; not adopted |
| Bonato (PhysioNet) | `bonato/` | Candidate exploratory addition; never formally decided |
| Fritschi et al. (2021) | `fritschi/` | Assessed and rejected as unsuitable |

These three are kept for now but contribute nothing to any current result — a decision on whether to remove them entirely is open (see `docs/TODO_GitHubReleasePrep_v002.md`).

## Contact

- Jennifer Cook, University of Birmingham — j.l.cook@bham.ac.uk
- Dagmar Fraser, University of Birmingham — d.s.fraser@bham.ac.uk
