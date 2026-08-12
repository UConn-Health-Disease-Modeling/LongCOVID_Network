# Network Drivers of Long COVID

Analysis code for:

> Jo Y., Jung J., Lee J. *Network Drivers of Long COVID: SARS-CoV-2
> Infection Reshapes Symptom Interdependencies.*

The study builds weighted symptom co-occurrence networks from a Korean
Long COVID cohort (n = 1,361; 787 infected, 574 uninfected; 29
symptoms) and compares infected and uninfected networks using weighted
density, node strength, weighted PageRank, Louvain modularity, and
strength assortativity, together with prevalence-preserving null
models, longitudinal (time-sliced) analyses, and an independent
validation with psychometric network methods (eLasso Ising estimation
and the Network Comparison Test).

## Repository contents

| File | Purpose |
|---|---|
| `LongCOVID_Network_Unified_v3_Revision.R` | Complete analysis pipeline (all figures and tables, including the revision analyses) |
| `run_synthetic_reproduction.R` | One-command reproduction of the full pipeline on the synthetic dataset |
| `data/synthetic_longcovid.xlsx` | Fully synthetic dataset mirroring the cohort's variable schema and per-group marginal symptom prevalences |
| `worked_example_public_data.R` | Worked example applying the framework to public data (NHANES 2017–2018 PHQ-9 depression symptoms) |

## Pipeline overview (`LongCOVID_Network_Unified_v3_Revision.R`)

- **4A** Baseline co-occurrence networks and differential (Δstrength /
  ΔPageRank) analysis
- **4B–4C** Propensity-score IPW and IPW-weighted networks
- **4D** Longitudinal time-sliced networks (0–3, 3–6, 6–12, 12–24 months)
- **4E** Bootstrap and permutation inference
- **4F–4G** Symptom × time-bin ΔPageRank quantification and slopes
- **5A** Prevalence-preserving null models (column permutation and
  fixed-margin/curveball randomization), density-gap decomposition,
  PageRank-reordering test, normalization sensitivity
- **5B** Node strength vs weighted PageRank concordance
- **5C** Longitudinal robustness: attrition table, individual-level
  (cluster) bootstrap, one-window-per-person sensitivity analysis,
  within-window null comparisons
- **5D** Psychometric validation: eLasso Ising networks (IsingFit) and
  the Network Comparison Test
- **5E** Synthetic dataset generation and verification

## Reproducibility

The original cohort data are governed by the Korea Disease Control and
Prevention Agency (KDCA) and cannot be redistributed. To allow the
pipeline to be inspected and tested end-to-end despite this
restriction, the repository provides:

**1. Synthetic reproduction.** `data/synthetic_longcovid.xlsx` mirrors
the variable schema and per-group marginal symptom prevalences of the
cohort; covariates are drawn independently from their per-group
marginal distributions and symptoms independently per visit, so no
real participant is contained or reconstructable. The main script
points at this file by default:

```r
source("run_synthetic_reproduction.R")
```

Outputs are written to `./output`. Results resemble the published
prevalence-driven patterns but are not the published estimates. For a
quick smoke test, reduce `n_boot`, `n_perm`, `n_null`, and
`n_boot_cluster` in the USER CONFIGURATION section first (the
bootstrap/permutation sections dominate the run time). Authorized
users of the original KDCA data can set `data_path` to the cohort file
to reproduce the published results.

**2. Worked example on public data.** `worked_example_public_data.R`
applies the same framework (co-occurrence adjacency, weighted density,
strength, weighted PageRank, modularity, assortativity, differential
analysis, permutation test) to the PHQ-9 depression items from NHANES
2017–2018, downloaded directly from the CDC (no registration
required):

```r
source("worked_example_public_data.R")
```

Outputs are written to `./output_worked_example`.

## Requirements

R (≥ 4.1) with: `dplyr`, `tidyr`, `purrr`, `stringr`, `readxl`,
`igraph`, `ggplot2`, `ggrepel`, `ineq`, `tableone`, `survey`,
`scales`. Optional: `IsingFit` and `NetworkComparisonTest` (Section
5D; skipped with a message if absent), `writexl` (xlsx export in 5E),
`haven` (worked example).

```r
install.packages(c("dplyr","tidyr","purrr","stringr","readxl","igraph",
                   "ggplot2","ggrepel","ineq","tableone","survey",
                   "scales","IsingFit","NetworkComparisonTest",
                   "writexl","haven"))
```

## Data availability

The cohort data used in the study are not publicly available in
accordance with KDCA policies; they may become available in the
future, subject to the authority of the KDCA.
