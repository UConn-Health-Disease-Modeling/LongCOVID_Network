# LongCOVID_Network

R code accompanying:

> **Jo Y., Jung J., Lee J.**
> *Network Drivers of Long COVID: SARS-CoV-2 Infection Reshapes Symptom Interdependencies.*

This repository contains the full analysis pipeline used to construct and
compare weighted symptom co-occurrence networks across SARS-CoV-2-infected
and uninfected individuals, with longitudinal time-bin and bootstrap-based
inferential extensions.

---

## Contents

| File | Description |
|------|-------------|
| `LongCOVID_Network_Unified.R` | End-to-end analysis pipeline. One self-contained R script that reproduces every figure and table reported in the manuscript and supplement. |
| `README.md` | This file. |
| `.gitignore` | Excludes R session state and auto-generated outputs from version control. |

---

## What the script produces

All outputs land under the user-defined `out_dir` (default
`Output/figures/`, `Output/tables/`, `Output/data/`). Mapping to the paper:

### Main manuscript

| Object | Output filename |
|---|---|
| **Table 1** — participant characteristics (unweighted + IPW-weighted) | `Table1_unweighted.csv`, `Table1_weighted.csv` |
| **Figure 1** — A. Uninfected / B. Infected symptom-influence lollipop | `overall_group_lollipop.png` |
| **Figure 2** — Δ Strength vs Δ PageRank scatter, coloured by 6 clinical domains | `baseline_delta_scatter.png` |
| **Figure 3** — Bootstrap forest plot of Δ PageRank with 95 % CI and FDR-adjusted significance | `bootstrap_forest_pagerank.png` |

### Supplementary

| Object | Output filename |
|---|---|
| Figure S1 — participant flow diagram | *not produced by R* (external diagram) |
| Table S1 — 6-domain symptom grouping | *manually constructed* (in the manuscript) |
| **Figure S2** — IPW-adjusted two-panel network plot | `ipw_network_two_panel.png` |
| **Table S2** — Differential network analysis (Δ Strength, Δ PageRank, FDR) | `bootstrap_node_results.csv` |
| **Figure S3** — Time-bin-specific symptom influence (8 panels) | `timebin_symptom_influence_forest_full.png` |
| **Figure S4** — Symptom × time-bin Δ PageRank heatmap (B = 300) | `delta_pagerank_heatmap.png` |

The full bootstrap and permutation distributions are also serialised to
`Output/data/*.rds` for downstream re-analysis without re-running the
inference loops.

---

## Pipeline structure

The script is divided into clearly labelled sections:

| Section | Purpose | Manuscript outputs |
|---|---|---|
| **0** | Setup, library imports, output directories | — |
| **1** | User configuration (paths, column names, PS covariates, bootstrap iterations) | — |
| **2** | Data loading and preprocessing (BMI / age / CCI grouping, factor levels) | — |
| **3** | Helper functions (`to01`, `clean_symptom_names`, `symptom_category`, `weighted_density`, `weighted_assortativity_strength`, `global_network_metrics`, `calc_node_metrics`) | — |
| **4A** | Baseline (unweighted) network construction and per-symptom centrality | Fig 1, Fig 2 |
| **4B** | Propensity score model, stabilised IPW (truncated 99 %), weighted Table 1 | Table 1 |
| **4C** | IPW-weighted symptom network for each group | Fig S2 |
| **4D** | Time-bin (0–3 / 3–6 / 6–12 / 12–24 mo) network analysis with bootstrap CIs for global metrics | Fig S3 inputs |
| **4E** | Bootstrap and permutation inference for Δ Strength and Δ PageRank, BH-FDR correction | Fig 3, Table S2 |
| **4F** | Symptom × time-bin Δ PageRank quantification (B = 300 bootstrap) | Fig S4 |

---

## Software requirements

- **R ≥ 4.3**
- CRAN packages used (loaded at the top of the script):

  ```
  dplyr, tidyr, purrr, stringr, readxl,
  igraph, ggplot2, ggrepel,
  ineq, tableone, survey, scales
  ```

  Install missing packages with:

  ```r
  install.packages(c("dplyr","tidyr","purrr","stringr","readxl",
                     "igraph","ggplot2","ggrepel",
                     "ineq","tableone","survey","scales"))
  ```

---

## Data

The Korean cohort dataset analysed in the manuscript is governed by KDCA
data-sharing policy and is not publicly redistributable. The script can
be adapted to other cohorts by editing the user-configuration block at
the top of the file (Section 1):

```r
data_path   <- "path/to/your/dataset.xlsx"
id_col      <- "SUBJ_ID"          # subject identifier column
infect_col  <- "infected"         # 0 = uninfected, 1 = infected
time_col    <- "months_since_latest"  # or set NULL to skip time-bin analysis

# Symptom columns: hard-code or auto-detect by position
symptom_cols_hardcoded <- c("symptom1", "symptom2", ...)

# Propensity-score covariates
ps_vars <- c("Age_group", "SEX", "CCI_group", "BMI_group",
             "smoking", "alcohol",
             "Accumulate_vac_group", "marital_status",
             "househole_member_group", "education_status", "JOB_PREINF")
```

Helper functions for category creation (BMI → group, CCI → group, age →
band, etc.) are in Section 2 and can be modified to match the variable
names in your dataset.

---

## Running

After updating the configuration block:

```r
source("LongCOVID_Network_Unified.R")
```

Approximate runtime on a standard laptop:

| Section | Iterations | Wall-clock |
|---|---|---|
| Baseline + IPW + Time-bin | — | < 1 min |
| Bootstrap inference | B = 500 | ~ 2–3 min |
| Permutation inference | B = 1000 | ~ 4–5 min |
| Time-bin Δ PageRank bootstrap | B = 300 | ~ 2–3 min |
| **Total** | | **~ 10–12 min** |

For faster iteration during development, lower `n_boot`, `n_perm`, and
`n_boot_pr` near the top of the script (set to 100 for testing) and
restore them before the final paper figures.

---

## Methodological notes

- **Symptom encoding** is binary (present/absent) per visit; severity
  and duration are not modelled.
- **Co-occurrence edges** reflect within-individual cross-sectional
  association at the visit level; they do not encode temporal
  precedence or causation between symptoms.
- **Inverse probability weighting (IPW)** addresses observed
  confounders (age, sex, BMI, CCI, smoking, alcohol, vaccination dose,
  marital status, household members, education, pre-infection job
  status). Weights are stabilised and truncated at the 99 th
  percentile.
- **Community labels (1, 2, 3 …)** are computed independently within
  each network and are not directly comparable across panels.
- **Bootstrap and permutation** confidence intervals are computed at
  the participant level within group × time-bin strata.
- All randomness is controlled by user-set seeds (`set_seed_boot`,
  `set_seed_perm`).

---

## Citation

If you use this code, please cite the accompanying manuscript:

```
Jo Y., Jung J., Lee J. Network Drivers of Long COVID:
SARS-CoV-2 Infection Reshapes Symptom Interdependencies.
[journal] [year];[vol]:[pages]. doi:[…]
```

---

## License

Released under the MIT License unless otherwise noted. See `LICENSE`
file (to be added) for terms.

---

## Contact

For questions about the analysis or implementation, please open a
GitHub issue or contact the corresponding authors:

- Youngji Jo — `jo@uchc.edu`
- Jaehun Jung — `eastside1st@gmail.com`
