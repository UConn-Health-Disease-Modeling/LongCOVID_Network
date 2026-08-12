########################################################################
# Reproducibility runner
# ---------------------------------------------------------------------
# Runs the complete analysis pipeline end-to-end on the fully synthetic
# dataset shipped with this repository (data/synthetic_longcovid.xlsx),
# regenerating all figures and tables under ./output.
#
# The synthetic dataset mirrors the variable schema and per-group
# marginal symptom prevalences of the original cohort but contains no
# real participant information; results will therefore resemble the
# published prevalence-driven patterns but are NOT the published
# estimates.
#
# Notes
# - Full run time is dominated by the bootstrap/permutation sections;
#   for a quick smoke test, reduce n_boot, n_perm, n_null, and
#   n_boot_cluster in the USER CONFIGURATION section of the main script
#   (e.g., to 50, 100, 50, 50).
# - Sections 5D (IsingFit / NetworkComparisonTest) and the xlsx export
#   in 5E are skipped automatically if the optional packages are not
#   installed.
########################################################################

required <- c("dplyr", "tidyr", "purrr", "stringr", "readxl", "igraph",
              "ggplot2", "ggrepel", "ineq", "tableone", "survey", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1),
                            quietly = TRUE)]
if (length(missing) > 0) {
  stop("Please install the required packages first:\n  install.packages(c(",
       paste0('"', missing, '"', collapse = ", "), "))")
}

message("Running the full pipeline on the synthetic dataset ...")
source("LongCOVID_Network_Unified_v3_Revision.R")
message("Done. All outputs are under ./output")
