########################################################################
# Worked example: symptom co-occurrence network analysis on PUBLIC data
# ---------------------------------------------------------------------
# Applies the analytical framework of
#   Jo Y., Jung J., Lee J. "Network Drivers of Long COVID:
#   SARS-CoV-2 Infection Reshapes Symptom Interdependencies."
# to a real, openly accessible symptom dataset: the PHQ-9 depression
# screener (DPQ_J) from NHANES 2017-2018 (US CDC), merged with
# demographics (DEMO_J).
#
# Nine depression symptoms are binarized (present = experienced on
# "several days" or more over the past two weeks) and weighted symptom
# co-occurrence networks are compared between two groups (here: male vs
# female, purely as a demonstration of the two-group differential
# framework). The example exercises the same machinery as the main
# pipeline: weighted co-occurrence adjacency, weighted density, node
# strength, weighted PageRank, Louvain modularity, strength
# assortativity, differential (delta) centrality, and a label-
# permutation test.
#
# Requirements: haven, dplyr, tidyr, igraph, ggplot2, ggrepel
# The two NHANES files (~3 MB total) are downloaded directly from the
# CDC; no registration is required.
########################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(igraph); library(ggplot2)
  library(ggrepel)
})
if (!requireNamespace("haven", quietly = TRUE)) {
  stop("Please install the 'haven' package: install.packages('haven')")
}

out_dir <- "output_worked_example"
dir.create(out_dir, showWarnings = FALSE)
set.seed(123)

########################################################################
# 1. Download and prepare the public data
########################################################################

base_url <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/"
tmp_dpq  <- file.path(tempdir(), "DPQ_J.xpt")
tmp_demo <- file.path(tempdir(), "DEMO_J.xpt")

download_robust <- function(url, dest) {
  ok <- tryCatch({
    download.file(url, dest, mode = "wb", quiet = TRUE); TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) {
    # fallback for networks with TLS revocation-check issues (Windows)
    ok <- tryCatch({
      download.file(url, dest, mode = "wb", quiet = TRUE,
                    method = "curl", extra = "--ssl-no-revoke -L"); TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
  }
  if (!ok) stop("Could not download ", url,
                "\nPlease download the file manually and place it at: ", dest)
  invisible(TRUE)
}

download_robust(paste0(base_url, "DPQ_J.xpt"),  tmp_dpq)
download_robust(paste0(base_url, "DEMO_J.xpt"), tmp_demo)
dpq  <- haven::read_xpt(tmp_dpq)
demo <- haven::read_xpt(tmp_demo)

symptom_labels <- c(
  DPQ010 = "Anhedonia",
  DPQ020 = "Depressed mood",
  DPQ030 = "Sleep problems",
  DPQ040 = "Fatigue",
  DPQ050 = "Appetite change",
  DPQ060 = "Feeling bad about self",
  DPQ070 = "Concentration problems",
  DPQ080 = "Psychomotor change",
  DPQ090 = "Self-harm thoughts")
symptom_cols <- names(symptom_labels)

df <- dpq %>%
  inner_join(demo %>% select(SEQN, RIAGENDR), by = "SEQN") %>%
  select(SEQN, RIAGENDR, all_of(symptom_cols)) %>%
  # item scores: 0 = not at all ... 3 = nearly every day; 7/9 = missing
  mutate(across(all_of(symptom_cols),
                ~ ifelse(.x %in% 0:3, as.integer(.x >= 1), NA_integer_))) %>%
  filter(if_all(all_of(symptom_cols), ~ !is.na(.x))) %>%
  mutate(group = as.integer(RIAGENDR == 2))   # 0 = male, 1 = female

message("Analytic sample: ", nrow(df), " participants (",
        sum(df$group == 0), " male / ", sum(df$group == 1), " female)")

########################################################################
# 2. Framework functions (compact versions of the main pipeline)
########################################################################

build_adj <- function(X) {
  O <- crossprod(as.matrix(X))
  diag(O) <- 0
  O
}

weighted_density <- function(A) {
  n <- nrow(A)
  sum(A[upper.tri(A)]) / choose(n, 2)
}

weighted_assortativity_strength <- function(g) {
  edf <- as_data_frame(g, what = "edges")
  if (nrow(edf) == 0) return(NA_real_)
  s  <- strength(g, weights = E(g)$weight)
  su <- s[edf$from]; sv <- s[edf$to]; w <- edf$weight
  mu_u <- weighted.mean(su, w); mu_v <- weighted.mean(sv, w)
  cov_w <- sum(w * (su - mu_u) * (sv - mu_v)) / sum(w)
  sd_u  <- sqrt(sum(w * (su - mu_u)^2) / sum(w))
  sd_v  <- sqrt(sum(w * (sv - mu_v)^2) / sum(w))
  if (sd_u == 0 || sd_v == 0) return(NA_real_)
  cov_w / (sd_u * sd_v)
}

network_metrics <- function(X) {
  A <- build_adj(X)
  g <- graph_from_adjacency_matrix(A, mode = "undirected",
                                   weighted = TRUE, diag = FALSE)
  comm <- cluster_louvain(g, weights = E(g)$weight)
  list(
    graph        = g,
    density      = weighted_density(A),
    modularity   = modularity(comm),
    assortativity = weighted_assortativity_strength(g),
    strength     = strength(g, weights = E(g)$weight)[colnames(X)],
    pagerank     = page_rank(g, weights = E(g)$weight)$vector[colnames(X)],
    community    = membership(comm)[colnames(X)]
  )
}

########################################################################
# 3. Two-group networks and differential analysis
########################################################################

X0 <- df %>% filter(group == 0) %>% select(all_of(symptom_cols))
X1 <- df %>% filter(group == 1) %>% select(all_of(symptom_cols))
m0 <- network_metrics(X0)
m1 <- network_metrics(X1)

global_tab <- tibble(
  group         = c("Male", "Female"),
  n             = c(nrow(X0), nrow(X1)),
  density       = c(m0$density, m1$density),
  modularity    = c(m0$modularity, m1$modularity),
  assortativity = c(m0$assortativity, m1$assortativity)
)
print(global_tab)
write.csv(global_tab, file.path(out_dir, "global_metrics.csv"),
          row.names = FALSE)

node_tab <- tibble(
  symptom        = unname(symptom_labels[symptom_cols]),
  strength_g0    = as.numeric(m0$strength),
  strength_g1    = as.numeric(m1$strength),
  pagerank_g0    = as.numeric(m0$pagerank),
  pagerank_g1    = as.numeric(m1$pagerank),
  delta_strength = strength_g1 - strength_g0,
  delta_pagerank = pagerank_g1 - pagerank_g0
) %>% arrange(desc(delta_pagerank))
print(node_tab, n = Inf)
write.csv(node_tab, file.path(out_dir, "node_metrics.csv"),
          row.names = FALSE)

########################################################################
# 4. Permutation test on the differential metrics
########################################################################

B <- 1000
perm <- replicate(B, {
  g_perm <- sample(df$group)
  Xp0 <- df[g_perm == 0, symptom_cols]
  Xp1 <- df[g_perm == 1, symptom_cols]
  mp0 <- network_metrics(Xp0)
  mp1 <- network_metrics(Xp1)
  c(delta_density    = mp1$density - mp0$density,
    delta_modularity = mp1$modularity - mp0$modularity)
})
obs <- c(delta_density    = m1$density - m0$density,
         delta_modularity = m1$modularity - m0$modularity)
pvals <- vapply(names(obs), function(k)
  mean(abs(perm[k, ]) >= abs(obs[[k]])), numeric(1))
perm_tab <- tibble(metric = names(obs), observed = unname(obs),
                   p_permutation = unname(pvals))
print(perm_tab)
write.csv(perm_tab, file.path(out_dir, "permutation_results.csv"),
          row.names = FALSE)

########################################################################
# 5. Figures
########################################################################

plot_df <- bind_rows(
  tibble(group = "A. Male",   symptom = unname(symptom_labels[symptom_cols]),
         strength = as.numeric(m0$strength),
         pagerank = as.numeric(m0$pagerank)),
  tibble(group = "B. Female", symptom = unname(symptom_labels[symptom_cols]),
         strength = as.numeric(m1$strength),
         pagerank = as.numeric(m1$pagerank))
)

p1 <- plot_df %>%
  group_by(group) %>%
  arrange(pagerank, .by_group = TRUE) %>%
  mutate(sym_f = factor(paste(group, symptom, sep = "___"),
                        levels = paste(group, symptom, sep = "___"))) %>%
  ggplot(aes(x = pagerank, y = sym_f)) +
  geom_segment(aes(x = 0, xend = pagerank, yend = sym_f),
               linewidth = 0.8, color = "#377EB8") +
  geom_point(aes(size = strength), color = "#377EB8") +
  facet_wrap(~ group, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  labs(x = "Weighted PageRank (global influence)", y = NULL,
       size = "Strength",
       title = "Worked example: PHQ-9 symptom influence by group (NHANES 2017-2018)") +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "pagerank_lollipop.png"), p1,
       width = 11, height = 5, dpi = 300)

p2 <- node_tab %>%
  mutate(z_s = as.numeric(scale(delta_strength)),
         z_p = as.numeric(scale(delta_pagerank))) %>%
  ggplot(aes(x = z_s, y = z_p, label = symptom)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(size = 2.5, color = "#E41A1C") +
  geom_text_repel(size = 3.5) +
  labs(x = "Delta strength (z)", y = "Delta PageRank (z)",
       title = "Worked example: differential symptom network centrality (Female - Male)") +
  theme_classic(base_size = 12)
ggsave(file.path(out_dir, "delta_scatter.png"), p2,
       width = 8, height = 6, dpi = 300)

message("Worked example complete. Outputs in ", normalizePath(out_dir))
