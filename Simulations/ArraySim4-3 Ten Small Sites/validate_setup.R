#!/usr/bin/env Rscript
# Lightweight preflight for ArraySim4-3 Ten Small Sites.

rm(list = ls())
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
  setwd(dirname(script_path))
}

source("Simulation.R")

tau <- 5L
t_star <- 8L
dat <- simulate_multisite_tv(
  K = 10, n_per_site = TEN_SITE_SIZES,
  init_intercepts = TEN_SITE_INIT_INTERCEPTS,
  tau = tau, t_star = t_star,
  beta_event = set_beta_trt(-0.7),
  beta_init = scale_confounding(DEFAULT_CONF_MULTS[["medium"]]),
  base_seed = 8675309
)

site_n <- table(dat$site)
stopifnot(
  length(site_n) == 10L,
  identical(as.integer(site_n), as.integer(TEN_SITE_SIZES))
)
init_rate <- tapply(dat$A_tau, dat$site, mean)
if (cor(seq_along(init_rate), init_rate, method = "spearman") < 0.8)
  stop("Site initiation prevalence does not follow the intended gradient.")
message(sprintf(
  "PASS: ten 300-patient sites generated; initiation by tau ranges %.3f to %.3f.",
  min(init_rate), max(init_rate)
))

site_data <- split(dat, dat$site)
hnum <- central_common_num_hazard(
  lapply(site_data, local_initiation_counts, t_star = t_star)
)
stopifnot(max(abs(hnum - .tv_common_num_hazard(dat, t_star))) < 1e-12)
message("PASS: common numerator uses only aggregated counts from ten sites.")

fed_pp <- run_fed_perprotocol(dat, tau, t_star)
pooled_dat <- dat
pooled_dat$site <- 1L
pooled_pp <- run_fed_perprotocol(pooled_dat, tau, t_star)
stopifnot(
  abs(fed_pp$RD - pooled_pp$RD) < 1e-12,
  abs(fed_pp$RMST_diff - pooled_pp$RMST_diff) < 1e-12
)
message("PASS: federated per-protocol matches pooled crude analysis.")

one_site <- dat[dat$site == 5L, , drop = FALSE]
fed1 <- run_fed_ccw_tvipcw(one_site, tau, t_star)
pool1 <- run_pooled_ccw_tvipcw(one_site, tau, t_star)
meta1 <- run_local_ccw_meta(one_site, tau, t_star)
stopifnot(
  abs(fed1$RD - pool1$RD) < 1e-12,
  abs(fed1$RMST_diff - pool1$RMST_diff) < 1e-12,
  abs(fed1$RD - meta1$RD) < 1e-12,
  abs(fed1$RMST_diff - meta1$RMST_diff) < 1e-12
)
message("PASS: all CCW aggregations agree when K=1.")

landmark <- run_fed_landmark_ipw(dat, tau, t_star)
stopifnot(
  all(landmark$S1[seq_len(tau + 1L)] == 1),
  all(landmark$S0[seq_len(tau + 1L)] == 1)
)
message("PASS: landmark follow-up begins after tau.")

truth <- compute_truth_ccw_tv(
  N = 50000, tau = tau, t_star = t_star,
  beta_event = set_beta_trt(-0.7),
  beta_init = scale_confounding(DEFAULT_CONF_MULTS[["medium"]]),
  site_sizes = TEN_SITE_SIZES,
  init_intercepts = TEN_SITE_INIT_INTERCEPTS,
  seed = 24680
)
out <- run_once_tv(
  K = 10, n_per_site = TEN_SITE_SIZES,
  init_intercepts = TEN_SITE_INIT_INTERCEPTS,
  tau = tau, t_star = t_star, beta_trt = -0.7,
  beta_init = scale_confounding(DEFAULT_CONF_MULTS[["medium"]]),
  base_seed = 13579, truth = truth, verbose = FALSE
)
expected <- c(
  "fed_ccw_tvipcw", "pooled_ccw_tvipcw",
  "fed_ipw_no_clone", "fed_perprotocol_naive",
  "fed_landmark_ipw", "local_ccw_meta"
)
observed <- sort(unique(out$results$method))
stopifnot(
  identical(observed, sort(expected)),
  nrow(out$results) == 6L * 8L,
  all(is.finite(out$results$estimate)),
  !any(grepl("gcomp", observed, ignore.case = TRUE))
)
message("PASS: end-to-end run returns six finite methods and no G-comp.")
message("All ArraySim4-3 Ten Small Sites preflight checks passed.")
