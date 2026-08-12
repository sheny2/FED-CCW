#!/usr/bin/env Rscript
# Lightweight structural and end-to-end preflight for ArraySim4.

rm(list = ls())
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
  setwd(dirname(script_path))
}

source("Simulation.R")

tau <- 3L
t_star <- 8L
dat <- simulate_multisite_tv(
  K = 3, n_per_site = 400,
  tau = tau, t_star = t_star,
  beta_event = set_beta_trt(-0.7),
  beta_init = scale_confounding(DEFAULT_CONF_MULTS[["strong"]]),
  base_seed = 8675309
)

# The numerator can be reconstructed using only site-level interval counts.
site_data <- split(dat, dat$site)
hnum <- central_common_num_hazard(
  lapply(site_data, local_initiation_counts, t_star = t_star)
)
stopifnot(max(abs(hnum - .tv_common_num_hazard(dat, t_star))) < 1e-12)
message("PASS: common numerator uses only aggregated site counts.")

# Federated per-protocol counts must reproduce the pooled crude calculation.
fed_pp <- run_fed_perprotocol(dat, tau, t_star)
pooled_dat <- dat
pooled_dat$site <- 1L
pooled_pp <- run_fed_perprotocol(pooled_dat, tau, t_star)
stopifnot(
  abs(fed_pp$RD - pooled_pp$RD) < 1e-12,
  abs(fed_pp$RMST_diff - pooled_pp$RMST_diff) < 1e-12
)
message("PASS: federated per-protocol exactly matches pooled crude analysis.")

# K=1 provides an identity check for the two CCW implementations.
one_site <- dat[dat$site == 1L, , drop = FALSE]
fed1 <- run_fed_ccw_tvipcw(one_site, tau, t_star)
pool1 <- run_pooled_ccw_tvipcw(one_site, tau, t_star)
stopifnot(
  abs(fed1$RD - pool1$RD) < 1e-12,
  abs(fed1$RMST_diff - pool1$RMST_diff) < 1e-12
)
message("PASS: K=1 federated and pooled CCW are identical.")

# End-to-end smoke test checks the requested method contract.
truth <- compute_truth_ccw_tv(
  N = 30000, tau = tau, t_star = t_star,
  beta_event = set_beta_trt(-0.7),
  beta_init = scale_confounding(DEFAULT_CONF_MULTS[["strong"]]),
  seed = 24680
)
out <- run_once_tv(
  K = 3, n_per_site = 200,
  tau = tau, t_star = t_star, beta_trt = -0.7,
  beta_init = scale_confounding(DEFAULT_CONF_MULTS[["strong"]]),
  base_seed = 13579, truth = truth, verbose = FALSE
)
expected <- c(
  "fed_ccw_tvipcw", "pooled_ccw_tvipcw",
  "fed_ipw_no_clone", "fed_perprotocol_naive"
)
observed <- sort(unique(out$results$method))
stopifnot(
  identical(observed, sort(expected)),
  nrow(out$results) == 4L * 8L,
  all(is.finite(out$results$estimate)),
  !any(grepl("gcomp", observed, ignore.case = TRUE))
)
message("PASS: end-to-end run returns exactly four methods and no G-comp.")
message("All ArraySim4 preflight checks passed.")
