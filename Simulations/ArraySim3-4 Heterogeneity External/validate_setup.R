#!/usr/bin/env Rscript
# Lightweight preflight for the external-reference one-shot estimator.

rm(list = ls())
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
  setwd(dirname(script_path))
}

source("external_reference.R")
source("Simulation.R")

ref <- load_external_reference(
  EXTERNAL_REFERENCE_FILE,
  t_star = EXTERNAL_REFERENCE_TSTAR
)
hnum <- ref$numerator_hazard
message(sprintf(
  "Reference loaded: N=%d, heterogeneity=%s, seed=%d, t*=%d.",
  ref$N_used, ref$heterogeneity, ref$seed, ref$t_star
))

tau <- 3L
t_star <- 8L
dat <- simulate_multisite_tv(
  K = 3, n_per_site = 400,
  tau = tau, t_star = t_star,
  heterogeneity = "high",
  beta_event = set_beta_trt(-0.7),
  base_seed = 8675309
)

# Every site must use the exact external vector without estimating a numerator.
for (d in split(dat, dat$site)) {
  P <- .tv_prep(d, t_star)
  H <- .tv_init_hazard(P, hnum = hnum[seq_len(t_star)])
  for (m in seq_len(t_star)) {
    active <- !is.na(H$Hden[, m])
    if (any(active))
      stopifnot(max(abs(H$Hnum[active, m] - hnum[m])) < 1e-12)
  }
}
message("PASS: every site uses the identical fixed external numerator.")

# With K=1, federated and pooled implementations must be identical.
one_site <- dat[dat$site == 1L, , drop = FALSE]
f1 <- run_fed_ccw_tvipcw(
  one_site, tau, t_star, hnum = hnum[seq_len(t_star)]
)
p1 <- run_pooled_ccw_tvipcw(
  one_site, tau, t_star, hnum = hnum[seq_len(t_star)]
)
stopifnot(
  abs(f1$RD - p1$RD) < 1e-12,
  abs(f1$RMST_diff - p1$RMST_diff) < 1e-12
)
message("PASS: K=1 federated and pooled estimates are identical.")

fed <- run_fed_ccw_tvipcw(
  dat, tau, t_star, hnum = hnum[seq_len(t_star)]
)
pooled <- run_pooled_ccw_tvipcw(
  dat, tau, t_star, hnum = hnum[seq_len(t_star)]
)
gap <- fed$RMST_diff - pooled$RMST_diff
if (!is.finite(gap)) stop("Non-finite heterogeneous Fed-pooled gap.")
if (abs(gap) > 1.5)
  warning(sprintf("Large smoke-test RMST gap: %.3f. Inspect weights.", gap))
message(sprintf(
  "PASS: heterogeneous smoke test completed (Fed-pooled RMST gap %.3f).",
  gap
))

truth <- compute_truth_multisite_tv(
  N = 30000, K = 3, heterogeneity = "high",
  tau = tau, t_star = t_star,
  beta_event = set_beta_trt(-0.7),
  numerator_hazard = hnum,
  seed = 24680
)
out <- run_once_tv(
  K = 3, n_per_site = 200, heterogeneity = "high",
  tau = tau, t_star = t_star, beta_trt = -0.7,
  numerator_hazard = hnum,
  base_seed = 13579, truth = truth, verbose = FALSE
)
expected <- c(
  "fed_ccw_tvipcw", "pooled_ccw_tvipcw",
  "tvipcw_no_clone", "perprotocol_naive"
)
observed <- sort(unique(out$results$method))
stopifnot(
  identical(observed, sort(expected)),
  nrow(out$results) == 4L * 8L,
  all(is.finite(out$results$estimate)),
  !any(grepl("gcomp", observed, ignore.case = TRUE))
)
message("PASS: end-to-end run returns exactly four methods and no G-comp.")
message(sprintf(
  "Oracle smoke test: RD=%.4f, RMST difference=%.4f, N=%d.",
  truth$RD, truth$RMST_diff, truth$N_used
))
message("All external-reference preflight checks passed.")
