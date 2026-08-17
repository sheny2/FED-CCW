#!/usr/bin/env Rscript
# Lightweight preflight for ArraySim5-1 Five Methods.

rm(list = ls())
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
  setwd(dirname(script_path))
}

source("Simulation.R")
validate_params()

tau <- SIM_TAU_VALUES[[1L]]
t_star <- max(tau + 3L, 8L)
beta_test <- SIM_BETA_TRT_VALUES[[min(2L, length(SIM_BETA_TRT_VALUES))]]
conf_test <- if ("medium" %in% names(DEFAULT_CONF_MULTS)) "medium" else
  names(DEFAULT_CONF_MULTS)[[1L]]
dat <- simulate_multisite_tv(
  K = length(DEFAULT_SITE_SIZES), n_per_site = DEFAULT_SITE_SIZES,
  init_intercepts = DEFAULT_SITE_INIT_INTERCEPTS,
  tau = tau, t_star = t_star,
  beta_event = set_beta_trt(beta_test),
  beta_init = scale_confounding(DEFAULT_CONF_MULTS[[conf_test]]),
  base_seed = 8675309
)

site_n <- table(dat$site)
stopifnot(
  length(site_n) == length(DEFAULT_SITE_SIZES),
  identical(as.integer(site_n), as.integer(DEFAULT_SITE_SIZES))
)
init_rate <- tapply(dat$A_tau, dat$site, mean)
if (cor(seq_along(init_rate), init_rate, method = "spearman") < 0.8)
  stop("Site initiation prevalence does not follow the intended gradient.")
message(sprintf(
  "PASS: %d sites (%s patients) generated; initiation by tau ranges %.3f to %.3f.",
  length(DEFAULT_SITE_SIZES), paste(DEFAULT_SITE_SIZES, collapse = "/"),
  min(init_rate), max(init_rate)
))

post_event_init <- dat$S <= t_star & dat$S > dat$T_event
stopifnot(!any(post_event_init))
L1 <- as.matrix(dat[, paste0("L1_", seq_len(t_star)), drop = FALSE])
L2 <- as.matrix(dat[, paste0("L2_", seq_len(t_star)), drop = FALSE])
after_event <- outer(dat$T_event, seq_len(t_star), `<`)
stopifnot(
  all(is.na(L1[after_event])),
  all(is.na(L2[after_event]))
)
message("PASS: no post-event initiation or covariate evolution occurs.")

site_data <- split(dat, dat$site)
for (d in site_data) {
  counts <- local_initiation_counts(d, t_star)
  expected_initiated <- vapply(seq_len(t_star), function(m) {
    sum(d$S == m & d$T_event >= m)
  }, numeric(1))
  expected_eligible <- vapply(seq_len(t_star), function(m) {
    sum(d$S >= m & d$T_event >= m)
  }, numeric(1))
  stopifnot(
    identical(counts$initiated, expected_initiated),
    identical(counts$eligible, expected_eligible)
  )
}
hnum <- central_common_num_hazard(
  lapply(site_data, local_initiation_counts, t_star = t_star)
)
stopifnot(max(abs(hnum - .tv_common_num_hazard(dat, t_star))) < 1e-12)
message(paste(
  "PASS: common numerator uses only aggregated event-free initiation-risk",
  sprintf("counts from %d sites.", length(DEFAULT_SITE_SIZES))
))

fed_pp <- run_fed_perprotocol(dat, tau, t_star)
pooled_dat <- dat
pooled_dat$site <- 1L
pooled_pp <- run_fed_perprotocol(pooled_dat, tau, t_star)
stopifnot(
  abs(fed_pp$RD - pooled_pp$RD) < 1e-12,
  abs(fed_pp$RMST_diff - pooled_pp$RMST_diff) < 1e-12
)
message("PASS: federated per-protocol matches pooled crude analysis.")

one_site_id <- ceiling(length(DEFAULT_SITE_SIZES) / 2)
one_site <- dat[dat$site == one_site_id, , drop = FALSE]
fed1 <- run_fed_ccw_tvipcw(one_site, tau, t_star)
pool1 <- run_pooled_ccw_tvipcw(one_site, tau, t_star)
stopifnot(
  abs(fed1$RD - pool1$RD) < 1e-12,
  abs(fed1$RMST_diff - pool1$RMST_diff) < 1e-12
)
message("PASS: federated and pooled CCW agree when K=1.")

landmark <- run_fed_landmark_ipw(dat, tau, t_star)
stopifnot(
  all(landmark$S1[seq_len(tau + 1L)] == 1),
  all(landmark$S0[seq_len(tau + 1L)] == 1)
)
message("PASS: landmark follow-up begins after tau.")

truth <- compute_truth_ccw_tv(
  N = 50000, tau = tau, t_star = t_star,
  beta_event = set_beta_trt(beta_test),
  beta_init = scale_confounding(DEFAULT_CONF_MULTS[[conf_test]]),
  site_sizes = DEFAULT_SITE_SIZES,
  init_intercepts = DEFAULT_SITE_INIT_INTERCEPTS,
  seed = 24680
)
out <- run_once_tv(
  K = length(DEFAULT_SITE_SIZES), n_per_site = DEFAULT_SITE_SIZES,
  init_intercepts = DEFAULT_SITE_INIT_INTERCEPTS,
  tau = tau, t_star = t_star, beta_trt = beta_test,
  beta_init = scale_confounding(DEFAULT_CONF_MULTS[[conf_test]]),
  base_seed = 13579, truth = truth, verbose = FALSE
)
expected <- c(
  "fed_ccw_tvipcw", "pooled_ccw_tvipcw",
  "fed_ipw_no_clone", "fed_perprotocol_naive",
  "fed_landmark_ipw"
)
observed <- sort(unique(out$results$method))
stopifnot(
  identical(observed, sort(expected)),
  nrow(out$results) == 5L * 8L,
  all(is.finite(out$results$estimate)),
  !any(grepl("gcomp", observed, ignore.case = TRUE))
)
message("PASS: end-to-end run returns five finite methods and no G-comp.")
message("All ArraySim5-1 Five Methods preflight checks passed.")
