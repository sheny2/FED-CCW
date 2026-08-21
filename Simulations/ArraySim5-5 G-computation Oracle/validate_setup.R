#!/usr/bin/env Rscript
# Lightweight preflight for ArraySim5-5 G-computation Oracle.

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
heterogeneity_test <- "high"
setting_test <- get_heterogeneity_setting(heterogeneity_test)
init_intercepts_test <- setting_test$init_intercepts
mix_test <- setting_test$patient_mix
dat <- simulate_multisite_tv(
  K = length(setting_test$site_sizes), n_per_site = setting_test$site_sizes,
  init_intercepts = init_intercepts_test,
  patient_mix = mix_test,
  tau = tau, t_star = t_star,
  beta_event = set_beta_trt(beta_test),
  beta_init = DEFAULT_BETA_INIT,
  base_seed = 8675309
)

site_n <- table(dat$site)
stopifnot(
  length(site_n) == length(DEFAULT_SITE_SIZES),
  identical(as.integer(site_n), as.integer(setting_test$site_sizes))
)
init_rate <- tapply(dat$A_tau, dat$site, mean)
if (cor(seq_along(init_rate), init_rate, method = "spearman") < 0.8)
  stop("Site initiation prevalence does not follow the intended gradient.")
message(sprintf(
  "PASS: %d sites (%s patients) generated; high-heterogeneity initiation by tau ranges %.3f to %.3f.",
  length(setting_test$site_sizes), paste(setting_test$site_sizes, collapse = "/"),
  min(init_rate), max(init_rate)
))

# Check that the generated high-level site populations match their configured
# X distributions despite the unequal site sizes.
empirical_mix <- do.call(rbind, lapply(split(dat, dat$site), function(d) {
  data.frame(x1_mean = mean(d$x1), x1_sd = sd(d$x1), x2_prob = mean(d$x2))
}))
stopifnot(
  max(abs(empirical_mix$x1_mean - mix_test$x1_mean)) < 0.20,
  max(abs(empirical_mix$x1_sd - mix_test$x1_sd)) < 0.18,
  max(abs(empirical_mix$x2_prob - mix_test$x2_prob)) < 0.07
)
size_ratio <- vapply(HETEROGENEITY_SETTINGS, function(s) {
  max(s$site_sizes) / min(s$site_sizes)
}, numeric(1))
x1_spread <- vapply(HETEROGENEITY_SETTINGS, function(s) {
  diff(range(s$patient_mix$x1_mean))
}, numeric(1))
x2_spread <- vapply(HETEROGENEITY_SETTINGS, function(s) {
  diff(range(s$patient_mix$x2_prob))
}, numeric(1))
practice_spread <- vapply(HETEROGENEITY_SETTINGS, function(s) {
  diff(range(s$init_intercepts))
}, numeric(1))
stopifnot(
  all(diff(size_ratio) > 0),
  all(diff(x1_spread) > 0),
  all(diff(x2_spread) > 0),
  all(diff(practice_spread) > 0),
  all(vapply(HETEROGENEITY_SETTINGS, function(s) sum(s$site_sizes),
             numeric(1)) == 5000),
  all(vapply(HETEROGENEITY_SETTINGS, function(s) mean(s$init_intercepts),
             numeric(1)) == -3)
)
message(paste(
  "PASS: sample-size imbalance, patient-mix separation, and practice",
  "intercept separation all increase from low to moderate to high."
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

fed_ccw <- run_fed_ccw_tvipcw(dat, tau, t_star)
pooled_ccw <- run_pooled_ccw_site_stratified(dat, tau, t_star)
aligned_fields <- c(
  "psi1", "psi0", "RD", "RR", "OR", "RMST1", "RMST0", "RMST_diff",
  "SE_RD", "SE_logRR", "SE_logOR", "SE_RMSTdiff", "SE_RMST1", "SE_RMST0"
)
alignment_error <- max(abs(
  unlist(fed_ccw[aligned_fields]) - unlist(pooled_ccw[aligned_fields])
))
stopifnot(alignment_error < 1e-12)
message(paste(
  "PASS: Fed-CCW and fully site-stratified pooled CCW agree across five",
  "jointly heterogeneous sites (point estimates and model-based SEs)."
))

pooled_fe <- run_pooled_ccw_site_fe(dat, tau, t_star)
stopifnot(all(is.finite(unlist(pooled_fe[c(
  "psi1", "psi0", "RD", "RR", "OR", "RMST1", "RMST0", "RMST_diff",
  "SE_RD", "SE_logRR", "SE_logOR", "SE_RMSTdiff"
)]))))
message(paste(
  "PASS: conventional pooled CCW with site fixed effects and shared slopes",
  "returns finite estimates under high combined heterogeneity."
))

landmark <- run_fed_landmark_ipw(dat, tau, t_star)
stopifnot(
  all(landmark$S1[seq_len(tau + 1L)] == 1),
  all(landmark$S0[seq_len(tau + 1L)] == 1)
)
message("PASS: landmark follow-up begins after tau.")

truth <- compute_truth_gcomp_tv(
  N = 50000, tau = tau, t_star = t_star,
  beta_event = set_beta_trt(beta_test),
  beta_init = DEFAULT_BETA_INIT,
  site_sizes = setting_test$site_sizes,
  init_intercepts = init_intercepts_test,
  patient_mix = mix_test,
  seed = 24680
)
stopifnot(
  identical(truth$oracle_type, "parametric Monte Carlo g-computation"),
  length(truth$S1) == t_star + 1L,
  length(truth$S0) == t_star + 1L,
  all(diff(truth$S1) <= 1e-12),
  all(diff(truth$S0) <= 1e-12),
  truth$RMST1 > truth$RMST0,
  truth$psi1 < truth$psi0,
  abs(sum(truth$initiation_time_distribution_g1) - 1) < 1e-12
)
message(paste(
  "PASS: g-computation oracle generates valid survival curves under forced",
  "initiation by tau versus withholding treatment."
))
out <- run_once_tv(
  K = length(setting_test$site_sizes), n_per_site = setting_test$site_sizes,
  init_intercepts = init_intercepts_test,
  patient_mix = mix_test,
  tau = tau, t_star = t_star, beta_trt = beta_test,
  beta_init = DEFAULT_BETA_INIT,
  base_seed = 13579, truth = truth, verbose = FALSE
)
expected <- c(
  "fed_ccw_tvipcw", "pooled_ccw_site_fe",
  "fed_ipw_no_clone", "fed_perprotocol_naive",
  "fed_landmark_ipw"
)
observed <- sort(unique(out$results$method))
stopifnot(
  identical(observed, sort(expected)),
  nrow(out$results) == 5L * 8L,
  all(is.finite(out$results$estimate)),
  !"pooled_ccw_site_stratified" %in% observed
)
message(paste(
  "PASS: end-to-end run returns five finite methods and omits the redundant",
  "site-stratified pooled CCW comparator."
))
grid <- simulation_grid()
stopifnot(
  nrow(grid) == 27L,
  identical(unique(grid$heterogeneity), names(HETEROGENEITY_SETTINGS))
)
message("PASS: array grid contains 27 cells and all three heterogeneity levels.")
message("All ArraySim5-5 G-computation Oracle preflight checks passed.")
