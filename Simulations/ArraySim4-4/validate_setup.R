#!/usr/bin/env Rscript
# Preflight checks for ArraySim4-4.

rm(list = ls())
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
  setwd(dirname(script_path))
}
source("Simulation.R")

grid <- make_study_grid()
stopifnot(
  nrow(grid) == 12L,
  identical(sort(unique(grid$tau)), as.integer(STUDY_TAUS)),
  identical(sort(unique(grid$sample_size_scenario)), c("large", "low")),
  identical(sort(unique(grid$initiation_scenario)), c("high", "low")),
  all(table(grid$tau) == 4L),
  all(table(grid$sample_size_scenario) == 6L),
  all(table(grid$initiation_scenario) == 6L),
  STUDY_BETA_TRT == -0.7,
  DEFAULT_CONF_MULT == 1.0
)
message("PASS: 12-cell tau x sample-size x initiation grid is complete.")

make_cell_data <- function(n_per_site, init_intercept, seed) {
  simulate_multisite_tv(
    K = STUDY_K,
    n_per_site = rep(n_per_site, STUDY_K),
    init_intercepts = rep(init_intercept, STUDY_K),
    tau = 5L,
    t_star = 8L,
    beta_event = set_beta_trt(0),
    beta_init = scale_confounding(DEFAULT_CONF_MULT),
    base_seed = seed
  )
}

low_rare <- make_cell_data(SAMPLE_SIZE_LEVELS[["low"]],
                           INITIATION_LEVELS[["low"]], 8675309)
large_high <- make_cell_data(SAMPLE_SIZE_LEVELS[["large"]],
                             INITIATION_LEVELS[["high"]], 8675309)
stopifnot(
  nrow(low_rare) == STUDY_K * SAMPLE_SIZE_LEVELS[["low"]],
  nrow(large_high) == STUDY_K * SAMPLE_SIZE_LEVELS[["large"]],
  mean(low_rare$A_tau) < mean(large_high$A_tau),
  abs(plogis(INITIATION_LEVELS[["low"]]) - 0.01098694) < 1e-7,
  abs(plogis(INITIATION_LEVELS[["high"]]) - 0.1824255) < 1e-7
)
message(sprintf(
  paste0("PASS: boundary cells generate N=%d and N=%d; initiation by tau ",
         "is %.3f versus %.3f."),
  nrow(low_rare), nrow(large_high), mean(low_rare$A_tau), mean(large_high$A_tau)
))

check_event_free_dgp <- function(dat, t_star) {
  stopifnot(!any(dat$S <= t_star & dat$S > dat$T_event))
  L1 <- as.matrix(dat[, paste0("L1_", seq_len(t_star)), drop = FALSE])
  L2 <- as.matrix(dat[, paste0("L2_", seq_len(t_star)), drop = FALSE])
  after_event <- outer(dat$T_event, seq_len(t_star), `<`)
  stopifnot(all(is.na(L1[after_event])), all(is.na(L2[after_event])))

  for (d in split(dat, dat$site)) {
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
}
check_event_free_dgp(low_rare, 8L)
check_event_free_dgp(large_high, 8L)
message("PASS: both boundary cells obey the event-free initiation-risk DGP.")

site_data <- split(low_rare, low_rare$site)
hnum <- central_common_num_hazard(
  lapply(site_data, local_initiation_counts, t_star = 8L)
)
stopifnot(max(abs(hnum - .tv_common_num_hazard(low_rare, 8L))) < 1e-12)
message("PASS: common numerator is reconstructed from local event-free counts.")

truth <- compute_truth_ccw_tv(
  N = 50000, tau = 5L, t_star = 8L,
  beta_event = set_beta_trt(STUDY_BETA_TRT),
  beta_init = scale_confounding(DEFAULT_CONF_MULT),
  site_sizes = rep(SAMPLE_SIZE_LEVELS[["low"]], STUDY_K),
  init_intercepts = rep(INITIATION_LEVELS[["low"]], STUDY_K),
  seed = 24680
)
fit_warnings <- character()
out <- withCallingHandlers(
  run_once_tv(
    K = STUDY_K,
    n_per_site = rep(SAMPLE_SIZE_LEVELS[["low"]], STUDY_K),
    init_intercepts = rep(INITIATION_LEVELS[["low"]], STUDY_K),
    tau = 5L, t_star = 8L,
    beta_trt = STUDY_BETA_TRT,
    beta_init = scale_confounding(DEFAULT_CONF_MULT),
    base_seed = 13579, truth = truth, verbose = FALSE
  ),
  warning = function(w) {
    fit_warnings <<- c(fit_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
unexpected_warnings <- fit_warnings[!grepl("^glm.fit:", fit_warnings)]
if (length(unexpected_warnings))
  stop("Unexpected warning(s): ", paste(unique(unexpected_warnings), collapse = "; "))
if (length(fit_warnings))
  message(sprintf(
    "INFO: rare/small-cell smoke fit produced %d expected glm separation warning(s).",
    length(fit_warnings)
  ))
expected_methods <- c(
  "fed_ccw_tvipcw", "pooled_ccw_tvipcw", "fed_ipw_no_clone",
  "fed_perprotocol_naive", "fed_landmark_ipw", "local_ccw_meta"
)
stopifnot(
  identical(sort(unique(out$results$method)), sort(expected_methods)),
  nrow(out$results) == 6L * 8L,
  all(is.finite(out$results$estimate)),
  !any(grepl("gcomp", out$results$method, ignore.case = TRUE))
)
message("PASS: low-sample/rare-initiation end-to-end run returns six methods.")
message("All ArraySim4-4 preflight checks passed.")
