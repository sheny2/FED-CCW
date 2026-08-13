#!/usr/bin/env Rscript
# Preflight checks for ArraySim4-5.

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
  identical(sort(unique(grid$outcome_scenario)), c("common", "rare")),
  all(table(grid$tau) == 4L),
  all(table(grid$sample_size_scenario) == 6L),
  all(table(grid$outcome_scenario) == 6L),
  STUDY_BETA_TRT == -0.7,
  STUDY_CONF_MULT == 1.0,
  STUDY_INIT_INTERCEPT == -3.0
)
message("PASS: 12-cell tau x sample-size x outcome-frequency grid is complete.")

make_cell_data <- function(n_per_site, event_intercept, seed) {
  beta_event <- set_beta_trt(
    STUDY_BETA_TRT,
    set_event_intercept(event_intercept)
  )
  simulate_multisite_tv(
    K = STUDY_K,
    n_per_site = rep(n_per_site, STUDY_K),
    init_intercepts = rep(STUDY_INIT_INTERCEPT, STUDY_K),
    tau = 4L,
    t_star = STUDY_TSTAR,
    beta_event = beta_event,
    beta_init = scale_confounding(STUDY_CONF_MULT),
    base_seed = seed
  )
}

low_rare <- make_cell_data(SAMPLE_SIZE_LEVELS[["low"]],
                           OUTCOME_LEVELS[["rare"]], 8675309)
large_common <- make_cell_data(SAMPLE_SIZE_LEVELS[["large"]],
                               OUTCOME_LEVELS[["common"]], 8675309)
rare_death <- mean(low_rare$delta)
common_death <- mean(large_common$delta)
stopifnot(
  nrow(low_rare) == STUDY_K * SAMPLE_SIZE_LEVELS[["low"]],
  nrow(large_common) == STUDY_K * SAMPLE_SIZE_LEVELS[["large"]],
  rare_death < 0.15,
  common_death > 0.60,
  rare_death < common_death
)
message(sprintf(
  paste0("PASS: boundary cells generate N=%d and N=%d; cumulative mortality ",
         "is %.3f versus %.3f."),
  nrow(low_rare), nrow(large_common), rare_death, common_death
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
check_event_free_dgp(low_rare, STUDY_TSTAR)
check_event_free_dgp(large_common, STUDY_TSTAR)
message("PASS: both boundary cells obey the event-free initiation-risk DGP.")

site_data <- split(low_rare, low_rare$site)
hnum <- central_common_num_hazard(
  lapply(site_data, local_initiation_counts, t_star = STUDY_TSTAR)
)
stopifnot(max(abs(hnum - .tv_common_num_hazard(low_rare, STUDY_TSTAR))) < 1e-12)
message("PASS: common numerator is reconstructed from local event-free counts.")

rare_beta_event <- set_beta_trt(
  STUDY_BETA_TRT,
  set_event_intercept(OUTCOME_LEVELS[["rare"]])
)
truth <- compute_truth_ccw_tv(
  N = 50000, tau = 4L, t_star = STUDY_TSTAR,
  beta_event = rare_beta_event,
  beta_init = scale_confounding(STUDY_CONF_MULT),
  site_sizes = rep(SAMPLE_SIZE_LEVELS[["low"]], STUDY_K),
  init_intercepts = rep(STUDY_INIT_INTERCEPT, STUDY_K),
  seed = 24680
)

fit_warnings <- character()
out <- withCallingHandlers(
  run_once_tv(
    K = STUDY_K,
    n_per_site = rep(SAMPLE_SIZE_LEVELS[["low"]], STUDY_K),
    init_intercepts = rep(STUDY_INIT_INTERCEPT, STUDY_K),
    tau = 4L, t_star = STUDY_TSTAR,
    beta_trt = STUDY_BETA_TRT,
    beta_event = rare_beta_event,
    beta_init = scale_confounding(STUDY_CONF_MULT),
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
    "INFO: rare-outcome smoke fit produced %d expected glm warning(s).",
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
message("PASS: low-sample/rare-outcome end-to-end run returns six finite methods.")
message("All ArraySim4-5 preflight checks passed.")
