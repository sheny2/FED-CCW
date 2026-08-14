#!/usr/bin/env Rscript
# Preflight checks for ArraySim4-6.

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
  identical(sort(unique(grid$size_scenario)), c("balanced", "unbalanced")),
  identical(sort(unique(grid$mix_scenario)), c("heterogeneous", "homogeneous")),
  all(table(grid$tau) == 4L),
  all(table(grid$size_scenario) == 6L),
  all(table(grid$mix_scenario) == 6L),
  all(vapply(SITE_SIZE_LEVELS, sum, numeric(1)) == 3000),
  STUDY_BETA_TRT == -0.7,
  STUDY_CONF_MULT == 1,
  STUDY_INIT_INTERCEPT == -3
)
message("PASS: 12-cell tau x site-size x patient-mix grid is complete.")

hom <- get_site_mix("homogeneous")
het <- get_site_mix("heterogeneous")
stopifnot(
  identical(unname(get_site_sizes("balanced")), c(1000L, 1000L, 1000L)),
  identical(unname(get_site_sizes("unbalanced")), c(200L, 800L, 2000L)),
  all(hom$x1_mean == 0), all(hom$x1_sd == 1), all(hom$x2_prob == 0.4),
  identical(het$x1_mean, c(-1.5, 0, 1.5)),
  identical(het$x1_sd, c(0.7, 1, 1.3)),
  identical(het$x2_prob, c(0.1, 0.4, 0.7))
)
weighted_target <- function(size_level, mix_level) {
  n <- get_site_sizes(size_level)
  m <- get_site_mix(mix_level)
  w <- n / sum(n)
  c(x1 = sum(w * m$x1_mean), x2 = sum(w * m$x2_prob))
}
stopifnot(
  max(abs(weighted_target("balanced", "heterogeneous") - c(x1 = 0, x2 = 0.4))) < 1e-12,
  max(abs(weighted_target("unbalanced", "heterogeneous") - c(x1 = 0.9, x2 = 0.58))) < 1e-12
)
message("PASS: site sizes, X distributions, and size-weighted targets match the design.")

make_cell_data <- function(size_level, mix_level, seed) {
  simulate_multisite_tv(
    K = STUDY_K,
    n_per_site = get_site_sizes(size_level),
    init_intercepts = rep(STUDY_INIT_INTERCEPT, STUDY_K),
    site_mix = get_site_mix(mix_level),
    tau = 4L, t_star = STUDY_TSTAR,
    beta_event = set_beta_trt(STUDY_BETA_TRT),
    beta_init = scale_confounding(STUDY_CONF_MULT), base_seed = seed
  )
}

balanced_hom <- make_cell_data("balanced", "homogeneous", 8675309)
unbalanced_het <- make_cell_data("unbalanced", "heterogeneous", 135790)
stopifnot(
  nrow(balanced_hom) == 3000L, nrow(unbalanced_het) == 3000L,
  identical(as.integer(table(balanced_hom$site)), c(1000L, 1000L, 1000L)),
  identical(as.integer(table(unbalanced_het$site)), c(200L, 800L, 2000L))
)
obs_means <- tapply(unbalanced_het$x1, unbalanced_het$site, mean)
obs_sds <- tapply(unbalanced_het$x1, unbalanced_het$site, sd)
obs_x2 <- tapply(unbalanced_het$x2, unbalanced_het$site, mean)
stopifnot(
  max(abs(obs_means - het$x1_mean)) < 0.20,
  max(abs(obs_sds - het$x1_sd)) < 0.15,
  max(abs(obs_x2 - het$x2_prob)) < 0.08,
  abs(mean(unbalanced_het$x1) - 0.9) < 0.10,
  abs(mean(unbalanced_het$x2) - 0.58) < 0.05
)
message("PASS: simulated baseline X values recover the site-specific distributions.")

check_event_free_dgp <- function(dat, t_star) {
  stopifnot(!any(dat$S <= t_star & dat$S > dat$T_event))
  L1 <- as.matrix(dat[, paste0("L1_", seq_len(t_star)), drop = FALSE])
  L2 <- as.matrix(dat[, paste0("L2_", seq_len(t_star)), drop = FALSE])
  after_event <- outer(dat$T_event, seq_len(t_star), `<`)
  stopifnot(all(is.na(L1[after_event])), all(is.na(L2[after_event])))
  for (d in split(dat, dat$site)) {
    counts <- local_initiation_counts(d, t_star)
    expected_initiated <- vapply(seq_len(t_star), function(m)
      sum(d$S == m & d$T_event >= m), numeric(1))
    expected_eligible <- vapply(seq_len(t_star), function(m)
      sum(d$S >= m & d$T_event >= m), numeric(1))
    stopifnot(identical(counts$initiated, expected_initiated),
              identical(counts$eligible, expected_eligible))
  }
}
check_event_free_dgp(balanced_hom, STUDY_TSTAR)
check_event_free_dgp(unbalanced_het, STUDY_TSTAR)
message("PASS: both boundary cells obey the event-free initiation-risk DGP.")

site_data <- split(unbalanced_het, unbalanced_het$site)
hnum <- central_common_num_hazard(
  lapply(site_data, local_initiation_counts, t_star = STUDY_TSTAR)
)
stopifnot(max(abs(hnum - .tv_common_num_hazard(unbalanced_het, STUDY_TSTAR))) < 1e-12)
message("PASS: the common numerator is reconstructed from local counts.")

truth <- compute_truth_ccw_tv(
  N = 50000, tau = 4L, t_star = STUDY_TSTAR,
  beta_event = set_beta_trt(STUDY_BETA_TRT),
  beta_init = scale_confounding(STUDY_CONF_MULT),
  site_sizes = get_site_sizes("unbalanced"),
  init_intercepts = rep(STUDY_INIT_INTERCEPT, STUDY_K),
  site_mix = get_site_mix("heterogeneous"), seed = 24680
)
stopifnot(
  all(unname(truth$site_N) == c(3333, 13333, 33334)),
  max(abs(truth$site_weights - c(1, 4, 10) / 15)) < 1e-12,
  identical(truth$site_mix, get_site_mix("heterogeneous"))
)
message("PASS: oracle simulation uses the exact site-size-weighted patient mix.")

fit_warnings <- character()
out <- withCallingHandlers(
  run_once_tv(
    K = STUDY_K, n_per_site = get_site_sizes("unbalanced"),
    init_intercepts = rep(STUDY_INIT_INTERCEPT, STUDY_K),
    site_mix = get_site_mix("heterogeneous"), tau = 4L,
    t_star = STUDY_TSTAR, beta_trt = STUDY_BETA_TRT,
    beta_event = set_beta_trt(STUDY_BETA_TRT),
    beta_init = scale_confounding(STUDY_CONF_MULT),
    base_seed = 97531, truth = truth, verbose = FALSE
  ),
  warning = function(w) {
    fit_warnings <<- c(fit_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
unexpected_warnings <- fit_warnings[!grepl("^glm.fit:", fit_warnings)]
if (length(unexpected_warnings))
  stop("Unexpected warning(s): ", paste(unique(unexpected_warnings), collapse = "; "))

expected_methods <- c(
  "fed_ccw_tvipcw", "pooled_ccw_tvipcw", "fed_ipw_no_clone",
  "fed_perprotocol_naive", "fed_landmark_ipw", "local_ccw_meta"
)
stopifnot(
  identical(sort(unique(out$results$method)), sort(expected_methods)),
  nrow(out$results) == 6L * 8L,
  all(is.finite(out$results$estimate))
)
message("PASS: unbalanced/heterogeneous end-to-end run returns six finite methods.")
message("All ArraySim4-6 preflight checks passed.")
