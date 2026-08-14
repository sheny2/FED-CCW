# =============================================================
# ArraySim4-7 one-replicate comparison.
#
# The targeted methods are federated risk-set aggregation and local CCW curve
# meta-analysis. Pooled CCW is retained as a centralized benchmark.
# =============================================================

source("params.R")
source("DGP_tv.R")
source("Fed_CCW_TVIPCW.R")

run_once_tv <- function(
    n_per_site = STUDY_SITE_SIZES,
    init_intercept = get_init_intercept("sparse"),
    tau = STUDY_TAU,
    t_star = STUDY_TSTAR,
    beta_trt = STUDY_BETA_TRT,
    beta_event = DEFAULT_BETA_EVENT,
    beta_init = DEFAULT_BETA_INIT,
    beta_L = DEFAULT_BETA_L,
    sd_L = DEFAULT_SD_L,
    trunc = STUDY_TRUNC,
    base_seed = 2029,
    truth_N = STUDY_TRUTH_N,
    truth_seed = STUDY_TRUTH_SEED,
    truth = NULL,
    verbose = TRUE) {

  K <- length(n_per_site)
  site_mix <- homogeneous_site_mix(K)
  init_intercepts <- rep(init_intercept, K)
  beta_event <- set_beta_trt(beta_trt, beta_event)
  beta_init["int"] <- init_intercept

  if (is.null(truth)) {
    # Site distributions and initiation models are identical, so a one-site
    # oracle is the exact common target for every fragmentation level.
    truth <- compute_truth_ccw_tv(
      N = truth_N, tau = tau, t_star = t_star,
      beta_event = beta_event, beta_init = beta_init,
      site_sizes = sum(n_per_site), init_intercepts = init_intercept,
      site_mix = homogeneous_site_mix(1),
      beta_L = beta_L, sd_L = sd_L, seed = truth_seed
    )
  }

  dat <- simulate_multisite_tv(
    K = K, n_per_site = n_per_site,
    init_intercepts = init_intercepts, site_mix = site_mix,
    tau = tau, t_star = t_star, base_seed = base_seed,
    beta_event = beta_event, beta_init = beta_init,
    beta_L = beta_L, sd_L = sd_L
  )

  pair <- run_fed_and_local_ccw(dat, tau = tau, t_star = t_star,
                                trunc = trunc)
  fits <- list(
    fed_ccw_tvipcw = pair$fed,
    pooled_ccw_tvipcw = run_pooled_ccw_tvipcw(
      dat, tau = tau, t_star = t_star, trunc = trunc
    ),
    local_ccw_meta = pair$local_meta
  )

  estimands <- c("psi1", "psi0", "RD", "RR", "OR",
                 "RMST1", "RMST0", "RMST_diff")
  se_name <- c(RD = "SE_RD", RR = "SE_logRR", OR = "SE_logOR",
               RMST_diff = "SE_RMSTdiff", RMST1 = "SE_RMST1",
               RMST0 = "SE_RMST0")
  ci_name <- c(RD = "CI_RD", RR = "CI_RR", OR = "CI_OR",
               RMST_diff = "CI_RMSTdiff")

  make_rows <- function(fit, method) {
    do.call(rbind, lapply(estimands, function(est) {
      se <- fit[[se_name[est]]]
      if (is.null(se)) se <- NA_real_
      ci <- fit[[ci_name[est]]]
      if (is.null(ci)) ci <- c(NA_real_, NA_real_)
      data.frame(
        method = method, estimand = est,
        estimate = fit[[est]], truth = truth[[est]],
        bias = fit[[est]] - truth[[est]],
        se = se, ci_lo = ci[1], ci_hi = ci[2],
        covered = if (is.na(ci[1])) NA_integer_ else
          as.integer(truth[[est]] >= ci[1] & truth[[est]] <= ci[2]),
        stringsAsFactors = FALSE
      )
    }))
  }
  results <- do.call(rbind, Map(make_rows, fits, names(fits)))
  rownames(results) <- NULL

  stats <- pair$stats
  initiated_by_site <- vapply(split(dat, dat$site), function(d) sum(d$S <= tau),
                                integer(1))
  support <- list(
    n_sites = K,
    n_per_site = unique(n_per_site),
    zero_initiation_sites = sum(initiated_by_site == 0),
    mean_initiators_per_site = mean(initiated_by_site),
    min_initiators_per_site = min(initiated_by_site),
    max_initiators_per_site = max(initiated_by_site),
    empty_g1_at_tau = sum(vapply(stats, function(s) s$rw1[tau] <= 0,
                                logical(1))),
    empty_g1_at_tstar = sum(vapply(stats, function(s) s$rw1[t_star] <= 0,
                                  logical(1))),
    empty_g0_at_tstar = sum(vapply(stats, function(s) s$rw0[t_star] <= 0,
                                  logical(1))),
    any_post_tau_empty_g1 = sum(vapply(stats, function(s)
      any(s$rw1[seq.int(tau + 1L, t_star)] <= 0), logical(1))),
    observed_init_rate = mean(dat$A_tau),
    observed_death_rate = mean(dat$delta)
  )

  if (verbose) {
    cat(sprintf(
      paste0("K=%d, n/site=%d, initiation=%.3f; empty g=1 sites at t*=%d\n"),
      K, unique(n_per_site), support$observed_init_rate,
      support$empty_g1_at_tstar
    ))
    show <- results[results$estimand %in% c("RD", "RMST_diff"),
                    c("method", "estimand", "estimate", "truth", "bias")]
    print(show, digits = 4, row.names = FALSE)
  }

  invisible(list(results = results, truth = truth, fits = fits,
                 data = dat, support = support, pair = pair))
}
