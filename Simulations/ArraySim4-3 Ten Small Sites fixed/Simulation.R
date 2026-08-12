# =============================================================
# ArraySim4-3: one replicate comparing six estimators across ten small sites
# with site-specific treatment-initiation prevalence:
#   common-numerator federated CCW, aligned pooled CCW,
#   federated IPW without cloning, federated unweighted PP,
#   federated landmark IPW, and local CCW curve meta-analysis.
# =============================================================

source("params.R")
source("DGP_tv.R")
source("Fed_CCW_TVIPCW.R")

run_once_tv <- function(K          = 10,
                        n_per_site = TEN_SITE_SIZES,
                        init_intercepts = TEN_SITE_INIT_INTERCEPTS,
                        tau        = DEFAULT_TAU,
                        t_star     = DEFAULT_TSTAR,
                        beta_trt   = DEFAULT_BETA_EVENT[["trt"]],
                        beta_event = DEFAULT_BETA_EVENT,
                        beta_init  = DEFAULT_BETA_INIT,
                        beta_L     = DEFAULT_BETA_L,
                        sd_L       = DEFAULT_SD_L,
                        trunc      = DEFAULT_TRUNC,
                        base_seed  = 2024,
                        truth_N    = 2e6,
                        truth_seed = 999,
                        truth      = NULL,
                        verbose    = TRUE) {

  beta_event <- set_beta_trt(beta_trt, beta_event)

  if (is.null(truth)) {
    truth <- compute_truth_ccw_tv(
      N = truth_N, tau = tau, t_star = t_star,
      beta_event = beta_event, beta_init = beta_init,
      site_sizes = n_per_site, init_intercepts = init_intercepts,
      beta_L = beta_L, sd_L = sd_L, seed = truth_seed
    )
  }

  dat <- simulate_multisite_tv(
    K = K, n_per_site = n_per_site,
    init_intercepts = init_intercepts,
    tau = tau, t_star = t_star, base_seed = base_seed,
    beta_event = beta_event, beta_init = beta_init,
    beta_L = beta_L, sd_L = sd_L
  )

  fits <- list(
    fed_ccw_tvipcw = run_fed_ccw_tvipcw(
      dat, tau = tau, t_star = t_star, trunc = trunc
    ),
    pooled_ccw_tvipcw = run_pooled_ccw_tvipcw(
      dat, tau = tau, t_star = t_star, trunc = trunc
    ),
    fed_ipw_no_clone = run_fed_ipw_nocloning(
      dat, tau = tau, t_star = t_star, trunc = trunc
    ),
    fed_perprotocol_naive = run_fed_perprotocol(
      dat, tau = tau, t_star = t_star
    ),
    fed_landmark_ipw = run_fed_landmark_ipw(
      dat, tau = tau, t_star = t_star, trunc = trunc
    ),
    local_ccw_meta = run_local_ccw_meta(
      dat, tau = tau, t_star = t_star, trunc = trunc
    )
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

  res <- do.call(rbind, Map(make_rows, fits, names(fits)))
  rownames(res) <- NULL

  if (verbose) {
    cat("\n==== Ground truth ====\n")
    cat(sprintf("psi1=%.4f psi0=%.4f RD=%.4f RR=%.4f RMSTd=%.4f\n",
                truth$psi1, truth$psi0, truth$RD, truth$RR,
                truth$RMST_diff))
    cat("\n==== Estimates vs truth ====\n")
    show <- res[
      res$estimand %in% c("RD", "RR", "OR", "RMST_diff"),
      c("method", "estimand", "estimate", "truth", "bias",
        "se", "ci_lo", "ci_hi", "covered")
    ]
    print(show, digits = 4, row.names = FALSE)
  }

  invisible(list(results = res, truth = truth, fits = fits, data = dat))
}
