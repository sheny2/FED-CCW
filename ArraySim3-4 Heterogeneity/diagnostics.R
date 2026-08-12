#!/usr/bin/env Rscript
# ============================================================================
# Diagnostics for the ArraySim3-4 federated-versus-pooled CCW discrepancy.
#
# Default (quick) run:
#   Rscript diagnostics.R
#
# Example with optional, more expensive diagnostics:
#   Rscript diagnostics.R --heterogeneity=high --tau=5 --beta-trt=-0.7 \
#     --run-convergence=true --convergence-reps=20 \
#     --run-oracle=true --oracle-N=300000
#
# All outputs are written to diagnostics/.
# ============================================================================

rm(list = ls())

# Run relative to this script, not the caller's submission directory.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
  setwd(dirname(script_path))
}

source("params.R")
source("DGP_tv.R")
source("Fed_CCW_TVIPCW.R")

# ---- command-line settings ----------------------------------------------
parse_cli <- function(x) {
  out <- list()
  for (arg in x) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) next
    z <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    out[[z[[1L]]]] <- paste(z[-1L], collapse = "=")
  }
  out
}

as_flag <- function(x) tolower(x) %in% c("1", "true", "yes", "y")
as_int_vector <- function(x) as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
CFG <- list(
  heterogeneity = if (!is.null(cli$heterogeneity)) cli$heterogeneity else "high",
  tau = if (!is.null(cli$tau)) as.integer(cli$tau) else 5L,
  t_star = if (!is.null(cli[["t-star"]])) as.integer(cli[["t-star"]]) else 25L,
  beta_trt = if (!is.null(cli[["beta-trt"]])) as.numeric(cli[["beta-trt"]]) else -0.7,
  n_per_site = if (!is.null(cli[["n-per-site"]])) {
    as.integer(cli[["n-per-site"]])
  } else 1000L,
  seed = if (!is.null(cli$seed)) as.integer(cli$seed) else 20260806L,
  run_convergence = if (!is.null(cli[["run-convergence"]])) {
    as_flag(cli[["run-convergence"]])
  } else FALSE,
  convergence_n = if (!is.null(cli[["convergence-n"]])) {
    as_int_vector(cli[["convergence-n"]])
  } else c(1000L, 5000L, 20000L),
  convergence_reps = if (!is.null(cli[["convergence-reps"]])) {
    as.integer(cli[["convergence-reps"]])
  } else 20L,
  run_oracle = if (!is.null(cli[["run-oracle"]])) {
    as_flag(cli[["run-oracle"]])
  } else FALSE,
  oracle_N = if (!is.null(cli[["oracle-N"]])) {
    as.integer(cli[["oracle-N"]])
  } else 300000L
)

if (!CFG$heterogeneity %in% names(DEFAULT_SITE_MIXES))
  stop("Unknown heterogeneity level: ", CFG$heterogeneity)
if (CFG$tau < 1L || CFG$tau > CFG$t_star)
  stop("tau must be between 1 and t_star.")
if (CFG$n_per_site < 20L) stop("n_per_site must be at least 20.")

out_dir <- "diagnostics"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_out <- function(x, filename) {
  path <- file.path(out_dir, filename)
  write.csv(x, path, row.names = FALSE)
  message("Wrote ", path)
}

message(sprintf(
  "Diagnostic cell: heterogeneity=%s, tau=%d, beta_trt=%g, n/site=%d, t*=%d",
  CFG$heterogeneity, CFG$tau, CFG$beta_trt, CFG$n_per_site, CFG$t_star
))

# ---- common data helpers -------------------------------------------------
interval_hazard <- function(S, M) {
  vapply(seq_len(M), function(m) {
    at_risk <- S >= m
    if (!any(at_risk)) return(NA_real_)
    .clamp_prob(sum(S == m) / sum(at_risk))
  }, numeric(1))
}

site_numerator_matrix <- function(dat, M) {
  sites <- sort(unique(dat$site))
  ans <- sapply(sites, function(k) interval_hazard(dat$S[dat$site == k], M))
  if (is.null(dim(ans))) ans <- matrix(ans, ncol = 1L)
  colnames(ans) <- paste0("site_", sites)
  ans
}

replace_numerator <- function(H, hnum) {
  for (m in seq_len(ncol(H$Hden))) {
    active <- !is.na(H$Hden[, m])
    H$Hnum[active, m] <- hnum[m]
  }
  H
}

true_hazard <- function(P, beta_init = DEFAULT_BETA_INIT) {
  n <- P$n
  M <- P$M
  Hden <- matrix(NA_real_, n, M)
  last_m <- pmin(P$S, M)

  for (m in seq_len(M)) {
    active <- m <= last_m
    lin <- beta_init["int"] +
      beta_init["x1"] * P$x1 +
      beta_init["x2"] * P$x2 +
      beta_init["L1"] * P$L1[, m] +
      beta_init["L2"] * P$L2[, m]
    Hden[active, m] <- .clamp_prob(plogis(lin[active]))
  }
  list(Hden = Hden, Hnum = matrix(NA_real_, n, M))
}

# Construct the same local sufficient statistics as local_ccw_tvipcw(), but
# permit diagnostic numerator/denominator hazards to be supplied explicitly.
local_stats_from_H <- function(site_data, H, tau, t_star,
                               trunc = DEFAULT_TRUNC) {
  P <- .tv_prep(site_data, t_star)
  M <- P$M
  Wt <- .tv_weights(P, H, tau, trunc)
  Cs <- .clone_censoring(P$S, tau, t_star)

  dw1 <- rw1 <- dw0 <- rw0 <- numeric(M)
  for (m in seq_len(M)) {
    tm1 <- m - 1L
    R1 <- P$Tcap > tm1 & Cs$C1 > tm1
    D1 <- R1 & P$Tev == m & P$Tev <= Cs$C1
    R0 <- P$Tcap > tm1 & Cs$C0 > tm1
    D0 <- R0 & P$Tev == m & P$Tev <= Cs$C0
    dw1[m] <- sum(Wt$SW1[, m] * D1)
    rw1[m] <- sum(Wt$SW1[, m] * R1)
    dw0[m] <- sum(Wt$SW0[, m] * D0)
    rw0[m] <- sum(Wt$SW0[, m] * R0)
  }

  list(
    M = M, dw1 = dw1, rw1 = rw1, dw0 = dw0, rw0 = rw0,
    local = list(
      SW1 = Wt$SW1, SW0 = Wt$SW0,
      C1 = Cs$C1, C0 = Cs$C0,
      Tev = P$Tev, Tcap = P$Tcap, n = P$n
    ),
    H = H
  )
}

run_local_den_common_num <- function(dat, tau, t_star, trunc = c(0, 1)) {
  hnum <- interval_hazard(dat$S, t_star)
  stats <- lapply(split(dat, dat$site), function(d) {
    P <- .tv_prep(d, t_star)
    H <- replace_numerator(.tv_init_hazard(P), hnum)
    local_stats_from_H(d, H, tau, t_star, trunc)
  })
  central_ccw_tvipcw(stats, tau, t_star)
}

run_true_den <- function(dat, tau, t_star, numerator = c("common", "local"),
                         trunc = c(0, 1)) {
  numerator <- match.arg(numerator)
  common_hnum <- interval_hazard(dat$S, t_star)
  local_hnum <- site_numerator_matrix(dat, t_star)
  sites <- sort(unique(dat$site))

  stats <- lapply(seq_along(sites), function(j) {
    d <- dat[dat$site == sites[j], , drop = FALSE]
    P <- .tv_prep(d, t_star)
    hnum <- if (numerator == "common") common_hnum else local_hnum[, j]
    H <- replace_numerator(true_hazard(P), hnum)
    local_stats_from_H(d, H, tau, t_star, trunc)
  })
  central_ccw_tvipcw(stats, tau, t_star)
}

# Pooled denominator with row-specific local numerator hazards. This is the
# complement of local-denominator/common-numerator and helps attribute a gap.
run_pooled_den_local_num <- function(dat, tau, t_star, trunc = c(0, 1)) {
  original_site <- dat$site
  local_hnum <- site_numerator_matrix(dat, t_star)
  sites <- sort(unique(original_site))

  pooled <- dat
  pooled$site <- 1L
  P <- .tv_prep(pooled, t_star)
  H <- .tv_init_hazard(P)

  for (j in seq_along(sites)) {
    rows <- original_site == sites[j]
    for (m in seq_len(t_star)) {
      active <- rows & !is.na(H$Hden[, m])
      H$Hnum[active, m] <- local_hnum[m, j]
    }
  }

  stats <- list(local_stats_from_H(pooled, H, tau, t_star, trunc))
  central_ccw_tvipcw(stats, tau, t_star)
}

extract_estimates <- function(fits) {
  estimands <- c("psi1", "psi0", "RD", "RR", "OR", "RMST1", "RMST0",
                 "RMST_diff")
  do.call(rbind, lapply(names(fits), function(nm) {
    data.frame(
      method = nm,
      estimand = estimands,
      estimate = vapply(estimands, function(e) fits[[nm]][[e]], numeric(1))
    )
  }))
}

# ---- generate the main diagnostic dataset -------------------------------
beta_event <- set_beta_trt(CFG$beta_trt)
dat <- simulate_multisite_tv(
  K = 3, n_per_site = CFG$n_per_site,
  tau = CFG$tau, t_star = CFG$t_star,
  heterogeneity = CFG$heterogeneity,
  beta_event = beta_event,
  base_seed = CFG$seed
)

# ---- 1. patient mix, uptake, and support --------------------------------
support <- do.call(rbind, lapply(split(dat, dat$site), function(d) {
  data.frame(
    site = unique(d$site),
    n = nrow(d),
    x1_mean = mean(d$x1),
    x1_sd = sd(d$x1),
    x2_prevalence = mean(d$x2),
    initiated_by_tau = mean(d$S <= CFG$tau),
    event_risk = mean(d$delta),
    n_initiated_by_tau = sum(d$S <= CFG$tau),
    n_not_initiated_by_tau = sum(d$S > CFG$tau)
  )
}))
write_out(support, "01_site_support.csv")

init_dist <- do.call(rbind, lapply(split(dat, dat$site), function(d) {
  tab <- table(factor(ifelse(d$S <= CFG$tau, d$S, CFG$tau + 1L),
                      levels = seq_len(CFG$tau + 1L)))
  data.frame(
    site = unique(d$site),
    initiation = c(paste0("m", seq_len(CFG$tau)), "after_tau_or_never"),
    n = as.integer(tab),
    proportion = as.numeric(prop.table(tab))
  )
}))
write_out(init_dist, "02_initiation_distribution.csv")

# ---- 2. local versus pooled numerator hazards ---------------------------
local_num <- site_numerator_matrix(dat, CFG$t_star)
common_num <- interval_hazard(dat$S, CFG$t_star)
num_df <- data.frame(
  interval = seq_len(CFG$t_star),
  common = common_num,
  local_num,
  check.names = FALSE
)
write_out(num_df, "03_numerator_hazards.csv")

png(file.path(out_dir, "03_numerator_hazards.png"), width = 1100, height = 700)
matplot(
  num_df$interval, as.matrix(num_df[, -1L]), type = "b", pch = 19,
  lty = 1, xlab = "Interval", ylab = "Numerator initiation hazard",
  main = "Local and common stabilizing numerator hazards"
)
legend("topright", legend = names(num_df)[-1L],
       col = seq_len(ncol(num_df) - 1L), lty = 1, pch = 19, bty = "n")
dev.off()

# ---- 3. fitted treatment models and probability diagnostics -------------
make_pi_df <- function(d, M) {
  P <- .tv_prep(d, M)
  last_m <- pmin(P$S, M)
  keep <- outer(seq_len(P$n), seq_len(M),
                function(i, m) m <= last_m[i])
  idx <- which(keep, arr.ind = TRUE)
  i <- idx[, "row"]
  m <- idx[, "col"]
  data.frame(
    y = as.integer(P$S[i] == m),
    m = factor(m),
    x1 = P$x1[i], x2 = P$x2[i],
    L1 = P$L1[keep], L2 = P$L2[keep],
    L1lag = .lag_matrix(P$L1)[keep],
    L2lag = .lag_matrix(P$L2)[keep]
  )
}

model_rows <- list()
prob_rows <- list()
fit_and_summarize <- function(d, label) {
  pi <- make_pi_df(d, CFG$t_star)
  fit <- glm(
    y ~ m + x1 + x2 + L1 + L2 + L1lag + L2lag,
    data = pi, family = binomial()
  )
  p <- .clamp_prob(predict(fit, type = "response"))
  model_rows[[length(model_rows) + 1L]] <<- data.frame(
    model = label, term = names(coef(fit)), estimate = unname(coef(fit))
  )
  prob_rows[[length(prob_rows) + 1L]] <<- data.frame(
    model = label,
    n_intervals = nrow(pi),
    event_rate = mean(pi$y),
    brier = mean((pi$y - p)^2),
    p_min = min(p), p_q01 = unname(quantile(p, 0.01)),
    p_q50 = median(p), p_q99 = unname(quantile(p, 0.99)),
    p_max = max(p),
    prop_p_lt_001 = mean(p < 0.01),
    prop_p_gt_099 = mean(p > 0.99)
  )
  invisible(fit)
}

invisible(lapply(split(dat, dat$site), function(d) {
  fit_and_summarize(d, paste0("site_", unique(d$site)))
}))
fit_and_summarize(dat, "pooled")
write_out(do.call(rbind, model_rows), "04_initiation_model_coefficients.csv")
write_out(do.call(rbind, prob_rows), "05_propensity_diagnostics.csv")

# ---- 4. method decomposition --------------------------------------------
fits <- list(
  fed_local_den_local_num = run_fed_ccw_tvipcw(
    dat, CFG$tau, CFG$t_star, trunc = c(0, 1)
  ),
  pooled_den_pooled_num = run_pooled_ccw_tvipcw(
    dat, CFG$tau, CFG$t_star, trunc = c(0, 1)
  ),
  fed_local_den_common_num = run_local_den_common_num(
    dat, CFG$tau, CFG$t_star, trunc = c(0, 1)
  ),
  pooled_den_local_num = run_pooled_den_local_num(
    dat, CFG$tau, CFG$t_star, trunc = c(0, 1)
  ),
  true_den_local_num = run_true_den(
    dat, CFG$tau, CFG$t_star, numerator = "local", trunc = c(0, 1)
  ),
  true_den_common_num = run_true_den(
    dat, CFG$tau, CFG$t_star, numerator = "common", trunc = c(0, 1)
  )
)
write_out(extract_estimates(fits), "06_method_decomposition.csv")

# ---- 5. weight, risk-set, and ESS diagnostics ----------------------------
weight_rows <- list()
prop_rows2 <- list()
for (num_type in c("local", "common")) {
  common_h <- interval_hazard(dat$S, CFG$t_star)
  for (d in split(dat, dat$site)) {
    site <- unique(d$site)
    P <- .tv_prep(d, CFG$t_star)
    H <- .tv_init_hazard(P)
    if (num_type == "common") H <- replace_numerator(H, common_h)
    st <- local_stats_from_H(d, H, CFG$tau, CFG$t_star, c(0, 1))
    L <- st$local

    for (m in seq_len(CFG$t_star)) {
      active <- !is.na(H$Hden[, m])
      if (any(active)) {
        pred <- H$Hden[active, m]
        prop_rows2[[length(prop_rows2) + 1L]] <- data.frame(
          numerator = num_type, site = site, interval = m,
          n_treatment_risk = sum(active),
          observed_initiations = sum(P$S == m),
          pred_min = min(pred), pred_q01 = unname(quantile(pred, 0.01)),
          pred_q50 = median(pred), pred_q99 = unname(quantile(pred, 0.99)),
          pred_max = max(pred),
          prop_extreme = mean(pred < 0.01 | pred > 0.99)
        )
      }

      tm1 <- m - 1L
      for (arm in c(1L, 0L)) {
        W <- if (arm == 1L) L$SW1[, m] else L$SW0[, m]
        Cg <- if (arm == 1L) L$C1 else L$C0
        risk <- L$Tcap > tm1 & Cg > tm1
        wr <- W[risk]
        ess <- if (length(wr) && sum(wr^2) > 0) {
          sum(wr)^2 / sum(wr^2)
        } else NA_real_
        weight_rows[[length(weight_rows) + 1L]] <- data.frame(
          numerator = num_type, site = site, arm = arm, interval = m,
          risk_n = length(wr), weight_sum = sum(wr),
          ess = ess,
          ess_fraction = if (length(wr)) ess / length(wr) else NA_real_,
          w_median = if (length(wr)) median(wr) else NA_real_,
          w_q95 = if (length(wr)) unname(quantile(wr, 0.95)) else NA_real_,
          w_q99 = if (length(wr)) unname(quantile(wr, 0.99)) else NA_real_,
          w_max = if (length(wr)) max(wr) else NA_real_
        )
      }
    }
  }
}
weight_diag <- do.call(rbind, weight_rows)
write_out(weight_diag, "07_weight_ess_diagnostics.csv")
write_out(do.call(rbind, prop_rows2), "08_interval_propensity_diagnostics.csv")

weight_flags <- subset(
  weight_diag,
  is.finite(ess_fraction) &
    (ess_fraction < 0.20 | w_max > 50 | risk_n < 30)
)
write_out(weight_flags, "09_weight_warnings.csv")

# ---- 6. truncation sensitivity ------------------------------------------
truncations <- list(
  none = c(0, 1),
  p01_p99 = c(0.01, 0.99),
  p025_p975 = c(0.025, 0.975)
)
trunc_fits <- list()
for (nm in names(truncations)) {
  tr <- truncations[[nm]]
  trunc_fits[[paste0("fed_", nm)]] <- run_fed_ccw_tvipcw(
    dat, CFG$tau, CFG$t_star, trunc = tr
  )
  trunc_fits[[paste0("pooled_", nm)]] <- run_pooled_ccw_tvipcw(
    dat, CFG$tau, CFG$t_star, trunc = tr
  )
  trunc_fits[[paste0("fed_common_num_", nm)]] <- run_local_den_common_num(
    dat, CFG$tau, CFG$t_star, trunc = tr
  )
}
write_out(extract_estimates(trunc_fits), "10_truncation_sensitivity.csv")

# ---- 7. required controls -----------------------------------------------
# K=1: fed and pooled must agree to floating-point tolerance.
one_site <- dat[dat$site == sort(unique(dat$site))[1L], , drop = FALSE]
k1_fed <- run_fed_ccw_tvipcw(one_site, CFG$tau, CFG$t_star, c(0, 1))
k1_pool <- run_pooled_ccw_tvipcw(one_site, CFG$tau, CFG$t_star, c(0, 1))

# Completely homogeneous three-site patient mix.
homogeneous <- do.call(rbind, lapply(1:3, function(k) {
  simulate_site_tv(
    n = CFG$n_per_site, tau = CFG$tau, t_star = CFG$t_star,
    beta_event = beta_event, site_id = k,
    x1_mean = 0, x1_sd = 1, x2_prob = 0.4,
    seed = CFG$seed + 1000L + k
  )
}))
hom_fed <- run_fed_ccw_tvipcw(homogeneous, CFG$tau, CFG$t_star, c(0, 1))
hom_pool <- run_pooled_ccw_tvipcw(homogeneous, CFG$tau, CFG$t_star, c(0, 1))

controls <- data.frame(
  check = c("K1_identity", "K1_identity",
            "homogeneous_three_sites", "homogeneous_three_sites"),
  estimand = rep(c("RD", "RMST_diff"), 2L),
  fed = c(k1_fed$RD, k1_fed$RMST_diff, hom_fed$RD, hom_fed$RMST_diff),
  pooled = c(k1_pool$RD, k1_pool$RMST_diff, hom_pool$RD, hom_pool$RMST_diff)
)
controls$difference <- controls$fed - controls$pooled
controls$passed <- ifelse(
  controls$check == "K1_identity",
  abs(controls$difference) < 1e-10,
  NA
)
write_out(controls, "11_identity_and_homogeneous_controls.csv")

# ---- 8. optional sample-size convergence --------------------------------
if (CFG$run_convergence) {
  conv_rows <- list()
  counter <- 0L
  for (n in CFG$convergence_n) {
    for (r in seq_len(CFG$convergence_reps)) {
      counter <- counter + 1L
      d <- simulate_multisite_tv(
        K = 3, n_per_site = n,
        tau = CFG$tau, t_star = CFG$t_star,
        heterogeneity = CFG$heterogeneity,
        beta_event = beta_event,
        base_seed = CFG$seed + 10000L + counter
      )
      f <- run_fed_ccw_tvipcw(d, CFG$tau, CFG$t_star, c(0, 1))
      p <- run_pooled_ccw_tvipcw(d, CFG$tau, CFG$t_star, c(0, 1))
      fc <- run_local_den_common_num(d, CFG$tau, CFG$t_star, c(0, 1))
      conv_rows[[counter]] <- data.frame(
        n_per_site = n, rep = r,
        fed_RD = f$RD, pooled_RD = p$RD, fed_common_RD = fc$RD,
        fed_RMST = f$RMST_diff, pooled_RMST = p$RMST_diff,
        fed_common_RMST = fc$RMST_diff,
        fed_minus_pooled_RMST = f$RMST_diff - p$RMST_diff,
        fed_common_minus_pooled_RMST = fc$RMST_diff - p$RMST_diff
      )
      message(sprintf("Convergence: n/site=%d replicate=%d/%d",
                      n, r, CFG$convergence_reps))
    }
  }
  conv <- do.call(rbind, conv_rows)
  write_out(conv, "12_sample_size_convergence.csv")

  conv_summary <- aggregate(
    cbind(fed_minus_pooled_RMST, fed_common_minus_pooled_RMST) ~ n_per_site,
    conv,
    function(x) c(mean = mean(x), sd = sd(x), rmse = sqrt(mean(x^2)))
  )
  write_out(conv_summary, "13_sample_size_convergence_summary.csv")
}

# ---- 9. optional large-sample oracle-alignment check --------------------
if (CFG$run_oracle) {
  n_oracle_site <- max(1000L, floor(CFG$oracle_N / 3L))
  message("Generating large diagnostic cohort: ", 3L * n_oracle_site,
          " patients.")
  oracle_dat <- simulate_multisite_tv(
    K = 3, n_per_site = n_oracle_site,
    tau = CFG$tau, t_star = CFG$t_star,
    heterogeneity = CFG$heterogeneity,
    beta_event = beta_event,
    base_seed = CFG$seed + 500000L
  )
  aggregated_true <- run_true_den(
    oracle_dat, CFG$tau, CFG$t_star,
    numerator = "common", trunc = c(0, 1)
  )
  existing_site_average <- compute_truth_multisite_tv(
    N = 3L * n_oracle_site, K = 3,
    heterogeneity = CFG$heterogeneity,
    tau = CFG$tau, t_star = CFG$t_star,
    beta_event = beta_event,
    seed = CFG$seed + 600000L
  )
  oracle_compare <- extract_estimates(list(
    aggregated_true_weights_common_num = aggregated_true,
    existing_site_average_oracle = existing_site_average
  ))
  write_out(oracle_compare, "14_oracle_alignment.csv")
}

# ---- concise console interpretation -------------------------------------
decomp <- subset(
  extract_estimates(fits),
  estimand %in% c("RD", "RMST_diff")
)
message("\nKey decomposition:")
print(decomp, row.names = FALSE)
message("\nK=1 identity differences:")
print(subset(controls, check == "K1_identity"), row.names = FALSE)
message(sprintf(
  "\nWeight warnings: %d rows (ESS fraction <0.20, max weight >50, or risk set <30).",
  nrow(weight_flags)
))
message("\nDiagnostics complete. See diagnostics/ for CSV files and plots.")

