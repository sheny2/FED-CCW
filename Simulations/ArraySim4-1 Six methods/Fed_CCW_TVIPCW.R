# =============================================================
# Fed_CCW_TVIPCW.R
#   Federated clone-censor-weight with time-varying IPCW for the
#   time-varying-confounding DGP (simulate_*_tv).
#
#   Requires wide L columns L1_1..L1_M, L2_1..L2_M alongside the baseline
#   columns (site, id, x1, x2, S, A_tau, T_event, T_obs, delta).
#
# WEIGHTING MODEL
#   Artificial censoring in CCW is deterministic given (S, tau), so a clone's
#   censoring hazard in interval m is the hazard of the deviation event, i.e.
#   treatment initiation among the not-yet-initiated. One pooled-over-intervals
#   logistic initiation-hazard model is fit per site on person-interval rows
#   using baseline X, current L_m and lagged L_{m-1}. The stabilizing
#   numerator is calculated once from federated interval-level counts and
#   broadcast to every site, ensuring a common target across sites.
#
#   Strategy weights built from the per-interval initiation hazard h_m:
#     g = 1 (initiate by tau): weight accrues over m = 1..min(S, tau).
#     g = 0 (never initiate by tau): weight accrues over m = 1..min(tau, S-1).
#   After tau no further artificial censoring occurs, so weights are carried
#   forward unchanged. Weights are stabilized and truncated at local
#   percentiles.
# =============================================================

source("params.R")

# ---- extract the per-site matrices needed downstream --------------------
.tv_prep <- function(d, t_star) {
  M <- t_star
  list(n  = nrow(d), M = M,
       L1 = as.matrix(d[, paste0("L1_", seq_len(M)), drop = FALSE]),
       L2 = as.matrix(d[, paste0("L2_", seq_len(M)), drop = FALSE]),
       x1 = d$x1, x2 = d$x2, S = d$S,
       Tev = d$T_event, Tcap = pmin(d$T_event, M))
}

# ---- lagged L matrix, with L_0 = 0 to match the DGP initialization ------
.lag_matrix <- function(Lmat) {
  cbind(0, Lmat[, -ncol(Lmat), drop = FALSE])
}

# ---- federated common numerator hazard ----------------------------------
# Each site need only share interval-level initiation and at-risk counts.
# Their sums define one common stabilizing numerator, so every site targets
# the same distribution of treatment-initiation times despite different
# patient mixes.
local_initiation_counts <- function(site_data, t_star) {
  list(
    initiated = vapply(seq_len(t_star), function(m) sum(site_data$S == m),
                       numeric(1)),
    eligible = vapply(seq_len(t_star), function(m) sum(site_data$S >= m),
                      numeric(1))
  )
}

central_common_num_hazard <- function(local_counts) {
  initiated <- Reduce(`+`, lapply(local_counts, `[[`, "initiated"))
  eligible <- Reduce(`+`, lapply(local_counts, `[[`, "eligible"))
  ifelse(eligible > 0, .clamp_prob(initiated / eligible), NA_real_)
}

.tv_common_num_hazard <- function(dat, t_star) {
  central_common_num_hazard(
    lapply(split(dat, dat$site), local_initiation_counts, t_star = t_star)
  )
}

# ---- per-site initiation-hazard models ---------------------------------
# Denominator conditions on X, current L_m and lagged L_{m-1}, matching the
# outcome hazard's dependence on current + lagged L. If `hnum` is supplied,
# it is the common central numerator vector broadcast to all sites.
.tv_init_hazard <- function(P, hnum = NULL) {
  n <- P$n; M <- P$M

  # person-interval rows: m = 1..min(S, M) for each person
  last_m <- pmin(P$S, M)
  keep   <- outer(seq_len(n), seq_len(M), function(i, m) m <= last_m[i])

  idx     <- which(keep, arr.ind = TRUE)
  rows_id <- idx[, "row"]
  rows_m  <- idx[, "col"]

  L1lag <- .lag_matrix(P$L1)
  L2lag <- .lag_matrix(P$L2)

  pi_df <- data.frame(
    y     = as.integer(P$S[rows_id] == rows_m),
    m     = factor(rows_m, levels = seq_len(M)),
    x1    = P$x1[rows_id],
    x2    = P$x2[rows_id],
    L1    = P$L1[keep],
    L2    = P$L2[keep],
    L1lag = L1lag[keep],
    L2lag = L2lag[keep],
    id    = rows_id
  )
  pi_df$m <- droplevels(pi_df$m)

  den_fit <- glm(y ~ m + x1 + x2 + L1 + L2 + L1lag + L2lag,
                 data = pi_df, family = binomial())
  Hden <- matrix(NA_real_, n, M)
  Hnum <- matrix(NA_real_, n, M)
  cell <- cbind(rows_id, rows_m)
  Hden[cell] <- .clamp_prob(predict(den_fit, type = "response"))
  if (is.null(hnum)) {
    num_fit <- glm(y ~ m, data = pi_df, family = binomial())
    Hnum[cell] <- .clamp_prob(predict(num_fit, type = "response"))
  } else {
    if (length(hnum) != M)
      stop("Common numerator must have one value per follow-up interval.")
    Hnum[cell] <- .clamp_prob(hnum[rows_m])
  }

  list(Hden = Hden, Hnum = Hnum)
}

# ---- stabilized, truncated strategy weights ----------------------------
# Returns n x M matrices SW1, SW0; entry [i, m] is the weight applicable
# DURING interval m.
.tv_weights <- function(P, H, tau, trunc = DEFAULT_TRUNC) {
  n <- P$n; M <- P$M
  Hden <- H$Hden; Hnum <- H$Hnum

  # per-interval factors, defaulting to 1 where no hazard is defined or m > tau
  ratio_init <- matrix(1, n, M)   # used in the interval where S == m
  ratio_stay <- matrix(1, n, M)   # used while still uninitiated

  active <- !is.na(Hden)
  active[, seq_len(M) > tau] <- FALSE

  ratio_init[active] <- Hnum[active] / Hden[active]
  ratio_stay[active] <- (1 - Hnum[active]) / (1 - Hden[active])

  init_at_m <- outer(P$S, seq_len(M), `==`)

  # g = 1: initiation ratio in the interval of initiation, stay-ratio before
  F1 <- ifelse(init_at_m, ratio_init, ratio_stay)
  F1[!active] <- 1

  # g = 0: stay-ratio only; once the person initiates the clone leaves the
  # risk set and the weight simply stops accruing
  F0 <- ratio_stay
  F0[!active | init_at_m] <- 1

  SW1 <- t(apply(F1, 1, cumprod))
  SW0 <- t(apply(F0, 1, cumprod))

  clip <- function(W) {
    q <- quantile(W[is.finite(W)], probs = trunc, names = FALSE)
    pmin(pmax(W, q[1]), q[2])
  }
  list(SW1 = clip(SW1), SW0 = clip(SW0))
}

# ---- clone censoring times ---------------------------------------------
# Deviation occurs IN the interval, so the clone is censored at the end of
# the previous one.
#   C1 (g = 1, "initiate by tau"): a clone that never initiates by tau
#       deviates AT tau, so it is protocol-consistent through tau-1.
#   C0 (g = 0, "never initiate by tau"): a clone that initiates at S <= tau
#       deviates IN interval S, so it is protocol-consistent through S-1.
# Cinf = t_star + 1 means "never artificially censored".
.clone_censoring <- function(S, tau, t_star) {
  Cinf <- t_star + 1L
  list(C1 = ifelse(S >  tau, tau - 1L, Cinf),
       C0 = ifelse(S <= tau, S   - 1L, Cinf))
}

# -------------------------------------------------------------
# LOCAL site computation
# -------------------------------------------------------------
local_ccw_tvipcw <- function(site_data, tau, t_star, trunc = DEFAULT_TRUNC,
                             hnum = NULL) {
  P <- .tv_prep(site_data, t_star)
  M <- P$M

  Wt  <- .tv_weights(P, .tv_init_hazard(P, hnum = hnum),
                     tau, trunc = trunc)
  SW1 <- Wt$SW1; SW0 <- Wt$SW0

  Cs <- .clone_censoring(P$S, tau, t_star)
  C1 <- Cs$C1; C0 <- Cs$C0

  dw1 <- numeric(M); rw1 <- numeric(M)
  dw0 <- numeric(M); rw0 <- numeric(M)
  Tev <- P$Tev; Tcap <- P$Tcap

  for (m in seq_len(M)) {
    tm1 <- m - 1
    w1 <- SW1[, m]; w0 <- SW0[, m]

    R1m <- as.integer(Tcap > tm1 & C1 > tm1)
    D1m <- as.integer(R1m == 1 & Tev > tm1 & Tev <= m & Tev <= C1)
    dw1[m] <- sum(w1 * D1m); rw1[m] <- sum(w1 * R1m)

    R0m <- as.integer(Tcap > tm1 & C0 > tm1)
    D0m <- as.integer(R0m == 1 & Tev > tm1 & Tev <= m & Tev <= C0)
    dw0[m] <- sum(w0 * D0m); rw0[m] <- sum(w0 * R0m)
  }

  list(
    M = M,
    dw1 = dw1, rw1 = rw1, dw0 = dw0, rw0 = rw0,
    local = list(SW1 = SW1, SW0 = SW0, C1 = C1, C0 = C0,
                 Tev = Tev, Tcap = Tcap, n = P$n)
  )
}

# -------------------------------------------------------------
# LOCAL second pass: influence functions.
#   Same algebra as the time-constant version with SW_{i,m} in place of W_i.
#   Weights are treated as known, so the variance is conservative.
# -------------------------------------------------------------
local_influence_tv <- function(stats_local, lambda1, lambda0,
                               S1_curve, S0_curve, RW1, RW0, M,
                               psi1, psi0, dt) {
  L <- stats_local; n <- L$n
  S1_tstar <- S1_curve[M + 1]; S0_tstar <- S0_curve[M + 1]

  cumU1 <- numeric(n); cumU0 <- numeric(n)
  IF_RMST1 <- numeric(n); IF_RMST0 <- numeric(n)

  for (m in seq_len(M)) {
    tm1 <- m - 1

    IF_RMST1 <- IF_RMST1 + dt[m] * (-S1_curve[m] * cumU1)
    IF_RMST0 <- IF_RMST0 + dt[m] * (-S0_curve[m] * cumU0)

    w1 <- L$SW1[, m]; w0 <- L$SW0[, m]

    R1m <- as.integer(L$Tcap > tm1 & L$C1 > tm1)
    D1m <- as.integer(R1m == 1 & L$Tev > tm1 & L$Tev <= m & L$Tev <= L$C1)
    if (RW1[m] > 0 && lambda1[m] < 1) {
      cumU1 <- cumU1 + (1 / (1 - lambda1[m])) * (w1 / RW1[m]) *
        (D1m - lambda1[m] * R1m)
    }

    R0m <- as.integer(L$Tcap > tm1 & L$C0 > tm1)
    D0m <- as.integer(R0m == 1 & L$Tev > tm1 & L$Tev <= m & L$Tev <= L$C0)
    if (RW0[m] > 0 && lambda0[m] < 1) {
      cumU0 <- cumU0 + (1 / (1 - lambda0[m])) * (w0 / RW0[m]) *
        (D0m - lambda0[m] * R0m)
    }
  }

  IF1 <- S1_tstar * cumU1
  IF0 <- S0_tstar * cumU0

  list(
    n = n,
    sum_RD2       = sum((IF1 - IF0)^2),
    sum_logRR2    = sum((IF1 / psi1 - IF0 / psi0)^2),
    sum_logOR2    = sum((IF1 / (psi1 * (1 - psi1)) -
                         IF0 / (psi0 * (1 - psi0)))^2),
    sum_RMSTdiff2 = sum((IF_RMST1 - IF_RMST0)^2),
    sum_RMST1_2   = sum(IF_RMST1^2),
    sum_RMST0_2   = sum(IF_RMST0^2)
  )
}

# ---- shared contrast constructor ---------------------------------------
.contrasts <- function(S1, S0, M, dt = rep(1, M)) {
  psi1 <- 1 - S1[M + 1]; psi0 <- 1 - S0[M + 1]
  RMST1 <- sum(S1[seq_len(M)] * dt); RMST0 <- sum(S0[seq_len(M)] * dt)
  list(S1 = S1, S0 = S0, psi1 = psi1, psi0 = psi0,
       RD = psi1 - psi0,
       RR = psi1 / psi0,
       OR = (psi1 / (1 - psi1)) / (psi0 / (1 - psi0)),
       RMST1 = RMST1, RMST0 = RMST0, RMST_diff = RMST1 - RMST0)
}

# -------------------------------------------------------------
# CENTRAL aggregation + estimation
# -------------------------------------------------------------
central_ccw_tvipcw <- function(site_stats_list, tau, t_star) {
  M <- t_star; dt <- rep(1, M)

  DW1 <- Reduce(`+`, lapply(site_stats_list, `[[`, "dw1"))
  RW1 <- Reduce(`+`, lapply(site_stats_list, `[[`, "rw1"))
  DW0 <- Reduce(`+`, lapply(site_stats_list, `[[`, "dw0"))
  RW0 <- Reduce(`+`, lapply(site_stats_list, `[[`, "rw0"))

  lambda1 <- ifelse(RW1 > 0, DW1 / RW1, 0)
  lambda0 <- ifelse(RW0 > 0, DW0 / RW0, 0)

  S1 <- c(1, cumprod(1 - lambda1))
  S0 <- c(1, cumprod(1 - lambda0))
  est <- .contrasts(S1, S0, M, dt)

  inf_list <- lapply(site_stats_list, function(s) {
    local_influence_tv(s$local, lambda1, lambda0, S1, S0,
                       RW1, RW0, M, est$psi1, est$psi0, dt)
  })
  agg <- function(f) sqrt(sum(vapply(inf_list, `[[`, numeric(1), f)))

  SE_RD       <- agg("sum_RD2")
  SE_logRR    <- agg("sum_logRR2")
  SE_logOR    <- agg("sum_logOR2")
  SE_RMSTdiff <- agg("sum_RMSTdiff2")

  z <- qnorm(0.975)
  c(est, list(
    lambda1 = lambda1, lambda0 = lambda0,
    SE_RD = SE_RD, SE_logRR = SE_logRR, SE_logOR = SE_logOR,
    SE_RMSTdiff = SE_RMSTdiff,
    SE_RMST1 = agg("sum_RMST1_2"), SE_RMST0 = agg("sum_RMST0_2"),
    CI_RD       = est$RD + c(-1, 1) * z * SE_RD,
    CI_RR       = exp(log(est$RR) + c(-1, 1) * z * SE_logRR),
    CI_OR       = exp(log(est$OR) + c(-1, 1) * z * SE_logOR),
    CI_RMSTdiff = est$RMST_diff + c(-1, 1) * z * SE_RMSTdiff
  ))
}

run_fed_ccw_tvipcw <- function(dat, tau, t_star, trunc = DEFAULT_TRUNC) {
  hnum <- .tv_common_num_hazard(dat, t_star)
  stats <- lapply(split(dat, dat$site), local_ccw_tvipcw,
                  tau = tau, t_star = t_star, trunc = trunc, hnum = hnum)
  central_ccw_tvipcw(stats, tau = tau, t_star = t_star)
}

# -------------------------------------------------------------
# COMPARATOR: pooled (oracle-federation) time-varying IPCW.
#   Same clone-censor + TV-IPCW pipeline, but all sites collapsed into one
#   dataset: a single initiation-hazard model and a single set of sufficient
#   statistics. A gap between fed and pooled is not causal bias, it is the
#   cost of fitting the weight model separately per site.
# -------------------------------------------------------------
run_pooled_ccw_tvipcw <- function(dat, tau, t_star, trunc = DEFAULT_TRUNC) {
  hnum <- .tv_common_num_hazard(dat, t_star)
  pooled <- dat
  pooled$site <- 1L
  stats <- list(local_ccw_tvipcw(pooled, tau = tau, t_star = t_star,
                                 trunc = trunc, hnum = hnum))
  central_ccw_tvipcw(stats, tau = tau, t_star = t_star)
}

# -------------------------------------------------------------
# FEDERATED NAIVE per-protocol comparator. Each site releases only event and
# risk-set counts by observed A_tau and interval. Central aggregation exactly
# reproduces the corresponding pooled crude survival estimator.
# -------------------------------------------------------------
local_perprotocol_counts <- function(site_data, tau, t_star) {
  count_arm <- function(arm) {
    sub <- site_data[site_data$A_tau == arm, , drop = FALSE]
    Tev <- sub$T_event
    Tcap <- pmin(Tev, t_star)
    event <- risk <- numeric(t_star)
    for (m in seq_len(t_star)) {
      tm1 <- m - 1
      Rm <- as.integer(Tcap > tm1)
      Dm <- as.integer(Rm == 1 & Tev > tm1 & Tev <= m)
      event[m] <- sum(Dm)
      risk[m] <- sum(Rm)
    }
    list(event = event, risk = risk)
  }
  list(arm1 = count_arm(1L), arm0 = count_arm(0L))
}

central_perprotocol <- function(site_counts, t_star) {
  sum_field <- function(arm, field) {
    Reduce(`+`, lapply(site_counts, function(x) x[[arm]][[field]]))
  }
  survival <- function(arm) {
    event <- sum_field(arm, "event")
    risk <- sum_field(arm, "risk")
    hazard <- ifelse(risk > 0, event / risk, 0)
    c(1, cumprod(1 - hazard))
  }
  .contrasts(survival("arm1"), survival("arm0"), t_star)
}

run_fed_perprotocol <- function(dat, tau, t_star) {
  site_counts <- lapply(split(dat, dat$site), local_perprotocol_counts,
                        tau = tau, t_star = t_star)
  central_perprotocol(site_counts, t_star = t_star)
}

# -------------------------------------------------------------
# COMPARATOR: time-varying IPCW WITHOUT cloning.
#   Same per-interval stabilized, truncated TV weights as fed-CCW, but no
#   cloning and no artificial censoring. The observed cohort is split by
#   realized A_tau and each group's crude discrete-time survival is estimated
#   with its own strategy weight (SW1 for initiators, SW0 for the rest).
#
#   Purpose: isolate the contribution of cloning. This is an as-treated
#   analysis, so it conditions on realized initiation among survivors and
#   carries immortal time bias that the weights do not remove. Point
#   estimates only.
# -------------------------------------------------------------
local_ipw_nocloning <- function(site_data, tau, t_star, hnum,
                                trunc = DEFAULT_TRUNC) {
  M <- t_star
  DW1 <- numeric(M); RW1 <- numeric(M)
  DW0 <- numeric(M); RW0 <- numeric(M)

  P  <- .tv_prep(site_data, t_star)
  Wt <- .tv_weights(P, .tv_init_hazard(P, hnum = hnum),
                    tau, trunc = trunc)

  Tev <- P$Tev; Tcap <- P$Tcap
  init <- P$S <= tau

  for (m in seq_len(M)) {
    tm1 <- m - 1

    Rm1 <- as.integer(Tcap[init] > tm1)
    Dm1 <- as.integer(Rm1 == 1 & Tev[init] > tm1 & Tev[init] <= m)
    w1  <- Wt$SW1[init, m]
    DW1[m] <- sum(w1 * Dm1); RW1[m] <- sum(w1 * Rm1)

    Rm0 <- as.integer(Tcap[!init] > tm1)
    Dm0 <- as.integer(Rm0 == 1 & Tev[!init] > tm1 & Tev[!init] <= m)
    w0  <- Wt$SW0[!init, m]
    DW0[m] <- sum(w0 * Dm0); RW0[m] <- sum(w0 * Rm0)
  }

  list(dw1 = DW1, rw1 = RW1, dw0 = DW0, rw0 = RW0)
}

central_ipw_nocloning <- function(site_stats, t_star) {
  add <- function(field) Reduce(`+`, lapply(site_stats, `[[`, field))
  DW1 <- add("dw1"); RW1 <- add("rw1")
  DW0 <- add("dw0"); RW0 <- add("rw0")

  lam1 <- ifelse(RW1 > 0, DW1 / RW1, 0)
  lam0 <- ifelse(RW0 > 0, DW0 / RW0, 0)
  .contrasts(c(1, cumprod(1 - lam1)),
             c(1, cumprod(1 - lam0)), t_star)
}

run_fed_ipw_nocloning <- function(dat, tau, t_star,
                                  trunc = DEFAULT_TRUNC) {
  site_data <- split(dat, dat$site)
  hnum <- central_common_num_hazard(
    lapply(site_data, local_initiation_counts, t_star = t_star)
  )
  site_stats <- lapply(site_data, local_ipw_nocloning,
                       tau = tau, t_star = t_star,
                       hnum = hnum, trunc = trunc)
  central_ipw_nocloning(site_stats, t_star = t_star)
}

# Compatibility aliases for the naming used in earlier simulation folders.
run_tvipcw_nocloning <- run_fed_ipw_nocloning
run_perprotocol_noiptw <- run_fed_perprotocol

# -------------------------------------------------------------
# FEDERATED LANDMARK IPW
#   Restrict to patients event-free through tau, classify them by observed
#   initiation by tau, and estimate post-landmark survival. Each site fits a
#   local landmark propensity model and releases weighted event/risk totals.
#   This deliberately targets survival conditional on reaching tau rather
#   than the original time-zero grace-period estimand.
# -------------------------------------------------------------
local_landmark_counts <- function(site_data, tau) {
  eligible <- site_data$T_event > tau
  c(treated = sum(site_data$A_tau[eligible] == 1L),
    eligible = sum(eligible))
}

central_landmark_numerator <- function(local_counts) {
  totals <- Reduce(`+`, local_counts)
  if (totals[["eligible"]] <= 0)
    stop("No patients are event-free at the landmark.")
  .clamp_prob(totals[["treated"]] / totals[["eligible"]])
}

local_landmark_ipw <- function(site_data, tau, t_star, p_num,
                               trunc = DEFAULT_TRUNC) {
  d <- site_data[site_data$T_event > tau, , drop = FALSE]
  if (!nrow(d)) stop("A site has no patients event-free through tau.")

  L1 <- as.matrix(d[, paste0("L1_", seq_len(t_star)), drop = FALSE])
  L2 <- as.matrix(d[, paste0("L2_", seq_len(t_star)), drop = FALSE])
  fit_df <- data.frame(
    A = d$A_tau,
    x1 = d$x1,
    x2 = d$x2,
    L1 = L1[, tau],
    L2 = L2[, tau],
    L1lag = L1[, tau - 1L],
    L2lag = L2[, tau - 1L]
  )
  ps_fit <- glm(A ~ x1 + x2 + L1 + L2 + L1lag + L2lag,
                data = fit_df, family = binomial())
  p_den <- .clamp_prob(predict(ps_fit, type = "response"))
  w <- ifelse(d$A_tau == 1L, p_num / p_den,
              (1 - p_num) / (1 - p_den))
  q <- quantile(w[is.finite(w)], probs = trunc, names = FALSE)
  w <- pmin(pmax(w, q[1]), q[2])

  arm_stats <- function(arm) {
    use <- d$A_tau == arm
    Tev <- d$T_event[use]
    Tcap <- pmin(Tev, t_star)
    wa <- w[use]
    event <- risk <- numeric(t_star)
    for (m in seq_len(t_star)) {
      if (m <= tau) next
      tm1 <- m - 1
      Rm <- as.integer(Tcap > tm1)
      Dm <- as.integer(Rm == 1 & Tev > tm1 & Tev <= m)
      event[m] <- sum(wa * Dm)
      risk[m] <- sum(wa * Rm)
    }
    list(event = event, risk = risk)
  }

  list(arm1 = arm_stats(1L), arm0 = arm_stats(0L), n = nrow(d))
}

central_landmark_ipw <- function(site_stats, t_star) {
  sum_field <- function(arm, field) {
    Reduce(`+`, lapply(site_stats, function(x) x[[arm]][[field]]))
  }
  survival <- function(arm) {
    event <- sum_field(arm, "event")
    risk <- sum_field(arm, "risk")
    hazard <- ifelse(risk > 0, event / risk, 0)
    c(1, cumprod(1 - hazard))
  }
  .contrasts(survival("arm1"), survival("arm0"), t_star)
}

run_fed_landmark_ipw <- function(dat, tau, t_star,
                                 trunc = DEFAULT_TRUNC) {
  site_data <- split(dat, dat$site)
  p_num <- central_landmark_numerator(
    lapply(site_data, local_landmark_counts, tau = tau)
  )
  site_stats <- lapply(site_data, local_landmark_ipw,
                       tau = tau, t_star = t_star,
                       p_num = p_num, trunc = trunc)
  central_landmark_ipw(site_stats, t_star = t_star)
}

# -------------------------------------------------------------
# LOCAL CCW + CURVE META-ANALYSIS
#   Sites use the common numerator but otherwise complete CCW locally. They
#   release only their two survival curves and sample size. The center takes
#   a sample-size-weighted average of the curves, then constructs contrasts.
#   This is intentionally different from Fed-CCW, which aggregates weighted
#   event/risk totals before constructing the survival curves.
# -------------------------------------------------------------
local_ccw_curve_summary <- function(site_data, tau, t_star, hnum,
                                    trunc = DEFAULT_TRUNC) {
  stats <- local_ccw_tvipcw(site_data, tau = tau, t_star = t_star,
                            trunc = trunc, hnum = hnum)
  lambda1 <- ifelse(stats$rw1 > 0, stats$dw1 / stats$rw1, 0)
  lambda0 <- ifelse(stats$rw0 > 0, stats$dw0 / stats$rw0, 0)
  list(
    n = nrow(site_data),
    S1 = c(1, cumprod(1 - lambda1)),
    S0 = c(1, cumprod(1 - lambda0))
  )
}

central_ccw_curve_meta <- function(site_summaries, t_star) {
  n <- vapply(site_summaries, `[[`, numeric(1), "n")
  weighted_curve <- function(arm) {
    curves <- lapply(site_summaries, `[[`, arm)
    Reduce(`+`, Map(`*`, curves, n)) / sum(n)
  }
  .contrasts(weighted_curve("S1"), weighted_curve("S0"), t_star)
}

run_local_ccw_meta <- function(dat, tau, t_star,
                               trunc = DEFAULT_TRUNC) {
  site_data <- split(dat, dat$site)
  hnum <- central_common_num_hazard(
    lapply(site_data, local_initiation_counts, t_star = t_star)
  )
  site_summaries <- lapply(site_data, local_ccw_curve_summary,
                           tau = tau, t_star = t_star,
                           hnum = hnum, trunc = trunc)
  central_ccw_curve_meta(site_summaries, t_star = t_star)
}
