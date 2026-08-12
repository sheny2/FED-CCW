# =============================================================
# Fed_CCW_TVIPCW.R  (trimmed two-method version)
#
# Estimators retained:
#   run_fed_ccw_tvipcw  -- federated clone-censor-weight with time-varying
#                          IPCW. Weight (initiation-hazard) model is subject
#                          to the shared misspecification axis.
#   run_ccw_gcomp       -- clone-censor + parametric g-computation. Outcome
#                          hazard model is subject to the same axis.
#
# The misspecification axis (correct / no_tv / coarse_L; see params.R)
# degrades how the time-varying covariates L enter each method's nuisance
# model. In the SAME cell it is applied to fed-CCW's weight denominator and
# to g-comp's outcome hazard, so the two methods are compared under a common
# nuisance error.
#
# Removed relative to the full study: pooled CCW, plain (no-clone) g-comp,
# IPTW-without-cloning, and the naive per-protocol comparator.
#
# Requires wide L columns L1_1..L1_M, L2_1..L2_M alongside the baseline
# columns (site, id, x1, x2, S, A_tau, T_event, T_obs, delta).
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

# ---- coarsened L, shared by both nuisance models ------------------------
# Binary indicator of being above a FIXED cut point. Cut points are the
# marginal medians of the fitted data, carried on the frame as an attribute
# so any prediction step reuses them rather than recomputing (which would
# make the coarsening depend on the data being predicted).
.add_coarse_L <- function(long, cuts = NULL) {
  if (is.null(cuts)) {
    cuts <- c(L1 = median(long$L1), L2 = median(long$L2))
  }
  long$L1_hi    <- as.integer(long$L1    > cuts[["L1"]])
  long$L2_hi    <- as.integer(long$L2    > cuts[["L2"]])
  long$L1lag_hi <- as.integer(long$L1lag > cuts[["L1"]])
  long$L2lag_hi <- as.integer(long$L2lag > cuts[["L2"]])
  attr(long, "L_cuts") <- cuts
  long
}


# =============================================================
# FEDERATED CCW with time-varying IPCW
# =============================================================

# ---- per-site initiation-hazard (weight) model -------------------------
# Denominator conditions on X and the time-varying covariates as dictated by
# `misspec`; numerator is the interval baseline only. The stabilized weight
# is prod_m (1 - Hnum)/(1 - Hden) while uninitiated, with the initiation
# ratio Hnum/Hden in the interval of initiation.
#
# misspec degrades the tv block of the DENOMINATOR only (the numerator has no
# L terms in any spec), which is what removes / coarsens the confounding
# adjustment carried by the weights.
.tv_init_hazard <- function(P, misspec = DEFAULT_MISSPEC) {
  n <- P$n; M <- P$M

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
    L2lag = L2lag[keep]
  )
  pi_df$m <- droplevels(pi_df$m)
  pi_df   <- .add_coarse_L(pi_df)

  den_rhs <- .build_rhs(misspec, base = "m + x1 + x2")
  den_fit <- glm(stats::as.formula(paste("y ~", den_rhs)),
                 data = pi_df, family = binomial())
  num_fit <- glm(y ~ m, data = pi_df, family = binomial())

  Hden <- matrix(NA_real_, n, M)
  Hnum <- matrix(NA_real_, n, M)
  cell <- cbind(rows_id, rows_m)
  Hden[cell] <- .clamp_prob(predict(den_fit, type = "response"))
  Hnum[cell] <- .clamp_prob(predict(num_fit, type = "response"))

  list(Hden = Hden, Hnum = Hnum)
}

# ---- stabilized, truncated strategy weights ----------------------------
# Returns n x M matrices SW1, SW0; entry [i, m] is the weight during m.
.tv_weights <- function(P, H, tau, trunc = DEFAULT_TRUNC) {
  n <- P$n; M <- P$M
  Hden <- H$Hden; Hnum <- H$Hnum

  ratio_init <- matrix(1, n, M)   # interval where S == m
  ratio_stay <- matrix(1, n, M)   # while still uninitiated

  active <- !is.na(Hden)
  active[, seq_len(M) > tau] <- FALSE

  ratio_init[active] <- Hnum[active] / Hden[active]
  ratio_stay[active] <- (1 - Hnum[active]) / (1 - Hden[active])

  init_at_m <- outer(P$S, seq_len(M), `==`)

  F1 <- ifelse(init_at_m, ratio_init, ratio_stay)
  F1[!active] <- 1

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
.clone_censoring <- function(S, tau, t_star) {
  Cinf <- t_star + 1L
  list(C1 = ifelse(S >  tau, tau - 1L, Cinf),
       C0 = ifelse(S <= tau, S   - 1L, Cinf))
}

# ---- LOCAL site computation --------------------------------------------
local_ccw_tvipcw <- function(site_data, tau, t_star, trunc = DEFAULT_TRUNC,
                             misspec = DEFAULT_MISSPEC) {
  P <- .tv_prep(site_data, t_star)
  M <- P$M

  Wt  <- .tv_weights(P, .tv_init_hazard(P, misspec = misspec), tau, trunc = trunc)
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

# ---- LOCAL second pass: influence functions ----------------------------
# Weights are treated as known, so the variance is conservative.
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

# ---- CENTRAL aggregation + estimation ----------------------------------
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

run_fed_ccw_tvipcw <- function(dat, tau, t_star, trunc = DEFAULT_TRUNC,
                               misspec = DEFAULT_MISSPEC) {
  stats <- lapply(split(dat, dat$site), local_ccw_tvipcw,
                  tau = tau, t_star = t_star, trunc = trunc, misspec = misspec)
  central_ccw_tvipcw(stats, tau = tau, t_star = t_star)
}


# =============================================================
# CLONE-CENSOR + G-COMPUTATION
# =============================================================
#   (1) fit pooled transition models for L1, L2 (L_m given L_{m-1}, X);
#   (2) fit a pooled discrete-time outcome hazard on the cloned/censored
#       person-interval data, with the tv block set by `misspec`;
#   (3) reconstruct each strategy's covariate trajectory by Monte Carlo and
#       plug into the outcome model for a person-specific survival curve;
#   (4) average and form contrasts.
#
#   NOT FEDERATED: nuisance models are fit on pooled person-interval data.
#   The L transition models are always correctly specified, so any change is
#   attributable to the outcome hazard model.

# ---- pooled person-interval frame (clone-censored) ---------------------
# Person-time is restricted to deviation-free follow-up: each person
# contributes rows up to pmin(Tcap, max(C1, C0)), the last interval in which
# they remain protocol-consistent under at least one strategy.
.gcomp_long <- function(dat, tau, t_star) {
  M <- t_star
  L1 <- as.matrix(dat[, paste0("L1_", seq_len(M)), drop = FALSE])
  L2 <- as.matrix(dat[, paste0("L2_", seq_len(M)), drop = FALSE])
  n  <- nrow(dat)

  Tcap <- pmin(dat$T_event, M)
  Cs   <- .clone_censoring(dat$S, tau, t_star)
  last <- pmin(Tcap, pmax(Cs$C1, Cs$C0))

  keep <- outer(seq_len(n), seq_len(M), function(i, m) m <= last[i])
  idx  <- which(keep, arr.ind = TRUE)
  i_id <- idx[, "row"]; m_id <- idx[, "col"]

  long <- data.frame(
    id    = i_id,
    m     = m_id,
    x1    = dat$x1[i_id],
    x2    = dat$x2[i_id],
    L1    = L1[keep],
    L2    = L2[keep],
    L1lag = .lag_matrix(L1)[keep],
    L2lag = .lag_matrix(L2)[keep],
    trt   = as.integer(dat$S[i_id] <= tau & m_id >= dat$S[i_id]),
    y     = as.integer(dat$T_event[i_id] == m_id)
  )
  .add_coarse_L(long)
}

# ---- nuisance models on pooled data ------------------------------------
.gcomp_fit <- function(dat, tau, t_star, misspec = DEFAULT_MISSPEC) {
  M    <- t_star
  long <- .gcomp_long(dat, tau, t_star)

  # outcome hazard: event in m given X, the tv block dictated by `misspec`,
  # and gated treatment
  out_rhs <- .build_rhs(misspec, base = "factor(m) + x1 + x2", tail = "trt")
  out_fit <- glm(stats::as.formula(paste("y ~", out_rhs)),
                 data = long, family = binomial())

  # L transition models (linear-Gaussian, matching the DGP); always correct.
  tr <- long[long$m >= 2, ]
  L1_fit <- lm(L1 ~ L1lag + x1 + x2, data = tr)
  L2_fit <- lm(L2 ~ L2lag + x1 + x2, data = tr)

  list(out_fit = out_fit,
       L1_fit = L1_fit, L2_fit = L2_fit,
       sd1 = sd(residuals(L1_fit)), sd2 = sd(residuals(L2_fit)),
       M = M, misspec = misspec,
       L_cuts = attr(long, "L_cuts"))
}

# ---- reconstruct trajectories under a strategy and predict survival ----
.gcomp_strategy_survival <- function(dat, fit, S_arm, t_star, mc = 20,
                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  M <- t_star; n <- nrow(dat)
  b_out <- fit$out_fit
  b1 <- coef(fit$L1_fit); b2 <- coef(fit$L2_fit)

  haz <- function(m, x1, x2, L1, L2, L1lag, L2lag, trt) {
    nd <- data.frame(m = m, x1 = x1, x2 = x2,
                     L1 = L1, L2 = L2, L1lag = L1lag, L2lag = L2lag,
                     trt = trt)
    # coarse_L needs the binarized columns; reuse the fitted cut points so
    # the coarsening does not shift by arm
    if (!is.null(fit$L_cuts)) nd <- .add_coarse_L(nd, cuts = fit$L_cuts)
    .clamp_prob(predict(b_out, newdata = nd, type = "response"))
  }

  Ssum <- numeric(M + 1)
  Ssum[1] <- n * mc

  for (i in seq_len(n)) {
    x1 <- dat$x1[i]; x2 <- dat$x2[i]
    Son <- S_arm[i]

    surv_i <- matrix(1, mc, M + 1)
    L1prev <- rep(0, mc); L2prev <- rep(0, mc)
    runS   <- rep(1, mc)

    for (m in seq_len(M)) {
      trt <- as.integer(m >= Son)

      mu1 <- b1["(Intercept)"] + b1["L1lag"] * L1prev +
        b1["x1"] * x1 + b1["x2"] * x2
      mu2 <- b2["(Intercept)"] + b2["L2lag"] * L2prev +
        b2["x1"] * x1 + b2["x2"] * x2
      L1m <- rnorm(mc, mu1, fit$sd1)
      L2m <- rnorm(mc, mu2, fit$sd2)

      runS <- runS * (1 - haz(m, x1, x2, L1m, L2m, L1prev, L2prev, trt))
      surv_i[, m + 1] <- runS

      L1prev <- L1m; L2prev <- L2m
    }
    Ssum <- Ssum + colSums(surv_i)
  }

  Ssum / (n * mc)
}

# ---- top level: clone-censor + g-computation ---------------------------
run_ccw_gcomp <- function(dat, tau, t_star, mc = 20, seed = 12345,
                          misspec = DEFAULT_MISSPEC) {
  M   <- t_star
  fit <- .gcomp_fit(dat, tau, t_star, misspec = misspec)

  S1_arm <- ifelse(dat$S <= tau, dat$S, tau)   # g = 1
  S0_arm <- rep(t_star + 1L, nrow(dat))        # g = 0

  S1 <- .gcomp_strategy_survival(dat, fit, S1_arm, t_star, mc = mc, seed = seed)
  S0 <- .gcomp_strategy_survival(dat, fit, S0_arm, t_star, mc = mc, seed = seed + 1)

  c(.contrasts(S1, S0, M), list(fit = fit))
}