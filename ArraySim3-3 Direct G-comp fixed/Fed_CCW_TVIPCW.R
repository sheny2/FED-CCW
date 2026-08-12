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
#   numerator uses interval effects only.
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

# ---- per-site initiation-hazard models ---------------------------------
# Denominator conditions on X, current L_m and lagged L_{m-1}, matching the
# outcome hazard's dependence on current + lagged L. Numerator is the
# interval baseline only.
.tv_init_hazard <- function(P) {
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
  num_fit <- glm(y ~ m, data = pi_df, family = binomial())

  Hden <- matrix(NA_real_, n, M)
  Hnum <- matrix(NA_real_, n, M)
  cell <- cbind(rows_id, rows_m)
  Hden[cell] <- .clamp_prob(predict(den_fit, type = "response"))
  Hnum[cell] <- .clamp_prob(predict(num_fit, type = "response"))

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
local_ccw_tvipcw <- function(site_data, tau, t_star, trunc = DEFAULT_TRUNC) {
  P <- .tv_prep(site_data, t_star)
  M <- P$M

  Wt  <- .tv_weights(P, .tv_init_hazard(P), tau, trunc = trunc)
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
  stats <- lapply(split(dat, dat$site), local_ccw_tvipcw,
                  tau = tau, t_star = t_star, trunc = trunc)
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
  pooled <- dat
  pooled$site <- 1L
  stats <- list(local_ccw_tvipcw(pooled, tau = tau, t_star = t_star,
                                 trunc = trunc))
  central_ccw_tvipcw(stats, tau = tau, t_star = t_star)
}

# -------------------------------------------------------------
# NAIVE per-protocol comparator: crude as-treated survival by observed
# A_tau, no weighting and no cloning. Expected to be biased here.
# -------------------------------------------------------------
run_perprotocol_noiptw <- function(dat, tau, t_star) {
  M <- t_star
  crude_survival <- function(sub) {
    Tev <- sub$T_event; Tcap <- pmin(Tev, M); lam <- numeric(M)
    for (m in seq_len(M)) {
      tm1 <- m - 1
      Rm <- as.integer(Tcap > tm1)
      Dm <- as.integer(Rm == 1 & Tev > tm1 & Tev <= m)
      r <- sum(Rm)
      lam[m] <- if (r > 0) sum(Dm) / r else 0
    }
    c(1, cumprod(1 - lam))
  }
  .contrasts(crude_survival(dat[dat$A_tau == 1, , drop = FALSE]),
             crude_survival(dat[dat$A_tau == 0, , drop = FALSE]), M)
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
run_tvipcw_nocloning <- function(dat, tau, t_star, trunc = DEFAULT_TRUNC) {
  M <- t_star

  DW1 <- numeric(M); RW1 <- numeric(M)
  DW0 <- numeric(M); RW0 <- numeric(M)

  # weights are built per site, matching how fed-CCW fits them, then the
  # per-interval weighted sums are pooled across sites
  for (site_dat in split(dat, dat$site)) {
    P  <- .tv_prep(site_dat, t_star)
    Wt <- .tv_weights(P, .tv_init_hazard(P), tau, trunc = trunc)

    Tev <- P$Tev; Tcap <- P$Tcap
    init <- P$S <= tau

    for (m in seq_len(M)) {
      tm1 <- m - 1

      Rm1 <- as.integer(Tcap[init] > tm1)
      Dm1 <- as.integer(Rm1 == 1 & Tev[init] > tm1 & Tev[init] <= m)
      w1  <- Wt$SW1[init, m]
      DW1[m] <- DW1[m] + sum(w1 * Dm1); RW1[m] <- RW1[m] + sum(w1 * Rm1)

      Rm0 <- as.integer(Tcap[!init] > tm1)
      Dm0 <- as.integer(Rm0 == 1 & Tev[!init] > tm1 & Tev[!init] <= m)
      w0  <- Wt$SW0[!init, m]
      DW0[m] <- DW0[m] + sum(w0 * Dm0); RW0[m] <- RW0[m] + sum(w0 * Rm0)
    }
  }

  lam1 <- ifelse(RW1 > 0, DW1 / RW1, 0)
  lam0 <- ifelse(RW0 > 0, DW0 / RW0, 0)
  .contrasts(c(1, cumprod(1 - lam1)), c(1, cumprod(1 - lam0)), M)
}


# =============================================================
# G-COMPUTATION (two variants)
# =============================================================
#   Both variants fit a pooled discrete-time outcome hazard and pooled L
#   transition models, then standardize by reconstructing each strategy's
#   covariate trajectory (Monte Carlo) and averaging person-specific survival
#   curves. They differ ONLY in the person-time the nuisance models are fit
#   on:
#
#     run_ccw_gcomp     -- CLONE-CENSOR g-computation. Each person is cloned
#         into a g=1 and a g=0 copy; each copy is artificially censored at its
#         OWN deviation time (C1 for g=1, C0 for g=0); the outcome and
#         transition models are fit on the STACKED cloned/censored rows. Thus
#         deviation-free early follow-up is contributed by BOTH clones, while
#         post-deviation person-time is dropped per arm. This is the paper's
#         Algorithm-2 input (fit on the cloned/censored dataset).
#
#     run_gcomp_noclone -- PLAIN g-computation. No cloning, no artificial
#         censoring: the models are fit ONCE on each person's observed
#         follow-up with the observed gated treatment status.
#
#   Both then standardize over the SAME strategy schedules via
#   .gcomp_standardize, so any difference between them is attributable purely
#   to the fitting data (cloned/censored vs observed).
#
#   Strategy schedules, matching the DGP's gated initiation:
#     g = 1 "initiate by tau": on from the natural S if S <= tau, otherwise
#           the deviator is reconstructed as initiating at tau.
#     g = 0 "never initiate by tau": never on.
#
#   NOT FEDERATED: transition + outcome models are fit on pooled data.
#   The L transition model has no treatment term, matching the DGP in which
#   treatment does not shift L.
# =============================================================

# ---- plain person-interval frame (observed follow-up) ------------------
# One row per (person, interval) for m = 1..Tcap_i, with the OBSERVED gated
# treatment status. This is the plain-g-computation input.
.gcomp_long_plain <- function(dat, t_star, tau) {
  M  <- t_star
  L1 <- as.matrix(dat[, paste0("L1_", seq_len(M)), drop = FALSE])
  L2 <- as.matrix(dat[, paste0("L2_", seq_len(M)), drop = FALSE])
  n  <- nrow(dat)

  Tcap <- pmin(dat$T_event, M)
  keep <- outer(seq_len(n), seq_len(M), function(i, m) m <= Tcap[i])
  idx  <- which(keep, arr.ind = TRUE)
  i_id <- idx[, "row"]; m_id <- idx[, "col"]

  data.frame(
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
}

# ---- cloned/censored person-interval frame -----------------------------
# Build ONE arm's clones: keep rows m = 1..min(Tcap_i, C_i) where C_i is that
# arm's clone-censoring time, and stamp `trt` with the strategy-consistent
# indicator. Within a clone's retained (deviation-free) rows the strategy-
# consistent indicator equals the observed gated status, because the clone is
# censored exactly when its observed treatment would first contradict its
# assignment -- but we set it explicitly so the frame is correct by
# construction rather than by that coincidence.
#
# The event indicator y is zeroed in the deviation interval and beyond by the
# censoring (those rows are dropped): an artificially censored clone is a
# non-event at its censoring time, exactly as in the IPCW risk sets.
.gcomp_clone_arm <- function(dat, t_star, tau, arm) {
  M  <- t_star
  L1 <- as.matrix(dat[, paste0("L1_", seq_len(M)), drop = FALSE])
  L2 <- as.matrix(dat[, paste0("L2_", seq_len(M)), drop = FALSE])
  n  <- nrow(dat)

  Tcap <- pmin(dat$T_event, M)
  Cs   <- .clone_censoring(dat$S, tau, t_star)
  Cg   <- if (arm == 1L) Cs$C1 else Cs$C0
  last <- pmin(Tcap, Cg)                 # rows m = 1..min(Tcap, C_arm)

  keep <- outer(seq_len(n), seq_len(M), function(i, m) m <= last[i])
  idx  <- which(keep, arr.ind = TRUE)
  i_id <- idx[, "row"]; m_id <- idx[, "col"]

  # strategy-consistent treatment schedule for this arm
  if (arm == 1L) {
    # g = 1: initiate by tau -> on from min(S, tau) onward
    Son <- pmin(dat$S, tau)
    trt <- as.integer(m_id >= Son[i_id])
  } else {
    # g = 0: never initiate by tau -> never on within retained (pre-tau) rows
    trt <- 0L
  }

  # an event counts only if it occurs within the retained (uncensored) rows,
  # i.e. at or before the clone's censoring time; rows are already truncated
  # at min(Tcap, C_arm), so a retained event row has T_event == m <= C_arm.
  y <- as.integer(dat$T_event[i_id] == m_id & dat$T_event[i_id] <= Cg[i_id])

  data.frame(
    id    = i_id,
    m     = m_id,
    x1    = dat$x1[i_id],
    x2    = dat$x2[i_id],
    L1    = L1[keep],
    L2    = L2[keep],
    L1lag = .lag_matrix(L1)[keep],
    L2lag = .lag_matrix(L2)[keep],
    trt   = trt,
    y     = y
  )
}

# Stack the two arms' clones into the cloned/censored fitting dataset.
.gcomp_long_cloned <- function(dat, t_star, tau) {
  rbind(.gcomp_clone_arm(dat, t_star, tau, arm = 1L),
        .gcomp_clone_arm(dat, t_star, tau, arm = 0L))
}

# ---- nuisance models on a given person-interval frame ------------------
.gcomp_fit_long <- function(long, t_star) {
  M <- t_star

  # outcome hazard: event in m given X, current + lagged L, gated treatment
  out_fit <- glm(y ~ factor(m) + x1 + x2 + L1 + L2 + L1lag + L2lag + trt,
                 data = long, family = binomial())

  # L transition models from consecutive-interval pairs (linear-Gaussian,
  # matching the DGP). Intervals m >= 2 only, since L_1 has no observed lag.
  tr <- long[long$m >= 2, ]
  L1_fit <- lm(L1 ~ L1lag + x1 + x2, data = tr)
  L2_fit <- lm(L2 ~ L2lag + x1 + x2, data = tr)

  list(out_fit = out_fit,
       L1_fit = L1_fit, L2_fit = L2_fit,
       sd1 = sd(residuals(L1_fit)), sd2 = sd(residuals(L2_fit)),
       M = M)
}

# ---- reconstruct trajectories under a strategy and predict survival ----
# S_arm[i] is the interval treatment turns on for person i under the
# strategy (t_star + 1 means never). mc = Monte Carlo covariate draws.
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
    predict(b_out, newdata = nd, type = "response")
  }

  Ssum <- numeric(M + 1)
  Ssum[1] <- n * mc   # S(t_0) = 1 for every (person, draw)

  for (i in seq_len(n)) {
    x1 <- dat$x1[i]; x2 <- dat$x2[i]
    Son <- S_arm[i]

    surv_i <- matrix(1, mc, M + 1)
    L1prev <- rep(0, mc); L2prev <- rep(0, mc)
    runS   <- rep(1, mc)

    for (m in seq_len(M)) {
      trt <- as.integer(m >= Son)   # gated, strategy-consistent

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

# ---- shared standardization step ---------------------------------------
# Given fitted nuisance models, reconstruct both strategy arms and form the
# contrasts. The strategy schedules are identical for both g-comp variants;
# only the data the nuisance models were fit on differs. Standardization is
# over the FULL observed cohort (`dat`), not the cloned frame -- cloning
# affects only how the nuisance models are estimated.
.gcomp_standardize <- function(dat, fit, tau, t_star, mc = 20, seed = 12345) {
  M <- t_star

  # g = 1: natural S when S <= tau, deviators reconstructed as initiating at tau
  S1_arm <- ifelse(dat$S <= tau, dat$S, tau)
  # g = 0: never on
  S0_arm <- rep(t_star + 1L, nrow(dat))

  S1 <- .gcomp_strategy_survival(dat, fit, S1_arm, t_star, mc = mc, seed = seed)
  S0 <- .gcomp_strategy_survival(dat, fit, S0_arm, t_star, mc = mc, seed = seed + 1)

  c(.contrasts(S1, S0, M), list(fit = fit))
}

# ---- top level: CLONE-CENSOR g-computation -----------------------------
# Nuisance models fit on the STACKED cloned/censored dataset.
run_ccw_gcomp <- function(dat, tau, t_star, mc = 20, seed = 12345) {
  long <- .gcomp_long_cloned(dat, t_star, tau)
  fit  <- .gcomp_fit_long(long, t_star)
  .gcomp_standardize(dat, fit, tau, t_star, mc = mc, seed = seed)
}

# -------------------------------------------------------------
# COMPARATOR: PLAIN G-COMPUTATION (no cloning, no censoring)
#   Nuisance models fit ONCE on each person's observed follow-up with the
#   observed gated treatment status. Nothing is cloned and no artificial
#   censoring is applied; the standardization step over the two strategy
#   schedules is identical to run_ccw_gcomp.
#
#   PURPOSE: isolate the contribution of CLONE-CENSORING to the g-formula
#   arm, exactly as run_tvipcw_nocloning does for the weighting arm. Any gap
#   between this and run_ccw_gcomp is attributable to fitting the hazard and
#   transition models on cloned/censored vs observed person-time.
# -------------------------------------------------------------
run_gcomp_noclone <- function(dat, tau, t_star, mc = 20, seed = 12345) {
  long <- .gcomp_long_plain(dat, t_star, tau)
  fit  <- .gcomp_fit_long(long, t_star)
  .gcomp_standardize(dat, fit, tau, t_star, mc = mc, seed = seed)
}