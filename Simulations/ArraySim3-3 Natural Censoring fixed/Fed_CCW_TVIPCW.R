# =============================================================
# Fed_CCW_TVIPCW.R  (natural-censoring variant)
#   Federated clone-censor-weight with time-varying IPCW, under a second,
#   stochastic (natural / loss-to-follow-up) censoring process.
#
#   Requires wide L columns L1_1..L1_M, L2_1..L2_M plus the baseline columns
#   (site, id, x1, x2, S, A_tau, tau_observed, T_event, T_cens, T_obs,
#   delta, cens).
#
# WEIGHTING MODEL
#   Artificial censoring in CCW is deterministic given (S, tau), so a clone's
#   artificial-censoring hazard in interval m is the hazard of the deviation
#   event = treatment initiation among the not-yet-initiated. One pooled
#   logistic initiation-hazard model is fit per site on person-interval rows
#   (baseline X, current L_m, lagged L_{m-1}); the stabilizing numerator uses
#   interval effects only.
#
#   Natural censoring is a SECOND, stochastic process. Because the two
#   mechanisms are conditionally independent given history, the total weight
#   factorizes:
#       SW^g_{i,m} = SW^{artificial,g}_{i,m} x SW^{natural,g}_{i,m}
#   The natural factor comes from one fitted censoring-hazard model per site,
#   evaluated SEPARATELY PER ARM (it conditions on the gated treatment
#   indicator, which differs between g=1 and g=0). With beta_cens["trt"] = 0
#   the two arm predictions coincide, but the code is general.
#
#   The artificial factor accrues only while m <= tau; the natural factor
#   accrues over all M intervals. The natural factor applicable DURING
#   interval m is accumulated through m-1 (standard IPCW timing). Weights are
#   stabilized and truncated on the PRODUCT at local percentiles.
# =============================================================

source("params.R")

# ---- per-site prep ------------------------------------------------------
.tv_prep <- function(d, tau, t_star) {
  M  <- t_star
  L1 <- as.matrix(d[, paste0("L1_", seq_len(M)), drop = FALSE])
  L2 <- as.matrix(d[, paste0("L2_", seq_len(M)), drop = FALSE])
  Tobs <- pmin(d$T_event, d$T_cens, t_star)
  tau_obs <- if (!is.null(d$tau_observed)) as.logical(d$tau_observed) else
    ((Tobs >= tau) | (d$S <= tau))

  list(n = nrow(d), M = M, L1 = L1, L2 = L2,
       x1 = d$x1, x2 = d$x2, S = d$S,
       Tev = d$T_event, Tcap = pmin(d$T_event, t_star),
       Tcens = d$T_cens, Tobs = Tobs, tau_obs = tau_obs)
}

# ---- lagged L matrix, with L_0 = 0 (matches DGP initialization) ---------
.lag_matrix <- function(Lmat) {
  cbind(0, Lmat[, -ncol(Lmat), drop = FALSE])
}

# ---- per-site initiation-hazard (artificial-censoring) model -----------
# Denominator conditions on X, current L_m, lagged L_{m-1}; numerator on the
# interval baseline only. Person-interval rows run to min(S, Tobs, M): a
# subject lost to follow-up contributes rows through Tobs (their L is
# observed) and nothing after.
.tv_init_hazard <- function(P) {
  n <- P$n; M <- P$M
  last_m <- pmin(P$S, P$Tobs, M)
  keep   <- outer(seq_len(n), seq_len(M), function(i, m) m <= last_m[i])

  idx     <- which(keep, arr.ind = TRUE)
  rows_id <- idx[, "row"]; rows_m <- idx[, "col"]
  L1lag <- .lag_matrix(P$L1); L2lag <- .lag_matrix(P$L2)

  pi_df <- data.frame(
    y     = as.integer(P$S[rows_id] == rows_m),
    m     = factor(rows_m, levels = seq_len(M)),
    x1    = P$x1[rows_id], x2 = P$x2[rows_id],
    L1    = P$L1[keep], L2 = P$L2[keep],
    L1lag = L1lag[keep], L2lag = L2lag[keep]
  )
  pi_df$m <- droplevels(pi_df$m)

  den_fit <- glm(y ~ m + x1 + x2 + L1 + L2 + L1lag + L2lag,
                 data = pi_df, family = binomial())
  num_fit <- glm(y ~ m, data = pi_df, family = binomial())

  Hden <- matrix(NA_real_, n, M); Hnum <- matrix(NA_real_, n, M)
  cell <- cbind(rows_id, rows_m)
  Hden[cell] <- .clamp_prob(predict(den_fit, type = "response"))
  Hnum[cell] <- .clamp_prob(predict(num_fit, type = "response"))
  list(Hden = Hden, Hnum = Hnum)
}

# ---- per-site natural-censoring hazard model ---------------------------
# Modelled event is loss-to-follow-up. Person-interval rows: every interval
# in which the subject was under observation AND event-free (m = 1..Tobs,
# excluding the interval where the event occurred, since no censoring draw
# was made there). Denominator conditions on X, current + lagged L, and the
# gated treatment indicator; numerator on interval baseline only.
#
# Returns three denominator matrices from the SAME fitted model, evaluated
# under different treatment schedules:
#   Cden1 : g = 1 (treated from min(S, tau) onward)
#   Cden0 : g = 0 (never treated)
#   Cnum  : marginal stabilizing numerator (arm-free)
.tv_cens_hazard <- function(P, tau) {
  n <- P$n; M <- P$M
  last_m  <- pmin(P$Tobs, M)
  in_grid <- outer(seq_len(n), seq_len(M), function(i, m) m <= last_m[i])
  # exclude the event interval (no censoring draw there)
  not_event <- outer(seq_len(n), seq_len(M), function(i, m) P$Tev[i] != m)
  keep <- in_grid & not_event

  idx     <- which(keep, arr.ind = TRUE)
  rows_id <- idx[, "row"]; rows_m <- idx[, "col"]
  L1lag <- .lag_matrix(P$L1); L2lag <- .lag_matrix(P$L2)
  trt_obs <- as.integer(P$S[rows_id] <= tau & rows_m >= P$S[rows_id])

  cd <- data.frame(
    y     = as.integer(P$Tcens[rows_id] == rows_m),
    m     = factor(rows_m, levels = seq_len(M)),
    x1    = P$x1[rows_id], x2 = P$x2[rows_id],
    L1    = P$L1[keep], L2 = P$L2[keep],
    L1lag = L1lag[keep], L2lag = L2lag[keep],
    trt   = trt_obs, id = rows_id
  )
  cd$m <- droplevels(cd$m)

  den_fit <- glm(y ~ m + x1 + x2 + L1 + L2 + L1lag + L2lag + trt,
                 data = cd, family = binomial())
  num_fit <- glm(y ~ m, data = cd, family = binomial())

  # arm-specific predictions from the same model: replace trt with the
  # strategy-consistent indicator, then re-predict.
  cd1 <- cd; cd1$trt <- as.integer(rows_m >= pmin(P$S[cd$id], tau))
  cd0 <- cd; cd0$trt <- 0L

  Cden1 <- matrix(NA_real_, n, M); Cden0 <- matrix(NA_real_, n, M)
  Cnum  <- matrix(NA_real_, n, M)
  cell  <- cbind(rows_id, rows_m)
  Cden1[cell] <- .clamp_prob(predict(den_fit, newdata = cd1, type = "response"))
  Cden0[cell] <- .clamp_prob(predict(den_fit, newdata = cd0, type = "response"))
  Cnum[cell]  <- .clamp_prob(predict(num_fit, type = "response"))
  list(Cden1 = Cden1, Cden0 = Cden0, Cnum = Cnum)
}

# ---- stabilized, truncated strategy weights ----------------------------
# SW1, SW0 are n x M; entry [i, m] is the weight DURING interval m, the
# product of an artificial factor (accrues while m <= tau) and a natural
# factor (accrues over all intervals, applied through m-1). Truncation is on
# the product.
.tv_weights <- function(P, H, CH, tau, trunc = DEFAULT_TRUNC) {
  n <- P$n; M <- P$M
  Hden <- H$Hden; Hnum <- H$Hnum
  Cden1 <- CH$Cden1; Cden0 <- CH$Cden0; Cnum <- CH$Cnum
  SW1 <- matrix(1, n, M); SW0 <- matrix(1, n, M)

  for (i in seq_len(n)) {
    aw1 <- 1; aw0 <- 1   # artificial factors
    nw1 <- 1; nw0 <- 1   # natural factors (through m-1 when applied)

    for (m in seq_len(M)) {
      # --- artificial component (m <= tau) ---
      hd <- Hden[i, m]; hn <- Hnum[i, m]
      if (m <= tau && !is.na(hd)) {
        if (P$S[i] == m) {
          aw1 <- aw1 * (hn / hd)
        } else {
          aw1 <- aw1 * ((1 - hn) / (1 - hd))
          aw0 <- aw0 * ((1 - hn) / (1 - hd))
        }
      }

      # --- total weight DURING interval m: natural factor through m-1 ---
      SW1[i, m] <- aw1 * nw1
      SW0[i, m] <- aw0 * nw0

      # --- update natural factor for interval m+1 onward ---
      cd1 <- Cden1[i, m]; cd0 <- Cden0[i, m]; cn <- Cnum[i, m]
      if (!is.na(cd1) && P$Tcens[i] > m) {
        nw1 <- nw1 * ((1 - cn) / (1 - cd1))
        nw0 <- nw0 * ((1 - cn) / (1 - cd0))
      }
    }
  }

  clip <- function(W) {
    q <- quantile(W[is.finite(W)], probs = trunc, names = FALSE)
    pmin(pmax(W, q[1]), q[2])
  }
  list(SW1 = clip(SW1), SW0 = clip(SW0))
}

# ---- clone artificial-censoring times ----------------------------------
# g=1: if still uninitiated when interval tau begins, deviate in interval
#      tau -> protocol-consistent through tau-1.
# g=0: deviate in interval S if initiated within grace -> censor at S-1. If
#      lost to follow-up before tau without confirmed non-initiation, censor
#      at Tobs (natural weights carry the correction).
.clone_censoring <- function(P, tau, t_star) {
  Cinf <- t_star + 1
  C1 <- ifelse(P$S >  tau, tau - 1L, Cinf)
  C0 <- ifelse(P$S <= tau, P$S - 1,
               ifelse(P$tau_obs, Cinf, P$Tobs))
  list(C1 = C1, C0 = C0)
}

# -------------------------------------------------------------
# LOCAL site computation
# -------------------------------------------------------------
local_ccw_tvipcw <- function(site_data, tau, t_star, trunc = DEFAULT_TRUNC) {
  P <- .tv_prep(site_data, tau, t_star)
  M <- P$M

  H  <- .tv_init_hazard(P)
  CH <- .tv_cens_hazard(P, tau)
  Wt <- .tv_weights(P, H, CH, tau, trunc = trunc)
  SW1 <- Wt$SW1; SW0 <- Wt$SW0

  Cs <- .clone_censoring(P, tau, t_star)
  C1 <- Cs$C1; C0 <- Cs$C0

  dw1 <- numeric(M); rw1 <- numeric(M)
  dw0 <- numeric(M); rw0 <- numeric(M)
  Tev <- P$Tev; Tcap <- P$Tcap; Tcens <- P$Tcens

  for (m in seq_len(M)) {
    tm1 <- m - 1
    w1 <- SW1[, m]; w0 <- SW0[, m]

    # risk sets additionally require being under observation (Tcens > m-1)
    R1m <- as.integer(Tcap > tm1 & C1 > tm1 & Tcens > tm1)
    D1m <- as.integer(R1m == 1 & Tev > tm1 & Tev <= m & Tev <= C1 & Tev <= Tcens)
    dw1[m] <- sum(w1 * D1m); rw1[m] <- sum(w1 * R1m)

    R0m <- as.integer(Tcap > tm1 & C0 > tm1 & Tcens > tm1)
    D0m <- as.integer(R0m == 1 & Tev > tm1 & Tev <= m & Tev <= C0 & Tev <= Tcens)
    dw0[m] <- sum(w0 * D0m); rw0[m] <- sum(w0 * R0m)
  }

  list(
    M = M,
    dw1 = dw1, rw1 = rw1, dw0 = dw0, rw0 = rw0,
    local = list(SW1 = SW1, SW0 = SW0, C1 = C1, C0 = C0,
                 Tev = Tev, Tcap = Tcap, Tcens = Tcens, n = P$n)
  )
}

# -------------------------------------------------------------
# LOCAL second pass: influence functions.
#   Weights (both nuisances) treated as known => conservative variance and
#   coverage that may run above nominal. Deriving the two-nuisance correction
#   is out of scope.
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

    R1m <- as.integer(L$Tcap > tm1 & L$C1 > tm1 & L$Tcens > tm1)
    D1m <- as.integer(R1m == 1 & L$Tev > tm1 & L$Tev <= m & L$Tev <= L$C1 &
                        L$Tev <= L$Tcens)
    if (RW1[m] > 0 && lambda1[m] < 1) {
      cumU1 <- cumU1 + (1 / (1 - lambda1[m])) * (w1 / RW1[m]) *
        (D1m - lambda1[m] * R1m)
    }

    R0m <- as.integer(L$Tcap > tm1 & L$C0 > tm1 & L$Tcens > tm1)
    D0m <- as.integer(R0m == 1 & L$Tev > tm1 & L$Tev <= m & L$Tev <= L$C0 &
                        L$Tev <= L$Tcens)
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
#   Same pipeline on one collapsed "site". A gap vs fed is the cost of
#   federating the weight models (per-site glm vs one pooled glm), not
#   causal bias.
# -------------------------------------------------------------
run_pooled_ccw_tvipcw <- function(dat, tau, t_star, trunc = DEFAULT_TRUNC) {
  pooled <- dat; pooled$site <- 1L
  stats <- list(local_ccw_tvipcw(pooled, tau = tau, t_star = t_star,
                                 trunc = trunc))
  central_ccw_tvipcw(stats, tau = tau, t_star = t_star)
}

# -------------------------------------------------------------
# NAIVE per-protocol comparator: crude as-treated survival by observed
# A_tau, no weighting, no cloning. Natural censoring handled only by removing
# censored subjects from the risk set (treated as independent). Expected to
# be biased -- shows what goes wrong when cloning AND censoring-weighting are
# both omitted.
# -------------------------------------------------------------
run_perprotocol_noiptw <- function(dat, tau, t_star) {
  M <- t_star
  crude_survival <- function(sub) {
    Tev <- sub$T_event; Tcap <- pmin(Tev, t_star); Tcens <- sub$T_cens
    lam <- numeric(M)
    for (m in seq_len(M)) {
      tm1 <- m - 1
      Rm <- as.integer(Tcap > tm1 & Tcens > tm1)
      Dm <- as.integer(Rm == 1 & Tev > tm1 & Tev <= m & Tev <= Tcens)
      r <- sum(Rm)
      lam[m] <- if (r > 0) sum(Dm) / r else 0
    }
    c(1, cumprod(1 - lam))
  }
  .contrasts(crude_survival(dat[dat$A_tau == 1, , drop = FALSE]),
             crude_survival(dat[dat$A_tau == 0, , drop = FALSE]), M)
}

# -------------------------------------------------------------
# NAIVE comparator: time-varying IPCW WITHOUT cloning.
#   Same per-interval stabilized, truncated TV weights (natural censoring IS
#   corrected), but no cloning and no artificial censoring. The observed
#   cohort is split by realized A_tau; each group's crude survival uses its
#   own strategy weight. Isolates the contribution of cloning; carries
#   immortal time bias the weights do not remove. Point estimates only.
# -------------------------------------------------------------
run_tvipcw_nocloning <- function(dat, tau, t_star, trunc = DEFAULT_TRUNC) {
  M <- t_star
  DW1 <- numeric(M); RW1 <- numeric(M)
  DW0 <- numeric(M); RW0 <- numeric(M)

  for (site_dat in split(dat, dat$site)) {
    P  <- .tv_prep(site_dat, tau, t_star)
    H  <- .tv_init_hazard(P)
    CH <- .tv_cens_hazard(P, tau)
    Wt <- .tv_weights(P, H, CH, tau, trunc = trunc)

    Tev <- P$Tev; Tcap <- P$Tcap; Tcens <- P$Tcens
    init <- P$S <= tau

    for (m in seq_len(M)) {
      tm1 <- m - 1
      Rm1 <- as.integer(Tcap[init] > tm1 & Tcens[init] > tm1)
      Dm1 <- as.integer(Rm1 == 1 & Tev[init] > tm1 & Tev[init] <= m &
                          Tev[init] <= Tcens[init])
      w1  <- Wt$SW1[init, m]
      DW1[m] <- DW1[m] + sum(w1 * Dm1); RW1[m] <- RW1[m] + sum(w1 * Rm1)

      Rm0 <- as.integer(Tcap[!init] > tm1 & Tcens[!init] > tm1)
      Dm0 <- as.integer(Rm0 == 1 & Tev[!init] > tm1 & Tev[!init] <= m &
                          Tev[!init] <= Tcens[!init])
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
#   covariate trajectory and averaging person-specific survival curves.
#
#     run_ccw_gcomp     -- each person is cloned into a g=1 and a g=0 copy.
#         Each copy is restricted to observed follow-up and artificially
#         censored at its OWN deviation time. The two arm-specific frames are
#         stacked before fitting the nuisance models.
#
#     run_gcomp_noclone -- each person contributes observed follow-up once,
#         with the observed gated treatment status and no artificial
#         censoring.
#
#   Natural censoring is handled by restricting both fitting frames to
#   observed follow-up. Under the correctly specified conditional outcome
#   and transition models, forward simulation then reconstructs follow-up
#   beyond loss to follow-up. Unlike the weighting methods, g-computation
#   does not fit a separate natural-censoring model.
#
#   Both variants standardize over the same full baseline cohort and use the
#   same strategy schedules. They differ only in the fitting data. Neither
#   g-computation implementation is federated.

# ---- plain person-interval frame (observed follow-up) ------------------
.gcomp_long_plain <- function(dat, t_star, tau) {
  M  <- t_star
  L1 <- as.matrix(dat[, paste0("L1_", seq_len(M)), drop = FALSE])
  L2 <- as.matrix(dat[, paste0("L2_", seq_len(M)), drop = FALSE])
  n  <- nrow(dat)

  Tobs <- pmin(dat$T_event, dat$T_cens, M)
  keep <- outer(seq_len(n), seq_len(M), function(i, m) m <= Tobs[i])
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
    y     = as.integer(dat$T_event[i_id] == m_id &
                       dat$T_event[i_id] <= dat$T_cens[i_id])
  )
}

# ---- cloned/censored person-interval frame -----------------------------
# Build one strategy arm at a time. Rows are restricted by BOTH observed
# follow-up and that arm's artificial-censoring time. This is the key fix
# relative to the earlier natural-censoring implementation, which did not
# create arm-specific clones.
.gcomp_clone_arm <- function(dat, t_star, tau, arm) {
  M  <- t_star
  L1 <- as.matrix(dat[, paste0("L1_", seq_len(M)), drop = FALSE])
  L2 <- as.matrix(dat[, paste0("L2_", seq_len(M)), drop = FALSE])
  n  <- nrow(dat)

  P    <- .tv_prep(dat, tau, t_star)
  Cs   <- .clone_censoring(P, tau, t_star)
  Cg   <- if (arm == 1L) Cs$C1 else Cs$C0
  Tobs <- pmin(dat$T_event, dat$T_cens, M)
  last <- pmin(Tobs, Cg)

  keep <- outer(seq_len(n), seq_len(M), function(i, m) m <= last[i])
  idx  <- which(keep, arr.ind = TRUE)
  i_id <- idx[, "row"]; m_id <- idx[, "col"]

  if (arm == 1L) {
    Son <- pmin(dat$S, tau)
    trt <- as.integer(m_id >= Son[i_id])
  } else {
    trt <- 0L
  }

  y <- as.integer(dat$T_event[i_id] == m_id &
                  dat$T_event[i_id] <= dat$T_cens[i_id] &
                  dat$T_event[i_id] <= Cg[i_id])

  data.frame(
    id    = i_id,
    arm   = arm,
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

.gcomp_long_cloned <- function(dat, t_star, tau) {
  rbind(.gcomp_clone_arm(dat, t_star, tau, arm = 1L),
        .gcomp_clone_arm(dat, t_star, tau, arm = 0L))
}

# ---- nuisance models on a supplied person-interval frame --------------
.gcomp_fit_long <- function(long, t_star) {
  M <- t_star

  out_fit <- glm(y ~ factor(m) + x1 + x2 + L1 + L2 + L1lag + L2lag + trt,
                 data = long, family = binomial())

  tr <- long[long$m >= 2, ]
  L1_fit <- lm(L1 ~ L1lag + x1 + x2, data = tr)
  L2_fit <- lm(L2 ~ L2lag + x1 + x2, data = tr)

  list(out_fit = out_fit,
       L1_fit = L1_fit, L2_fit = L2_fit,
       sd1 = sd(residuals(L1_fit)), sd2 = sd(residuals(L2_fit)),
       M = M)
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

# ---- shared standardization step ---------------------------------------
.gcomp_standardize <- function(dat, fit, tau, t_star, mc = 20, seed = 12345) {
  M <- t_star

  S1_arm <- ifelse(dat$S <= tau, dat$S, tau)
  S0_arm <- rep(t_star + 1L, nrow(dat))

  S1 <- .gcomp_strategy_survival(dat, fit, S1_arm, t_star,
                                 mc = mc, seed = seed)
  S0 <- .gcomp_strategy_survival(dat, fit, S0_arm, t_star,
                                 mc = mc, seed = seed + 1)

  c(.contrasts(S1, S0, M), list(fit = fit))
}

# ---- clone-censor g-computation ----------------------------------------
run_ccw_gcomp <- function(dat, tau, t_star, mc = 20, seed = 12345) {
  long <- .gcomp_long_cloned(dat, t_star, tau)
  fit  <- .gcomp_fit_long(long, t_star)
  .gcomp_standardize(dat, fit, tau, t_star, mc = mc, seed = seed)
}

# ---- ordinary g-computation without cloning ----------------------------
run_gcomp_noclone <- function(dat, tau, t_star, mc = 20, seed = 12345) {
  long <- .gcomp_long_plain(dat, t_star, tau)
  fit  <- .gcomp_fit_long(long, t_star)
  .gcomp_standardize(dat, fit, tau, t_star, mc = mc, seed = seed)
}
