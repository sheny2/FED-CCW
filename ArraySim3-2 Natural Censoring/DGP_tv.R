# =============================================================
# Discrete-time DGP with time-varying confounding + NATURAL CENSORING
# (gated-treatment variant, time-varying L).
#
# Interval m = 1..t_star.
#   x1, x2            baseline covariates (continuous, binary)
#   L1_m, L2_m        time-varying confounders, AR(1) in own lag + baseline X;
#                     NOT affected by prior treatment. They drive initiation,
#                     the event hazard (current + one lag), AND the natural
#                     censoring hazard.
#   treatment         absorbing; initiated at S, gated so the treated hazard
#                     applies iff A_tau == 1 and m >= S
#   natural censoring stochastic loss-to-follow-up; hazard depends on X,
#                     current + lagged L, and (optionally) gated treatment.
#                     beta_cens["trt"] = 0 => censoring independent of
#                     treatment (tier 2, the default).
#
# Order within interval m:
#   (a) draw L1_m, L2_m
#   (b) if not yet initiated, draw initiation using X and L_m
#   (c) draw event using X, L_m, L_{m-1}, gated treatment
#   (d) if NO event, draw natural censoring using X, L_m, L_{m-1}, gated trt
#
# Censoring is drawn LAST, so event and censor in the same interval resolve
# in favour of the event (delta is unambiguous). A subject censored in m HAS
# observed L_m and contributes a person-interval row for m. Wide L columns
# are NA for m beyond observed follow-up, so any code assuming a complete L
# history fails loudly rather than using unobserved covariates.
#
# Observed-data columns:
#   T_event      true event interval (t_star+1 if none)
#   T_cens       natural censoring interval (t_star+1 if never)
#   T_obs        min(T_event, T_cens, t_star)
#   delta        1 iff event observed (T_event <= min(T_cens, t_star))
#   cens         1 iff natural censoring observed before event & horizon
#   S            initiation interval (t_star+1 if never); observable, since
#                initiation is only drawn among those under observation
#   A_tau        1 iff S <= tau
#   tau_observed 1 iff the grace period was fully resolved (initiated by tau,
#                or followed at least to tau). If 0, A_tau == 0 is NOT a
#                confirmed non-initiation and the g=0 clone must be
#                artificially censored at T_obs.
# =============================================================

source("params.R")

# L transition mean for one interval (shared by DGP and oracle).
.mu_L <- function(beta_L, L_prev, x1, x2) {
  beta_L["int"] + beta_L["ar"] * L_prev + beta_L["x1"] * x1 + beta_L["x2"] * x2
}

simulate_site_tv <- function(n,
                             tau        = DEFAULT_TAU,
                             t_star     = DEFAULT_TSTAR,
                             beta_event = DEFAULT_BETA_EVENT,
                             beta_init  = DEFAULT_BETA_INIT,
                             beta_L     = DEFAULT_BETA_L,
                             beta_cens  = DEFAULT_BETA_CENS,
                             sd_L       = DEFAULT_SD_L,
                             site_id    = 1,
                             seed       = NULL) {

  if (!is.null(seed)) set.seed(seed)

  x1 <- rnorm(n, 0, 1)
  x2 <- rbinom(n, 1, 0.4)

  L1 <- matrix(NA_real_, n, t_star)
  L2 <- matrix(NA_real_, n, t_star)
  S         <- rep(t_star + 1L, n)
  T_event   <- rep(t_star + 1L, n)
  T_cens    <- rep(t_star + 1L, n)
  initiated <- rep(FALSE, n)
  alive     <- rep(TRUE,  n)   # event-free entering m
  uncens    <- rep(TRUE,  n)   # under observation entering m

  L1_prev <- rep(0, n)
  L2_prev <- rep(0, n)

  for (m in seq_len(t_star)) {

    at_risk <- alive & uncens
    if (!any(at_risk)) break

    # (a) time-varying covariates (only for those under observation)
    L1[at_risk, m] <- rnorm(sum(at_risk),
                            .mu_L(beta_L, L1_prev, x1, x2)[at_risk], sd_L)
    L2[at_risk, m] <- rnorm(sum(at_risk),
                            .mu_L(beta_L, L2_prev, x1, x2)[at_risk], sd_L)
    L1_m <- L1[, m]; L2_m <- L2[, m]

    # (b) initiation among the not-yet-initiated, driven by X and current L
    not_yet <- at_risk & !initiated
    if (any(not_yet)) {
      lin_init <- beta_init["int"] + beta_init["x1"] * x1 + beta_init["x2"] * x2 +
        beta_init["L1"] * L1_m + beta_init["L2"] * L2_m
      fire <- not_yet & (runif(n) < plogis(lin_init))
      S[fire]         <- m
      initiated[fire] <- TRUE
    }

    # gated treatment indicator, shared by event and censoring hazards
    on_treatment <- (S <= tau) & (m >= S)

    # (c) event hazard: current L, lagged L, baseline X, gated treatment
    lin_ev <- beta_event["int"] +
      beta_event["x1"] * x1 + beta_event["x2"] * x2 +
      beta_event["L1"] * L1_m + beta_event["L2"] * L2_m +
      beta_event["L1lag"] * L1_prev + beta_event["L2lag"] * L2_prev +
      ifelse(on_treatment, beta_event["trt"], 0)
    die_now <- at_risk & (runif(n) < plogis(lin_ev))
    T_event[die_now] <- m
    alive[die_now]   <- FALSE

    # (d) natural censoring, only among the event-free this interval
    still_here <- at_risk & !die_now
    if (any(still_here)) {
      lin_c <- beta_cens["int"] +
        beta_cens["x1"] * x1 + beta_cens["x2"] * x2 +
        beta_cens["L1"] * L1_m + beta_cens["L2"] * L2_m +
        beta_cens["L1lag"] * L1_prev + beta_cens["L2lag"] * L2_prev +
        ifelse(on_treatment, beta_cens["trt"], 0)
      cens_now <- still_here & (runif(n) < plogis(lin_c))
      T_cens[cens_now] <- m
      uncens[cens_now] <- FALSE
    }

    # (e) roll lags forward (irrelevant for those who left observation)
    L1_prev <- ifelse(at_risk, L1_m, L1_prev)
    L2_prev <- ifelse(at_risk, L2_m, L2_prev)
  }

  # ---- observed-data summaries -------------------------------------------
  T_obs <- pmin(T_event, T_cens, t_star)
  delta <- as.integer(T_event <= t_star & T_event <= T_cens)
  cens  <- as.integer(T_cens  <= t_star & T_cens  <  T_event)
  A_tau <- as.integer(S <= tau)
  # grace period fully resolved? if FALSE, A_tau == 0 is unconfirmed
  tau_observed <- as.integer(T_obs >= tau | S <= tau)

  out <- data.frame(
    site         = site_id,
    id           = seq_len(n),
    x1           = x1,
    x2           = x2,
    S            = S,
    A_tau        = A_tau,
    tau_observed = tau_observed,
    T_event      = T_event,
    T_cens       = T_cens,
    T_obs        = T_obs,
    delta        = delta,
    cens         = cens
  )
  L1df <- as.data.frame(L1); names(L1df) <- paste0("L1_", seq_len(t_star))
  L2df <- as.data.frame(L2); names(L2df) <- paste0("L2_", seq_len(t_star))
  cbind(out, L1df, L2df)
}


simulate_multisite_tv <- function(K          = 3,
                                  n_per_site = 1000,
                                  tau        = DEFAULT_TAU,
                                  t_star     = DEFAULT_TSTAR,
                                  base_seed  = 2024,
                                  ...) {
  sites <- lapply(seq_len(K), function(k) {
    simulate_site_tv(n       = n_per_site,
                     tau     = tau,
                     t_star  = t_star,
                     site_id = k,
                     seed    = base_seed + k,
                     ...)
  })
  do.call(rbind, sites)
}


# =============================================================
# compute_truth_ccw_tv
#   Ground truth for the time-varying CCW estimand UNDER natural censoring.
#   Simulates a large natural cohort (including natural censoring, so the
#   cohort matches the observed-data DGP), then applies the same artificial
#   censoring rules and stabilized weighting as Fed_CCW_TVIPCW.R, using the
#   TRUE data-generating probabilities instead of fitted ones.
#
#   Total weight = SW_artificial x SW_natural. The natural factor is
#   evaluated per arm because the censoring hazard depends on the gated
#   treatment indicator (which differs between g=1 and g=0). With
#   beta_cens["trt"] = 0 the two arms coincide, but the code is general.
# =============================================================
compute_truth_ccw_tv <- function(N          = 2e6,
                                 tau        = DEFAULT_TAU,
                                 t_star     = DEFAULT_TSTAR,
                                 beta_event = DEFAULT_BETA_EVENT,
                                 beta_init  = DEFAULT_BETA_INIT,
                                 beta_L     = DEFAULT_BETA_L,
                                 beta_cens  = DEFAULT_BETA_CENS,
                                 sd_L       = DEFAULT_SD_L,
                                 seed       = 999) {

  set.seed(seed)
  M <- t_star

  x1 <- rnorm(N, 0, 1)
  x2 <- rbinom(N, 1, 0.4)

  S         <- rep(M + 1L, N)
  T_event   <- rep(M + 1L, N)
  T_cens    <- rep(M + 1L, N)
  initiated <- rep(FALSE, N)
  alive     <- rep(TRUE,  N)
  uncens    <- rep(TRUE,  N)
  L1_prev   <- rep(0, N)
  L2_prev   <- rep(0, N)

  Hden  <- matrix(NA_real_, N, tau)   # true initiation hazard, m <= tau
  Cden1 <- matrix(NA_real_, N, M)     # true censoring hazard under g = 1
  Cden0 <- matrix(NA_real_, N, M)     # true censoring hazard under g = 0
  Cobs  <- matrix(NA_real_, N, M)     # censoring hazard actually realized

  for (m in seq_len(M)) {
    at_risk <- alive & uncens
    if (!any(at_risk)) break

    L1_m <- rnorm(N, .mu_L(beta_L, L1_prev, x1, x2), sd_L)
    L2_m <- rnorm(N, .mu_L(beta_L, L2_prev, x1, x2), sd_L)

    not_yet <- at_risk & !initiated
    if (any(not_yet)) {
      lin_init <- beta_init["int"] + beta_init["x1"] * x1 + beta_init["x2"] * x2 +
        beta_init["L1"] * L1_m + beta_init["L2"] * L2_m
      p_init <- plogis(lin_init)
      if (m <= tau) Hden[not_yet, m] <- p_init[not_yet]
      fire <- not_yet & (runif(N) < p_init)
      S[fire]         <- m
      initiated[fire] <- TRUE
    }

    on_treatment <- (S <= tau) & (m >= S)
    lin_ev <- beta_event["int"] +
      beta_event["x1"] * x1 + beta_event["x2"] * x2 +
      beta_event["L1"] * L1_m + beta_event["L2"] * L2_m +
      beta_event["L1lag"] * L1_prev + beta_event["L2lag"] * L2_prev +
      ifelse(on_treatment, beta_event["trt"], 0)
    die_now <- at_risk & (runif(N) < plogis(lin_ev))
    T_event[die_now] <- m
    alive[die_now]   <- FALSE

    still_here <- at_risk & !die_now
    lin_c_base <- beta_cens["int"] +
      beta_cens["x1"] * x1 + beta_cens["x2"] * x2 +
      beta_cens["L1"] * L1_m + beta_cens["L2"] * L2_m +
      beta_cens["L1lag"] * L1_prev + beta_cens["L2lag"] * L2_prev

    # strategy-consistent treatment for the censoring hazard
    trt1_m <- as.integer(m >= pmin(S, tau))   # g=1: treated from min(S,tau)
    Cden1[, m] <- plogis(lin_c_base + beta_cens["trt"] * trt1_m)
    Cden0[, m] <- plogis(lin_c_base)          # g=0: never treated
    Cobs[, m]  <- plogis(lin_c_base + ifelse(on_treatment, beta_cens["trt"], 0))

    if (any(still_here)) {
      cens_now <- still_here & (runif(N) < Cobs[, m])
      T_cens[cens_now] <- m
      uncens[cens_now] <- FALSE
    }

    L1_prev <- L1_m
    L2_prev <- L2_m
  }

  Tcap    <- pmin(T_event, M)
  T_obs_o <- pmin(T_event, T_cens, M)
  tau_obs <- (T_obs_o >= tau) | (S <= tau)

  # ---- artificial censoring (CCW step 2) ---------------------------------
  # g=1: censor at tau if still uninitiated at tau.
  # g=0: censor at S-1 if initiation occurs within the grace period; if lost
  #      to follow-up before tau without confirmed non-initiation, censor at
  #      the last observed interval (natural weights carry the correction).
  Cinf <- M + 1L
  C1 <- ifelse(S >  tau, tau, Cinf)
  C0 <- ifelse(S <= tau, S - 1L,
               ifelse(tau_obs, Cinf, T_obs_o))

  # ---- true stabilized weights: artificial x natural ---------------------
  SW1 <- matrix(1, N, M)
  SW0 <- matrix(1, N, M)
  aw1 <- rep(1, N); aw0 <- rep(1, N)   # artificial factors
  nw1 <- rep(1, N); nw0 <- rep(1, N)   # natural factors

  for (m in seq_len(M)) {

    # (i) artificial factor, accrues while m <= tau
    if (m <= tau) {
      at_risk_m <- (S >= m)
      hn_m <- .clamp_prob(sum(S == m) / sum(at_risk_m))
      hd_m <- .clamp_prob(Hden[, m])

      idx_init_m <- at_risk_m & (S == m)
      idx_wait_m <- at_risk_m & (S >  m)

      aw1[idx_init_m] <- aw1[idx_init_m] * (hn_m / hd_m[idx_init_m])
      aw1[idx_wait_m] <- aw1[idx_wait_m] * ((1 - hn_m) / (1 - hd_m[idx_wait_m]))
      aw0[idx_wait_m] <- aw0[idx_wait_m] * ((1 - hn_m) / (1 - hd_m[idx_wait_m]))
    }

    # (ii) total weight DURING interval m uses natural factor through m-1
    SW1[, m] <- aw1 * nw1
    SW0[, m] <- aw0 * nw0

    # (iii) update natural factor for use in interval m+1 onward
    at_obs_m <- (T_cens >= m) & (pmin(T_event, M + 1L) >= m)
    n_at_obs <- sum(at_obs_m)
    cn_m <- if (n_at_obs > 0) .clamp_prob(sum(T_cens == m) / n_at_obs) else PROB_EPS
    cd1_m <- .clamp_prob(Cden1[, m])
    cd0_m <- .clamp_prob(Cden0[, m])

    not_cens_m <- (T_cens > m)
    nw1[not_cens_m] <- nw1[not_cens_m] * ((1 - cn_m) / (1 - cd1_m[not_cens_m]))
    nw0[not_cens_m] <- nw0[not_cens_m] * ((1 - cn_m) / (1 - cd0_m[not_cens_m]))
  }

  # ---- population-standardized discrete hazards --------------------------
  strategy_survival <- function(Cg, SWg) {
    lam <- numeric(M)
    for (m in seq_len(M)) {
      tm1 <- m - 1
      Rm <- (Tcap > tm1) & (Cg > tm1) & (T_cens > tm1)
      Dm <- Rm & (T_event == m) & (T_event <= Cg) & (T_event <= T_cens)
      w_m <- SWg[, m]
      rw <- sum(w_m[Rm])   # boolean subset avoids NA * FALSE propagation
      dw <- sum(w_m[Dm])
      lam[m] <- if (rw > 0) dw / rw else 0
    }
    c(1, cumprod(1 - lam))
  }

  S1 <- strategy_survival(C1, SW1)
  S0 <- strategy_survival(C0, SW0)

  psi1 <- 1 - S1[M + 1]
  psi0 <- 1 - S0[M + 1]

  list(
    S1 = S1, S0 = S0,
    psi1 = psi1, psi0 = psi0,
    RD = psi1 - psi0,
    RR = psi1 / psi0,
    OR = (psi1 / (1 - psi1)) / (psi0 / (1 - psi0)),
    RMST1 = sum(S1[seq_len(M)]), RMST0 = sum(S0[seq_len(M)]),
    RMST_diff = sum(S1[seq_len(M)]) - sum(S0[seq_len(M)])
  )
}