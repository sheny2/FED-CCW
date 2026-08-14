# =============================================================
# ArraySim4-4 discrete-time DGP with time-varying confounding.
#
# Interval m = 1..t_star.
#   x1, x2            baseline covariates (continuous, binary)
#   L1_m, L2_m        time-varying covariates, AR(1) in own lag + baseline X
#   treatment         absorbing; initiated at S, gated so the treated hazard
#                     applies iff A_tau == 1 and m >= S
#   no independent censoring; everyone followed to t_star unless event
#   treatment eligibility and covariate evolution stop at the terminal event
#
# Order within interval m:
#   (a) draw L1_m, L2_m
#   (b) if not yet initiated, draw initiation using X and L_m
#   (c) draw event using X, L_m, L_{m-1} and gated treatment status
#
# Output is wide: baseline columns plus L1_1..L1_M, L2_1..L2_M.
# =============================================================

source("params.R")

# ---- L transition mean for one interval --------------------------------
.mu_L <- function(beta_L, L_prev, x1, x2) {
  beta_L["int"] + beta_L["ar"] * L_prev + beta_L["x1"] * x1 + beta_L["x2"] * x2
}

simulate_site_tv <- function(n,
                             tau        = DEFAULT_TAU,
                             t_star     = DEFAULT_TSTAR,
                             beta_event = DEFAULT_BETA_EVENT,
                             beta_init  = DEFAULT_BETA_INIT,
                             beta_L     = DEFAULT_BETA_L,
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
  initiated <- rep(FALSE, n)
  alive     <- rep(TRUE,  n)

  L1_prev <- rep(0, n)
  L2_prev <- rep(0, n)

  for (m in seq_len(t_star)) {

    # (a) time-varying covariates among those alive at interval start.
    # Values after a terminal event remain NA and are never used downstream.
    live_idx <- which(alive)
    if (length(live_idx)) {
      L1[live_idx, m] <- rnorm(
        length(live_idx),
        .mu_L(beta_L, L1_prev[live_idx], x1[live_idx], x2[live_idx]), sd_L
      )
      L2[live_idx, m] <- rnorm(
        length(live_idx),
        .mu_L(beta_L, L2_prev[live_idx], x1[live_idx], x2[live_idx]), sd_L
      )
    }

    # (b) initiation among the not-yet-initiated, driven by X and current L
    not_yet <- alive & !initiated
    if (any(not_yet)) {
      lin_init <- beta_init["int"] + beta_init["x1"] * x1 + beta_init["x2"] * x2 +
        beta_init["L1"] * L1[, m] + beta_init["L2"] * L2[, m]
      fire <- not_yet & (runif(n) < plogis(lin_init))
      S[fire]         <- m
      initiated[fire] <- TRUE
    }

    # (c) event hazard: current L, lagged L, baseline X, gated treatment
    on_treatment <- (S <= tau) & (m >= S)
    lin_ev <- beta_event["int"] +
      beta_event["x1"] * x1 + beta_event["x2"] * x2 +
      beta_event["L1"] * L1[, m] + beta_event["L2"] * L2[, m] +
      beta_event["L1lag"] * L1_prev + beta_event["L2lag"] * L2_prev +
      ifelse(on_treatment, beta_event["trt"], 0)

    die_now <- alive & (runif(n) < plogis(lin_ev))
    T_event[die_now] <- m
    alive[die_now]   <- FALSE

    # Retain the last observed value for people who have just died; already
    # dead participants are never evaluated again.
    L1_prev[live_idx] <- L1[live_idx, m]
    L2_prev[live_idx] <- L2[live_idx, m]
  }

  out <- data.frame(
    site    = site_id,
    id      = seq_len(n),
    x1      = x1,
    x2      = x2,
    S       = S,
    A_tau   = as.integer(S <= tau),
    T_event = T_event,
    T_obs   = pmin(T_event, t_star),
    delta   = as.integer(T_event <= t_star)
  )
  L1df <- as.data.frame(L1); names(L1df) <- paste0("L1_", seq_len(t_star))
  L2df <- as.data.frame(L2); names(L2df) <- paste0("L2_", seq_len(t_star))
  cbind(out, L1df, L2df)
}


simulate_multisite_tv <- function(K          = 3,
                                  n_per_site = 1000,
                                  init_intercepts = NULL,
                                  tau        = DEFAULT_TAU,
                                  t_star     = DEFAULT_TSTAR,
                                  base_seed  = 2024,
                                  ...) {
  if (length(n_per_site) == 1L) n_per_site <- rep(n_per_site, K)
  if (length(n_per_site) != K)
    stop("n_per_site must have length 1 or K.")
  if (is.null(init_intercepts)) init_intercepts <- rep(NA_real_, K)
  if (length(init_intercepts) != K)
    stop("init_intercepts must be NULL or have length K.")

  sites <- lapply(seq_len(K), function(k) {
    site_beta_init <- list(...)$beta_init
    if (is.null(site_beta_init)) site_beta_init <- DEFAULT_BETA_INIT
    if (is.finite(init_intercepts[k]))
      site_beta_init["int"] <- init_intercepts[k]
    extra <- list(...)
    extra$beta_init <- NULL
    do.call(simulate_site_tv, c(list(
                     n       = n_per_site[k],
                     tau     = tau,
                     t_star  = t_star,
                     site_id = k,
                     seed    = base_seed + k,
                     beta_init = site_beta_init), extra))
  })
  do.call(rbind, sites)
}


# =============================================================
# compute_truth_ccw_tv
#   Ground truth for the time-varying CCW estimand. Simulates a large
#   natural cohort, applies the same artificial censoring and stabilized
#   weighting logic as Fed_CCW_TVIPCW.R, but uses the TRUE data-generating
#   initiation probabilities instead of fitted ones.
# =============================================================
compute_truth_ccw_tv <- function(N          = 2e6,
                                 tau        = DEFAULT_TAU,
                                 t_star     = DEFAULT_TSTAR,
                                 beta_event = DEFAULT_BETA_EVENT,
                                 beta_init  = DEFAULT_BETA_INIT,
                                 site_sizes = TEN_SITE_SIZES,
                                 init_intercepts = TEN_SITE_INIT_INTERCEPTS,
                                 beta_L     = DEFAULT_BETA_L,
                                 sd_L       = DEFAULT_SD_L,
                                 seed       = 999) {

  set.seed(seed)
  M <- t_star

  if (length(site_sizes) != length(init_intercepts))
    stop("site_sizes and init_intercepts must have the same length.")
  site_prob <- site_sizes / sum(site_sizes)
  site_N <- as.integer(floor(N * site_prob))
  site_N[length(site_N)] <- site_N[length(site_N)] + (N - sum(site_N))
  site <- rep(seq_along(site_N), times = site_N)
  init_intercept_i <- unname(init_intercepts[site])

  x1 <- rnorm(N, 0, 1)
  x2 <- rbinom(N, 1, 0.4)

  S         <- rep(M + 1L, N)
  T_event   <- rep(M + 1L, N)
  initiated <- rep(FALSE, N)
  alive     <- rep(TRUE,  N)
  L1_prev   <- rep(0, N)
  L2_prev   <- rep(0, N)

  # true initiation hazard per person-interval, needed only for m <= tau
  Hden <- matrix(NA_real_, N, tau)

  # ---- forward simulation of the natural cohort --------------------------
  for (m in seq_len(M)) {

    L1_m <- rep(NA_real_, N)
    L2_m <- rep(NA_real_, N)
    live_idx <- which(alive)
    if (length(live_idx)) {
      L1_m[live_idx] <- rnorm(
        length(live_idx),
        .mu_L(beta_L, L1_prev[live_idx], x1[live_idx], x2[live_idx]), sd_L
      )
      L2_m[live_idx] <- rnorm(
        length(live_idx),
        .mu_L(beta_L, L2_prev[live_idx], x1[live_idx], x2[live_idx]), sd_L
      )
    }

    not_yet <- alive & !initiated
    if (any(not_yet)) {
      lin_init <- init_intercept_i + beta_init["x1"] * x1 + beta_init["x2"] * x2 +
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

    die_now <- alive & (runif(N) < plogis(lin_ev))
    T_event[die_now] <- m
    alive[die_now]   <- FALSE

    L1_prev[live_idx] <- L1_m[live_idx]
    L2_prev[live_idx] <- L2_m[live_idx]
  }

  Tcap <- pmin(T_event, M)

  # ---- artificial censoring (CCW step 2) ---------------------------------
  # Non-initiators by tau deviate IN interval tau  -> censor at tau - 1.
  # Initiators at S <= tau deviate from g=0 IN S   -> censor at S - 1.
  Cinf <- M + 1L
  C1 <- ifelse(S >  tau, tau - 1L, Cinf)
  C0 <- ifelse(S <= tau, S   - 1L, Cinf)

  # ---- true stabilized time-varying IPCW ---------------------------------
  SW1 <- matrix(1, N, M)
  SW0 <- matrix(1, N, M)
  cw1 <- rep(1, N)
  cw0 <- rep(1, N)

  for (m in seq_len(M)) {
    if (m <= tau) {
      # Alive and not yet initiated at the start of interval m. Events occur
      # after the initiation decision within an interval, so T_event == m is
      # still eligible at m.
      at_risk_m <- (S >= m) & (T_event >= m)
      hn_m <- .clamp_prob(sum(S == m) / sum(at_risk_m))
      hd_m <- .clamp_prob(Hden[, m])

      idx_init_m <- at_risk_m & (S == m)
      idx_wait_m <- at_risk_m & (S >  m)

      # g = 1 (initiate by tau)
      cw1[idx_init_m] <- cw1[idx_init_m] * (hn_m / hd_m[idx_init_m])
      cw1[idx_wait_m] <- cw1[idx_wait_m] * ((1 - hn_m) / (1 - hd_m[idx_wait_m]))

      # g = 0 (never initiate by tau)
      cw0[idx_wait_m] <- cw0[idx_wait_m] * ((1 - hn_m) / (1 - hd_m[idx_wait_m]))
    }
    SW1[, m] <- cw1
    SW0[, m] <- cw0
  }

  # ---- population-standardized discrete hazards --------------------------
  strategy_survival <- function(Cg, SWg) {
    lam <- numeric(M)
    for (m in seq_len(M)) {
      tm1 <- m - 1
      Rm <- (Tcap > tm1) & (Cg > tm1)
      Dm <- Rm & (T_event == m) & (T_event <= Cg)
      w_m <- SWg[, m]
      rw <- sum(w_m * Rm)
      lam[m] <- if (rw > 0) sum(w_m * Dm) / rw else 0
    }
    c(1, cumprod(1 - lam))
  }

  S1 <- strategy_survival(C1, SW1)
  S0 <- strategy_survival(C0, SW0)

  psi1 <- 1 - S1[M + 1]
  psi0 <- 1 - S0[M + 1]
  RMST1 <- sum(S1[seq_len(M)])
  RMST0 <- sum(S0[seq_len(M)])

  list(
    S1 = S1, S0 = S0,
    psi1 = psi1, psi0 = psi0,
    RD = psi1 - psi0,
    RR = psi1 / psi0,
    OR = (psi1 / (1 - psi1)) / (psi0 / (1 - psi0)),
    RMST1 = RMST1, RMST0 = RMST0, RMST_diff = RMST1 - RMST0,
    N_used = N, site_N = site_N,
    site_weights = site_prob,
    init_intercepts = init_intercepts
  )
}
