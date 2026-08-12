# =============================================================
# Independent external-reference numerator construction.
# No outcome data are needed: only baseline/time-varying covariates and
# natural treatment initiation are simulated.
# =============================================================

source("params.R")

simulate_external_reference_numerator <- function(
    N = EXTERNAL_REFERENCE_N,
    t_star = EXTERNAL_REFERENCE_TSTAR,
    heterogeneity = EXTERNAL_REFERENCE_HETEROGENEITY,
    seed = EXTERNAL_REFERENCE_SEED,
    beta_init = DEFAULT_BETA_INIT,
    beta_L = DEFAULT_BETA_L,
    sd_L = DEFAULT_SD_L) {

  mix <- get_site_mix(heterogeneity, K = 3L)
  n_site <- as.integer(floor(N / 3L))
  if (n_site < 1L) stop("External reference N must be at least 3.")

  site_counts <- lapply(seq_len(3L), function(k) {
    set.seed(seed + k)
    x1 <- rnorm(n_site, mix$x1_mean[k], mix$x1_sd[k])
    x2 <- rbinom(n_site, 1, mix$x2_prob[k])
    initiated <- rep(FALSE, n_site)
    L1_prev <- rep(0, n_site)
    L2_prev <- rep(0, n_site)
    eligible_n <- initiated_n <- numeric(t_star)

    for (m in seq_len(t_star)) {
      L1_m <- rnorm(
        n_site,
        beta_L["int"] + beta_L["ar"] * L1_prev +
          beta_L["x1"] * x1 + beta_L["x2"] * x2,
        sd_L
      )
      L2_m <- rnorm(
        n_site,
        beta_L["int"] + beta_L["ar"] * L2_prev +
          beta_L["x1"] * x1 + beta_L["x2"] * x2,
        sd_L
      )

      eligible <- !initiated
      eligible_n[m] <- sum(eligible)
      lin_init <- beta_init["int"] +
        beta_init["x1"] * x1 + beta_init["x2"] * x2 +
        beta_init["L1"] * L1_m + beta_init["L2"] * L2_m
      fire <- eligible & runif(n_site) < plogis(lin_init)
      initiated_n[m] <- sum(fire)
      initiated[fire] <- TRUE

      L1_prev <- L1_m
      L2_prev <- L2_m
    }
    list(eligible = eligible_n, initiated = initiated_n)
  })

  eligible <- Reduce(`+`, lapply(site_counts, `[[`, "eligible"))
  initiated <- Reduce(`+`, lapply(site_counts, `[[`, "initiated"))
  hazard <- ifelse(
    eligible > 0,
    .clamp_prob(initiated / eligible),
    NA_real_
  )

  list(
    numerator_hazard = hazard,
    initiated = initiated,
    eligible = eligible,
    N_requested = N,
    N_used = n_site * 3L,
    t_star = as.integer(t_star),
    heterogeneity = heterogeneity,
    seed = as.integer(seed),
    beta_init = beta_init,
    beta_L = beta_L,
    sd_L = sd_L,
    site_mix = mix,
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
}

validate_external_reference <- function(ref, t_star = NULL, tau = NULL) {
  required <- c(
    "numerator_hazard", "N_used", "t_star", "heterogeneity", "seed",
    "beta_init", "beta_L", "sd_L"
  )
  missing <- setdiff(required, names(ref))
  if (length(missing))
    stop("External reference is missing: ", paste(missing, collapse = ", "))
  if (any(!is.finite(ref$numerator_hazard)) ||
      any(ref$numerator_hazard <= 0 | ref$numerator_hazard >= 1))
    stop("External numerator hazards must be finite and strictly between 0 and 1.")
  if (!is.null(t_star) && length(ref$numerator_hazard) < t_star)
    stop("External numerator does not cover t_star=", t_star, ".")
  if (!is.null(tau) && length(ref$numerator_hazard) < tau)
    stop("External numerator does not cover tau=", tau, ".")
  invisible(TRUE)
}

load_external_reference <- function(
    path = EXTERNAL_REFERENCE_FILE,
    t_star = NULL,
    tau = NULL) {
  if (!file.exists(path))
    stop(
      "Missing external numerator file: ", path,
      ". Run: Rscript prepare_external_reference.R"
    )
  ref <- readRDS(path)
  validate_external_reference(ref, t_star = t_star, tau = tau)
  ref
}
