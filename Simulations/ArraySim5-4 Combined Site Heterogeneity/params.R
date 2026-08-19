# =============================================================
# ArraySim5-4 Combined Site Heterogeneity params.R -- DGP parameters.
# Every other file pulls its defaults from here; nothing is
# duplicated or hard-coded downstream.
# =============================================================

# ---- user-configurable study design ------------------------------------
# Each level jointly changes site size, baseline patient mix, and the residual
# site-practice intercept in the treatment-initiation model. Total N remains
# 5,000 at every level. The gradients are aligned by site: later/larger sites
# have larger X values and higher residual initiation propensity.
.site_mix <- function(x1_mean, x1_sd, x2_prob) {
  data.frame(site = 1:5, x1_mean = x1_mean,
             x1_sd = x1_sd, x2_prob = x2_prob)
}

HETEROGENEITY_SETTINGS <- list(
  low = list(
    site_sizes = setNames(c(900L, 950L, 1000L, 1050L, 1100L), 1:5),
    patient_mix = .site_mix(
      c(-0.25, -0.125, 0, 0.125, 0.25),
      c(0.95, 0.975, 1, 1.025, 1.05),
      c(0.35, 0.375, 0.40, 0.425, 0.45)
    ),
    init_intercepts = setNames(c(-3.25, -3.125, -3, -2.875, -2.75), 1:5)
  ),
  moderate = list(
    site_sizes = setNames(c(600L, 800L, 1000L, 1200L, 1400L), 1:5),
    patient_mix = .site_mix(
      c(-0.75, -0.375, 0, 0.375, 0.75),
      c(0.85, 0.925, 1, 1.075, 1.15),
      c(0.20, 0.30, 0.40, 0.50, 0.60)
    ),
    init_intercepts = setNames(c(-3.75, -3.375, -3, -2.625, -2.25), 1:5)
  ),
  high = list(
    site_sizes = setNames(c(400L, 600L, 800L, 1200L, 2000L), 1:5),
    patient_mix = .site_mix(
      c(-1.50, -0.75, 0, 0.75, 1.50),
      c(0.70, 0.85, 1, 1.15, 1.30),
      c(0.10, 0.25, 0.40, 0.55, 0.70)
    ),
    init_intercepts = setNames(c(-4.50, -3.75, -3, -2.25, -1.50), 1:5)
  )
)

get_heterogeneity_setting <- function(level) {
  level <- as.character(level)
  if (length(level) != 1L || !level %in% names(HETEROGENEITY_SETTINGS))
    stop("Unknown combined heterogeneity level: ", paste(level, collapse = ", "))
  HETEROGENEITY_SETTINGS[[level]]
}

DEFAULT_HETEROGENEITY_LEVEL <- "moderate"
.default_setting <- get_heterogeneity_setting(DEFAULT_HETEROGENEITY_LEVEL)
DEFAULT_SITE_SIZES <- .default_setting$site_sizes
DEFAULT_SITE_INIT_INTERCEPTS <- .default_setting$init_intercepts
DEFAULT_PATIENT_MIX <- .default_setting$patient_mix

SIM_TAU_VALUES <- c(5L, 10L, 15L)
SIM_BETA_TRT_VALUES <- c(-1.0, -0.7, -0.5)
SIM_N_REPS <- 500L
SIM_TSTAR <- 25L
SIM_TRUTH_N <- 3e6
SIM_TRUTH_SEED <- 321L
SIM_BASE_SEED <- 2026L
SIM_TASK_SEED_STRIDE <- 10000L

# Function defaults used for interactive calls.
DEFAULT_TAU <- SIM_TAU_VALUES[[1L]]
DEFAULT_TSTAR <- SIM_TSTAR

DEFAULT_BETA_EVENT <- c(int   = -2.5,
                        x1    =  0.5,
                        x2    =  0.4,
                        L1    =  0.4,
                        L2    =  0.3,
                        L1lag =  0.2,
                        L2lag =  0.15,
                        trt   = -0.7)

DEFAULT_BETA_INIT <- c(int = -3,
                       x1  =  0.8,
                       x2  = -0.3,
                       L1  =  0.5,
                       L2  =  0.4)

# L transition: AR(1) in own lag + baseline covariates. No treatment effect
# on L (beta_L_trt was zero and has been removed from the model entirely).
DEFAULT_BETA_L <- c(int = 0.0,
                    ar  = 0.6,
                    x1  = 0.3,
                    x2  = 0.2)

DEFAULT_SD_L <- 1.0

# DEFAULT_TRUNC <- c(0.001, 0.999)
DEFAULT_TRUNC <- c(0, 1)

# Numerical guard for probabilities entering weight ratios.
PROB_EPS <- 1e-6

.clamp_prob <- function(p, eps = PROB_EPS) pmin(pmax(p, eps), 1 - eps)

# Replace the treatment coefficient of the event model without restating
# the rest of the vector.
set_beta_trt <- function(beta_trt, beta_event = DEFAULT_BETA_EVENT) {
  beta_event["trt"] <- beta_trt
  beta_event
}

# Construct the full array-job grid from the settings above. Keeping this in
# params.R prevents the runner, aggregation script, and submission command
# from drifting apart when a user changes a scenario vector.
simulation_grid <- function() {
  grid <- expand.grid(
    tau = SIM_TAU_VALUES,
    beta_trt = SIM_BETA_TRT_VALUES,
    heterogeneity = names(HETEROGENEITY_SETTINGS),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid
}

validate_params <- function() {
  n_sites <- length(DEFAULT_SITE_SIZES)
  if (!n_sites || any(!is.finite(DEFAULT_SITE_SIZES)) ||
      any(DEFAULT_SITE_SIZES <= 0) || any(DEFAULT_SITE_SIZES %% 1 != 0))
    stop("DEFAULT_SITE_SIZES must contain positive integer site sizes.")
  if (length(DEFAULT_SITE_INIT_INTERCEPTS) != n_sites ||
      any(!is.finite(DEFAULT_SITE_INIT_INTERCEPTS)))
    stop("DEFAULT_SITE_INIT_INTERCEPTS must contain one finite value per site.")
  if (!length(SIM_TAU_VALUES) || any(SIM_TAU_VALUES <= 0) ||
      any(SIM_TAU_VALUES >= SIM_TSTAR))
    stop("SIM_TAU_VALUES must be positive and strictly less than SIM_TSTAR.")
  if (!length(SIM_BETA_TRT_VALUES) || any(!is.finite(SIM_BETA_TRT_VALUES)))
    stop("SIM_BETA_TRT_VALUES must contain finite values.")
  if (n_sites != 5L)
    stop("ArraySim5-4 requires exactly five sites.")
  if (!length(HETEROGENEITY_SETTINGS) ||
      is.null(names(HETEROGENEITY_SETTINGS)))
    stop("HETEROGENEITY_SETTINGS must be a named list.")
  required_mix_cols <- c("site", "x1_mean", "x1_sd", "x2_prob")
  for (level in names(HETEROGENEITY_SETTINGS)) {
    setting <- HETEROGENEITY_SETTINGS[[level]]
    sizes <- setting$site_sizes
    ints <- setting$init_intercepts
    mix <- setting$patient_mix
    if (length(sizes) != n_sites || any(!is.finite(sizes)) ||
        any(sizes <= 0) || any(sizes %% 1 != 0) || sum(sizes) != 5000L)
      stop("Invalid site sizes for level: ", level)
    if (length(ints) != n_sites || any(!is.finite(ints)))
      stop("Invalid initiation intercepts for level: ", level)
    if (!identical(names(mix), required_mix_cols) || nrow(mix) != n_sites ||
        !identical(as.integer(mix$site), seq_len(n_sites)) ||
        any(!is.finite(as.matrix(mix))) || any(mix$x1_sd <= 0) ||
        any(mix$x2_prob < 0 | mix$x2_prob > 1))
      stop("Invalid patient mix for level: ", level)
  }
  if (length(DEFAULT_TRUNC) != 2L || any(!is.finite(DEFAULT_TRUNC)) ||
      DEFAULT_TRUNC[[1L]] < 0 || DEFAULT_TRUNC[[2L]] > 1 ||
      DEFAULT_TRUNC[[1L]] >= DEFAULT_TRUNC[[2L]])
    stop("DEFAULT_TRUNC must be an increasing two-value interval in [0, 1].")
  positive_scalars <- c(
    SIM_N_REPS = SIM_N_REPS, SIM_TSTAR = SIM_TSTAR,
    SIM_TRUTH_N = SIM_TRUTH_N, SIM_TASK_SEED_STRIDE = SIM_TASK_SEED_STRIDE
  )
  if (any(!is.finite(positive_scalars)) || any(positive_scalars <= 0))
    stop("Replicate, horizon, truth-size, and seed-stride settings must be positive.")
  invisible(TRUE)
}
