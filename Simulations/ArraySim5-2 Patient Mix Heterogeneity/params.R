# =============================================================
# ArraySim5-2 Patient Mix Heterogeneity params.R -- DGP parameters.
# Every other file pulls its defaults from here; nothing is
# duplicated or hard-coded downstream.
# =============================================================

# ---- user-configurable study design ------------------------------------
# Edit this section to change the simulation. The number of sites is derived
# from DEFAULT_SITE_SIZES, so unequal site sizes are also supported.
DEFAULT_SITE_SIZES <- setNames(rep(1000L, 3L), seq_len(3L))
DEFAULT_SITE_INIT_INTERCEPTS <- setNames(
  seq(-4.5, -1.5, length.out = length(DEFAULT_SITE_SIZES)),
  names(DEFAULT_SITE_SIZES)
)

SIM_TAU_VALUES <- c(5L, 10L, 15L)
SIM_BETA_TRT_VALUES <- c(-1.0, -0.7, -0.5)
SIM_N_REPS <- 100L
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

# ---- site-specific patient-mix heterogeneity ----------------------------
# Each row gives one site's baseline X distribution. Confounding-model
# coefficients are fixed at DEFAULT_BETA_INIT in every scenario; only the
# across-site distribution of baseline covariates changes.
DEFAULT_PATIENT_MIX <- list(
  low = data.frame(
    site = 1:3,
    x1_mean = c(-0.25, 0.00, 0.25),
    x1_sd = c(0.95, 1.00, 1.05),
    x2_prob = c(0.35, 0.40, 0.45)
  ),
  moderate = data.frame(
    site = 1:3,
    x1_mean = c(-0.75, 0.00, 0.75),
    x1_sd = c(0.85, 1.00, 1.15),
    x2_prob = c(0.20, 0.40, 0.60)
  ),
  high = data.frame(
    site = 1:3,
    x1_mean = c(-1.50, 0.00, 1.50),
    x1_sd = c(0.70, 1.00, 1.30),
    x2_prob = c(0.10, 0.40, 0.70)
  )
)

get_patient_mix <- function(level) {
  level <- as.character(level)
  if (length(level) != 1L || !level %in% names(DEFAULT_PATIENT_MIX))
    stop("Unknown patient-mix heterogeneity level: ", paste(level, collapse = ", "))
  DEFAULT_PATIENT_MIX[[level]]
}

# Construct the full array-job grid from the settings above. Keeping this in
# params.R prevents the runner, aggregation script, and submission command
# from drifting apart when a user changes a scenario vector.
simulation_grid <- function() {
  grid <- expand.grid(
    tau = SIM_TAU_VALUES,
    beta_trt = SIM_BETA_TRT_VALUES,
    heterogeneity = names(DEFAULT_PATIENT_MIX),
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
  if (!length(DEFAULT_PATIENT_MIX) || is.null(names(DEFAULT_PATIENT_MIX)))
    stop("DEFAULT_PATIENT_MIX must be a named list.")
  required_mix_cols <- c("site", "x1_mean", "x1_sd", "x2_prob")
  for (level in names(DEFAULT_PATIENT_MIX)) {
    mix <- DEFAULT_PATIENT_MIX[[level]]
    if (!identical(names(mix), required_mix_cols) || nrow(mix) != n_sites ||
        !identical(as.integer(mix$site), seq_len(n_sites)) ||
        any(!is.finite(as.matrix(mix))) || any(mix$x1_sd <= 0) ||
        any(mix$x2_prob < 0 | mix$x2_prob > 1))
      stop("Invalid patient-mix specification for level: ", level)
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
