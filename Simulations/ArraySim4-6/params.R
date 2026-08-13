# =============================================================
# ArraySim4-6 parameters: site-size balance x patient-mix study.
# =============================================================

DEFAULT_TAU <- 4L
DEFAULT_TSTAR <- 25L

DEFAULT_BETA_EVENT <- c(
  int = -2.5, x1 = 0.5, x2 = 0.4, L1 = 0.4, L2 = 0.3,
  L1lag = 0.2, L2lag = 0.15, trt = -0.7
)
DEFAULT_BETA_INIT <- c(int = -3.0, x1 = 0.8, x2 = -0.3, L1 = 0.5, L2 = 0.4)
DEFAULT_BETA_L <- c(int = 0.0, ar = 0.6, x1 = 0.3, x2 = 0.2)
DEFAULT_SD_L <- 1.0
DEFAULT_TRUNC <- c(0, 1)
PROB_EPS <- 1e-6

.clamp_prob <- function(p, eps = PROB_EPS) pmin(pmax(p, eps), 1 - eps)

set_beta_trt <- function(beta_trt, beta_event = DEFAULT_BETA_EVENT) {
  beta_event["trt"] <- beta_trt
  beta_event
}

CONFOUNDING_SLOPES <- c("x1", "x2", "L1", "L2")
scale_confounding <- function(mult, beta_init = DEFAULT_BETA_INIT) {
  beta_init[CONFOUNDING_SLOPES] <- beta_init[CONFOUNDING_SLOPES] * mult
  beta_init
}

STUDY_K <- 3L
STUDY_TAUS <- c(4L, 6L, 8L)
STUDY_BETA_TRT <- -0.7
STUDY_TSTAR <- 25L
STUDY_N_REPS <- 100L
STUDY_TRUTH_N <- 3e6
STUDY_CONF_MULT <- 1.0
STUDY_INIT_INTERCEPT <- -4
STUDY_TRUNC <- c(0, 1)

SITE_SIZE_LEVELS <- list(
  balanced = c(`1` = 500L, `2` = 500L, `3` = 500L),
  unbalanced = c(`1` = 50L, `2` = 100L, `3` = 1350L)
)

homogeneous_site_mix <- function(K = STUDY_K) {
  data.frame(
    site = seq_len(K),
    x1_mean = rep(0, K),
    x1_sd = rep(1, K),
    x2_prob = rep(0.4, K)
  )
}

# High-heterogeneity patient mix from the supplied design screenshot.
HETEROGENEOUS_SITE_MIX <- data.frame(
  site = 1:3,
  x1_mean = c(-1.50, 0.00, 1.50),
  x1_sd = c(0.70, 1.00, 1.30),
  x2_prob = c(0.10, 0.40, 0.70)
)

PATIENT_MIX_LEVELS <- list(
  homogeneous = homogeneous_site_mix(STUDY_K),
  heterogeneous = HETEROGENEOUS_SITE_MIX
)

validate_site_mix <- function(site_mix, K = STUDY_K) {
  needed <- c("site", "x1_mean", "x1_sd", "x2_prob")
  if (!is.data.frame(site_mix) || !all(needed %in% names(site_mix)))
    stop("site_mix must contain: ", paste(needed, collapse = ", "))
  if (nrow(site_mix) != K || !identical(as.integer(site_mix$site), seq_len(K)))
    stop("site_mix must have one ordered row per site.")
  if (any(!is.finite(site_mix$x1_mean)) || any(site_mix$x1_sd <= 0) ||
      any(site_mix$x2_prob <= 0 | site_mix$x2_prob >= 1))
    stop("Invalid site-mix parameters.")
  invisible(TRUE)
}

get_site_sizes <- function(level) {
  if (!level %in% names(SITE_SIZE_LEVELS)) stop("Unknown size level: ", level)
  SITE_SIZE_LEVELS[[level]]
}

get_site_mix <- function(level) {
  if (!level %in% names(PATIENT_MIX_LEVELS)) stop("Unknown mix level: ", level)
  PATIENT_MIX_LEVELS[[level]]
}

make_study_grid <- function() {
  expand.grid(
    tau = STUDY_TAUS,
    size_scenario = names(SITE_SIZE_LEVELS),
    mix_scenario = names(PATIENT_MIX_LEVELS),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}

STUDY_SITE_SIZES <- get_site_sizes("balanced")
STUDY_INIT_INTERCEPTS <- setNames(rep(STUDY_INIT_INTERCEPT, STUDY_K),
                                 seq_len(STUDY_K))
