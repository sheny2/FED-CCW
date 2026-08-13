# =============================================================
# ArraySim4-7 parameters: fragmentation x treatment-support study.
# =============================================================

DEFAULT_TAU <- 4L
DEFAULT_TSTAR <- 25L

DEFAULT_BETA_EVENT <- c(
  int = -2.5, x1 = 0.5, x2 = 0.4, L1 = 0.4, L2 = 0.3,
  L1lag = 0.2, L2lag = 0.15, trt = -0.7
)
DEFAULT_BETA_INIT <- c(int = -4.5, x1 = 0, x2 = 0, L1 = 0, L2 = 0)
DEFAULT_BETA_L <- c(int = 0, ar = 0.6, x1 = 0.3, x2 = 0.2)
DEFAULT_SD_L <- 1
DEFAULT_TRUNC <- c(0.01, 0.99)
PROB_EPS <- 1e-6

.clamp_prob <- function(p, eps = PROB_EPS) pmin(pmax(p, eps), 1 - eps)

set_beta_trt <- function(beta_trt, beta_event = DEFAULT_BETA_EVENT) {
  beta_event["trt"] <- beta_trt
  beta_event
}

CONFOUNDING_SLOPES <- c("x1", "x2", "L1", "L2")
scale_confounding <- function(mult, beta_init = c(
  int = -4.5, x1 = 0.8, x2 = -0.3, L1 = 0.5, L2 = 0.4
)) {
  beta_init[CONFOUNDING_SLOPES] <- beta_init[CONFOUNDING_SLOPES] * mult
  beta_init
}

STUDY_TAU <- 4L
STUDY_TSTAR <- 25L
STUDY_TOTAL_N <- 5000L
STUDY_BETA_TRT <- -0.7
STUDY_CONF_MULT <- 0
STUDY_N_REPS <- 300L
STUDY_TRUTH_N <- 5e6
STUDY_TRUTH_SEED <- 47001L
STUDY_TRUNC <- c(0.01, 0.99)

# Holding total N fixed isolates data fragmentation across sites.
FRAGMENTATION_LEVELS <- list(
  `10_sites_x_500` = rep(500L, 10L),
  `20_sites_x_250` = rep(250L, 20L),
  `50_sites_x_100` = rep(100L, 50L),
  `100_sites_x_50` = rep(50L, 100L)
)

INITIATION_LEVELS <- c(
  moderate = -3.5,
  sparse = -4.5
)

homogeneous_site_mix <- function(K) {
  data.frame(
    site = seq_len(K),
    x1_mean = rep(0, K),
    x1_sd = rep(1, K),
    x2_prob = rep(0.4, K)
  )
}

validate_site_mix <- function(site_mix, K) {
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
  if (!level %in% names(FRAGMENTATION_LEVELS))
    stop("Unknown fragmentation level: ", level)
  FRAGMENTATION_LEVELS[[level]]
}

get_init_intercept <- function(level) {
  if (!level %in% names(INITIATION_LEVELS))
    stop("Unknown initiation level: ", level)
  unname(INITIATION_LEVELS[[level]])
}

make_study_grid <- function() {
  expand.grid(
    fragmentation = names(FRAGMENTATION_LEVELS),
    initiation = names(INITIATION_LEVELS),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}

# Defaults for interactive one-replicate calls.
STUDY_SITE_SIZES <- get_site_sizes("10_sites_x_500")
STUDY_INIT_INTERCEPTS <- rep(get_init_intercept("sparse"),
                             length(STUDY_SITE_SIZES))
