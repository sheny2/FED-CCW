# =============================================================
# ArraySim4-5 parameters: sample size x outcome-frequency study.
# =============================================================

DEFAULT_TAU <- 4L
DEFAULT_TSTAR <- 25L

DEFAULT_BETA_EVENT <- c(
  int = -3.0,
  x1 = 0.5,
  x2 = 0.4,
  L1 = 0.4,
  L2 = 0.3,
  L1lag = 0.2,
  L2lag = 0.15,
  trt = -0.7
)

DEFAULT_BETA_INIT <- c(
  int = -3.0,
  x1 = 0.8,
  x2 = -0.3,
  L1 = 0.5,
  L2 = 0.4
)

DEFAULT_BETA_L <- c(int = 0.0, ar = 0.6, x1 = 0.3, x2 = 0.2)
DEFAULT_SD_L <- 1.0
DEFAULT_TRUNC <- c(0, 1)
PROB_EPS <- 1e-6

.clamp_prob <- function(p, eps = PROB_EPS) pmin(pmax(p, eps), 1 - eps)

set_beta_trt <- function(beta_trt, beta_event = DEFAULT_BETA_EVENT) {
  beta_event["trt"] <- beta_trt
  beta_event
}

set_event_intercept <- function(intercept, beta_event = DEFAULT_BETA_EVENT) {
  beta_event["int"] <- intercept
  beta_event
}

CONFOUNDING_SLOPES <- c("x1", "x2", "L1", "L2")

scale_confounding <- function(mult, beta_init = DEFAULT_BETA_INIT) {
  beta_init[CONFOUNDING_SLOPES] <- beta_init[CONFOUNDING_SLOPES] * mult
  beta_init
}

# ---- ArraySim4-5 factorial design ----------------------------------------
STUDY_K <- 10L
SAMPLE_SIZE_LEVELS <- c(low = 100L, large = 500L)
OUTCOME_LEVELS <- c(rare = -7.0, common = -3.0)
STUDY_TAUS <- c(4L, 6L, 8L)
STUDY_BETA_TRT <- -0.7
STUDY_TSTAR <- 25L
STUDY_N_REPS <- 300
STUDY_TRUTH_N <- 3e6
STUDY_CONF_MULT <- 1.0
STUDY_INIT_INTERCEPT <- -3.0

make_study_grid <- function() {
  expand.grid(
    tau = STUDY_TAUS,
    sample_size_scenario = names(SAMPLE_SIZE_LEVELS),
    outcome_scenario = names(OUTCOME_LEVELS),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}

# Defaults for the shared one-replicate wrapper. Array jobs replace sample
# size and event intercept with the selected factorial cell.
TEN_SITE_SIZES <- setNames(rep(SAMPLE_SIZE_LEVELS[["low"]], STUDY_K),
                           seq_len(STUDY_K))
TEN_SITE_INIT_INTERCEPTS <- setNames(rep(STUDY_INIT_INTERCEPT, STUDY_K),
                                    seq_len(STUDY_K))
