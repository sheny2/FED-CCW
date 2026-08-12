# =============================================================
# ArraySim4-3 Ten Small Sites params.R -- DGP parameters.
# Every other file pulls its defaults from here; nothing is
# duplicated or hard-coded downstream.
# =============================================================

DEFAULT_TAU    <- 3
DEFAULT_TSTAR  <- 12

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

# ---- confounding strength ----------------------------------------------
# Confounding is governed by beta_init, the covariate coefficients in the
# initiation model. Those same covariates (x1, x2 at baseline; L1, L2
# time-varying) also drive the event hazard through beta_event, so scaling
# the beta_init slopes scales the strength of confounding.
#
# The intercept is held fixed because it mainly sets the marginal initiation
# rate rather than the covariate-outcome association; scaling it too would
# confound "more confounding" with "more/less treatment uptake". All four
# slopes are scaled together so baseline and time-varying confounding move
# in step. Note x2 is negative, so a common multiplier scales magnitudes
# without flipping any sign.
CONFOUNDING_SLOPES <- c("x1", "x2", "L1", "L2")

scale_confounding <- function(mult, beta_init = DEFAULT_BETA_INIT) {
  beta_init[CONFOUNDING_SLOPES] <- beta_init[CONFOUNDING_SLOPES] * mult
  beta_init
}

# Multipliers defining the small / medium / strong scenarios. medium = 1
# reproduces DEFAULT_BETA_INIT exactly.
DEFAULT_CONF_MULTS <- c(small = 0.4, medium = 1.0, strong = 2.0)
DEFAULT_CONF_MULTS <- c(small = 0.5, medium = 1.0, strong = 1.5)

# ---- ten noisy local sites ----------------------------------------------
# Total N remains 3,000, but it is divided into ten independently estimated
# site nuisance models and ten completed local CCW curves. Treatment uptake
# ranges smoothly from rare at site 1 to common at site 10.
TEN_SITE_SIZES <- setNames(rep(300L, 10L), seq_len(10L))
TEN_SITE_SIZES <- setNames(rep(100L, 10L), seq_len(10L))
TEN_SITE_INIT_INTERCEPTS <- setNames(
  seq(-4.5, -1.5, length.out = 10L), seq_len(10L)
)
