# =============================================================
# params.R -- single source of truth for the NATURAL-CENSORING DGP.
# Every other file pulls its defaults from here; nothing is
# duplicated or hard-coded downstream.
#
# This study compares SIX correctly-specified methods (no
# misspecification axis) across the tau x beta_trt x confounding grid:
#   fed-CCW-TVIPCW, pooled-CCW-TVIPCW, CC-gcomp,
#   plain g-comp, per-protocol (naive), IPCW-without-cloning (naive).
#
# Natural (loss-to-follow-up) censoring is present and, by default,
# INDEPENDENT of treatment (tier 2: beta_cens["trt"] = 0).
# =============================================================

DEFAULT_TAU    <- 3
DEFAULT_TSTAR  <- 12

DEFAULT_BETA_EVENT <- c(int   = -2.5,
                        x1    =  0.5,
                        x2    =  0.4,
                        L1    =  0.4,   # current L1
                        L2    =  0.3,   # current L2
                        L1lag =  0.2,   # lagged L1
                        L2lag =  0.15,  # lagged L2
                        trt   = -0.7)

DEFAULT_BETA_INIT <- c(int = -3,
                       x1  =  0.8,
                       x2  = -0.3,
                       L1  =  0.5,      # L drives initiation
                       L2  =  0.4)

# L transition: AR(1) in own lag + baseline covariates. Treatment does not
# affect L (pure time-varying confounder).
DEFAULT_BETA_L <- c(int = 0.0,
                    ar  = 0.6,
                    x1  = 0.3,
                    x2  = 0.2)

# Natural-censoring hazard. trt = 0 => censoring independent of treatment
# (tier 2). The intercept sets the marginal dropout rate.
DEFAULT_BETA_CENS <- c(int   = -4.0,
                       x1    =  0.3,
                       x2    =  0.2,
                       L1    =  0.3,
                       L2    =  0.2,
                       L1lag =  0.1,
                       L2lag =  0.1,
                       trt   =  0.0)

DEFAULT_SD_L <- 1.0

DEFAULT_TRUNC <- c(0.01, 0.99)

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
# time-varying) also drive the event hazard, so scaling the beta_init slopes
# scales the strength of confounding. The intercept is held fixed so that
# "more confounding" does not also mean "more/less treatment uptake". x2 is
# negative, so a common multiplier scales magnitudes without flipping signs.
CONFOUNDING_SLOPES <- c("x1", "x2", "L1", "L2")

scale_confounding <- function(mult, beta_init = DEFAULT_BETA_INIT) {
  beta_init[CONFOUNDING_SLOPES] <- beta_init[CONFOUNDING_SLOPES] * mult
  beta_init
}

# Multipliers defining the small / medium / strong scenarios. medium = 1
# reproduces DEFAULT_BETA_INIT exactly.
DEFAULT_CONF_MULTS <- c(small = 0.4, medium = 1.0, strong = 2.0)
