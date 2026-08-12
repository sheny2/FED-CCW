# =============================================================
# params.R -- single source of truth for DGP parameters.
# Every other file pulls its defaults from here; nothing is
# duplicated or hard-coded downstream.
#
# This is the trimmed two-method study: only G-computation and
# federated CCW are estimated, each in a correct and a misspecified
# form. The misspecification axis (no_tv / coarse_L) is applied to
# G-comp's OUTCOME model and to fed-CCW's WEIGHT model in the same
# cell, so a single axis degrades both nuisance models together.
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

# ---- nuisance-model misspecification -----------------------------------
# A single axis degrades the time-varying covariates in BOTH nuisance
# models: G-comp's outcome hazard and fed-CCW's weight (initiation-hazard)
# denominator. The linear predictor is otherwise the correct one, so the
# effect is attributable to how L enters.
#
#   correct  -- L1, L2, L1lag, L2lag enter linearly (matches the DGP).
#   no_tv    -- drop the time-varying covariates entirely, leaving only
#               baseline X. L is the time-varying confounder, so this is
#               the most damaging error for both methods.
#   coarse_L -- keep L but coarsen each to a binary indicator of being above
#               its pooled median, discarding within-category variation
#               (measurement-error-style misspecification).
#
# The specs are expressed as the covariate block that replaces
#   L1 + L2 + L1lag + L2lag
# in each nuisance model's linear predictor. factor(m) / interval baseline,
# x1, x2, and (for the outcome model) trt are always retained.
GCOMP_MISSPEC_SPECS <- c("correct", "no_tv", "coarse_L")
DEFAULT_MISSPEC     <- "correct"

# Covariate block for a given spec. "" means no L terms at all (no_tv).
.tv_block <- function(misspec = DEFAULT_MISSPEC) {
  switch(misspec,
    correct  = "L1 + L2 + L1lag + L2lag",
    no_tv    = "",
    coarse_L = "L1_hi + L2_hi + L1lag_hi + L2lag_hi",
    stop("Unknown misspecification '", misspec, "'. Options: ",
         paste(GCOMP_MISSPEC_SPECS, collapse = ", "))
  )
}

# Assemble a full RHS by inserting the tv block into a template that already
# carries the always-present terms. `base` is the leading terms
# (e.g. "factor(m) + x1 + x2" for the outcome model, "m + x1 + x2" for the
# weight model); `tail` is appended after the tv block (e.g. "+ trt").
.build_rhs <- function(misspec, base, tail = "") {
  blk <- .tv_block(misspec)
  parts <- c(base, if (nzchar(blk)) blk else NULL, if (nzchar(tail)) tail else NULL)
  paste(parts, collapse = " + ")
}