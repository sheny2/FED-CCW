# =============================================================
# params.R -- single source of truth for DGP parameters.
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

# ---- g-computation outcome-model misspecification -----------------------
# The DGP hazard is linear in (x1, x2, L1, L2, L1lag, L2lag) plus a gated
# treatment term, with no interactions and no nonlinearity. The "correct"
# outcome model therefore matches the DGP exactly, and misspecification has
# to be induced deliberately.
#
# Each option below is the RIGHT-HAND SIDE of the outcome hazard glm. The
# interval baseline factor(m), the treatment term trt, and x1/x2 are kept in
# every variant so the treatment contrast remains estimable and the arms stay
# comparable; what varies is how the time-varying covariates enter.
#
#   correct     -- matches the DGP. Reference case.
#   no_lag      -- drops the lagged L terms. The DGP hazard depends on
#                  L_{m-1} through beta_event[c("L1lag","L2lag")], so this
#                  omits a genuine cause of the outcome that is also
#                  correlated with treatment timing.
#   no_tv       -- drops the time-varying covariates entirely, leaving only
#                  baseline X. This is the classic "adjust for baseline only"
#                  error and should be the most damaging, since L is the
#                  time-varying confounder.
#   coarse_L    -- keeps L but coarsens it to a binary indicator of being
#                  above the interval median, discarding within-category
#                  variation. Measurement-error-style misspecification.
#   nonlinear_L -- fits a quadratic in L when the truth is linear. This is
#                  OVER-specification rather than under-specification and is
#                  included as a control: it should cost variance but little
#                  bias.
#
# Note that only the OUTCOME model is varied here. The L transition models
# are left correctly specified so the effect is attributable to the hazard
# model alone; see GCOMP_MISSPEC_TRANSITION below to vary them together.
GCOMP_OUTCOME_RHS <- list(
  correct     = "factor(m) + x1 + x2 + L1 + L2 + L1lag + L2lag + trt",
  no_lag      = "factor(m) + x1 + x2 + L1 + L2 + trt",
  no_tv       = "factor(m) + x1 + x2 + trt",
  coarse_L    = "factor(m) + x1 + x2 + L1_hi + L2_hi + L1lag_hi + L2lag_hi + trt",
  nonlinear_L = paste("factor(m) + x1 + x2 + L1 + L2 + I(L1^2) + I(L2^2)",
                      "+ L1lag + L2lag + trt"),
  # Link-misspecification variants. The linear predictor is the CORRECT one;
  # only the link/family is wrong, isolating link misspecification from
  # covariate misspecification.
  #   cloglog -- asymmetric complementary log-log link instead of logit. A
  #              genuinely different hazard shape (not a near-twin of logit
  #              the way probit is), so it is the link most likely to bite.
  #   lpm     -- linear probability model: identity link, Gaussian family.
  #              Fits hazards on the probability scale and can predict values
  #              outside [0,1]; those are clamped at prediction time (see
  #              PROB_EPS), which is part of why the model is misspecified.
  cloglog     = "factor(m) + x1 + x2 + L1 + L2 + L1lag + L2lag + trt",
  lpm         = "factor(m) + x1 + x2 + L1 + L2 + L1lag + L2lag + trt"
)

# Family/link used to FIT each spec. Anything not listed here defaults to
# binomial() (logit), so the five original specs are unchanged.
GCOMP_OUTCOME_FAMILY <- list(
  cloglog = binomial(link = "cloglog"),
  lpm     = gaussian(link = "identity")
)

DEFAULT_GCOMP_MISSPEC <- "correct"

# When TRUE, the L transition models are degraded alongside the outcome model
# (their lag term is dropped for no_lag/no_tv). Left FALSE by default so the
# outcome model is the only thing changing.
GCOMP_MISSPEC_TRANSITION <- FALSE

.gcomp_outcome_formula <- function(misspec = DEFAULT_GCOMP_MISSPEC) {
  if (!misspec %in% names(GCOMP_OUTCOME_RHS))
    stop("Unknown gcomp misspecification '", misspec, "'. Options: ",
         paste(names(GCOMP_OUTCOME_RHS), collapse = ", "))
  stats::as.formula(paste("y ~", GCOMP_OUTCOME_RHS[[misspec]]))
}

# Family for a spec; defaults to logistic when not overridden above.
.gcomp_outcome_family <- function(misspec = DEFAULT_GCOMP_MISSPEC) {
  if (!is.null(GCOMP_OUTCOME_FAMILY[[misspec]])) GCOMP_OUTCOME_FAMILY[[misspec]]
  else binomial()
}