# =============================================================
# params.R -- single source of truth for the heterogeneous-site DGP.
# Every other file pulls its defaults from here.
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

# ---- between-site patient-mix heterogeneity -----------------------------
# Each level has three equally sized sites.  The conditional treatment,
# covariate-transition, and event models stay common across sites; only the
# baseline patient mix changes.  This isolates the impact of non-IID sites
# while preserving the pooled marginal means E[x1] = 0 and E[x2] = 0.4.
#
# Separation increases through:
#   * wider differences in the site-specific mean of continuous x1;
#   * increasingly different x1 variances; and
#   * wider differences in the prevalence of binary x2.
#
# "low" is intentionally not perfectly homogeneous, so every grid cell is a
# multisite-heterogeneity experiment.
DEFAULT_SITE_MIXES <- list(
  low = data.frame(
    site    = 1:3,
    x1_mean = c(-0.25, 0.00, 0.25),
    x1_sd   = c( 0.95, 1.00, 1.05),
    x2_prob = c( 0.35, 0.40, 0.45)
  ),
  moderate = data.frame(
    site    = 1:3,
    x1_mean = c(-0.75, 0.00, 0.75),
    x1_sd   = c( 0.85, 1.00, 1.15),
    x2_prob = c( 0.20, 0.40, 0.60)
  ),
  high = data.frame(
    site    = 1:3,
    x1_mean = c(-1.50, 0.00, 1.50),
    x1_sd   = c( 0.70, 1.00, 1.30),
    x2_prob = c( 0.10, 0.40, 0.70)
  )
)

get_site_mix <- function(level, K = 3L) {
  if (!level %in% names(DEFAULT_SITE_MIXES))
    stop("Unknown heterogeneity level: ", level)
  mix <- DEFAULT_SITE_MIXES[[level]]
  if (K != nrow(mix))
    stop(sprintf("Heterogeneity profiles require K=%d sites; received K=%d.",
                 nrow(mix), K))
  mix
}
