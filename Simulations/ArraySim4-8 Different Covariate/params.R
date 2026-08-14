# =============================================================
# ArraySim4-8: different covariate availability across sites.
# =============================================================

DEFAULT_TAU <- 4L
DEFAULT_TSTAR <- 25L
DEFAULT_BETA_EVENT <- c(
  int=-2.5,x1=.5,x2=.4,L1=.4,L2=.3,L1lag=.2,L2lag=.15,trt=-.7
)
DEFAULT_BETA_INIT <- c(int=-3,x1=.8,x2=-.3,L1=.5,L2=.4)
DEFAULT_BETA_L <- c(int=0,ar=.6,x1=.3,x2=.2)
DEFAULT_SD_L <- 1
DEFAULT_TRUNC <- c(0,1)
PROB_EPS <- 1e-6

.clamp_prob <- function(p,eps=PROB_EPS) pmin(pmax(p,eps),1-eps)
set_beta_trt <- function(beta_trt,beta_event=DEFAULT_BETA_EVENT) {
  beta_event["trt"]<-beta_trt; beta_event
}
CONFOUNDING_SLOPES <- c("x1","x2","L1","L2")
scale_confounding <- function(mult,beta_init=DEFAULT_BETA_INIT) {
  beta_init[CONFOUNDING_SLOPES]<-beta_init[CONFOUNDING_SLOPES]*mult
  beta_init
}

STUDY_K <- 3L
STUDY_TAU <- 4L
STUDY_TSTAR <- 25L
STUDY_BETA_TRT <- -.7
STUDY_CONF_MULT <- 1
STUDY_INIT_INTERCEPT <- -3
STUDY_N_REPS <- 300L
STUDY_TRUTH_N <- 4e6
STUDY_TRUTH_SEED <- 48001L
STUDY_TRUNC <- c(0,1)

SAMPLE_SIZE_LEVELS <- c(small=100L,medium=500L,large=1000L)
WEIGHT_COVARIATES <- c("x1","x2","L1","L2","L1lag","L2lag")

homogeneous_site_mix <- function(K=STUDY_K) data.frame(
  site=seq_len(K),x1_mean=rep(0,K),x1_sd=rep(1,K),x2_prob=rep(.4,K)
)
validate_site_mix <- function(site_mix,K=STUDY_K) {
  needed<-c("site","x1_mean","x1_sd","x2_prob")
  if(!is.data.frame(site_mix)||!all(needed%in%names(site_mix)))
    stop("site_mix has invalid columns.")
  if(nrow(site_mix)!=K||!identical(as.integer(site_mix$site),seq_len(K)))
    stop("site_mix must have one ordered row per site.")
  if(any(site_mix$x1_sd<=0)||any(site_mix$x2_prob<=0|site_mix$x2_prob>=1))
    stop("Invalid site-mix parameters.")
  invisible(TRUE)
}

.site_specific_dgp <- function() {
  base<-scale_confounding(STUDY_CONF_MULT)
  base["int"]<-STUDY_INIT_INTERCEPT
  site1<-base; site1[c("x2","L2")]<-0
  site2<-base; site2[c("x1","L1")]<-0
  list(site1,site2,base)
}

full_covariates <- rep(list(WEIGHT_COVARIATES),STUDY_K)

COVARIATE_SCENARIOS <- list(
  common_full=list(
    description="All sites observe and use the full sufficient history.",
    beta_init_by_site=.site_specific_dgp(),
    analysis_covariates=full_covariates
  ),
  different_sufficient=list(
    description="Different local variable sets; every local set remains sufficient.",
    beta_init_by_site=.site_specific_dgp(),
    analysis_covariates=list(
      c("x1","L1","L1lag"),
      c("x2","L2","L2lag"),
      WEIGHT_COVARIATES
    )
  ),
  different_insufficient=list(
    description="Different local sets omit active time-varying confounders.",
    beta_init_by_site=.site_specific_dgp(),
    analysis_covariates=list(c("x1"),c("x2"),c("x1","x2"))
  )
)

validate_covariate_config <- function(config,K=STUDY_K) {
  if(length(config$beta_init_by_site)!=K||length(config$analysis_covariates)!=K)
    stop("Covariate configuration must have one entry per site.")
  ok_beta<-vapply(config$beta_init_by_site,function(x)
    all(names(DEFAULT_BETA_INIT)%in%names(x)),logical(1))
  ok_cov<-vapply(config$analysis_covariates,function(x)
    all(x%in%WEIGHT_COVARIATES),logical(1))
  if(!all(ok_beta)||!all(ok_cov)) stop("Invalid covariate configuration.")
  invisible(TRUE)
}

get_covariate_config <- function(level) {
  if(!level%in%names(COVARIATE_SCENARIOS)) stop("Unknown covariate scenario: ",level)
  z<-COVARIATE_SCENARIOS[[level]]; validate_covariate_config(z); z
}

common_covariates <- function(site_covariates) {
  if(!length(site_covariates)) return(character())
  Reduce(intersect,site_covariates)
}

make_study_grid <- function() expand.grid(
  sample_size=names(SAMPLE_SIZE_LEVELS),
  covariate_scenario=names(COVARIATE_SCENARIOS),
  KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE
)

STUDY_SITE_SIZES <- rep(SAMPLE_SIZE_LEVELS[["medium"]],STUDY_K)
STUDY_INIT_INTERCEPTS <- rep(STUDY_INIT_INTERCEPT,STUDY_K)
