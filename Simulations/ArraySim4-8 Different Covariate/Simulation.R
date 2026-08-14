# =============================================================
# ArraySim4-8 one-replicate wrapper.
# =============================================================
source("params.R")
source("DGP_tv.R")
source("Fed_CCW_TVIPCW.R")

run_once_tv <- function(n_per_site=STUDY_SITE_SIZES,
                        covariate_scenario="common_full",
                        tau=STUDY_TAU,t_star=STUDY_TSTAR,
                        beta_trt=STUDY_BETA_TRT,
                        beta_event=DEFAULT_BETA_EVENT,
                        beta_L=DEFAULT_BETA_L,sd_L=DEFAULT_SD_L,
                        trunc=STUDY_TRUNC,base_seed=2048,
                        truth_N=STUDY_TRUTH_N,
                        truth_seed=STUDY_TRUTH_SEED,
                        truth=NULL,verbose=TRUE) {
  if(length(n_per_site)==1L) n_per_site<-rep(n_per_site,STUDY_K)
  if(length(n_per_site)!=STUDY_K) stop("n_per_site must have length 1 or 3.")
  config<-get_covariate_config(covariate_scenario)
  beta_event<-set_beta_trt(beta_trt,beta_event)
  common_covs<-common_covariates(config$analysis_covariates)

  if(is.null(truth)) truth<-compute_truth_ccw_tv(
    N=truth_N,tau=tau,t_star=t_star,beta_event=beta_event,
    beta_init=DEFAULT_BETA_INIT,beta_init_by_site=config$beta_init_by_site,
    site_sizes=n_per_site,
    init_intercepts=rep(STUDY_INIT_INTERCEPT,STUDY_K),
    site_mix=homogeneous_site_mix(STUDY_K),beta_L=beta_L,sd_L=sd_L,
    seed=truth_seed)

  dat<-simulate_multisite_tv(
    K=STUDY_K,n_per_site=n_per_site,
    init_intercepts=rep(STUDY_INIT_INTERCEPT,STUDY_K),
    beta_init_by_site=config$beta_init_by_site,
    site_mix=homogeneous_site_mix(STUDY_K),tau=tau,t_star=t_star,
    base_seed=base_seed,beta_event=beta_event,beta_init=DEFAULT_BETA_INIT,
    beta_L=beta_L,sd_L=sd_L)

  fits<-list(
    fed_ccw_site_specific=run_fed_ccw_tvipcw(
      dat,tau,t_star,trunc,site_covariates=config$analysis_covariates),
    pooled_ccw_common_set=run_pooled_ccw_tvipcw(
      dat,tau,t_star,trunc,covariates=common_covs),
    fed_ipw_no_clone=run_fed_ipw_nocloning(
      dat,tau,t_star,trunc,site_covariates=config$analysis_covariates),
    fed_perprotocol_naive=run_fed_perprotocol(dat,tau,t_star),
    fed_landmark_site_specific=run_fed_landmark_ipw(
      dat,tau,t_star,trunc,site_covariates=config$analysis_covariates),
    local_ccw_meta_site_specific=run_local_ccw_meta(
      dat,tau,t_star,trunc,site_covariates=config$analysis_covariates)
  )

  estimands<-c("psi1","psi0","RD","RR","OR","RMST1","RMST0","RMST_diff")
  se_name<-c(RD="SE_RD",RR="SE_logRR",OR="SE_logOR",
             RMST_diff="SE_RMSTdiff",RMST1="SE_RMST1",RMST0="SE_RMST0")
  ci_name<-c(RD="CI_RD",RR="CI_RR",OR="CI_OR",RMST_diff="CI_RMSTdiff")
  make_rows<-function(fit,method) do.call(rbind,lapply(estimands,function(est){
    se<-fit[[se_name[est]]];if(is.null(se))se<-NA_real_
    ci<-fit[[ci_name[est]]];if(is.null(ci))ci<-c(NA_real_,NA_real_)
    data.frame(method=method,estimand=est,estimate=fit[[est]],truth=truth[[est]],
      bias=fit[[est]]-truth[[est]],se=se,ci_lo=ci[1],ci_hi=ci[2],
      covered=if(is.na(ci[1]))NA_integer_ else
        as.integer(truth[[est]]>=ci[1]&truth[[est]]<=ci[2]),
      stringsAsFactors=FALSE)
  }))
  results<-do.call(rbind,Map(make_rows,fits,names(fits)));rownames(results)<-NULL
  if(verbose) {
    cat("Scenario:",covariate_scenario,"; common pooled set:",
        if(length(common_covs))paste(common_covs,collapse="+") else "intercept only","\n")
    print(results[results$estimand%in%c("RD","RMST_diff"),
      c("method","estimand","estimate","truth","bias")],digits=4,row.names=FALSE)
  }
  invisible(list(results=results,truth=truth,fits=fits,data=dat,
    config=config,common_covariates=common_covs))
}
