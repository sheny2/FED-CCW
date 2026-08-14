#!/usr/bin/env Rscript
# Preflight checks for ArraySim4-8.
rm(list=ls())
file_arg<-grep("^--file=",commandArgs(FALSE),value=TRUE)
if(length(file_arg)){script_path<-normalizePath(sub("^--file=","",file_arg[[1L]]));setwd(dirname(script_path))}
source("Simulation.R")
g<-make_study_grid()
stopifnot(nrow(g)==9L,all(table(g$sample_size)==3L),all(table(g$covariate_scenario)==3L),
 identical(unname(SAMPLE_SIZE_LEVELS),c(100L,500L,1000L)),STUDY_K==3L,
 STUDY_TAU==4L,STUDY_BETA_TRT==-.7,STUDY_CONF_MULT==1,
 STUDY_INIT_INTERCEPT==-3,STUDY_N_REPS==100L)
message("PASS: nine-cell sample-size x covariate-availability grid is complete.")

full<-get_covariate_config("common_full")
suf<-get_covariate_config("different_sufficient")
ins<-get_covariate_config("different_insufficient")
stopifnot(all(vapply(full$analysis_covariates,identical,logical(1),WEIGHT_COVARIATES)),
 identical(suf$analysis_covariates[[1]],c("x1","L1","L1lag")),
 identical(suf$analysis_covariates[[2]],c("x2","L2","L2lag")),
 length(common_covariates(suf$analysis_covariates))==0L,
 identical(full$beta_init_by_site,suf$beta_init_by_site),
 identical(suf$beta_init_by_site,ins$beta_init_by_site),
 suf$beta_init_by_site[[1]]["L1"]!=0,suf$beta_init_by_site[[1]]["L2"]==0,
 suf$beta_init_by_site[[2]]["L1"]==0,suf$beta_init_by_site[[2]]["L2"]!=0,
 !"L1"%in%ins$analysis_covariates[[1]],!"L2"%in%ins$analysis_covariates[[2]])
message("PASS: all analyses share one DGP; sufficient and insufficient sets are correctly encoded.")

make_dat<-function(scenario,seed){z<-get_covariate_config(scenario);simulate_multisite_tv(
 K=3,n_per_site=rep(100L,3),init_intercepts=rep(-3,3),
 beta_init_by_site=z$beta_init_by_site,site_mix=homogeneous_site_mix(3),
 tau=4,t_star=25,base_seed=seed,beta_event=set_beta_trt(-.7))}
d<-make_dat("different_sufficient",12345)
stopifnot(nrow(d)==300L,!any(d$S<=25&d$S>d$T_event))
L1<-as.matrix(d[,paste0("L1_",1:25)]);after<-outer(d$T_event,1:25,`<`)
stopifnot(all(is.na(L1[after])))
message("PASS: site-specific DGP retains the corrected event-free risk process.")

truth_s<-compute_truth_ccw_tv(N=50000,tau=4,t_star=25,
 beta_event=set_beta_trt(-.7),beta_init_by_site=suf$beta_init_by_site,
 site_sizes=rep(100,3),init_intercepts=rep(-3,3),
 site_mix=homogeneous_site_mix(3),seed=888)
truth_i<-compute_truth_ccw_tv(N=50000,tau=4,t_star=25,
 beta_event=set_beta_trt(-.7),beta_init_by_site=ins$beta_init_by_site,
 site_sizes=rep(100,3),init_intercepts=rep(-3,3),
 site_mix=homogeneous_site_mix(3),seed=888)
truth_f<-compute_truth_ccw_tv(N=50000,tau=4,t_star=25,
 beta_event=set_beta_trt(-.7),beta_init_by_site=full$beta_init_by_site,
 site_sizes=rep(100,3),init_intercepts=rep(-3,3),
 site_mix=homogeneous_site_mix(3),seed=888)
stopifnot(identical(truth_s[c("RD","RMST_diff")],truth_i[c("RD","RMST_diff")]),
          identical(truth_s[c("RD","RMST_diff")],truth_f[c("RD","RMST_diff")]))
message("PASS: all covariate-availability analyses have the exact same DGP truth.")

o<-suppressWarnings(run_once_tv(rep(100L,3),"different_sufficient",
 truth=truth_s,base_seed=54321,verbose=FALSE))
expected<-c("fed_ccw_site_specific","pooled_ccw_common_set","fed_ipw_no_clone",
 "fed_perprotocol_naive","fed_landmark_site_specific","local_ccw_meta_site_specific")
stopifnot(identical(sort(unique(o$results$method)),sort(expected)),
 nrow(o$results)==6L*8L,all(is.finite(o$results$estimate)),
 identical(o$common_covariates,character()))
message("PASS: end-to-end heterogeneous-covariate run returns six finite methods.")

of<-suppressWarnings(run_once_tv(rep(100L,3),"common_full",
 truth_N=30000,base_seed=111,verbose=FALSE))
legacy<-suppressWarnings(run_fed_ccw_tvipcw(of$data,4,25,STUDY_TRUNC))
stopifnot(abs(legacy$RD-of$fits$fed_ccw_site_specific$RD)<1e-12)
message("PASS: site-specific implementation reproduces the legacy full-set fit.")
message("All ArraySim4-8 preflight checks passed.")
