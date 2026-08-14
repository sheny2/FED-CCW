#!/usr/bin/env Rscript
# ArraySim4-8 runner: sample size x covariate-availability scenario.
rm(list=ls())
file_arg<-grep("^--file=",commandArgs(FALSE),value=TRUE)
if(length(file_arg)){script_path<-normalizePath(sub("^--file=","",file_arg[[1L]]));setwd(dirname(script_path))}
suppressPackageStartupMessages(library(parallel))
SRC_FILES<-c("params.R","DGP_tv.R","Fed_CCW_TVIPCW.R","Simulation.R")
invisible(lapply(SRC_FILES,source))
args<-commandArgs(trailingOnly=TRUE)
task_id<-if(length(args))as.integer(args[[1L]]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID",unset="1"))
grid<-make_study_grid()
if(is.na(task_id)||task_id<1L||task_id>nrow(grid))stop("task_id must be 1-",nrow(grid),".")
cell<-grid[task_id,];n_site<-unname(SAMPLE_SIZE_LEVELS[[cell$sample_size]])
sizes<-rep(n_site,STUDY_K);config<-get_covariate_config(cell$covariate_scenario)
message(sprintf("[task %d/%d] n/site=%d; covariates=%s",task_id,nrow(grid),n_site,cell$covariate_scenario))

truth_shared<-compute_truth_ccw_tv(N=STUDY_TRUTH_N,tau=STUDY_TAU,t_star=STUDY_TSTAR,
  beta_event=set_beta_trt(STUDY_BETA_TRT),beta_init=DEFAULT_BETA_INIT,
  beta_init_by_site=config$beta_init_by_site,site_sizes=sizes,
  init_intercepts=rep(STUDY_INIT_INTERCEPT,STUDY_K),
  site_mix=homogeneous_site_mix(STUDY_K),seed=STUDY_TRUTH_SEED)
message(sprintf("[task %d] truth RD=%.4f; RMSTd=%.4f",task_id,truth_shared$RD,truth_shared$RMST_diff))

one_rep<-function(r,sizes,scenario,truth_shared,task_id){
  warns<-character()
  out<-withCallingHandlers(try(run_once_tv(n_per_site=sizes,
    covariate_scenario=scenario,base_seed=480000+task_id*10000+r,
    truth=truth_shared,verbose=FALSE),silent=TRUE),warning=function(w){
      warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")})
  if(inherits(out,"try-error"))return(out)
  z<-out$results;d<-out$data
  z$rep<-r;z$fit_warning_count<-length(warns)
  z$observed_init_rate<-mean(d$A_tau);z$observed_death_rate<-mean(d$delta)
  z$min_site_init_rate<-min(tapply(d$A_tau,d$site,mean))
  z$max_site_init_rate<-max(tapply(d$A_tau,d$site,mean))
  z
}
n_cores<-suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",unset=as.character(min(6L,max(1L,detectCores()))))))
if(!is.finite(n_cores)||n_cores<1L)n_cores<-1L;n_cores<-min(n_cores,STUDY_N_REPS)
message(sprintf("[task %d] running %d reps on %d cores",task_id,STUDY_N_REPS,n_cores))
if(.Platform$OS.type=="windows"){
  cl<-makeCluster(n_cores);on.exit(stopCluster(cl),add=TRUE)
  clusterExport(cl,"SRC_FILES",envir=environment());clusterEvalQ(cl,invisible(lapply(SRC_FILES,source)))
  clusterExport(cl,c("sizes","cell","truth_shared","task_id","one_rep"),envir=environment())
  ans<-parLapply(cl,seq_len(STUDY_N_REPS),one_rep,sizes,cell$covariate_scenario,truth_shared,task_id)
}else ans<-mclapply(seq_len(STUDY_N_REPS),one_rep,sizes=sizes,
  scenario=cell$covariate_scenario,truth_shared=truth_shared,task_id=task_id,mc.cores=n_cores)
ok<-!vapply(ans,inherits,logical(1),what="try-error")
if(any(!ok))warning(sum(!ok)," replicate(s) failed and were dropped.")
if(!any(ok))stop("All replicates failed.")
results<-do.call(rbind,ans[ok]);rownames(results)<-NULL
results$sample_size<-cell$sample_size;results$covariate_scenario<-cell$covariate_scenario
results$n_per_site<-n_site;results$total_n<-sum(sizes);results$n_sites<-STUDY_K
results$tau<-STUDY_TAU;results$t_star<-STUDY_TSTAR;results$beta_trt<-STUDY_BETA_TRT
results$conf_mult<-STUDY_CONF_MULT;results$init_intercept<-STUDY_INIT_INTERCEPT
results$site_covariates<-paste(vapply(config$analysis_covariates,function(x)paste(x,collapse="+"),character(1)),collapse=" | ")
cc<-common_covariates(config$analysis_covariates)
results$pooled_common_covariates<-if(length(cc))paste(cc,collapse="+") else "intercept_only"
results$task_id<-task_id
dir.create("results",showWarnings=FALSE)
f<-file.path("results",sprintf("res_task_%03d.rds",task_id));saveRDS(results,f)
message(sprintf("[task %d] wrote %d rows -> %s",task_id,nrow(results),f))
