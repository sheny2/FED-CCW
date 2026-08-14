#!/usr/bin/env Rscript
# Aggregate ArraySim4-8 results.
file_arg<-grep("^--file=",commandArgs(FALSE),value=TRUE)
if(length(file_arg)){script_path<-normalizePath(sub("^--file=","",file_arg[[1L]]));setwd(dirname(script_path))}
suppressPackageStartupMessages({library(dplyr);library(ggplot2)})
source("params.R")
files<-list.files("results","^res_task_[0-9]+[.]rds$",full.names=TRUE)
if(!length(files))stop("No result files in results/.")
x<-do.call(rbind,lapply(files,readRDS));rownames(x)<-NULL
message(sprintf("Loaded %d rows from %d task files.",nrow(x),length(files)))
req<-c("method","estimand","estimate","truth","bias","rep","task_id",
 "sample_size","covariate_scenario","n_per_site","total_n",
 "site_covariates","pooled_common_covariates","fit_warning_count")
miss<-setdiff(req,names(x));if(length(miss))stop("Missing: ",paste(miss,collapse=", "))
if(!identical(sort(unique(x$task_id)),seq_len(nrow(make_study_grid()))))stop("Not all nine tasks are present.")
methods<-c("fed_ccw_site_specific","pooled_ccw_common_set","fed_ipw_no_clone",
 "fed_perprotocol_naive","fed_landmark_site_specific","local_ccw_meta_site_specific")
if(!identical(sort(unique(x$method)),sort(methods)))stop("Unexpected method set.")
counts<-x%>%count(task_id,method,estimand)
if(any(counts$n!=STUDY_N_REPS))warning("Unexpected replicate count in at least one cell.")
x<-x%>%mutate(squared_error=bias^2)
summ<-x%>%group_by(method,estimand,sample_size,covariate_scenario,n_per_site,
 total_n,n_sites,tau,t_star,beta_trt,conf_mult,init_intercept,
 site_covariates,pooled_common_covariates)%>%summarise(
 mean_est=mean(estimate),median_est=median(estimate),truth=first(truth),
 mean_bias=mean(bias),median_bias=median(bias),emp_sd=sd(estimate),
 rmse=sqrt(mean(squared_error)),
 coverage=if(all(is.na(covered)))NA_real_ else mean(covered,na.rm=TRUE),
 mean_init_rate=mean(observed_init_rate),mean_death_rate=mean(observed_death_rate),
 mean_fit_warning_count=mean(fit_warning_count),n_reps=n(),.groups="drop")%>%as.data.frame()
write.csv(summ,"summary_all_scenarios.csv",row.names=FALSE)

diag<-x%>%distinct(task_id,rep,sample_size,covariate_scenario,n_per_site,total_n,
 site_covariates,pooled_common_covariates,observed_init_rate,min_site_init_rate,
 max_site_init_rate,observed_death_rate,fit_warning_count)%>%arrange(task_id,rep)
write.csv(diag,"design_diagnostics.csv",row.names=FALSE)

fed_boundary<-summ%>%filter(method=="fed_ccw_site_specific",
 estimand%in%c("RD","RMST_diff"))%>%
 select(estimand,sample_size,covariate_scenario,n_per_site,truth,mean_est,
        mean_bias,emp_sd,rmse,coverage,n_reps)
write.csv(fed_boundary,"fed_ccw_boundary_summary.csv",row.names=FALSE)

method_labels<-c(fed_ccw_site_specific="Federated CCW (site-specific sets)",
 pooled_ccw_common_set="Pooled CCW (common set)",
 fed_ipw_no_clone="Federated IPW (no cloning)",
 fed_perprotocol_naive="Federated per-protocol",
 fed_landmark_site_specific="Federated landmark IPW",
 local_ccw_meta_site_specific="Local CCW + curve meta-analysis")
colors<-c("Federated CCW (site-specific sets)"="#2c7fb8",
 "Pooled CCW (common set)"="#41b6c4","Federated IPW (no cloning)"="#f5a623",
 "Federated per-protocol"="#e8482c","Federated landmark IPW"="#7a5195",
 "Local CCW + curve meta-analysis"="#2ca25f")
scenario_labels<-c(common_full="Common full set",
 different_sufficient="Different sets, locally sufficient",
 different_insufficient="Different sets, confounders missing")
pd<-x%>%filter(estimand%in%c("RD","RMST_diff"))%>%mutate(
 method=factor(method,levels=names(method_labels),labels=method_labels),
 sample_size=factor(sample_size,levels=names(SAMPLE_SIZE_LEVELS),
  labels=paste0(names(SAMPLE_SIZE_LEVELS)," (n=",SAMPLE_SIZE_LEVELS,"/site)")),
 covariate_scenario=factor(covariate_scenario,levels=names(scenario_labels),
  labels=scenario_labels))
for(est in c("RD","RMST_diff")){
 d<-pd%>%filter(estimand==est) %>% filter(method != "Pooled CCW (common set)")
 p<-ggplot(d,aes(sample_size,bias,fill=method))+geom_boxplot(outlier.size=.35,
  position=position_dodge(.8),linewidth=.25)+geom_hline(yintercept=0,linetype="dashed")+
  facet_wrap(~covariate_scenario,nrow=1)+scale_fill_manual(values=colors)+
  labs(title=paste("Bias of",est,"with different covariates across sites"),
   subtitle="3 sites | tau=4 | beta=-0.7 | medium confounding | initiation intercept=-3",
   x="Sample size",y="Estimate - truth",fill="Method")+theme_bw(base_size=11)+
  theme(legend.position="bottom",axis.text.x=element_text(angle=20,hjust=1),panel.grid.minor=element_blank())
 ggsave(paste0("bias_boxplot_",est,".png"),p,width=15,height=7,dpi=150)
}
message("Wrote summaries, diagnostics, and plots.")
