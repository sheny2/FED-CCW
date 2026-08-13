#!/usr/bin/env Rscript
# Aggregate ArraySim4-7 fragmentation x treatment-support results.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
  setwd(dirname(script_path))
}
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
source("params.R")

files <- list.files("results", "^res_task_[0-9]+[.]rds$", full.names = TRUE)
if (!length(files)) stop("No result files found in results/.")
x <- do.call(rbind, lapply(files, readRDS)); rownames(x) <- NULL
message(sprintf("Loaded %d rows from %d task files.", nrow(x), length(files)))

required <- c("method","estimand","estimate","truth","bias","rep","task_id",
  "fragmentation","initiation","init_intercept","n_sites","n_per_site",
  "total_n","empty_g1_at_tstar","zero_initiation_sites",
  "observed_init_rate","observed_death_rate","fit_warning_count")
miss <- setdiff(required, names(x))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
if (!identical(sort(unique(x$task_id)), seq_len(nrow(make_study_grid()))))
  stop("Expected all task IDs 1:", nrow(make_study_grid()), ".")
expected_methods <- c("fed_ccw_tvipcw","pooled_ccw_tvipcw","local_ccw_meta")
if (!identical(sort(unique(x$method)), sort(expected_methods)))
  stop("Unexpected method set.")
counts <- x %>% count(task_id,method,estimand)
if (any(counts$n != STUDY_N_REPS))
  warning("At least one cell has fewer than the planned replicates.")

x <- x %>% mutate(squared_error=bias^2,
                  empty_g1_fraction=empty_g1_at_tstar/n_sites)
summary_all <- x %>%
  group_by(method,estimand,fragmentation,initiation,init_intercept,
           n_sites,n_per_site,total_n,tau,t_star,beta_trt,conf_mult,truncation) %>%
  summarise(mean_est=mean(estimate),median_est=median(estimate),truth=first(truth),
            mean_bias=mean(bias),median_bias=median(bias),emp_sd=sd(estimate),
            rmse=sqrt(mean(squared_error)),
            coverage=if(all(is.na(covered))) NA_real_ else mean(covered,na.rm=TRUE),
            mean_init_rate=mean(observed_init_rate),
            mean_death_rate=mean(observed_death_rate),
            mean_zero_initiation_sites=mean(zero_initiation_sites),
            mean_empty_g1_at_tstar=mean(empty_g1_at_tstar),
            mean_empty_g1_fraction=mean(empty_g1_fraction),
            mean_fit_warning_count=mean(fit_warning_count),
            n_reps=n(),.groups="drop") %>% as.data.frame()
write.csv(summary_all,"summary_all_scenarios.csv",row.names=FALSE)

diagnostics <- x %>% distinct(task_id,rep,fragmentation,initiation,init_intercept,
  n_sites,n_per_site,total_n,observed_init_rate,observed_death_rate,
  zero_initiation_sites,mean_initiators_per_site,min_initiators_per_site,
  max_initiators_per_site,empty_g1_at_tau,empty_g1_at_tstar,
  empty_g0_at_tstar,any_post_tau_empty_g1,fit_warning_count) %>%
  arrange(task_id,rep)
write.csv(diagnostics,"support_diagnostics.csv",row.names=FALSE)

# Paired fed-vs-local errors for the two primary contrasts.
primary <- x %>% filter(method %in% c("fed_ccw_tvipcw","local_ccw_meta"),
                         estimand %in% c("RD","RMST_diff"))
fed <- primary %>% filter(method=="fed_ccw_tvipcw") %>%
  select(task_id,rep,estimand,fragmentation,initiation,n_sites,n_per_site,
         empty_g1_at_tstar,truth,fed_estimate=estimate,fed_bias=bias)
loc <- primary %>% filter(method=="local_ccw_meta") %>%
  select(task_id,rep,estimand,local_estimate=estimate,local_bias=bias)
paired <- inner_join(fed,loc,by=c("task_id","rep","estimand")) %>%
  mutate(local_minus_fed=local_estimate-fed_estimate,
         absolute_bias_gap=abs(local_bias)-abs(fed_bias),
         empty_g1_fraction=empty_g1_at_tstar/n_sites)
write.csv(paired,"fed_vs_local_paired.csv",row.names=FALSE)

method_labels <- c(fed_ccw_tvipcw="Federated CCW",
  pooled_ccw_tvipcw="Pooled CCW",local_ccw_meta="Local CCW + curve meta-analysis")
colors <- c("Federated CCW"="#2c7fb8","Pooled CCW"="#41b6c4",
            "Local CCW + curve meta-analysis"="#d95f0e")
fragment_levels <- names(FRAGMENTATION_LEVELS)
fragment_labels <- c("10 sites × 500","20 sites × 250",
                     "50 sites × 100","100 sites × 50")
plot_data <- x %>% filter(estimand %in% c("RD","RMST_diff")) %>%
  mutate(method=factor(method,levels=names(method_labels),labels=method_labels),
         fragmentation=factor(fragmentation,levels=fragment_levels,
                              labels=fragment_labels),
         initiation=factor(initiation,levels=c("moderate","sparse"),
                           labels=c("Moderate initiation (int. -3.5)",
                                    "Sparse initiation (int. -4.5)")))

for (est in c("RD","RMST_diff")) {
  d <- plot_data %>% filter(estimand==est)
  p <- ggplot(d,aes(fragmentation,bias,fill=method))+
    geom_boxplot(outlier.size=.35,position=position_dodge(.8),linewidth=.25)+
    geom_hline(yintercept=0,linetype="dashed")+facet_wrap(~initiation,nrow=1)+
    scale_fill_manual(values=colors)+
    labs(title=paste("Bias of",est,"under site fragmentation"),
      subtitle="Fixed total N=5,000; homogeneous sites; tau=4; beta=-0.7; no treatment confounding",
      x="Fragmentation (number of sites × patients/site)",
      y="Estimate - truth",fill="Method")+theme_bw(base_size=11)+
    theme(legend.position="bottom",axis.text.x=element_text(angle=20,hjust=1),
          panel.grid.minor=element_blank())
  ggsave(paste0("bias_boxplot_",est,".png"),p,width=13,height=7,dpi=150)
}

scatter <- paired %>%
  select(task_id,rep,estimand,fragmentation,initiation,n_sites,
         empty_g1_fraction,fed_bias,local_bias) %>%
  tidyr::pivot_longer(c(fed_bias,local_bias),names_to="method",values_to="bias") %>%
  mutate(method=factor(method,levels=c("fed_bias","local_bias"),
                       labels=c("Federated CCW","Local CCW + curve meta-analysis")),
         initiation=factor(initiation,levels=c("moderate","sparse")))
p <- ggplot(scatter,aes(empty_g1_fraction,bias,color=method))+
  geom_point(alpha=.35,size=1)+geom_smooth(method="lm",se=FALSE,linewidth=.8)+
  geom_hline(yintercept=0,linetype="dashed")+
  facet_grid(estimand~initiation,scales="free_y")+scale_color_manual(values=colors)+
  labs(title="Bias versus loss of local treated-strategy support",
       x="Fraction of sites with zero treated-strategy risk set at t*",
       y="Estimate - truth",color="Method")+theme_bw(base_size=11)+
  theme(legend.position="bottom",panel.grid.minor=element_blank())
ggsave("bias_vs_empty_treated_support.png",p,width=12,height=8,dpi=150)
message("Wrote summaries, diagnostics, paired comparison, and plots.")
