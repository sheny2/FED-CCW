#!/usr/bin/env Rscript
# Preflight checks for ArraySim4-7.

rm(list=ls())
file_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(length(file_arg)) {
  script_path <- normalizePath(sub("^--file=","",file_arg[[1L]]))
  setwd(dirname(script_path))
}
source("Simulation.R")

grid <- make_study_grid()
stopifnot(
  nrow(grid)==8L,
  identical(unique(grid$fragmentation),names(FRAGMENTATION_LEVELS)),
  identical(sort(unique(grid$initiation)),sort(names(INITIATION_LEVELS))),
  all(vapply(FRAGMENTATION_LEVELS,sum,numeric(1))==STUDY_TOTAL_N),
  identical(unname(vapply(FRAGMENTATION_LEVELS,length,integer(1))),c(10L,20L,50L,100L)),
  identical(unname(vapply(FRAGMENTATION_LEVELS,function(x)unique(x),integer(1))),
            c(500L,250L,100L,50L)),
  STUDY_TAU==4L,STUDY_BETA_TRT==-0.7,STUDY_CONF_MULT==0,
  identical(STUDY_TRUNC,c(.01,.99))
)
message("PASS: eight-cell fragmentation x initiation-support grid is complete.")

b <- scale_confounding(STUDY_CONF_MULT)
stopifnot(all(b[CONFOUNDING_SLOPES]==0))
message("PASS: initiation is independent of X/L, isolating aggregation and support.")

make_data <- function(level,init_level,seed) {
  sizes <- get_site_sizes(level); int <- get_init_intercept(init_level)
  bi <- scale_confounding(0); bi["int"] <- int
  simulate_multisite_tv(K=length(sizes),n_per_site=sizes,
    init_intercepts=rep(int,length(sizes)),site_mix=homogeneous_site_mix(length(sizes)),
    tau=STUDY_TAU,t_star=STUDY_TSTAR,base_seed=seed,
    beta_event=set_beta_trt(STUDY_BETA_TRT),beta_init=bi)
}
control <- make_data("10_sites_x_500","sparse",13579)
stress <- make_data("100_sites_x_50","sparse",24680)
stopifnot(nrow(control)==5000L,nrow(stress)==5000L,
  nlevels(factor(control$site))==10L,nlevels(factor(stress$site))==100L,
  abs(mean(control$A_tau)-mean(stress$A_tau))<.02)
message(sprintf("PASS: control and stress cells both have N=5,000 and initiation rates %.3f/%.3f.",
                mean(control$A_tau),mean(stress$A_tau)))

check_event_free <- function(dat) {
  stopifnot(!any(dat$S<=STUDY_TSTAR & dat$S>dat$T_event))
  L1<-as.matrix(dat[,paste0("L1_",seq_len(STUDY_TSTAR)),drop=FALSE])
  after<-outer(dat$T_event,seq_len(STUDY_TSTAR),`<`)
  stopifnot(all(is.na(L1[after])))
}
check_event_free(control); check_event_free(stress)
message("PASS: corrected event-free DGP invariant holds.")

support_count <- function(dat) {
  ds<-split(dat,dat$site)
  hn<-central_common_num_hazard(lapply(ds,local_initiation_counts,t_star=STUDY_TSTAR))
  st<-suppressWarnings(lapply(ds,local_ccw_tvipcw,tau=STUDY_TAU,
    t_star=STUDY_TSTAR,trunc=STUDY_TRUNC,hnum=hn))
  c(zero_init=sum(vapply(ds,function(d)sum(d$S<=STUDY_TAU)==0,logical(1))),
    empty_end=sum(vapply(st,function(s)s$rw1[STUDY_TSTAR]<=0,logical(1))))
}
control_support<-support_count(control); stress_support<-support_count(stress)
stopifnot(stress_support[["empty_end"]]>20,
          stress_support[["empty_end"]]>control_support[["empty_end"]])
message(sprintf("PASS: stress cell has %d empty endpoint treated risk sets versus %d in control.",
  stress_support[["empty_end"]],control_support[["empty_end"]]))

truth <- compute_truth_ccw_tv(N=50000,tau=STUDY_TAU,t_star=STUDY_TSTAR,
  beta_event=set_beta_trt(STUDY_BETA_TRT),
  beta_init={z<-scale_confounding(0);z["int"]<-get_init_intercept("sparse");z},
  site_sizes=STUDY_TOTAL_N,init_intercepts=get_init_intercept("sparse"),
  site_mix=homogeneous_site_mix(1),seed=97531)
out <- suppressWarnings(run_once_tv(n_per_site=get_site_sizes("10_sites_x_500"),
  init_intercept=get_init_intercept("sparse"),base_seed=86420,truth=truth,verbose=FALSE))
stopifnot(
  identical(sort(unique(out$results$method)),
    sort(c("fed_ccw_tvipcw","pooled_ccw_tvipcw","local_ccw_meta"))),
  nrow(out$results)==3L*8L,all(is.finite(out$results$estimate))
)

ind_fed <- suppressWarnings(run_fed_ccw_tvipcw(out$data,STUDY_TAU,STUDY_TSTAR,STUDY_TRUNC))
ind_local <- suppressWarnings(run_local_ccw_meta(out$data,STUDY_TAU,STUDY_TSTAR,STUDY_TRUNC))
stopifnot(abs(ind_fed$RD-out$fits$fed_ccw_tvipcw$RD)<1e-12,
          abs(ind_local$RD-out$fits$local_ccw_meta$RD)<1e-12)
message("PASS: joint implementation exactly reproduces independent fed and local fits.")

pilot<-read.csv("pilot_search_summary.csv")
extreme<-pilot[pilot$K==100 & pilot$n_per_site==50 & pilot$init_intercept==-4.5,]
stopifnot(nrow(extreme)==1L,abs(extreme$fed_RMST_bias)<.10,
          abs(extreme$local_RMST_bias)>1,
          extreme$mean_empty_treated_sites>50)
message("PASS: recorded pilot identifies the intended high-separation stress cell.")
message("All ArraySim4-7 preflight checks passed.")
