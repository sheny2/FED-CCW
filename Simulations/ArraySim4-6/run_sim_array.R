#!/usr/bin/env Rscript
# ArraySim4-6 array runner: tau x site-size balance x patient-mix heterogeneity.

rm(list = ls())
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
  setwd(dirname(script_path))
}

suppressPackageStartupMessages(library(parallel))
SRC_FILES <- c("params.R", "DGP_tv.R", "Fed_CCW_TVIPCW.R", "Simulation.R")
invisible(lapply(SRC_FILES, source))

args <- commandArgs(trailingOnly = TRUE)
task_id <- if (length(args)) as.integer(args[[1L]]) else
  as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
if (is.na(task_id) || task_id < 1L) stop("Invalid task_id: ", task_id)

grid <- make_study_grid()
n_combos <- nrow(grid)
if (task_id > n_combos)
  stop(sprintf("task_id %d exceeds number of combinations (%d).",
               task_id, n_combos))

this <- grid[task_id, ]
site_sizes <- get_site_sizes(this$size_scenario)
site_mix <- get_site_mix(this$mix_scenario)
site_weights <- site_sizes / sum(site_sizes)

SIM <- list(
  n_reps = STUDY_N_REPS,
  K = STUDY_K,
  n_per_site = site_sizes,
  init_intercepts = setNames(rep(STUDY_INIT_INTERCEPT, STUDY_K),
                             seq_len(STUDY_K)),
  site_mix = site_mix,
  tau = this$tau,
  t_star = STUDY_TSTAR,
  beta_trt = STUDY_BETA_TRT,
  beta_event = set_beta_trt(STUDY_BETA_TRT),
  beta_init = scale_confounding(STUDY_CONF_MULT),
  trunc = STUDY_TRUNC,
  base_seed = 2028 + task_id * 1000,
  truth_N = STUDY_TRUTH_N,
  truth_seed = 654,
  size_scenario = this$size_scenario,
  mix_scenario = this$mix_scenario,
  total_n = sum(site_sizes),
  target_x1_mean = sum(site_weights * site_mix$x1_mean),
  target_x2_prob = sum(site_weights * site_mix$x2_prob)
)

message(sprintf(
  paste0("[task %d/%d] tau=%g; size=%s (%s; total N=%d); ",
         "patient mix=%s"),
  task_id, n_combos, SIM$tau, SIM$size_scenario,
  paste(SIM$n_per_site, collapse = "/"), SIM$total_n, SIM$mix_scenario
))
message(sprintf(
  paste0("  target E[X1]=%.3f; target P(X2=1)=%.3f; beta_trt=%.1f; ",
         "medium confounding; initiation intercept=%.1f"),
  SIM$target_x1_mean, SIM$target_x2_prob, SIM$beta_trt,
  STUDY_INIT_INTERCEPT
))

truth_shared <- compute_truth_ccw_tv(
  N = SIM$truth_N,
  tau = SIM$tau,
  t_star = SIM$t_star,
  beta_event = SIM$beta_event,
  beta_init = SIM$beta_init,
  site_sizes = SIM$n_per_site,
  init_intercepts = SIM$init_intercepts,
  site_mix = SIM$site_mix,
  seed = SIM$truth_seed
)
message(sprintf("[task %d] truth: RD=%.4f RR=%.4f OR=%.4f RMSTd=%.4f",
                task_id, truth_shared$RD, truth_shared$RR,
                truth_shared$OR, truth_shared$RMST_diff))

one_rep <- function(r, SIM, truth_shared) {
  out <- try(run_once_tv(
    K = SIM$K,
    n_per_site = SIM$n_per_site,
    init_intercepts = SIM$init_intercepts,
    site_mix = SIM$site_mix,
    tau = SIM$tau,
    t_star = SIM$t_star,
    beta_trt = SIM$beta_trt,
    beta_event = SIM$beta_event,
    beta_init = SIM$beta_init,
    trunc = SIM$trunc,
    base_seed = SIM$base_seed + r,
    truth = truth_shared,
    verbose = FALSE
  ), silent = TRUE)
  if (inherits(out, "try-error")) return(out)

  res <- out$results
  dat <- out$data
  site_death <- tapply(dat$delta, dat$site, mean)
  site_init <- tapply(dat$A_tau, dat$site, mean)
  site_x1_mean <- tapply(dat$x1, dat$site, mean)
  site_x1_sd <- tapply(dat$x1, dat$site, sd)
  site_x2_prob <- tapply(dat$x2, dat$site, mean)
  landmark_eligible <- dat$T_event > SIM$tau

  res$rep <- r
  res$observed_x1_mean <- mean(dat$x1)
  res$observed_x1_sd <- sd(dat$x1)
  res$observed_x2_prob <- mean(dat$x2)
  for (k in seq_len(SIM$K)) {
    res[[paste0("site", k, "_x1_mean")]] <- site_x1_mean[k]
    res[[paste0("site", k, "_x1_sd")]] <- site_x1_sd[k]
    res[[paste0("site", k, "_x2_prob")]] <- site_x2_prob[k]
  }
  res$observed_death_rate <- mean(dat$delta)
  res$min_site_death_rate <- min(site_death)
  res$max_site_death_rate <- max(site_death)
  res$observed_init_rate <- mean(dat$A_tau)
  res$min_site_init_rate <- min(site_init)
  res$max_site_init_rate <- max(site_init)
  res$landmark_n <- sum(landmark_eligible)
  res$landmark_event_n <- sum(landmark_eligible & dat$T_event <= SIM$t_star)
  res$landmark_event_rate <- if (any(landmark_eligible)) {
    mean(dat$T_event[landmark_eligible] <= SIM$t_star)
  } else NA_real_
  res$landmark_treated_rate <- if (any(landmark_eligible)) {
    mean(dat$A_tau[landmark_eligible])
  } else NA_real_
  res
}

n_cores <- suppressWarnings(as.integer(Sys.getenv(
  "SLURM_CPUS_PER_TASK", unset = as.character(min(10L, max(1L, detectCores())))
)))
if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L
n_cores <- min(n_cores, SIM$n_reps)
message(sprintf("[task %d] running %d reps on %d cores",
                task_id, SIM$n_reps, n_cores))

if (.Platform$OS.type == "windows") {
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, "SRC_FILES", envir = environment())
  clusterEvalQ(cl, invisible(lapply(SRC_FILES, source)))
  clusterExport(cl, c("SIM", "truth_shared", "one_rep"), envir = environment())
  res_list <- parLapply(cl, seq_len(SIM$n_reps), one_rep,
                       SIM = SIM, truth_shared = truth_shared)
} else {
  res_list <- mclapply(seq_len(SIM$n_reps), one_rep,
                       SIM = SIM, truth_shared = truth_shared,
                       mc.cores = n_cores)
}

ok <- !vapply(res_list, inherits, logical(1), what = "try-error")
if (any(!ok))
  warning(sprintf("[task %d] %d replicate(s) failed and were dropped.",
                  task_id, sum(!ok)))
if (!any(ok)) stop("Every replicate failed for task ", task_id, ".")

results <- do.call(rbind, res_list[ok])
rownames(results) <- NULL
results$tau <- SIM$tau
results$t_star <- SIM$t_star
results$beta_trt <- SIM$beta_trt
results$scenario <- "medium"
results$conf_mult <- STUDY_CONF_MULT
results$size_scenario <- SIM$size_scenario
results$mix_scenario <- SIM$mix_scenario
results$n_sites <- SIM$K
results$total_n <- SIM$total_n
results$site_sizes <- paste(SIM$n_per_site, collapse = "/")
results$site_weights <- paste(round(SIM$n_per_site / SIM$total_n, 6),
                              collapse = "/")
results$init_intercept <- STUDY_INIT_INTERCEPT
results$target_x1_mean <- SIM$target_x1_mean
results$target_x2_prob <- SIM$target_x2_prob
results$site_x1_means <- paste(SIM$site_mix$x1_mean, collapse = "/")
results$site_x1_sds <- paste(SIM$site_mix$x1_sd, collapse = "/")
results$site_x2_probs <- paste(SIM$site_mix$x2_prob, collapse = "/")
results$task_id <- task_id

dir.create("results", showWarnings = FALSE)
out_file <- file.path("results", sprintf("res_task_%03d.rds", task_id))
saveRDS(results, out_file)
message(sprintf("[task %d] wrote %d rows -> %s",
                task_id, nrow(results), out_file))
