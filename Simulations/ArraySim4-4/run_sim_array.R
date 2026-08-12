#!/usr/bin/env Rscript
# ArraySim4-4 array runner: tau x sample size x initiation prevalence.

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
n_per_site_value <- unname(SAMPLE_SIZE_LEVELS[this$sample_size_scenario])
init_intercept_value <- unname(INITIATION_LEVELS[this$initiation_scenario])

SIM <- list(
  n_reps = STUDY_N_REPS,
  K = STUDY_K,
  n_per_site = setNames(rep(n_per_site_value, STUDY_K), seq_len(STUDY_K)),
  init_intercepts = setNames(rep(init_intercept_value, STUDY_K),
                             seq_len(STUDY_K)),
  tau = this$tau,
  t_star = STUDY_TSTAR,
  beta_trt = STUDY_BETA_TRT,
  scenario = "medium",
  conf_mult = DEFAULT_CONF_MULT,
  beta_init = scale_confounding(DEFAULT_CONF_MULT),
  trunc = c(0, 1),
  base_seed = 2026 + task_id * 10000,
  truth_N = STUDY_TRUTH_N,
  truth_seed = 321,
  sample_size_scenario = this$sample_size_scenario,
  initiation_scenario = this$initiation_scenario,
  n_per_site_value = n_per_site_value,
  total_n = STUDY_K * n_per_site_value,
  init_intercept_value = init_intercept_value,
  zero_cov_init_prob = plogis(init_intercept_value)
)

message(sprintf(
  paste0("[task %d/%d] tau=%g; sample=%s (%d/site, total N=%d); ",
         "initiation=%s (intercept=%.1f, zero-covariate p=%.4f)"),
  task_id, n_combos, SIM$tau, SIM$sample_size_scenario,
  SIM$n_per_site_value, SIM$total_n, SIM$initiation_scenario,
  SIM$init_intercept_value, SIM$zero_cov_init_prob
))
message(sprintf("  beta_trt=%.1f; confounding=medium (x%.1f)",
                SIM$beta_trt, SIM$conf_mult))

truth_shared <- compute_truth_ccw_tv(
  N = SIM$truth_N,
  tau = SIM$tau,
  t_star = SIM$t_star,
  beta_event = set_beta_trt(SIM$beta_trt),
  beta_init = SIM$beta_init,
  site_sizes = SIM$n_per_site,
  init_intercepts = SIM$init_intercepts,
  seed = SIM$truth_seed
)
message(sprintf("[task %d] truth: RD=%.4f RR=%.4f OR=%.4f RMSTd=%.4f",
                task_id, truth_shared$RD, truth_shared$RR,
                truth_shared$OR, truth_shared$RMST_diff))

one_rep <- function(r, SIM, truth_shared) {
  out <- run_once_tv(
    K = SIM$K,
    n_per_site = SIM$n_per_site,
    init_intercepts = SIM$init_intercepts,
    tau = SIM$tau,
    t_star = SIM$t_star,
    beta_trt = SIM$beta_trt,
    beta_init = SIM$beta_init,
    trunc = SIM$trunc,
    base_seed = SIM$base_seed + r,
    truth = truth_shared,
    verbose = FALSE
  )
  res <- out$results
  site_init_rate <- tapply(out$data$A_tau, out$data$site, mean)
  landmark_eligible <- out$data$T_event > SIM$tau
  res$rep <- r
  res$observed_init_rate <- mean(out$data$A_tau)
  res$min_site_init_rate <- min(site_init_rate)
  res$max_site_init_rate <- max(site_init_rate)
  res$landmark_n <- sum(landmark_eligible)
  res$landmark_treated_rate <- if (any(landmark_eligible)) {
    mean(out$data$A_tau[landmark_eligible])
  } else {
    NA_real_
  }
  res
}

n_cores <- as.integer(Sys.getenv(
  "SLURM_CPUS_PER_TASK", unset = as.character(max(1L, detectCores()))
))
n_cores <- max(1L, n_cores)
message(sprintf("[task %d] running %d reps on %d cores",
                task_id, SIM$n_reps, n_cores))

if (.Platform$OS.type == "windows") {
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, "SRC_FILES", envir = environment())
  clusterEvalQ(cl, {
    invisible(lapply(SRC_FILES, source))
    library(stats)
  })
  clusterExport(cl, c("SIM", "truth_shared", "one_rep"),
                envir = environment())
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
results$scenario <- SIM$scenario
results$conf_mult <- SIM$conf_mult
results$sample_size_scenario <- SIM$sample_size_scenario
results$initiation_scenario <- SIM$initiation_scenario
results$n_sites <- SIM$K
results$n_per_site <- SIM$n_per_site_value
results$total_n <- SIM$total_n
results$init_intercept <- SIM$init_intercept_value
results$zero_cov_init_prob <- SIM$zero_cov_init_prob
results$site_sizes <- paste(SIM$n_per_site, collapse = "/")
results$init_intercepts <- paste(SIM$init_intercepts, collapse = "/")
results$task_id <- task_id

dir.create("results", showWarnings = FALSE)
out_file <- file.path("results", sprintf("res_task_%03d.rds", task_id))
saveRDS(results, out_file)
message(sprintf("[task %d] wrote %d rows -> %s",
                task_id, nrow(results), out_file))
