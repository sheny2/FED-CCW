#!/usr/bin/env Rscript
# ArraySim4-5 array runner: tau x sample size x outcome frequency.

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
event_intercept_value <- unname(OUTCOME_LEVELS[this$outcome_scenario])
beta_event_value <- set_beta_trt(
  STUDY_BETA_TRT,
  set_event_intercept(event_intercept_value)
)

SIM <- list(
  n_reps = STUDY_N_REPS,
  K = STUDY_K,
  n_per_site = setNames(rep(n_per_site_value, STUDY_K), seq_len(STUDY_K)),
  init_intercepts = setNames(rep(STUDY_INIT_INTERCEPT, STUDY_K),
                             seq_len(STUDY_K)),
  tau = this$tau,
  t_star = STUDY_TSTAR,
  beta_trt = STUDY_BETA_TRT,
  beta_event = beta_event_value,
  beta_init = scale_confounding(STUDY_CONF_MULT),
  trunc = c(0, 1),
  base_seed = 2027 + task_id * 1000,
  truth_N = STUDY_TRUTH_N,
  truth_seed = 456,
  sample_size_scenario = this$sample_size_scenario,
  outcome_scenario = this$outcome_scenario,
  n_per_site_value = n_per_site_value,
  total_n = STUDY_K * n_per_site_value,
  event_intercept = event_intercept_value
)

message(sprintf(
  paste0("[task %d/%d] tau=%g; sample=%s (%d/site, total N=%d); ",
         "outcome=%s (event intercept=%.1f)"),
  task_id, n_combos, SIM$tau, SIM$sample_size_scenario,
  SIM$n_per_site_value, SIM$total_n, SIM$outcome_scenario,
  SIM$event_intercept
))
message(sprintf(
  "  beta_trt=%.1f; confounding=medium (x%.1f); initiation intercept=%.1f",
  SIM$beta_trt, STUDY_CONF_MULT, STUDY_INIT_INTERCEPT
))

truth_shared <- compute_truth_ccw_tv(
  N = SIM$truth_N,
  tau = SIM$tau,
  t_star = SIM$t_star,
  beta_event = SIM$beta_event,
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
    beta_event = SIM$beta_event,
    beta_init = SIM$beta_init,
    trunc = SIM$trunc,
    base_seed = SIM$base_seed + r,
    truth = truth_shared,
    verbose = FALSE
  )
  res <- out$results
  site_death_rate <- tapply(out$data$delta, out$data$site, mean)
  site_init_rate <- tapply(out$data$A_tau, out$data$site, mean)
  landmark_eligible <- out$data$T_event > SIM$tau
  landmark_event <- landmark_eligible & out$data$T_event <= SIM$t_star

  res$rep <- r
  res$observed_death_rate <- mean(out$data$delta)
  res$min_site_death_rate <- min(site_death_rate)
  res$max_site_death_rate <- max(site_death_rate)
  res$observed_init_rate <- mean(out$data$A_tau)
  res$min_site_init_rate <- min(site_init_rate)
  res$max_site_init_rate <- max(site_init_rate)
  res$landmark_n <- sum(landmark_eligible)
  res$landmark_event_n <- sum(landmark_event)
  res$landmark_event_rate <- if (any(landmark_eligible)) {
    mean(out$data$T_event[landmark_eligible] <= SIM$t_star)
  } else {
    NA_real_
  }
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
results$scenario <- "medium"
results$conf_mult <- STUDY_CONF_MULT
results$sample_size_scenario <- SIM$sample_size_scenario
results$outcome_scenario <- SIM$outcome_scenario
results$n_sites <- SIM$K
results$n_per_site <- SIM$n_per_site_value
results$total_n <- SIM$total_n
results$event_intercept <- SIM$event_intercept
results$init_intercept <- STUDY_INIT_INTERCEPT
results$site_sizes <- paste(SIM$n_per_site, collapse = "/")
results$init_intercepts <- paste(SIM$init_intercepts, collapse = "/")
results$task_id <- task_id

dir.create("results", showWarnings = FALSE)
out_file <- file.path("results", sprintf("res_task_%03d.rds", task_id))
saveRDS(results, out_file)
message(sprintf("[task %d] wrote %d rows -> %s",
                task_id, nrow(results), out_file))
