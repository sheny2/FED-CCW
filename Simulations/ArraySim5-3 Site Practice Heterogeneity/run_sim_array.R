#!/usr/bin/env Rscript
# ============================================================================
# ArraySim5-3 Site Practice Heterogeneity array-job simulation over the
# (tau, beta_trt, site-practice heterogeneity) grid.
#
# Patient mix, covariate effects, and the outcome model are fixed. Across the
# low, moderate and high scenarios, only the site-specific intercepts in the
# treatment-initiation DGP change. The oracle is recomputed for every cell.
#
# t_star is fixed; tau, beta_trt and site-practice heterogeneity vary.
#
# Each SLURM array task handles one grid cell and runs all replicates for it
# in parallel across the cores it was given.
#
# Usage (invoked by job.sh):
#   Rscript run_sim_array.R <task_id>
# where <task_id> is 1..(#combos). If omitted, falls back to
# SLURM_ARRAY_TASK_ID, then to 1.
# ============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(parallel)
})

SRC_FILES <- c("params.R", "DGP_tv.R", "Fed_CCW_TVIPCW.R", "Simulation.R")
invisible(lapply(SRC_FILES, source))
validate_params()

# ---- 0. resolve which array task this is -------------------------------
args <- commandArgs(trailingOnly = TRUE)
task_id <- if (length(args) >= 1) {
  as.integer(args[[1]])
} else {
  as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
}
if (is.na(task_id) || task_id < 1) stop("Invalid task_id: ", task_id)

# ---- 1. parameter grid --------------------------------------------------
# Practice scenarios come from SITE_INIT_INTERCEPTS in params.R.
grid <- simulation_grid()
n_combos <- nrow(grid)

if (task_id > n_combos)
  stop(sprintf("task_id %d exceeds number of combos (%d).", task_id, n_combos))

this <- grid[task_id, ]

# ---- 2. settings shared across all replicates --------------------------
SIM <- list(
  n_reps     = SIM_N_REPS,
  K          = length(DEFAULT_SITE_SIZES),
  n_per_site = DEFAULT_SITE_SIZES,
  init_intercepts = get_init_intercepts(this$practice_heterogeneity),
  tau        = this$tau,
  t_star     = SIM_TSTAR,
  beta_trt   = this$beta_trt,
  practice_heterogeneity = this$practice_heterogeneity,
  patient_mix = DEFAULT_PATIENT_MIX,
  beta_init  = DEFAULT_BETA_INIT,
  trunc      = DEFAULT_TRUNC,
  base_seed  = SIM_BASE_SEED + task_id * SIM_TASK_SEED_STRIDE,
  truth_N    = SIM_TRUTH_N,
  truth_seed = SIM_TRUTH_SEED
)

message(sprintf("[task %d/%d] tau=%g  beta_trt=%g  site-practice heterogeneity=%s",
                task_id, n_combos, SIM$tau, SIM$beta_trt,
                SIM$practice_heterogeneity))
message("  beta_init = ",
        paste(sprintf("%s=%.3f", names(SIM$beta_init), SIM$beta_init),
              collapse = ", "))
message("  site sizes = ", paste(SIM$n_per_site, collapse = "/"),
        "; initiation intercepts = ",
        paste(SIM$init_intercepts, collapse = "/"))
message("  patient mix is homogeneous: X1~N(0,1), P(X2=1)=0.4")

# ---- 3. compute the oracle once per grid cell --------------------------
truth_shared <- compute_truth_ccw_tv(
  N          = SIM$truth_N,
  tau        = SIM$tau,
  t_star     = SIM$t_star,
  beta_event = set_beta_trt(SIM$beta_trt),
  beta_init  = SIM$beta_init,
  site_sizes = SIM$n_per_site,
  init_intercepts = SIM$init_intercepts,
  patient_mix = SIM$patient_mix,
  seed       = SIM$truth_seed
)
message(sprintf("[task %d] truth: RD=%.4f RR=%.4f OR=%.4f RMSTd=%.4f",
                task_id, truth_shared$RD, truth_shared$RR,
                truth_shared$OR, truth_shared$RMST_diff))

# ---- 4. one replicate ---------------------------------------------------
one_rep <- function(r, SIM, truth_shared) {
  out <- run_once_tv(
    K          = SIM$K,
    n_per_site = SIM$n_per_site,
    init_intercepts = SIM$init_intercepts,
    patient_mix = SIM$patient_mix,
    tau        = SIM$tau,
    t_star     = SIM$t_star,
    beta_trt   = SIM$beta_trt,
    beta_init  = SIM$beta_init,
    trunc      = SIM$trunc,
    base_seed  = SIM$base_seed + r,
    truth      = truth_shared,
    verbose    = FALSE
  )
  res <- out$results
  res$rep <- r
  res
}

# ---- 5. run replicates in parallel across this task's cores ------------
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
                                 unset = as.character(max(1, detectCores()))))
n_cores <- max(1, n_cores)
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
  clusterExport(cl, c("SIM", "truth_shared", "one_rep"), envir = environment())
  res_list <- parLapply(cl, seq_len(SIM$n_reps), one_rep,
                        SIM = SIM, truth_shared = truth_shared)
} else {
  res_list <- mclapply(seq_len(SIM$n_reps), one_rep,
                       SIM = SIM, truth_shared = truth_shared,
                       mc.cores = n_cores)
}

# Guard against silent worker failures.
ok <- !vapply(res_list, inherits, logical(1), what = "try-error")
if (any(!ok))
  warning(sprintf("[task %d] %d replicate(s) failed and were dropped.",
                  task_id, sum(!ok)))
results <- do.call(rbind, res_list[ok])
rownames(results) <- NULL

# ---- 6. tag with the grid cell and save --------------------------------
results$tau       <- SIM$tau
results$t_star    <- SIM$t_star
results$beta_trt  <- SIM$beta_trt
results$practice_heterogeneity <- SIM$practice_heterogeneity
results$site_sizes <- paste(SIM$n_per_site, collapse = "/")
results$init_intercepts <- paste(SIM$init_intercepts, collapse = "/")
results$task_id   <- task_id

dir.create("results", showWarnings = FALSE)
out_file <- file.path("results", sprintf("res_task_%03d.rds", task_id))
saveRDS(results, out_file)
message(sprintf("[task %d] wrote %d rows -> %s",
                task_id, nrow(results), out_file))
