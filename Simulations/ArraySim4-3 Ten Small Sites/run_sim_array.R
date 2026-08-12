#!/usr/bin/env Rscript
# ============================================================================
# ArraySim4-3 Ten Small Sites array-job simulation over the
# (tau, beta_trt, confounding strength) grid.
#
# Confounding is governed by beta_init, the covariate coefficients in the
# initiation model. Because those covariates (x1, x2 at baseline; L1, L2
# time-varying) also drive the event hazard, scaling the beta_init slopes
# scales the strength of confounding. See scale_confounding() in params.R
# for the exact construction and for why the intercept is held fixed.
#
# The truth depends on beta_init -- changing it changes who initiates and
# when -- so the oracle is recomputed for every grid cell, not just for
# every beta_trt.
#
# t_star is FIXED across the grid; tau, beta_trt and confounding vary.
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

# ---- 0. resolve which array task this is -------------------------------
args <- commandArgs(trailingOnly = TRUE)
task_id <- if (length(args) >= 1) {
  as.integer(args[[1]])
} else {
  as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
}
if (is.na(task_id) || task_id < 1) stop("Invalid task_id: ", task_id)

# ---- 1. parameter grid --------------------------------------------------
# Confounding scenarios come from DEFAULT_CONF_MULTS (params.R); the grid
# carries the scenario NAME so the multiplier and its label cannot drift
# apart.
grid <- expand.grid(
  tau      = c(3, 6, 9),
  beta_trt = c(-1.0, -0.7, -0.5),
  scenario = names(DEFAULT_CONF_MULTS),
  KEEP.OUT.ATTRS   = FALSE,
  stringsAsFactors = FALSE
)
grid$conf_mult <- unname(DEFAULT_CONF_MULTS[grid$scenario])
n_combos <- nrow(grid)

if (task_id > n_combos)
  stop(sprintf("task_id %d exceeds number of combos (%d).", task_id, n_combos))

this <- grid[task_id, ]

# ---- 2. settings shared across all replicates --------------------------
SIM <- list(
  n_reps     = 500,
  K          = 10,
  n_per_site = TEN_SITE_SIZES,
  init_intercepts = TEN_SITE_INIT_INTERCEPTS,
  tau        = this$tau,
  t_star     = 25,                 # fixed follow-up horizon
  beta_trt   = this$beta_trt,
  scenario   = this$scenario,
  conf_mult  = this$conf_mult,
  beta_init  = scale_confounding(this$conf_mult),
  trunc      = c(0, 1),
  base_seed  = 2026 + task_id * 10000,  # disjoint replicate seeds per task
  truth_N    = 3e6,
  truth_seed = 321
)

message(sprintf("[task %d/%d] tau=%g  beta_trt=%g  confounding=%s (x%g)",
                task_id, n_combos, SIM$tau, SIM$beta_trt,
                SIM$scenario, SIM$conf_mult))
message("  beta_init = ",
        paste(sprintf("%s=%.3f", names(SIM$beta_init), SIM$beta_init),
              collapse = ", "))
message("  site sizes = ", paste(SIM$n_per_site, collapse = "/"),
        "; initiation intercepts = ",
        paste(SIM$init_intercepts, collapse = "/"))

# ---- 3. compute the oracle once per grid cell --------------------------
truth_shared <- compute_truth_ccw_tv(
  N          = SIM$truth_N,
  tau        = SIM$tau,
  t_star     = SIM$t_star,
  beta_event = set_beta_trt(SIM$beta_trt),
  beta_init  = SIM$beta_init,
  site_sizes = SIM$n_per_site,
  init_intercepts = SIM$init_intercepts,
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
results$scenario  <- SIM$scenario
results$conf_mult <- SIM$conf_mult
results$site_sizes <- paste(SIM$n_per_site, collapse = "/")
results$init_intercepts <- paste(SIM$init_intercepts, collapse = "/")
results$task_id   <- task_id

dir.create("results", showWarnings = FALSE)
out_file <- file.path("results", sprintf("res_task_%03d.rds", task_id))
saveRDS(results, out_file)
message(sprintf("[task %d] wrote %d rows -> %s",
                task_id, nrow(results), out_file))
