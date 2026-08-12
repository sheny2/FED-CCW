#!/usr/bin/env Rscript
# ============================================================================
# Array-job simulation over the (tau, beta_trt, confounding strength, misspec)
# grid, for the trimmed TWO-METHOD study: fed-CCW and CC-gcomp only.
#
# The `misspec` axis (correct / no_tv / coarse_L) degrades the time-varying
# covariates in BOTH nuisance models in the same cell: fed-CCW's weight
# (initiation-hazard) denominator and g-comp's outcome hazard. So each cell
# compares the two methods under a common nuisance error.
#
# Confounding strength is set by scaling beta_init (see scale_confounding()
# in params.R). The truth depends on beta_init, so the oracle is recomputed
# per DGP cell. It does NOT depend on misspec, so seeds are keyed on the DGP
# cell and every misspec variant sees identical datasets and oracle.
#
# t_star is FIXED across the grid.
#
# Grid: 3 tau x 3 beta_trt x 3 scenario x 3 misspec = 81 tasks.
#
# Usage (invoked by job.sh):
#   Rscript run_sim_array.R <task_id>    (task_id in 1..81)
# falls back to SLURM_ARRAY_TASK_ID, then 1.
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
grid <- expand.grid(
  tau      = c(2, 5, 8),
  beta_trt = c(-1.0, -0.7, -0.5),
  scenario = names(DEFAULT_CONF_MULTS),
  misspec  = GCOMP_MISSPEC_SPECS,
  KEEP.OUT.ATTRS   = FALSE,
  stringsAsFactors = FALSE
)
grid$conf_mult <- unname(DEFAULT_CONF_MULTS[grid$scenario])
n_combos <- nrow(grid)

message(sprintf("[grid] %d combos = %d tau x %d beta x %d scenario x %d misspec; array must be 1-%d",
                n_combos, length(unique(grid$tau)), length(unique(grid$beta_trt)),
                length(unique(grid$scenario)), length(unique(grid$misspec)),
                n_combos))

if (task_id > n_combos)
  stop(sprintf("task_id %d exceeds number of combos (%d).", task_id, n_combos))

this <- grid[task_id, ]

# Index of the DGP cell, ignoring misspec. Tasks sharing a DGP cell share
# their data and oracle, so all misspec variants are compared on identical
# datasets.
dgp_cells <- unique(grid[, c("tau", "beta_trt", "scenario")])
dgp_cell_id <- which(
  dgp_cells$tau      == this$tau &
  dgp_cells$beta_trt == this$beta_trt &
  dgp_cells$scenario == this$scenario
)

# ---- 2. settings shared across all replicates --------------------------
SIM <- list(
  n_reps     = 100,
  K          = 3,
  n_per_site = 1000,
  tau        = this$tau,
  t_star     = 25,
  beta_trt   = this$beta_trt,
  scenario   = this$scenario,
  conf_mult  = this$conf_mult,
  beta_init  = scale_confounding(this$conf_mult),
  misspec    = this$misspec,
  trunc      = c(0, 1),
  base_seed  = 2026 + dgp_cell_id * 1000,   # keyed on DGP cell, not misspec
  truth_N    = 3e6,
  truth_seed = 321
)

message(sprintf("[task %d/%d] tau=%g  beta_trt=%g  confounding=%s (x%g)  misspec=%s",
                task_id, n_combos, SIM$tau, SIM$beta_trt,
                SIM$scenario, SIM$conf_mult, SIM$misspec))
message("  beta_init = ",
        paste(sprintf("%s=%.3f", names(SIM$beta_init), SIM$beta_init),
              collapse = ", "))

# ---- 3. compute the oracle once per DGP cell ---------------------------
truth_shared <- compute_truth_ccw_tv(
  N          = SIM$truth_N,
  tau        = SIM$tau,
  t_star     = SIM$t_star,
  beta_event = set_beta_trt(SIM$beta_trt),
  beta_init  = SIM$beta_init,
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
    tau        = SIM$tau,
    t_star     = SIM$t_star,
    beta_trt   = SIM$beta_trt,
    beta_init  = SIM$beta_init,
    trunc      = SIM$trunc,
    misspec    = SIM$misspec,
    base_seed  = SIM$base_seed + r,
    truth      = truth_shared,
    verbose    = FALSE
  )
  res <- out$results
  res$rep <- r
  res
}

# ---- 5. run replicates in parallel -------------------------------------
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
results$misspec   <- SIM$misspec
results$task_id   <- task_id

dir.create("results", showWarnings = FALSE)
out_file <- file.path("results", sprintf("res_task_%03d.rds", task_id))
saveRDS(results, out_file)
message(sprintf("[task %d] wrote %d rows -> %s",
                task_id, nrow(results), out_file))