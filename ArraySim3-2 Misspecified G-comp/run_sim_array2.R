#!/usr/bin/env Rscript
# ============================================================================
# Array-job simulation over the (tau, beta_trt, confounding strength) grid.
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
# A fourth axis, `misspec`, varies the outcome-hazard specification used by
# the two g-computation methods only. The DGP, the data and the oracle are
# unchanged by it, so seeds are keyed on the DGP cell and every misspec
# variant sees identical datasets.
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
#
# GCOMP_MISSPEC_GRID selects which outcome-model specifications to run for
# the two g-computation methods (see GCOMP_OUTCOME_RHS in params.R). It
# multiplies the number of tasks. The IPCW and per-protocol methods do not
# use an outcome model and are unaffected by this axis -- their results are
# re-estimated in each misspec cell and should agree up to Monte Carlo noise,
# which doubles as a sanity check.
#
# APPEND RUN: this grid holds only the two link-misspecification specs, to be
# added on top of an existing completed run of the other five. Because
# dgp_cell_id (below) is computed from the tau/beta/scenario crossing only --
# which is invariant to the misspec set -- the seeds here match the earlier
# run exactly, so cloglog/lpm see the SAME datasets as the completed specs
# and are directly comparable. Output filenames are shifted by
# TASK_ID_OFFSET so they do not overwrite the existing res_task_*.rds.
GCOMP_MISSPEC_GRID <- c("cloglog", "lpm")

# Number of tasks already on disk from the previous run (5 specs x 27 DGP
# cells). New files are written as res_task_(OFFSET + task_id).rds.
TASK_ID_OFFSET <- 135L

grid <- expand.grid(
  tau      = c(2, 5, 8),
  beta_trt = c(-1.0, -0.7, -0.5),
  scenario = names(DEFAULT_CONF_MULTS),
  misspec  = GCOMP_MISSPEC_GRID,
  KEEP.OUT.ATTRS   = FALSE,
  stringsAsFactors = FALSE
)
grid$conf_mult <- unname(DEFAULT_CONF_MULTS[grid$scenario])
n_combos <- nrow(grid)

# The SLURM array for THIS append run must cover 1..n_combos (the offset is
# applied to filenames, not to the array index). Make the requirement loud.
message(sprintf("[grid] %d combos = %d tau x %d beta x %d scenario x %d misspec; array must be 1-%d (files offset by %d)",
                n_combos, length(unique(grid$tau)), length(unique(grid$beta_trt)),
                length(unique(grid$scenario)), length(unique(grid$misspec)),
                n_combos, TASK_ID_OFFSET))

if (task_id > n_combos)
  stop(sprintf("task_id %d exceeds number of combos (%d).", task_id, n_combos))

this <- grid[task_id, ]

# Index of the underlying DGP cell, ignoring the misspec axis. Tasks that
# share a DGP cell share their data and their oracle.
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
  t_star     = 25,                 # fixed follow-up horizon
  beta_trt   = this$beta_trt,
  scenario   = this$scenario,
  conf_mult  = this$conf_mult,
  beta_init  = scale_confounding(this$conf_mult),
  misspec    = this$misspec,
  trunc      = c(0, 1),
  # Seeds key on the DGP cell (tau, beta_trt, scenario) and NOT on misspec,
  # so every misspecification variant sees the SAME simulated datasets and
  # the same oracle. Differences between misspec cells are then attributable
  # to the outcome model alone rather than to sampling noise.
  base_seed  = 2026 + dgp_cell_id * 1000,
  truth_N    = 3e6,
  truth_seed = 321
)

message(sprintf("[task %d/%d] tau=%g  beta_trt=%g  confounding=%s (x%g)  gcomp=%s",
                task_id, n_combos, SIM$tau, SIM$beta_trt,
                SIM$scenario, SIM$conf_mult, SIM$misspec))
message("  beta_init = ",
        paste(sprintf("%s=%.3f", names(SIM$beta_init), SIM$beta_init),
              collapse = ", "))

# ---- 3. compute the oracle once per grid cell --------------------------
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
    gcomp_misspec = SIM$misspec,
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
# The stored task_id and the filename carry the offset so this append run's
# outputs sit after the existing files and never collide with them.
file_task_id <- TASK_ID_OFFSET + task_id

results$tau       <- SIM$tau
results$t_star    <- SIM$t_star
results$beta_trt  <- SIM$beta_trt
results$scenario  <- SIM$scenario
results$conf_mult <- SIM$conf_mult
results$misspec   <- SIM$misspec
results$task_id   <- file_task_id

dir.create("results", showWarnings = FALSE)
out_file <- file.path("results", sprintf("res_task_%03d.rds", file_task_id))
saveRDS(results, out_file)
message(sprintf("[task %d -> file %d] wrote %d rows -> %s",
                task_id, file_task_id, nrow(results), out_file))