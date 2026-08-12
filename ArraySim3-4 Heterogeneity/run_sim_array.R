#!/usr/bin/env Rscript
# ============================================================================
# Array simulation over (tau, treatment effect, site heterogeneity).
#
# Heterogeneity changes the baseline patient mix across three sites through
# site-specific x1 means/SDs and x2 prevalences.  See DEFAULT_SITE_MIXES in
# params.R.  The pooled marginal means remain fixed, allowing cleaner
# comparisons as the sites become less exchangeable.
#
# Each SLURM array task handles one grid cell and runs all replicates in
# parallel across the cores allocated to it.
#
# Usage:
#   Rscript run_sim_array.R <task_id>
# ============================================================================

rm(list = ls())
suppressPackageStartupMessages(library(parallel))

SRC_FILES <- c("params.R", "DGP_tv.R", "Fed_CCW_TVIPCW.R", "Simulation.R")
invisible(lapply(SRC_FILES, source))

# ---- 0. array task -------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
task_id <- if (length(args) >= 1L) {
  as.integer(args[[1L]])
} else {
  as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
}
if (is.na(task_id) || task_id < 1L) stop("Invalid task_id: ", task_id)

# ---- 1. parameter grid: 3 x 3 x 3 = 27 tasks ----------------------------
grid <- expand.grid(
  tau           = c(2, 5, 8),
  beta_trt      = c(-1.0, -0.7, -0.5),
  heterogeneity = names(DEFAULT_SITE_MIXES),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
n_combos <- nrow(grid)
if (task_id > n_combos)
  stop(sprintf("task_id %d exceeds number of combinations (%d).",
               task_id, n_combos))
this <- grid[task_id, ]

# ---- 2. settings shared by the task's replicates ------------------------
SIM <- list(
  n_reps       = 100,
  K            = 3,
  n_per_site   = 1000,
  tau          = this$tau,
  t_star       = 25,
  beta_trt     = this$beta_trt,
  heterogeneity = this$heterogeneity,
  trunc        = c(0, 1),
  mc_gcomp     = 20,
  base_seed    = 2026 + task_id * 1000,
  truth_N      = 3e6,
  truth_seed   = 321
)

mix <- get_site_mix(SIM$heterogeneity, SIM$K)
message(sprintf("[task %d/%d] tau=%g beta_trt=%g heterogeneity=%s",
                task_id, n_combos, SIM$tau, SIM$beta_trt,
                SIM$heterogeneity))
message("  site mix: ",
        paste(sprintf("site%d[x1~N(%.2f,%.2f^2), x2~Bern(%.2f)]",
                      mix$site, mix$x1_mean, mix$x1_sd, mix$x2_prob),
              collapse = "; "))

# ---- 3. heterogeneous-population oracle, once per grid cell ------------
truth_shared <- compute_truth_multisite_tv(
  N             = SIM$truth_N,
  K             = SIM$K,
  heterogeneity = SIM$heterogeneity,
  tau           = SIM$tau,
  t_star        = SIM$t_star,
  beta_event    = set_beta_trt(SIM$beta_trt),
  seed          = SIM$truth_seed
)
message(sprintf("[task %d] truth: RD=%.4f RR=%.4f OR=%.4f RMSTd=%.4f",
                task_id, truth_shared$RD, truth_shared$RR,
                truth_shared$OR, truth_shared$RMST_diff))

# ---- 4. one replicate ----------------------------------------------------
one_rep <- function(r, SIM, truth_shared) {
  out <- run_once_tv(
    K             = SIM$K,
    n_per_site    = SIM$n_per_site,
    heterogeneity = SIM$heterogeneity,
    tau           = SIM$tau,
    t_star        = SIM$t_star,
    beta_trt      = SIM$beta_trt,
    trunc         = SIM$trunc,
    mc_gcomp      = SIM$mc_gcomp,
    base_seed     = SIM$base_seed + r,
    truth         = truth_shared,
    verbose       = FALSE
  )
  res <- out$results
  res$rep <- r
  res
}

# ---- 5. parallel replicates ---------------------------------------------
n_cores <- as.integer(Sys.getenv(
  "SLURM_CPUS_PER_TASK",
  unset = as.character(max(1L, detectCores()))
))
n_cores <- max(1L, min(n_cores, SIM$n_reps))
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
  res_list <- mclapply(
    seq_len(SIM$n_reps), one_rep,
    SIM = SIM, truth_shared = truth_shared,
    mc.cores = n_cores
  )
}

ok <- !vapply(res_list, inherits, logical(1), what = "try-error")
if (any(!ok))
  warning(sprintf("[task %d] %d replicate(s) failed and were dropped.",
                  task_id, sum(!ok)))
if (!any(ok)) stop("All replicates failed for task ", task_id, ".")

results <- do.call(rbind, res_list[ok])
rownames(results) <- NULL

# ---- 6. tag and save -----------------------------------------------------
results$tau           <- SIM$tau
results$t_star        <- SIM$t_star
results$beta_trt      <- SIM$beta_trt
results$heterogeneity <- SIM$heterogeneity
results$task_id       <- task_id

dir.create("results", showWarnings = FALSE, recursive = TRUE)
out_file <- file.path("results", sprintf("res_task_%03d.rds", task_id))
saveRDS(results, out_file)
message(sprintf("[task %d] wrote %d rows -> %s",
                task_id, nrow(results), out_file))
