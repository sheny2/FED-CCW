#!/usr/bin/env Rscript
# ArraySim4-7 runner: fragmentation x initiation support.

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
grid <- make_study_grid()
if (is.na(task_id) || task_id < 1L || task_id > nrow(grid))
  stop("task_id must be between 1 and ", nrow(grid), ".")

cell <- grid[task_id, ]
sizes <- get_site_sizes(cell$fragmentation)
init_intercept <- get_init_intercept(cell$initiation)
K <- length(sizes)
beta_init <- scale_confounding(STUDY_CONF_MULT)
beta_init["int"] <- init_intercept

message(sprintf(
  paste0("[task %d/%d] %s: K=%d, n/site=%d, total N=%d; ",
         "initiation=%s (intercept %.1f)"),
  task_id, nrow(grid), cell$fragmentation, K, unique(sizes), sum(sizes),
  cell$initiation, init_intercept
))

# Identical site populations imply a common one-site oracle for every
# fragmentation level. The same seed makes truth exactly comparable.
truth_shared <- compute_truth_ccw_tv(
  N = STUDY_TRUTH_N, tau = STUDY_TAU, t_star = STUDY_TSTAR,
  beta_event = set_beta_trt(STUDY_BETA_TRT), beta_init = beta_init,
  site_sizes = STUDY_TOTAL_N, init_intercepts = init_intercept,
  site_mix = homogeneous_site_mix(1), seed = STUDY_TRUTH_SEED
)
message(sprintf("[task %d] truth: RD=%.4f, RMST difference=%.4f",
                task_id, truth_shared$RD, truth_shared$RMST_diff))

one_rep <- function(r, sizes, init_intercept, truth_shared, task_id) {
  fit_warnings <- character()
  out <- withCallingHandlers(
    try(run_once_tv(
      n_per_site = sizes, init_intercept = init_intercept,
      tau = STUDY_TAU, t_star = STUDY_TSTAR,
      beta_trt = STUDY_BETA_TRT,
      beta_event = set_beta_trt(STUDY_BETA_TRT),
      beta_init = scale_confounding(STUDY_CONF_MULT),
      trunc = STUDY_TRUNC, base_seed = 470000 + task_id * 10000 + r,
      truth = truth_shared, verbose = FALSE
    ), silent = TRUE),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (inherits(out, "try-error")) return(out)
  z <- out$results
  z$rep <- r
  z$fit_warning_count <- length(fit_warnings)
  for (nm in names(out$support)) z[[nm]] <- out$support[[nm]]
  z
}

n_cores <- suppressWarnings(as.integer(Sys.getenv(
  "SLURM_CPUS_PER_TASK", unset = as.character(min(10L, max(1L, detectCores())))
)))
if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L
n_cores <- min(n_cores, STUDY_N_REPS)
message(sprintf("[task %d] running %d reps on %d cores",
                task_id, STUDY_N_REPS, n_cores))

if (.Platform$OS.type == "windows") {
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, "SRC_FILES", envir = environment())
  clusterEvalQ(cl, invisible(lapply(SRC_FILES, source)))
  clusterExport(cl, c("sizes", "init_intercept", "truth_shared", "task_id",
                      "one_rep"), envir = environment())
  ans <- parLapply(cl, seq_len(STUDY_N_REPS), one_rep,
                  sizes, init_intercept, truth_shared, task_id)
} else {
  ans <- mclapply(seq_len(STUDY_N_REPS), one_rep,
                  sizes = sizes, init_intercept = init_intercept,
                  truth_shared = truth_shared, task_id = task_id,
                  mc.cores = n_cores)
}

ok <- !vapply(ans, inherits, logical(1), what = "try-error")
if (any(!ok)) warning(sum(!ok), " replicate(s) failed and were dropped.")
if (!any(ok)) stop("Every replicate failed for task ", task_id, ".")
results <- do.call(rbind, ans[ok])
rownames(results) <- NULL
results$fragmentation <- cell$fragmentation
results$initiation <- cell$initiation
results$init_intercept <- init_intercept
results$tau <- STUDY_TAU
results$t_star <- STUDY_TSTAR
results$total_n <- STUDY_TOTAL_N
results$site_sizes <- paste(sizes, collapse = "/")
results$beta_trt <- STUDY_BETA_TRT
results$conf_mult <- STUDY_CONF_MULT
results$truncation <- paste(STUDY_TRUNC, collapse = "/")
results$task_id <- task_id

dir.create("results", showWarnings = FALSE)
out_file <- file.path("results", sprintf("res_task_%03d.rds", task_id))
saveRDS(results, out_file)
message(sprintf("[task %d] wrote %d rows -> %s",
                task_id, nrow(results), out_file))
