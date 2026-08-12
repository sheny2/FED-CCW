#!/usr/bin/env Rscript
# ============================================================================
# Combine per-task result files and produce grid-faceted bias boxplots.
#
# For each estimand, one figure:
#   rows    = confounding strength (small / medium / strong)
#   columns = true effect beta_trt
#   x-axis  = grace period tau
#   fill    = method (dodged within each tau)
#
# t_star is fixed by the runner, so it is reported in the subtitle rather
# than used as a facet.
#
# Run AFTER the array job finishes:  Rscript aggregate.R
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

source("params.R")   # DEFAULT_CONF_MULTS, for scenario ordering

files <- list.files("results", pattern = "^res_task_.*\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No result files found in results/.")

results <- do.call(rbind, lapply(files, readRDS))
rownames(results) <- NULL
message(sprintf("Loaded %d rows from %d task files.", nrow(results), length(files)))

req_cols <- c("scenario", "conf_mult", "tau", "beta_trt", "t_star")
missing_cols <- setdiff(req_cols, names(results))
if (length(missing_cols))
  stop("Result files are missing column(s): ", paste(missing_cols, collapse = ", "),
       ". They were probably produced by an older run_sim_array.R.")

# ---- tidy --------------------------------------------------------------
plot_estimands <- c("RD", "RR", "OR", "RMST_diff")

# Order controls the dodge order and the legend.
# The two g-comp variants sit next to each other so the effect of
# clone-censoring on the g-formula arm reads directly off the plot.
method_levels <- c("fed_ccw_tvipcw", "pooled_ccw_tvipcw",
                   "ccw_gcomp", "gcomp_no_clone",
                   "tvipcw_no_clone", "perprotocol_naive")
method_labels <- c("Federated CCW", "Pooled CCW",
                   "CC-Gcomp", "G-comp (no cloning)",
                   "IPTW (no cloning)", "Per-protocol (no IPTW)")
method_colors <- c("Federated CCW"          = "#2c7fb8",
                   "Pooled CCW"             = "#41b6c4",
                   "CC-Gcomp"               = "#568203",
                   "G-comp (no cloning)"    = "#a3c853",
                   "IPTW (no cloning)"      = "#f5a623",
                   "Per-protocol (no IPTW)" = "#e8482c")

# Scenarios ordered by confounding strength, not alphabetically.
scen_levels <- names(sort(DEFAULT_CONF_MULTS))
scen_levels <- scen_levels[scen_levels %in% unique(results$scenario)]
scen_labels <- sprintf("%s confounding (x%g)",
                       scen_levels, DEFAULT_CONF_MULTS[scen_levels])

pdat <- results %>%
  filter(estimand %in% plot_estimands, method %in% method_levels) %>%
  mutate(
    estimand = factor(estimand, levels = plot_estimands),
    method   = factor(method, levels = method_levels, labels = method_labels),
    tau_f    = factor(tau, levels = sort(unique(tau)),
                      labels = paste0("tau=", sort(unique(tau)))),
    scen_f   = factor(scenario, levels = scen_levels, labels = scen_labels),
    beta_f   = factor(beta_trt, levels = sort(unique(beta_trt)),
                      labels = paste0("beta=", sort(unique(beta_trt))))
  )

n_reps_cell <- pdat %>%
  count(estimand, method, tau, scenario, beta_trt) %>%
  summarise(m = max(n)) %>% pull(m)

tstar_used <- sort(unique(pdat$t_star))
tstar_note <- if (length(tstar_used) == 1) sprintf("t*=%g", tstar_used) else
  sprintf("t* in {%s}", paste(tstar_used, collapse = ", "))

# ---- one figure per estimand -------------------------------------------
for (est in plot_estimands) {
  d <- pdat %>% filter(estimand == est)

  p <- ggplot(d, aes(x = tau_f, y = bias, fill = method)) +
    geom_boxplot(outlier.size = 0.4, alpha = 0.9,
                 position = position_dodge(width = 0.8),
                 linewidth = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    facet_grid(scen_f ~ beta_f) +
    scale_fill_manual(values = method_colors) +
    labs(
      title    = sprintf("Bias of %s by method (%d reps/cell)", est, n_reps_cell),
      subtitle = sprintf(paste("rows: confounding strength  |  cols: true effect beta",
                               "|  x: grace period tau  |  %s"), tstar_note),
      x = "Grace period (tau)",
      y = "Bias (estimate - truth)",
      fill = "Method"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey85", color = NA)
    )

  fname <- sprintf("bias_boxplot_%s.png", est)
  ggsave(fname, p, width = 16, height = 10, dpi = 150)
  message("Wrote ", fname)
}

# ---- numeric summary ---------------------------------------------------
summ <- pdat %>%
  group_by(method, estimand, tau, beta_trt, scenario, conf_mult, t_star) %>%
  summarise(
    mean_est = mean(estimate),
    truth    = first(truth),
    rmse     = sqrt(mean(bias^2)),   # before `bias` is collapsed to its mean
    bias     = mean(bias),
    emp_sd   = sd(estimate),
    coverage = if (all(is.na(covered))) NA_real_ else mean(covered, na.rm = TRUE),
    n_reps   = n(),
    .groups  = "drop"
  ) %>%
  arrange(estimand, method, conf_mult, beta_trt, tau) %>%
  select(method, estimand, tau, beta_trt, scenario, conf_mult, t_star,
         mean_est, truth, bias, emp_sd, rmse, coverage, n_reps) %>%
  as.data.frame()

write.csv(summ, "summary_all_scenarios.csv", row.names = FALSE)
message("Wrote summary_all_scenarios.csv")