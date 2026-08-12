#!/usr/bin/env Rscript
# ============================================================================
# Combine per-task result files and produce grid-faceted bias boxplots for
# the four-method misspecification study.
#
# For each estimand and each misspecification spec, one figure:
#   rows    = confounding strength (small / medium / strong)
#   columns = true effect beta_trt
#   x-axis  = grace period tau
#   fill    = method
#
# Plus, for each estimand, a misspecification-focused figure with the spec on
# the x-axis so the cost of each nuisance-model error reads directly.
#
# Run AFTER the array job finishes:  Rscript aggregate.R
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

source("params.R")   # DEFAULT_CONF_MULTS, MISSPEC_SPECS

files <- list.files("results", pattern = "^res_task_.*\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No result files found in results/.")

results <- do.call(rbind, lapply(files, readRDS))
rownames(results) <- NULL
message(sprintf("Loaded %d rows from %d task files.", nrow(results), length(files)))

req_cols <- c("scenario", "conf_mult", "tau", "beta_trt", "t_star", "misspec")
missing_cols <- setdiff(req_cols, names(results))
if (length(missing_cols))
  stop("Result files are missing column(s): ", paste(missing_cols, collapse = ", "),
       ". They were probably produced by an older run_sim_array.R.")

# ---- tidy --------------------------------------------------------------
plot_estimands <- c("RD", "RR", "OR", "RMST_diff")

method_levels <- c("fed_ccw_tvipcw", "pooled_ccw_tvipcw",
                   "gcomp_no_clone", "ccw_gcomp")
method_labels <- c("Federated CCW", "Pooled CCW",
                   "Direct G-comp", "CC-Gcomp")
method_colors <- c("Federated CCW" = "#2c7fb8",
                   "Pooled CCW"    = "#41b6c4",
                   "Direct G-comp" = "#a3c853",
                   "CC-Gcomp"      = "#568203")

# Scenarios ordered by confounding strength.
scen_levels <- names(sort(DEFAULT_CONF_MULTS))
scen_levels <- scen_levels[scen_levels %in% unique(results$scenario)]
scen_labels <- sprintf("%s confounding (x%g)",
                       scen_levels, DEFAULT_CONF_MULTS[scen_levels])

# Misspec ordered as declared in params.R.
misspec_levels <- intersect(MISSPEC_SPECS, unique(results$misspec))
message(sprintf("Specifications present: %s",
                paste(misspec_levels, collapse = ", ")))
missing_spec <- setdiff(MISSPEC_SPECS, unique(results$misspec))
if (length(missing_spec))
  message(sprintf("Specifications declared but absent from results (not plotted): %s",
                  paste(missing_spec, collapse = ", ")))

pdat <- results %>%
  filter(estimand %in% plot_estimands, method %in% method_levels) %>%
  mutate(
    estimand = factor(estimand, levels = plot_estimands),
    method   = factor(method, levels = method_levels, labels = method_labels),
    tau_f    = factor(tau, levels = sort(unique(tau)),
                      labels = paste0("tau=", sort(unique(tau)))),
    scen_f   = factor(scenario, levels = scen_levels, labels = scen_labels),
    beta_f   = factor(beta_trt, levels = sort(unique(beta_trt)),
                      labels = paste0("beta=", sort(unique(beta_trt)))),
    misspec  = factor(misspec, levels = misspec_levels)
  )

n_reps_cell <- pdat %>%
  count(estimand, method, tau, scenario, beta_trt, misspec) %>%
  summarise(m = max(n)) %>% pull(m)

tstar_used <- sort(unique(pdat$t_star))
tstar_note <- if (length(tstar_used) == 1) sprintf("t*=%g", tstar_used) else
  sprintf("t* in {%s}", paste(tstar_used, collapse = ", "))

# ---- one figure per estimand x misspecification ------------------------
for (est in plot_estimands) {
  for (ms in misspec_levels) {
    d <- pdat %>% filter(estimand == est, misspec == ms)
    if (nrow(d) == 0) next

    p <- ggplot(d, aes(x = tau_f, y = bias, fill = method)) +
      geom_boxplot(outlier.size = 0.4, alpha = 0.9,
                   position = position_dodge(width = 0.8),
                   linewidth = 0.25) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
      facet_grid(scen_f ~ beta_f) +
      scale_fill_manual(values = method_colors) +
      labs(
        title    = sprintf("Bias of %s -- nuisance model: %s (%d reps/cell)",
                           est, ms, n_reps_cell),
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

    fname <- sprintf("bias_boxplot_%s_%s.png", est, ms)
    ggsave(fname, p, width = 14, height = 10, dpi = 150)
    message("Wrote ", fname)
  }
}

# ---- misspecification-focused figure: spec on the x-axis ---------------
if (length(misspec_levels) > 1) {
  for (est in plot_estimands) {
    d <- pdat %>% filter(estimand == est)

    p <- ggplot(d, aes(x = misspec, y = bias, fill = method)) +
      geom_boxplot(outlier.size = 0.4, alpha = 0.9,
                   position = position_dodge(width = 0.8),
                   linewidth = 0.25) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
      facet_grid(scen_f ~ beta_f) +
      scale_fill_manual(values = method_colors) +
      labs(
        title    = sprintf("Bias of %s under nuisance-model misspecification (%d reps/cell)",
                           est, n_reps_cell),
        subtitle = paste("rows: confounding strength  |  cols: true effect beta",
                         "|  x: nuisance-model specification  |  tau pooled"),
        x = "Nuisance-model specification",
        y = "Bias (estimate - truth)",
        fill = "Method"
      ) +
      theme_bw(base_size = 11) +
      theme(
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_rect(fill = "grey85", color = NA)
      )

    fname <- sprintf("misspec_%s.png", est)
    ggsave(fname, p, width = 12, height = 10, dpi = 150)
    message("Wrote ", fname)
  }
}

# ---- numeric summary ---------------------------------------------------
summ <- pdat %>%
  group_by(method, estimand, tau, beta_trt, scenario, conf_mult, t_star,
           misspec) %>%
  summarise(
    mean_est = mean(estimate),
    truth    = first(truth),
    rmse     = sqrt(mean(bias^2)),
    bias     = mean(bias),
    emp_sd   = sd(estimate),
    coverage = if (all(is.na(covered))) NA_real_ else mean(covered, na.rm = TRUE),
    n_reps   = n(),
    .groups  = "drop"
  ) %>%
  arrange(estimand, method, conf_mult, beta_trt, tau, misspec) %>%
  select(method, estimand, tau, beta_trt, scenario, conf_mult, t_star, misspec,
         mean_est, truth, bias, emp_sd, rmse, coverage, n_reps) %>%
  as.data.frame()

write.csv(summ, "summary_all_scenarios.csv", row.names = FALSE)
message("Wrote summary_all_scenarios.csv")
