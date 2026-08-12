#!/usr/bin/env Rscript
# Combine all ArraySim4-2 Unbalance results and summarize performance.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

files <- list.files("results", pattern = "^res_task_.*\\.rds$",
                    full.names = TRUE)
if (!length(files)) stop("No result files found in results/.")

results <- do.call(rbind, lapply(files, readRDS))
rownames(results) <- NULL
message(sprintf("Loaded %d rows from %d task files.",
                nrow(results), length(files)))

req_cols <- c("method", "estimand", "estimate", "truth", "bias", "rep",
              "tau", "beta_trt", "scenario", "conf_mult", "t_star",
              "site_sizes", "init_intercepts", "task_id")
missing_cols <- setdiff(req_cols, names(results))
if (length(missing_cols))
  stop("Result files are missing columns: ",
       paste(missing_cols, collapse = ", "))

expected_tasks <- seq_len(27L)
observed_tasks <- sort(unique(results$task_id))
if (!identical(observed_tasks, expected_tasks)) {
  stop("Expected task IDs 1:27; observed: ",
       paste(observed_tasks, collapse = ", "))
}

plot_estimands <- c("RD", "RR", "OR", "RMST_diff")
method_levels <- c(
  "fed_ccw_tvipcw", "pooled_ccw_tvipcw",
  "fed_ipw_no_clone", "fed_perprotocol_naive",
  "fed_landmark_ipw", "local_ccw_meta"
)
method_labels <- c(
  "Federated CCW", "Pooled CCW",
  "Federated IPW (no cloning)", "Federated per-protocol (unweighted)",
  "Federated landmark IPW", "Local CCW + curve meta-analysis"
)
method_colors <- c(
  "Federated CCW"                         = "#2c7fb8",
  "Pooled CCW"                            = "#41b6c4",
  "Federated IPW (no cloning)"            = "#f5a623",
  "Federated per-protocol (unweighted)"   = "#e8482c",
  "Federated landmark IPW"                = "#7a5195",
  "Local CCW + curve meta-analysis"       = "#2ca25f"
)

observed_methods <- sort(unique(results$method))
if (!identical(observed_methods, sort(method_levels))) {
  stop("Unexpected method set: ", paste(observed_methods, collapse = ", "))
}

scenario_levels <- c("small", "medium", "strong")
pdat <- results %>%
  filter(estimand %in% plot_estimands) %>%
  mutate(
    estimand = factor(estimand, levels = plot_estimands),
    method = factor(method, levels = method_levels, labels = method_labels),
    tau_f = factor(tau, levels = sort(unique(tau)),
                   labels = paste0("tau=", sort(unique(tau)))),
    scenario_f = factor(
      scenario, levels = scenario_levels,
      labels = paste(tools::toTitleCase(scenario_levels), "confounding")
    ),
    beta_f = factor(beta_trt, levels = sort(unique(beta_trt)),
                    labels = paste0("beta=", sort(unique(beta_trt))))
  )

rep_counts <- pdat %>%
  count(estimand, method, tau, scenario, beta_trt)
if (any(rep_counts$n != 100L)) {
  warning("At least one method/estimand/scenario cell does not have 100 reps.")
}
rep_note <- sprintf("%d-%d reps/cell", min(rep_counts$n), max(rep_counts$n))
if (length(unique(rep_counts$n)) == 1L)
  rep_note <- sprintf("%d reps/cell", unique(rep_counts$n))

for (est in plot_estimands) {
  d <- pdat %>% filter(estimand == est)
  p <- ggplot(d, aes(x = tau_f, y = bias, fill = method)) +
    geom_boxplot(
      outlier.size = 0.4, alpha = 0.9,
      position = position_dodge(width = 0.8), linewidth = 0.25
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    facet_grid(scenario_f ~ beta_f) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    labs(
      title = sprintf("Bias of %s by confounding strength (%s)", est, rep_note),
      subtitle = paste(
        "rows: confounding | columns: true treatment effect | t*=25 |",
        "site n=200/800/2000; initiation intercepts=-4.5/-3/-1.5"
      ),
      x = "Grace period (tau)", y = "Bias (estimate - truth)", fill = "Method"
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

summ <- pdat %>%
  group_by(method, estimand, tau, beta_trt, scenario, conf_mult, t_star,
           site_sizes, init_intercepts) %>%
  summarise(
    mean_est = mean(estimate),
    truth = first(truth),
    bias = mean(bias),
    emp_sd = sd(estimate),
    rmse = sqrt(mean(bias^2)),
    coverage = if (all(is.na(covered))) NA_real_ else
      mean(covered, na.rm = TRUE),
    n_reps = n(),
    .groups = "drop"
  ) %>%
  arrange(estimand, method, scenario, beta_trt, tau) %>%
  as.data.frame()

write.csv(summ, "summary_all_scenarios.csv", row.names = FALSE)
message("Wrote summary_all_scenarios.csv")
