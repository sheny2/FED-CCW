#!/usr/bin/env Rscript
# Aggregate ArraySim4-5 outcome-frequency factorial results.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
  setwd(dirname(script_path))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("params.R")

files <- list.files("results", pattern = "^res_task_.*[.]rds$",
                    full.names = TRUE)
if (!length(files)) stop("No result files found in results/.")
results <- do.call(rbind, lapply(files, readRDS))
rownames(results) <- NULL
message(sprintf("Loaded %d rows from %d task files.",
                nrow(results), length(files)))

req_cols <- c(
  "method", "estimand", "estimate", "truth", "bias", "rep", "tau",
  "beta_trt", "sample_size_scenario", "outcome_scenario", "n_sites",
  "n_per_site", "total_n", "event_intercept", "init_intercept",
  "observed_death_rate", "min_site_death_rate", "max_site_death_rate",
  "observed_init_rate", "landmark_n", "landmark_event_n",
  "landmark_event_rate", "landmark_treated_rate", "task_id"
)
missing_cols <- setdiff(req_cols, names(results))
if (length(missing_cols))
  stop("Result files are missing columns: ", paste(missing_cols, collapse = ", "))

expected_tasks <- seq_len(12L)
observed_tasks <- sort(unique(results$task_id))
if (!identical(observed_tasks, expected_tasks))
  stop("Expected task IDs 1:12; observed: ",
       paste(observed_tasks, collapse = ", "))

method_levels <- c(
  "fed_ccw_tvipcw", "pooled_ccw_tvipcw", "fed_landmark_ipw",
  "local_ccw_meta", "fed_ipw_no_clone", "fed_perprotocol_naive"
)
method_labels <- c(
  "Federated CCW", "Pooled CCW", "Federated landmark IPW",
  "Local CCW + curve meta-analysis", "Federated IPW (no cloning)",
  "Federated per-protocol (unweighted)"
)
method_colors <- c(
  "Federated CCW" = "#2c7fb8",
  "Pooled CCW" = "#41b6c4",
  "Federated landmark IPW" = "#7a5195",
  "Local CCW + curve meta-analysis" = "#2ca25f",
  "Federated IPW (no cloning)" = "#f5a623",
  "Federated per-protocol (unweighted)" = "#e8482c"
)
if (!identical(sort(unique(results$method)), sort(method_levels)))
  stop("Unexpected method set: ", paste(sort(unique(results$method)),
                                        collapse = ", "))

cell_counts <- results %>% count(task_id, method, estimand)
if (any(cell_counts$n != STUDY_N_REPS))
  warning("At least one task/method/estimand cell has an unexpected replicate count.")

# Precompute squared error under a distinct name. This prevents the sequential
# dplyr summarise bug in ArraySim4-4, where overwriting `bias` before calculating
# RMSE made every reported RMSE equal abs(mean bias).
results <- results %>% mutate(squared_error = bias^2)

summ <- results %>%
  group_by(
    method, estimand, tau, beta_trt, scenario, conf_mult,
    sample_size_scenario, outcome_scenario, n_sites, n_per_site, total_n,
    event_intercept, init_intercept
  ) %>%
  summarise(
    mean_est = mean(estimate),
    median_est = median(estimate),
    truth = first(truth),
    mean_bias = mean(bias),
    median_bias = median(bias),
    emp_sd = sd(estimate),
    rmse = sqrt(mean(squared_error)),
    coverage = if (all(is.na(covered))) NA_real_ else mean(covered, na.rm = TRUE),
    mean_observed_death_rate = mean(observed_death_rate),
    mean_min_site_death_rate = mean(min_site_death_rate),
    mean_max_site_death_rate = mean(max_site_death_rate),
    mean_observed_init_rate = mean(observed_init_rate),
    mean_landmark_n = mean(landmark_n),
    mean_landmark_event_n = mean(landmark_event_n),
    mean_landmark_event_rate = mean(landmark_event_rate, na.rm = TRUE),
    mean_landmark_treated_rate = mean(landmark_treated_rate, na.rm = TRUE),
    n_reps = n(),
    .groups = "drop"
  ) %>%
  arrange(estimand, method, outcome_scenario, sample_size_scenario, tau) %>%
  as.data.frame()
write.csv(summ, "summary_all_scenarios.csv", row.names = FALSE)
message("Wrote summary_all_scenarios.csv")

design_diagnostics <- results %>%
  distinct(
    task_id, rep, tau, sample_size_scenario, outcome_scenario,
    n_per_site, total_n, event_intercept, init_intercept,
    observed_death_rate, min_site_death_rate, max_site_death_rate,
    observed_init_rate, min_site_init_rate, max_site_init_rate,
    landmark_n, landmark_event_n, landmark_event_rate, landmark_treated_rate
  ) %>%
  arrange(task_id, rep) %>%
  as.data.frame()
write.csv(design_diagnostics, "design_diagnostics.csv", row.names = FALSE)
message("Wrote design_diagnostics.csv")

ratio_extremes <- results %>%
  filter(estimand %in% c("RR", "OR")) %>%
  group_by(method, estimand, tau, sample_size_scenario, outcome_scenario) %>%
  summarise(
    n_reps = n(),
    n_gt_2 = sum(estimate > 2),
    n_gt_10 = sum(estimate > 10),
    n_gt_100 = sum(estimate > 100),
    maximum = max(estimate),
    .groups = "drop"
  ) %>%
  as.data.frame()
write.csv(ratio_extremes, "ratio_extreme_counts.csv", row.names = FALSE)
message("Wrote ratio_extreme_counts.csv")

sample_labels <- sprintf(
  "%s patients/site",
  format(unname(SAMPLE_SIZE_LEVELS[c("low", "large")]), big.mark = ",")
)
outcome_labels <- sprintf(
  "%s death (event intercept %g)",
  c("Rare", "Common"), unname(OUTCOME_LEVELS[c("rare", "common")])
)

plot_data <- results %>%
  filter(estimand %in% c("RD", "RR", "OR", "RMST_diff")) %>%
  mutate(
    method = factor(method, levels = method_levels, labels = method_labels),
    tau_f = factor(tau, levels = STUDY_TAUS,
                   labels = paste0("tau=", STUDY_TAUS)),
    sample_f = factor(sample_size_scenario, levels = c("low", "large"),
                      labels = sample_labels),
    outcome_f = factor(outcome_scenario, levels = c("rare", "common"),
                       labels = outcome_labels)
  )

for (est in c("RD", "RMST_diff")) {
  d <- plot_data %>% filter(estimand == est)
  p <- ggplot(d, aes(x = tau_f, y = bias, fill = method)) +
    geom_boxplot(outlier.size = 0.4, alpha = 0.9,
                 position = position_dodge(width = 0.8), linewidth = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    facet_grid(outcome_f ~ sample_f) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    labs(
      title = sprintf("Bias of %s: sample size x outcome frequency", est),
      subtitle = sprintf(
        "%d reps/cell | medium confounding | beta=-0.7 | t*=25 | 10 sites | initiation intercept=-3",
        STUDY_N_REPS
      ),
      x = "Grace period (tau)", y = "Bias (estimate - truth)", fill = "Method"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey85", color = NA))
  filename <- sprintf("bias_boxplot_%s.png", est)
  ggsave(filename, p, width = 14, height = 9, dpi = 150)
  message("Wrote ", filename)
}

for (est in c("RR", "OR")) {
  d <- plot_data %>%
    filter(estimand == est, estimate > 0, truth > 0) %>%
    mutate(log_ratio_error = log(estimate / truth))
  p <- ggplot(d, aes(x = tau_f, y = log_ratio_error, fill = method)) +
    geom_boxplot(outlier.size = 0.4, alpha = 0.9,
                 position = position_dodge(width = 0.8), linewidth = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    facet_grid(outcome_f ~ sample_f) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    labs(
      title = sprintf("Log-ratio error of %s: sample size x outcome frequency", est),
      subtitle = sprintf(
        "log(estimate/truth); zero is ideal | %d reps/cell | medium confounding | beta=-0.7 | t*=25",
        STUDY_N_REPS
      ),
      x = "Grace period (tau)", y = "log(estimate / truth)", fill = "Method"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey85", color = NA))
  filename <- sprintf("log_ratio_error_boxplot_%s.png", est)
  ggsave(filename, p, width = 14, height = 9, dpi = 150)
  message("Wrote ", filename)
}
