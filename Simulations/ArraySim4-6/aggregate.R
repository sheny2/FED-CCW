#!/usr/bin/env Rscript
# Aggregate ArraySim4-6 site-size balance x patient-mix results.

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

files <- list.files("results", pattern = "^res_task_[0-9]+[.]rds$",
                    full.names = TRUE)
if (!length(files)) stop("No result files found in results/.")
results <- do.call(rbind, lapply(files, readRDS))
rownames(results) <- NULL
message(sprintf("Loaded %d rows from %d task files.", nrow(results), length(files)))

site_diag_cols <- unlist(lapply(seq_len(STUDY_K), function(k) {
  paste0("site", k, c("_x1_mean", "_x1_sd", "_x2_prob"))
}))
req_cols <- c(
  "method", "estimand", "estimate", "truth", "bias", "rep", "tau",
  "beta_trt", "size_scenario", "mix_scenario", "n_sites", "total_n",
  "site_sizes", "site_weights", "target_x1_mean", "target_x2_prob",
  "observed_x1_mean", "observed_x1_sd", "observed_x2_prob",
  "observed_death_rate", "observed_init_rate", "landmark_n",
  "landmark_event_n", "landmark_event_rate", "landmark_treated_rate",
  "task_id", site_diag_cols
)
missing_cols <- setdiff(req_cols, names(results))
if (length(missing_cols))
  stop("Result files are missing columns: ", paste(missing_cols, collapse = ", "))

expected_tasks <- seq_len(nrow(make_study_grid()))
observed_tasks <- sort(unique(results$task_id))
if (!identical(observed_tasks, expected_tasks))
  stop("Expected task IDs ", paste(expected_tasks, collapse = ", "),
       "; observed: ", paste(observed_tasks, collapse = ", "))

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
  "Federated CCW" = "#2c7fb8", "Pooled CCW" = "#41b6c4",
  "Federated landmark IPW" = "#7a5195",
  "Local CCW + curve meta-analysis" = "#2ca25f",
  "Federated IPW (no cloning)" = "#f5a623",
  "Federated per-protocol (unweighted)" = "#e8482c"
)
if (!identical(sort(unique(results$method)), sort(method_levels)))
  stop("Unexpected method set: ", paste(sort(unique(results$method)), collapse = ", "))

cell_counts <- results %>% count(task_id, method, estimand)
if (any(cell_counts$n != STUDY_N_REPS))
  warning("At least one task/method/estimand cell has an unexpected replicate count.")

# Keep squared error separate from mean bias. Overwriting `bias` inside
# summarise before RMSE is computed would incorrectly make RMSE = |mean bias|.
results <- results %>% mutate(squared_error = bias^2)

summ <- results %>%
  group_by(
    method, estimand, tau, beta_trt, scenario, conf_mult,
    size_scenario, mix_scenario, n_sites, total_n, site_sizes, site_weights,
    target_x1_mean, target_x2_prob, init_intercept
  ) %>%
  summarise(
    mean_est = mean(estimate), median_est = median(estimate),
    truth = first(truth), mean_bias = mean(bias), median_bias = median(bias),
    emp_sd = sd(estimate), rmse = sqrt(mean(squared_error)),
    coverage = if (all(is.na(covered))) NA_real_ else mean(covered, na.rm = TRUE),
    mean_observed_x1 = mean(observed_x1_mean),
    mean_observed_x2_prob = mean(observed_x2_prob),
    mean_observed_death_rate = mean(observed_death_rate),
    mean_observed_init_rate = mean(observed_init_rate),
    mean_landmark_n = mean(landmark_n),
    mean_landmark_event_n = mean(landmark_event_n),
    mean_landmark_event_rate = mean(landmark_event_rate, na.rm = TRUE),
    mean_landmark_treated_rate = mean(landmark_treated_rate, na.rm = TRUE),
    n_reps = n(), .groups = "drop"
  ) %>%
  arrange(estimand, method, mix_scenario, size_scenario, tau) %>%
  as.data.frame()
write.csv(summ, "summary_all_scenarios.csv", row.names = FALSE)
message("Wrote summary_all_scenarios.csv")

diagnostic_cols <- c(
  "task_id", "rep", "tau", "size_scenario", "mix_scenario", "total_n",
  "site_sizes", "site_weights", "target_x1_mean", "target_x2_prob",
  "observed_x1_mean", "observed_x1_sd", "observed_x2_prob", site_diag_cols,
  "observed_death_rate", "min_site_death_rate", "max_site_death_rate",
  "observed_init_rate", "min_site_init_rate", "max_site_init_rate",
  "landmark_n", "landmark_event_n", "landmark_event_rate",
  "landmark_treated_rate"
)
design_diagnostics <- results %>%
  distinct(across(all_of(diagnostic_cols))) %>%
  arrange(task_id, rep) %>% as.data.frame()
write.csv(design_diagnostics, "design_diagnostics.csv", row.names = FALSE)
message("Wrote design_diagnostics.csv")

design_targets <- do.call(rbind, lapply(seq_len(nrow(make_study_grid())), function(i) {
  cell <- make_study_grid()[i, ]
  sizes <- get_site_sizes(cell$size_scenario)
  mix <- get_site_mix(cell$mix_scenario)
  weights <- sizes / sum(sizes)
  data.frame(
    task_id = i, tau = cell$tau, size_scenario = cell$size_scenario,
    mix_scenario = cell$mix_scenario, site = seq_len(STUDY_K),
    site_n = unname(sizes), site_weight = unname(weights),
    x1_mean = mix$x1_mean, x1_sd = mix$x1_sd, x2_prob = mix$x2_prob,
    pooled_x1_mean = sum(weights * mix$x1_mean),
    pooled_x2_prob = sum(weights * mix$x2_prob)
  )
}))
write.csv(design_targets, "design_targets.csv", row.names = FALSE)
message("Wrote design_targets.csv")

ratio_extremes <- results %>%
  filter(estimand %in% c("RR", "OR")) %>%
  group_by(method, estimand, tau, size_scenario, mix_scenario) %>%
  summarise(n_reps = n(), n_gt_2 = sum(estimate > 2),
            n_gt_10 = sum(estimate > 10), n_gt_100 = sum(estimate > 100),
            maximum = max(estimate), .groups = "drop") %>% as.data.frame()
write.csv(ratio_extremes, "ratio_extreme_counts.csv", row.names = FALSE)
message("Wrote ratio_extreme_counts.csv")

size_labels <- c(
  balanced = "Balanced sites (1,000 / 1,000 / 1,000)",
  unbalanced = "Unbalanced sites (200 / 800 / 2,000)"
)
mix_labels <- c(
  homogeneous = "Homogeneous patient mix",
  heterogeneous = "Heterogeneous patient mix"
)
plot_data <- results %>%
  filter(estimand %in% c("RD", "RR", "OR", "RMST_diff")) %>%
  mutate(
    method = factor(method, levels = method_levels, labels = method_labels),
    tau_f = factor(tau, levels = STUDY_TAUS, labels = paste0("tau=", STUDY_TAUS)),
    size_f = factor(size_scenario, levels = names(size_labels),
                    labels = unname(size_labels)),
    mix_f = factor(mix_scenario, levels = names(mix_labels),
                   labels = unname(mix_labels))
  )

common_theme <- theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey85", color = NA))

for (est in c("RD", "RMST_diff")) {
  d <- plot_data %>% filter(estimand == est)
  p <- ggplot(d, aes(x = tau_f, y = bias, fill = method)) +
    geom_boxplot(outlier.size = 0.4, alpha = 0.9,
                 position = position_dodge(width = 0.8), linewidth = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    facet_grid(mix_f ~ size_f) + scale_fill_manual(values = method_colors) +
    labs(title = sprintf("Bias of %s: site size x patient mix", est),
         subtitle = sprintf(paste0("%d reps/cell | 3 sites | total N=3,000 | ",
                                   "medium confounding | beta=-0.7 | t*=25"),
                            STUDY_N_REPS),
         x = "Grace period (tau)", y = "Bias (estimate - truth)", fill = "Method") +
    common_theme
  filename <- sprintf("bias_boxplot_%s.png", est)
  ggsave(filename, p, width = 14, height = 9, dpi = 150)
  message("Wrote ", filename)
}

for (est in c("RR", "OR")) {
  d <- plot_data %>% filter(estimand == est, estimate > 0, truth > 0) %>%
    mutate(log_ratio_error = log(estimate / truth))
  p <- ggplot(d, aes(x = tau_f, y = log_ratio_error, fill = method)) +
    geom_boxplot(outlier.size = 0.4, alpha = 0.9,
                 position = position_dodge(width = 0.8), linewidth = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    facet_grid(mix_f ~ size_f) + scale_fill_manual(values = method_colors) +
    labs(title = sprintf("Log-ratio error of %s: site size x patient mix", est),
         subtitle = sprintf(paste0("log(estimate/truth); zero is ideal | ",
                                   "%d reps/cell | 3 sites | total N=3,000"),
                            STUDY_N_REPS),
         x = "Grace period (tau)", y = "log(estimate / truth)", fill = "Method") +
    common_theme
  filename <- sprintf("log_ratio_error_boxplot_%s.png", est)
  ggsave(filename, p, width = 14, height = 9, dpi = 150)
  message("Wrote ", filename)
}
