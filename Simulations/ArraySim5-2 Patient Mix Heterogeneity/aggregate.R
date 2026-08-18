#!/usr/bin/env Rscript
# Combine all ArraySim5-2 Patient Mix Heterogeneity results.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("params.R")
validate_params()

files <- list.files("results", pattern = "^res_task_.*\\.rds$",
                    full.names = TRUE)
if (!length(files)) stop("No result files found in results/.")

results <- do.call(rbind, lapply(files, readRDS))
rownames(results) <- NULL
message(sprintf("Loaded %d rows from %d task files.",
                nrow(results), length(files)))

req_cols <- c("method", "estimand", "estimate", "truth", "bias", "rep",
              "tau", "beta_trt", "heterogeneity", "patient_mix", "t_star",
              "site_sizes", "init_intercepts", "task_id")
missing_cols <- setdiff(req_cols, names(results))
if (length(missing_cols))
  stop("Result files are missing columns: ",
       paste(missing_cols, collapse = ", "))

expected_tasks <- seq_len(nrow(simulation_grid()))
observed_tasks <- sort(unique(results$task_id))
if (!identical(observed_tasks, expected_tasks)) {
  stop("Expected task IDs 1:", length(expected_tasks), "; observed: ",
       paste(observed_tasks, collapse = ", "))
}

plot_estimands <- c("RD", "RR", "OR", "RMST_diff")
method_levels <- c(
  "fed_ccw_tvipcw", "pooled_ccw_tvipcw",
  "fed_landmark_ipw",
  "fed_ipw_no_clone", "fed_perprotocol_naive"
)
method_labels <- c(
  "Federated CCW", "Pooled CCW",
  "Federated landmark IPW",
  "Federated IPW (no cloning)", "Federated per-protocol (unweighted)"
)
method_colors <- c(
  "Federated CCW"                       = "#0072B2",
  "Pooled CCW"                          = "#56B4E9",
  "Federated IPW (no cloning)"          = "#E69F00",
  "Federated per-protocol (unweighted)" = "#D55E00",
  "Federated landmark IPW"              = "#7A5195"
)

observed_methods <- sort(unique(results$method))
if (!identical(observed_methods, sort(method_levels))) {
  stop("Unexpected method set: ", paste(observed_methods, collapse = ", "))
}

heterogeneity_levels <- names(DEFAULT_PATIENT_MIX)
tau_values <- sort(unique(results$tau))
tau_labels <- setNames(
  as.expression(lapply(tau_values, function(x) bquote(tau == .(x)))),
  as.character(tau_values)
)
pdat <- results %>%
  filter(estimand %in% plot_estimands) %>%
  mutate(
    estimand = factor(estimand, levels = plot_estimands),
    method = factor(method, levels = method_levels, labels = method_labels),
    tau_f = factor(tau, levels = tau_values),
    heterogeneity_f = factor(
      heterogeneity, levels = heterogeneity_levels,
      labels = paste(tools::toTitleCase(heterogeneity_levels),
                     "patient heterogeneity")
    ),
    beta_f = factor(
      beta_trt,
      levels = sort(unique(beta_trt)),
      labels = sprintf("beta[trt] == %.1f", sort(unique(beta_trt)))
    )
  )

rep_counts <- pdat %>%
  count(estimand, method, tau, heterogeneity, beta_trt)
# if (any(rep_counts$n != 100L)) {
#   warning("At least one method/estimand/scenario cell does not have 100 reps.")
# }
rep_note <- sprintf("%d-%d reps/cell", min(rep_counts$n), max(rep_counts$n))
if (length(unique(rep_counts$n)) == 1L)
  rep_note <- sprintf("%d reps/cell", unique(rep_counts$n))

for (est in plot_estimands) {
  d <- pdat %>% filter(estimand == est) %>% mutate(etimand = ifelse(est == "RMST_diff", "RMST difference", est))
  p <- ggplot(d, aes(x = tau_f, y = bias, fill = method)) +
    geom_boxplot(
      outlier.size = 0.45, outlier.alpha = 0.45, alpha = 0.88,
      position = position_dodge(width = 0.82), linewidth = 0.3,
      width = 0.72
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45,
               color = "grey30") +
    facet_grid(
      heterogeneity_f ~ beta_f,
      labeller = labeller(beta_f = label_parsed)
    ) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    scale_x_discrete(labels = tau_labels) +
    labs(
      title = paste0("Estimation performance: ", est),
      subtitle = paste0(
        length(DEFAULT_SITE_SIZES), " sites × ",
        if (length(unique(DEFAULT_SITE_SIZES)) == 1L) {
          paste0(unique(DEFAULT_SITE_SIZES), " patients")
        } else {
          paste0("site sizes ", paste(DEFAULT_SITE_SIZES, collapse = "/"))
        },
        "; initiation intercepts from ", min(DEFAULT_SITE_INIT_INTERCEPTS),
        " to ", max(DEFAULT_SITE_INIT_INTERCEPTS),
        "; follow-up horizon = ", SIM_TSTAR
      ),
      x = expression("Grace period"~(tau)),
      y = "Bias", fill = "Method"
      # caption = expression("Columns show treatment-effect coefficient "~beta[trt]~"; rows show patient-mix heterogeneity.")
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey92", color = "grey65"),
      strip.text = element_text(face = "bold", size = 10.5),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "grey25"),
      plot.caption = element_text(hjust = 0, color = "grey35"),
      plot.margin = margin(10, 14, 10, 10)
    )

  fname <- sprintf("bias_boxplot_%s.png", est)
  ggsave(fname, p, width = 15, height = 9, dpi = 300, bg = "white")
  message("Wrote ", fname)
}

summ <- pdat %>%
  group_by(method, estimand, tau, beta_trt, heterogeneity, patient_mix,
           t_star, site_sizes, init_intercepts) %>%
  summarise(
    mean_est = mean(estimate),
    truth = first(truth),
    rmse = sqrt(mean(bias^2)),
    bias = mean(bias),
    emp_sd = sd(estimate),
    coverage = if (all(is.na(covered))) NA_real_ else
      mean(covered, na.rm = TRUE),
    n_reps = n(),
    .groups = "drop"
  ) %>%
  arrange(estimand, method, heterogeneity, beta_trt, tau) %>%
  as.data.frame()

write.csv(summ, "summary_all_scenarios.csv", row.names = FALSE)
message("Wrote summary_all_scenarios.csv")
