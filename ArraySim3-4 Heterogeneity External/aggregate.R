#!/usr/bin/env Rscript
# Combine all array-task results, make heterogeneity-faceted bias plots, and
# write a numeric simulation summary.  Run from this directory after SLURM
# finishes: Rscript aggregate.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

source("params.R")

files <- list.files("results", pattern = "^res_task_.*\\.rds$",
                    full.names = TRUE)
if (!length(files)) stop("No result files found in results/.")

results <- do.call(rbind, lapply(files, readRDS))
rownames(results) <- NULL
message(sprintf("Loaded %d rows from %d task files.",
                nrow(results), length(files)))

req_cols <- c("heterogeneity", "tau", "beta_trt", "t_star",
              "numerator", "task_id")
missing_cols <- setdiff(req_cols, names(results))
if (length(missing_cols))
  stop("Result files are missing: ", paste(missing_cols, collapse = ", "))
if (length(unique(results$numerator)) != 1L ||
    !all(startsWith(results$numerator, "external_")))
  stop("Results include an incompatible numerator specification.")

plot_estimands <- c("RD", "RR", "OR", "RMST_diff")
method_levels <- c(
  "fed_ccw_tvipcw", "pooled_ccw_tvipcw",
  "tvipcw_no_clone", "perprotocol_naive"
)
method_labels <- c(
  "Federated CCW (external numerator)", "Pooled CCW",
  "TV-IPCW (no cloning)", "Per-protocol (naive)"
)
method_colors <- c(
  "Federated CCW (external numerator)" = "#2c7fb8",
  "Pooled CCW"                      = "#41b6c4",
  "TV-IPCW (no cloning)"            = "#f5a623",
  "Per-protocol (naive)"            = "#e8482c"
)

hetero_levels <- names(DEFAULT_SITE_MIXES)
pdat <- results %>%
  filter(estimand %in% plot_estimands, method %in% method_levels) %>%
  mutate(
    estimand = factor(estimand, levels = plot_estimands),
    method = factor(method, levels = method_levels, labels = method_labels),
    tau_f = factor(tau, levels = sort(unique(tau)),
                   labels = paste0("tau=", sort(unique(tau)))),
    hetero_f = factor(
      heterogeneity,
      levels = hetero_levels,
      labels = paste(tools::toTitleCase(hetero_levels), "heterogeneity")
    ),
    beta_f = factor(beta_trt, levels = sort(unique(beta_trt)),
                    labels = paste0("beta=", sort(unique(beta_trt))))
  )

if (!nrow(pdat)) stop("No rows remain after method/estimand filtering.")

rep_counts <- pdat %>%
  count(estimand, method, tau, heterogeneity, beta_trt)
rep_note <- if (length(unique(rep_counts$n)) == 1L) {
  sprintf("%d reps/cell", unique(rep_counts$n))
} else {
  sprintf("%d-%d reps/cell", min(rep_counts$n), max(rep_counts$n))
}
tstar_used <- sort(unique(pdat$t_star))
tstar_note <- if (length(tstar_used) == 1L) {
  sprintf("t*=%g", tstar_used)
} else {
  sprintf("t* in {%s}", paste(tstar_used, collapse = ", "))
}

for (est in plot_estimands) {
  d <- pdat %>% filter(estimand == est)
  p <- ggplot(d, aes(x = tau_f, y = bias, fill = method)) +
    geom_boxplot(
      outlier.size = 0.4, alpha = 0.9,
      position = position_dodge(width = 0.8), linewidth = 0.25
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    facet_grid(hetero_f ~ beta_f) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    labs(
      title = sprintf("Bias of %s by site heterogeneity (%s)", est, rep_note),
      subtitle = sprintf(
        "rows: site heterogeneity | columns: true treatment effect | %s",
        tstar_note
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
  group_by(method, estimand, tau, beta_trt, heterogeneity,
           numerator, t_star) %>%
  summarise(
    mean_est = mean(estimate),
    truth = first(truth),
    rmse = sqrt(mean(bias^2)),
    bias = mean(bias),
    emp_sd = sd(estimate),
    coverage = if (all(is.na(covered))) NA_real_ else mean(covered, na.rm = TRUE),
    n_reps = n(),
    .groups = "drop"
  ) %>%
  arrange(estimand, method, heterogeneity, beta_trt, tau) %>%
  as.data.frame()

write.csv(summ, "summary_all_scenarios.csv", row.names = FALSE)
message("Wrote summary_all_scenarios.csv")
