#!/usr/bin/env Rscript

# Weight diagnostics for ArraySim5-4 Combined Site Heterogeneity.
#
# Default design:
#   * 10 replicates per cell
#   * low, moderate, and high combined heterogeneity
#   * tau = 5, 10, and 15
#   * beta_trt = -0.7 (weights do not directly model beta_trt, although the
#     event process changes who remains at risk)
#   * no truncation, matching the production ArraySim5-4 analysis
#   * both the Fed/site-specific nuisance model and the pooled site-FE model
#     with shared slopes

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(file_arg)) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
script_dir <- dirname(script_path)
sim_dir <- dirname(script_dir)

parse_args <- function(x) {
  out <- list(
    n_reps = 10L,
    taus = c(5L, 10L, 15L),
    heterogeneity = c("low", "moderate", "high"),
    beta_trt = -0.7,
    base_seed = 54001L,
    include_pooled = TRUE,
    output_dir = file.path(script_dir, "outputs")
  )
  for (arg in x) {
    if (!grepl("^--[^=]+=", arg)) stop("Arguments must use --name=value: ", arg)
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", pieces[[1L]], fixed = TRUE)
    value <- paste(pieces[-1L], collapse = "=")
    if (key == "n_reps") out$n_reps <- as.integer(value)
    else if (key == "taus") out$taus <- as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])
    else if (key == "heterogeneity") out$heterogeneity <- strsplit(value, ",", fixed = TRUE)[[1L]]
    else if (key == "beta_trt") out$beta_trt <- as.numeric(value)
    else if (key == "base_seed") out$base_seed <- as.integer(value)
    else if (key == "include_pooled") out$include_pooled <- tolower(value) %in% c("true", "t", "1", "yes")
    else if (key == "output_dir") out$output_dir <- value
    else stop("Unknown argument: --", pieces[[1L]])
  }
  out
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
if (!is.finite(cfg$n_reps) || cfg$n_reps < 1L) stop("n_reps must be positive.")
if (!is.finite(cfg$beta_trt)) stop("beta_trt must be finite.")
if (!grepl("^(/|[A-Za-z]:[/\\\\])", cfg$output_dir))
  cfg$output_dir <- file.path(script_dir, cfg$output_dir)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(sim_dir)
source("params.R")
source("DGP_tv.R")
source("Fed_CCW_TVIPCW.R")
validate_params()

if (!all(cfg$taus %in% SIM_TAU_VALUES))
  stop("taus must be selected from: ", paste(SIM_TAU_VALUES, collapse = ", "))
if (!all(cfg$heterogeneity %in% names(HETEROGENEITY_SETTINGS)))
  stop("heterogeneity must be selected from: ",
       paste(names(HETEROGENEITY_SETTINGS), collapse = ", "))

dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

q_safe <- function(x, p) {
  if (!length(x) || !any(is.finite(x))) return(NA_real_)
  unname(quantile(x[is.finite(x)], p, names = FALSE, type = 8))
}

summarize_weights <- function(P, Wt, tau, model, heterogeneity, replicate) {
  Cs <- .clone_censoring(P$S, tau, P$M)
  groups <- c("Overall", sort(unique(as.character(P$site))))
  out <- vector("list", length = 2L * tau * length(groups))
  z <- 0L

  for (arm in c("Initiate by tau", "Do not initiate by tau")) {
    W <- if (arm == "Initiate by tau") Wt$SW1 else Wt$SW0
    Cg <- if (arm == "Initiate by tau") Cs$C1 else Cs$C0
    for (m in seq_len(tau)) {
      risk <- P$Tcap > (m - 1L) & Cg > (m - 1L)
      for (g in groups) {
        use <- risk & (g == "Overall" | as.character(P$site) == g)
        w <- W[use, m]
        n <- length(w)
        ess <- if (n && sum(w^2) > 0) sum(w)^2 / sum(w^2) else NA_real_
        z <- z + 1L
        out[[z]] <- data.frame(
          model = model, heterogeneity = heterogeneity,
          replicate = replicate, tau = tau, site = g,
          strategy = arm, interval = m, risk_n = n,
          weight_sum = sum(w), weight_mean = if (n) mean(w) else NA_real_,
          weight_sd = if (n > 1L) sd(w) else NA_real_,
          weight_cv = if (n > 1L && mean(w) > 0) sd(w) / mean(w) else NA_real_,
          weight_min = if (n) min(w) else NA_real_,
          weight_q01 = q_safe(w, 0.01), weight_q05 = q_safe(w, 0.05),
          weight_q50 = q_safe(w, 0.50), weight_q95 = q_safe(w, 0.95),
          weight_q99 = q_safe(w, 0.99), weight_max = if (n) max(w) else NA_real_,
          prop_weight_gt_5 = if (n) mean(w > 5) else NA_real_,
          prop_weight_gt_10 = if (n) mean(w > 10) else NA_real_,
          prop_weight_gt_20 = if (n) mean(w > 20) else NA_real_,
          ess = ess, ess_fraction = if (n) ess / n else NA_real_,
          design_effect = if (is.finite(ess) && ess > 0) n / ess else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  bind_rows(out[seq_len(z)])
}

summarize_propensity <- function(P, H, tau, model, heterogeneity, replicate,
                                 hnum) {
  groups <- c("Overall", sort(unique(as.character(P$site))))
  out <- vector("list", tau * length(groups))
  z <- 0L
  for (m in seq_len(tau)) {
    eligible <- !is.na(H$Hden[, m])
    y <- as.integer(P$S == m)
    for (g in groups) {
      use <- eligible & (g == "Overall" | as.character(P$site) == g)
      p <- H$Hden[use, m]
      yy <- y[use]
      p_observed <- ifelse(yy == 1L, p, 1 - p)
      n <- length(p)
      z <- z + 1L
      out[[z]] <- data.frame(
        model = model, heterogeneity = heterogeneity,
        replicate = replicate, tau = tau, site = g, interval = m,
        eligible_n = n, initiations = sum(yy),
        observed_initiation_rate = if (n) mean(yy) else NA_real_,
        common_numerator_hazard = hnum[[m]],
        denominator_mean = if (n) mean(p) else NA_real_,
        denominator_min = if (n) min(p) else NA_real_,
        denominator_q01 = q_safe(p, 0.01),
        denominator_q05 = q_safe(p, 0.05),
        denominator_q50 = q_safe(p, 0.50),
        denominator_q95 = q_safe(p, 0.95),
        denominator_q99 = q_safe(p, 0.99),
        denominator_max = if (n) max(p) else NA_real_,
        prop_den_lt_001 = if (n) mean(p < 0.01) else NA_real_,
        prop_den_gt_099 = if (n) mean(p > 0.99) else NA_real_,
        observed_prob_min = if (n) min(p_observed) else NA_real_,
        observed_prob_q01 = q_safe(p_observed, 0.01),
        observed_prob_q05 = q_safe(p_observed, 0.05),
        prop_observed_prob_lt_001 = if (n) mean(p_observed < 0.01) else NA_real_,
        prop_observed_prob_lt_005 = if (n) mean(p_observed < 0.05) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(out[seq_len(z)])
}

fit_with_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(list(message = conditionMessage(e)), class = "diag_error")),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings))
}

weight_rows <- list()
propensity_rows <- list()
prevalence_rows <- list()
log_rows <- list()
row_id <- 0L
log_id <- 0L

grid <- expand.grid(
  heterogeneity = cfg$heterogeneity,
  tau = cfg$taus,
  replicate = seq_len(cfg$n_reps),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)

for (i in seq_len(nrow(grid))) {
  level <- grid$heterogeneity[[i]]
  tau <- grid$tau[[i]]
  rep_id <- grid$replicate[[i]]
  setting <- get_heterogeneity_setting(level)
  cell_seed <- cfg$base_seed + match(level, names(HETEROGENEITY_SETTINGS)) * 100000L +
    tau * 1000L + rep_id * 10L

  message(sprintf("[%d/%d] heterogeneity=%s tau=%d replicate=%d",
                  i, nrow(grid), level, tau, rep_id))

  dat <- simulate_multisite_tv(
    K = length(setting$site_sizes), n_per_site = setting$site_sizes,
    init_intercepts = setting$init_intercepts,
    patient_mix = setting$patient_mix,
    tau = tau, t_star = SIM_TSTAR, base_seed = cell_seed,
    beta_event = set_beta_trt(cfg$beta_trt),
    beta_init = DEFAULT_BETA_INIT
  )
  hnum <- .tv_common_num_hazard(dat, SIM_TSTAR)

  prev_site <- dat %>%
    group_by(site) %>%
    summarise(n = n(), initiated_by_tau = sum(A_tau),
              initiation_prevalence = mean(A_tau), .groups = "drop") %>%
    mutate(heterogeneity = level, tau = tau, replicate = rep_id,
           beta_trt = cfg$beta_trt, .before = 1)
  prevalence_rows[[length(prevalence_rows) + 1L]] <- prev_site

  site_blocks <- split(dat, dat$site)
  fed_parts <- list()
  for (site_name in names(site_blocks)) {
    d <- site_blocks[[site_name]]
    P <- .tv_prep(d, SIM_TSTAR)
    fitted <- fit_with_warnings(.tv_init_hazard(P, hnum = hnum))
    log_id <- log_id + 1L
    failed <- inherits(fitted$value, "diag_error")
    log_rows[[log_id]] <- data.frame(
      heterogeneity = level, tau = tau, replicate = rep_id,
      model = "Federated: site-specific slopes", site = site_name,
      status = if (failed) "error" else if (length(fitted$warnings)) "warning" else "ok",
      message = if (failed) fitted$value$message else paste(fitted$warnings, collapse = " | "),
      stringsAsFactors = FALSE
    )
    if (failed) next
    H <- fitted$value
    Wt <- .tv_weights(P, H, tau = tau, trunc = DEFAULT_TRUNC)
    fed_parts[[site_name]] <- list(P = P, H = H, Wt = Wt)
    row_id <- row_id + 1L
    weight_rows[[row_id]] <- summarize_weights(
      P, Wt, tau, "Federated: site-specific slopes", level, rep_id
    ) %>% filter(site != "Overall")
    propensity_rows[[length(propensity_rows) + 1L]] <- summarize_propensity(
      P, H, tau, "Federated: site-specific slopes", level, rep_id, hnum
    ) %>% filter(site != "Overall")
  }

  # Reconstruct overall Fed summaries by combining the already calculated
  # site-specific weight matrices, preserving the exact local fits.
  fed_P <- .tv_prep(dat, SIM_TSTAR)
  fed_Hden <- matrix(NA_real_, nrow(dat), SIM_TSTAR)
  fed_Hnum <- matrix(NA_real_, nrow(dat), SIM_TSTAR)
  fed_SW1 <- matrix(NA_real_, nrow(dat), SIM_TSTAR)
  fed_SW0 <- matrix(NA_real_, nrow(dat), SIM_TSTAR)
  start <- 1L
  fed_complete <- length(fed_parts) == length(site_blocks)
  for (site_name in names(site_blocks)) {
    d <- site_blocks[[site_name]]
    idx <- start:(start + nrow(d) - 1L)
    start <- max(idx) + 1L
    if (!fed_complete) break
    part <- fed_parts[[site_name]]
    fed_Hden[idx, ] <- part$H$Hden; fed_Hnum[idx, ] <- part$H$Hnum
    fed_SW1[idx, ] <- part$Wt$SW1; fed_SW0[idx, ] <- part$Wt$SW0
  }
  if (fed_complete) {
    weight_rows[[length(weight_rows) + 1L]] <- summarize_weights(
      fed_P, list(SW1 = fed_SW1, SW0 = fed_SW0), tau,
      "Federated: site-specific slopes", level, rep_id
    ) %>% filter(site == "Overall")
    propensity_rows[[length(propensity_rows) + 1L]] <- summarize_propensity(
      fed_P, list(Hden = fed_Hden, Hnum = fed_Hnum), tau,
      "Federated: site-specific slopes", level, rep_id, hnum
    ) %>% filter(site == "Overall")
  }

  if (cfg$include_pooled) {
    P <- .tv_prep(dat, SIM_TSTAR)
    fitted <- fit_with_warnings(.tv_init_hazard(P, hnum = hnum))
    log_id <- log_id + 1L
    failed <- inherits(fitted$value, "diag_error")
    log_rows[[log_id]] <- data.frame(
      heterogeneity = level, tau = tau, replicate = rep_id,
      model = "Pooled site FE: shared slopes", site = "Overall",
      status = if (failed) "error" else if (length(fitted$warnings)) "warning" else "ok",
      message = if (failed) fitted$value$message else paste(fitted$warnings, collapse = " | "),
      stringsAsFactors = FALSE
    )
    if (!failed) {
      H <- fitted$value
      Wt <- .tv_weights(P, H, tau = tau, trunc = DEFAULT_TRUNC)
      weight_rows[[length(weight_rows) + 1L]] <- summarize_weights(
        P, Wt, tau, "Pooled site FE: shared slopes", level, rep_id
      )
      propensity_rows[[length(propensity_rows) + 1L]] <- summarize_propensity(
        P, H, tau, "Pooled site FE: shared slopes", level, rep_id, hnum
      )
    }
  }
}

weights <- bind_rows(weight_rows)
propensity <- bind_rows(propensity_rows)
prevalence <- bind_rows(prevalence_rows)
model_log <- bind_rows(log_rows)

level_order <- c("low", "moderate", "high")
weights$heterogeneity <- factor(weights$heterogeneity, level_order)
propensity$heterogeneity <- factor(propensity$heterogeneity, level_order)
prevalence$heterogeneity <- factor(prevalence$heterogeneity, level_order)

weight_summary <- weights %>%
  group_by(model, heterogeneity, tau, site, strategy, interval) %>%
  summarise(
    replicates = n_distinct(replicate), mean_risk_n = mean(risk_n),
    mean_weight = mean(weight_mean, na.rm = TRUE),
    mean_q95 = mean(weight_q95, na.rm = TRUE),
    mean_q99 = mean(weight_q99, na.rm = TRUE),
    median_max = median(weight_max, na.rm = TRUE),
    max_across_reps = max(weight_max, na.rm = TRUE),
    mean_ess = mean(ess, na.rm = TRUE),
    mean_ess_fraction = mean(ess_fraction, na.rm = TRUE),
    p10_ess_fraction = q_safe(ess_fraction, 0.10),
    mean_design_effect = mean(design_effect, na.rm = TRUE),
    mean_prop_weight_gt_10 = mean(prop_weight_gt_10, na.rm = TRUE),
    .groups = "drop"
  )

propensity_summary <- propensity %>%
  group_by(model, heterogeneity, tau, site, interval) %>%
  summarise(
    replicates = n_distinct(replicate), mean_eligible_n = mean(eligible_n),
    mean_observed_rate = mean(observed_initiation_rate, na.rm = TRUE),
    mean_numerator_hazard = mean(common_numerator_hazard, na.rm = TRUE),
    mean_denominator = mean(denominator_mean, na.rm = TRUE),
    mean_den_q01 = mean(denominator_q01, na.rm = TRUE),
    mean_den_q99 = mean(denominator_q99, na.rm = TRUE),
    mean_observed_prob_q01 = mean(observed_prob_q01, na.rm = TRUE),
    min_observed_prob = min(observed_prob_min, na.rm = TRUE),
    mean_prop_observed_prob_lt_001 = mean(prop_observed_prob_lt_001, na.rm = TRUE),
    mean_prop_observed_prob_lt_005 = mean(prop_observed_prob_lt_005, na.rm = TRUE),
    .groups = "drop"
  )

grace_summary <- weight_summary %>% filter(interval == tau)

write.csv(weights, file.path(cfg$output_dir, "weight_by_rep_site_interval.csv"), row.names = FALSE)
write.csv(weight_summary, file.path(cfg$output_dir, "weight_summary_across_reps.csv"), row.names = FALSE)
write.csv(grace_summary, file.path(cfg$output_dir, "weight_summary_at_tau.csv"), row.names = FALSE)
write.csv(propensity, file.path(cfg$output_dir, "propensity_by_rep_site_interval.csv"), row.names = FALSE)
write.csv(propensity_summary, file.path(cfg$output_dir, "propensity_summary_across_reps.csv"), row.names = FALSE)
write.csv(prevalence, file.path(cfg$output_dir, "initiation_prevalence_by_rep_site.csv"), row.names = FALSE)
write.csv(model_log, file.path(cfg$output_dir, "model_fit_log.csv"), row.names = FALSE)

plot_theme <- theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

overall_w <- weights %>% filter(site == "Overall")

p_ess <- overall_w %>%
  group_by(model, heterogeneity, tau, strategy, interval) %>%
  summarise(mean = mean(ess_fraction, na.rm = TRUE),
            lo = q_safe(ess_fraction, 0.10), hi = q_safe(ess_fraction, 0.90),
            .groups = "drop") %>%
  ggplot(aes(interval, mean, color = model, linetype = strategy,
             group = interaction(model, strategy))) +
  geom_hline(yintercept = 0.20, color = "grey50", linetype = "dotted") +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model, linetype = NULL),
              alpha = 0.10, color = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.1) +
  facet_grid(heterogeneity ~ tau, labeller = label_both) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Effective sample-size fraction over the grace period",
       subtitle = sprintf("Lines are means across %d replicates; ribbons are 10th-90th percentiles",
                          cfg$n_reps),
       x = "Interval", y = "ESS / unweighted risk-set size",
       color = "Weight model", fill = "Weight model", linetype = "Strategy") +
  plot_theme
ggsave(file.path(cfg$output_dir, "ess_fraction_trajectory.png"), p_ess,
       width = 14, height = 9, dpi = 180)

p_q99 <- overall_w %>%
  group_by(model, heterogeneity, tau, strategy, interval) %>%
  summarise(q99 = mean(weight_q99, na.rm = TRUE),
            max_median = median(weight_max, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(interval, q99, color = model, linetype = strategy,
             group = interaction(model, strategy))) +
  geom_hline(yintercept = 1, color = "grey65") +
  geom_line(linewidth = 0.8) + geom_point(size = 1.1) +
  facet_grid(heterogeneity ~ tau, labeller = label_both, scales = "free_y") +
  scale_y_log10() +
  labs(title = "Upper tail of the stabilized IPCW distribution",
       subtitle = "Mean 99th percentile across replicates; logarithmic vertical scale",
       x = "Interval", y = "99th percentile of weight",
       color = "Weight model", linetype = "Strategy") + plot_theme
ggsave(file.path(cfg$output_dir, "weight_q99_trajectory.png"), p_q99,
       width = 14, height = 9, dpi = 180)

p_site_ess <- weights %>%
  filter(site != "Overall", interval == tau) %>%
  ggplot(aes(factor(site), ess_fraction, fill = model)) +
  geom_hline(yintercept = 0.20, color = "grey50", linetype = "dotted") +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.alpha = 0.35) +
  facet_grid(interaction(heterogeneity, strategy, sep = ": ") ~ tau,
             labeller = label_both) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Site-level ESS fraction at the end of the grace period",
       x = "Site", y = "ESS / unweighted risk-set size", fill = "Weight model") +
  plot_theme
ggsave(file.path(cfg$output_dir, "site_ess_fraction_at_tau.png"), p_site_ess,
       width = 14, height = 11, dpi = 180)

p_support <- propensity %>%
  filter(site == "Overall") %>%
  group_by(model, heterogeneity, tau, interval) %>%
  summarise(q01 = mean(observed_prob_q01, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(interval, q01, color = model)) +
  geom_hline(yintercept = 0.01, color = "grey50", linetype = "dotted") +
  geom_line(linewidth = 0.8) + geom_point(size = 1.1) +
  facet_grid(heterogeneity ~ tau, labeller = label_both, scales = "free_y") +
  scale_y_log10() +
  labs(title = "Support for the observed initiation decision",
       subtitle = "Mean first percentile of P(observed initiation decision | history)",
       x = "Interval", y = "First percentile (log scale)", color = "Weight model") +
  plot_theme
ggsave(file.path(cfg$output_dir, "observed_probability_support.png"), p_support,
       width = 14, height = 9, dpi = 180)

p_prev <- prevalence %>%
  ggplot(aes(factor(site), initiation_prevalence,
             group = interaction(heterogeneity, tau, site))) +
  geom_boxplot(aes(fill = factor(site)), outlier.alpha = 0.3) +
  facet_grid(heterogeneity ~ tau, labeller = label_both) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Observed treatment initiation by the end of the grace period",
       x = "Site", y = "Initiation prevalence", fill = "Site") +
  plot_theme
ggsave(file.path(cfg$output_dir, "initiation_prevalence.png"), p_prev,
       width = 13, height = 8, dpi = 180)

run_info <- c(
  paste0("Run date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Replicates per cell: ", cfg$n_reps),
  paste0("Heterogeneity: ", paste(cfg$heterogeneity, collapse = ", ")),
  paste0("Tau: ", paste(cfg$taus, collapse = ", ")),
  paste0("beta_trt: ", cfg$beta_trt),
  paste0("Base seed: ", cfg$base_seed),
  paste0("Weight truncation: ", paste(DEFAULT_TRUNC, collapse = ", ")),
  paste0("Included pooled shared-slope model: ", cfg$include_pooled),
  paste0("Output directory: ", normalizePath(cfg$output_dir))
)
writeLines(run_info, file.path(cfg$output_dir, "run_info.txt"))

cat("\nWeight diagnostics complete.\n")
cat("Outputs:", normalizePath(cfg$output_dir), "\n")
cat("Model fit status:\n")
print(model_log %>% count(model, status), row.names = FALSE)
