#!/usr/bin/env Rscript
# ============================================================================
# Repeated-simulation Fed-CCW covariate-balance diagnostic.
#
# For each configured confounding strength, simulate multiple three-site
# datasets, fit the same site-local initiation models used by Fed-CCW, and
# calculate standardized mean differences (SMDs) between the two protocol-
# compatible groups at the end of the grace period, before and after weighting.
#
# Only federatable summaries are combined centrally: n, sum(w), sum(w^2),
# sum(w*x), and sum(w*x^2) by site, protocol group, and covariate.
# ============================================================================

rm(list = ls())
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
  setwd(dirname(script_path))
}

suppressPackageStartupMessages(library(tidyverse))
source("params.R")
source("DGP_tv.R")
source("Fed_CCW_TVIPCW.R")

parse_args <- function(x) {
  out <- list()
  for (arg in x) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) next
    z <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    out[[z[[1L]]]] <- paste(z[-1L], collapse = "=")
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
get_num <- function(name, default) {
  if (is.null(args[[name]])) default else as.numeric(args[[name]])
}
get_chr <- function(name, default) {
  if (is.null(args[[name]])) default else args[[name]]
}

CFG <- list(
  K = as.integer(get_num("K", 3)),
  n_per_site = as.integer(get_num("n-per-site", 1000)),
  tau = as.integer(get_num("tau", 5)),
  t_star = as.integer(get_num("t-star", 25)),
  beta_trt = get_num("beta-trt", -0.7),
  seed = as.integer(get_num("seed", 2026)),
  n_reps = as.integer(get_num("n-reps", 20)),
  trunc = DEFAULT_TRUNC,
  output_dir = get_chr("output-dir", "balance_diagnostics_repeated")
)

if (CFG$tau < 2L)
  stop("This diagnostic requires tau >= 2 to report lagged L at tau-1.")
if (CFG$tau >= CFG$t_star)
  stop("tau must be smaller than t_star.")
if (CFG$K < 1L || CFG$n_per_site < 10L)
  stop("K and n-per-site must define non-empty analysis sites.")
if (CFG$n_reps < 1L)
  stop("n-reps must be at least 1.")

dir.create(CFG$output_dir, recursive = TRUE, showWarnings = FALSE)

covariate_labels <- c(
  x1 = "Baseline x1",
  L1 = sprintf("L1 at tau=%d", CFG$tau),
  L2 = sprintf("L2 at tau=%d", CFG$tau),
  L1lag = sprintf("L1 at tau-1=%d", CFG$tau - 1L),
  L2lag = sprintf("L2 at tau-1=%d", CFG$tau - 1L)
)

# Site-level balance summaries. At the endpoint, arm 1 contains survivors who
# initiated by tau and arm 0 contains survivors who did not initiate by tau.
# These are exactly the clones still protocol-compatible through the grace
# period. Arm-specific weights are evaluated in interval tau.
local_balance_moments <- function(site_data, tau, t_star, hnum, trunc) {
  P <- .tv_prep(site_data, t_star)
  W <- .tv_weights(P, .tv_init_hazard(P, hnum = hnum),
                   tau = tau, trunc = trunc)

  survived_tau <- P$Tev > tau
  arm <- as.integer(P$S <= tau)
  endpoint_weight <- W$SW0[, tau]
  endpoint_weight[arm == 1L] <- W$SW1[arm == 1L, tau]

  x <- data.frame(
    x1 = P$x1,
    L1 = P$L1[, tau],
    L2 = P$L2[, tau],
    L1lag = P$L1[, tau - 1L],
    L2lag = P$L2[, tau - 1L]
  )

  keep <- survived_tau & is.finite(endpoint_weight)
  x <- x[keep, , drop = FALSE]
  arm <- arm[keep]
  endpoint_weight <- endpoint_weight[keep]

  if (!all(c(0L, 1L) %in% unique(arm)))
    stop("A site has no end-of-grace-period survivors in one protocol group.")

  rows <- list()
  idx <- 0L
  for (analysis in c("Before weighting", "After weighting")) {
    w_all <- if (analysis == "Before weighting") {
      rep(1, nrow(x))
    } else {
      endpoint_weight
    }
    for (g in c(0L, 1L)) {
      use <- arm == g
      wg <- w_all[use]
      for (v in names(x)) {
        xv <- x[[v]][use]
        ok <- is.finite(xv) & is.finite(wg)
        idx <- idx + 1L
        rows[[idx]] <- data.frame(
          analysis = analysis,
          arm = g,
          covariate = v,
          n = sum(ok),
          sw = sum(wg[ok]),
          sw2 = sum(wg[ok]^2),
          swx = sum(wg[ok] * xv[ok]),
          swx2 = sum(wg[ok] * xv[ok]^2),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

central_balance_smd <- function(local_moments) {
  keys <- c("analysis", "arm", "covariate")
  sums <- aggregate(
    local_moments[c("n", "sw", "sw2", "swx", "swx2")],
    local_moments[keys], sum
  )
  sums$mean <- sums$swx / sums$sw
  sums$variance <- pmax(sums$swx2 / sums$sw - sums$mean^2, 0)
  sums$ess <- ifelse(sums$sw2 > 0, sums$sw^2 / sums$sw2, NA_real_)

  rows <- list()
  idx <- 0L
  for (analysis in unique(sums$analysis)) {
    for (v in unique(sums$covariate)) {
      z <- sums[sums$analysis == analysis & sums$covariate == v, ]
      z0 <- z[z$arm == 0L, ]
      z1 <- z[z$arm == 1L, ]
      pooled_sd <- sqrt((z1$variance + z0$variance) / 2)
      smd <- if (is.finite(pooled_sd) && pooled_sd > 0) {
        (z1$mean - z0$mean) / pooled_sd
      } else {
        NA_real_
      }
      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        analysis = analysis,
        covariate = v,
        mean_arm1 = z1$mean,
        mean_arm0 = z0$mean,
        variance_arm1 = z1$variance,
        variance_arm0 = z0$variance,
        smd = smd,
        abs_smd = abs(smd),
        n_arm1 = z1$n,
        n_arm0 = z0$n,
        weight_sum_arm1 = z1$sw,
        weight_sum_arm0 = z0$sw,
        ess_arm1 = z1$ess,
        ess_arm0 = z0$ess,
        stringsAsFactors = FALSE
      )
    }
  }
  list(smd = do.call(rbind, rows), moments = sums)
}

all_smd <- list()
all_groups <- list()
all_hnum <- list()
out_idx <- 0L

scenario_names <- names(DEFAULT_CONF_MULTS)
for (s in seq_along(scenario_names)) {
  scenario <- scenario_names[[s]]
  conf_mult <- unname(DEFAULT_CONF_MULTS[[scenario]])

  for (iteration in seq_len(CFG$n_reps)) {
    # Widely separated deterministic seed blocks make every scenario/iteration
    # reproducible while avoiding accidental overlap in site-specific streams.
    scenario_seed <- CFG$seed + (s - 1L) * 100000L + iteration * 1000L
    message(sprintf(
      "Simulating %s confounding, iteration %d/%d (x%.2f; seed %d).",
      scenario, iteration, CFG$n_reps, conf_mult, scenario_seed
    ))

    dat <- simulate_multisite_tv(
      K = CFG$K, n_per_site = CFG$n_per_site,
      tau = CFG$tau, t_star = CFG$t_star,
      beta_event = set_beta_trt(CFG$beta_trt),
      beta_init = scale_confounding(conf_mult),
      base_seed = scenario_seed
    )
    site_data <- split(dat, dat$site)
    hnum <- central_common_num_hazard(
      lapply(site_data, local_initiation_counts, t_star = CFG$t_star)
    )

    local_moments <- do.call(
      rbind,
      lapply(site_data, local_balance_moments,
             tau = CFG$tau, t_star = CFG$t_star,
             hnum = hnum, trunc = CFG$trunc)
    )
    balance <- central_balance_smd(local_moments)

    out_idx <- out_idx + 1L
    balance$smd$scenario <- scenario
    balance$smd$conf_mult <- conf_mult
    balance$smd$iteration <- iteration
    balance$smd$seed <- scenario_seed
    all_smd[[out_idx]] <- balance$smd

    # One row per weighting status is enough for group size and ESS because
    # these values do not vary across the continuous covariates.
    g <- balance$smd[
      balance$smd$covariate == names(covariate_labels)[1L],
      c("analysis", "n_arm1", "n_arm0",
        "weight_sum_arm1", "weight_sum_arm0",
        "ess_arm1", "ess_arm0")
    ]
    g$scenario <- scenario
    g$conf_mult <- conf_mult
    g$iteration <- iteration
    g$seed <- scenario_seed
    all_groups[[out_idx]] <- g

    all_hnum[[out_idx]] <- data.frame(
      scenario = scenario,
      conf_mult = conf_mult,
      iteration = iteration,
      seed = scenario_seed,
      interval = seq_len(CFG$t_star),
      common_numerator_hazard = hnum,
      stringsAsFactors = FALSE
    )
  }
}

smd_long <- bind_rows(all_smd) |>
  mutate(
    covariate_label = unname(covariate_labels[covariate]),
    scenario = factor(scenario, levels = scenario_names),
    analysis = factor(
      analysis,
      levels = c("Before weighting", "After weighting")
    )
  ) |>
  arrange(scenario, iteration, analysis,
          match(covariate, names(covariate_labels)))

# Retain a paired before/after table for every iteration.
smd_by_iteration <- smd_long |>
  select(
    scenario, conf_mult, iteration, seed, covariate, covariate_label,
    analysis, smd, abs_smd, n_arm1, n_arm0,
    weight_sum_arm1, weight_sum_arm0, ess_arm1, ess_arm0
  ) |>
  mutate(analysis_key = if_else(
    analysis == "Before weighting", "before", "after"
  )) |>
  select(-analysis) |>
  pivot_wider(
    names_from = analysis_key,
    values_from = c(
      smd, abs_smd, weight_sum_arm1, weight_sum_arm0, ess_arm1, ess_arm0
    ),
    names_glue = "{.value}_{analysis_key}"
  ) |>
  arrange(scenario, iteration, match(covariate, names(covariate_labels)))

# Summaries across independent iterations. Mean absolute SMD is the primary
# balance measure; signed mean SMD is retained to show direction.
smd_average_long <- smd_long |>
  group_by(scenario, conf_mult, covariate, covariate_label, analysis) |>
  summarise(
    n_iterations = n(),
    mean_smd = mean(smd),
    sd_smd = sd(smd),
    mean_abs_smd = mean(abs_smd),
    sd_abs_smd = sd(abs_smd),
    median_abs_smd = median(abs_smd),
    q025_abs_smd = quantile(abs_smd, 0.025),
    q975_abs_smd = quantile(abs_smd, 0.975),
    proportion_abs_smd_below_0_10 = mean(abs_smd < 0.10),
    .groups = "drop"
  )

smd_average <- smd_average_long |>
  mutate(analysis_key = if_else(
    analysis == "Before weighting", "before", "after"
  )) |>
  select(-analysis) |>
  pivot_wider(
    names_from = analysis_key,
    values_from = c(
      mean_smd, sd_smd, mean_abs_smd, sd_abs_smd,
      median_abs_smd, q025_abs_smd, q975_abs_smd,
      proportion_abs_smd_below_0_10
    ),
    names_glue = "{.value}_{analysis_key}"
  ) |>
  mutate(
    mean_abs_smd_reduction =
      mean_abs_smd_before - mean_abs_smd_after
  ) |>
  arrange(scenario, match(covariate, names(covariate_labels)))

group_by_iteration <- bind_rows(all_groups) |>
  mutate(
    scenario = factor(scenario, levels = scenario_names),
    analysis = factor(
      analysis,
      levels = c("Before weighting", "After weighting")
    )
  ) |>
  arrange(scenario, iteration, analysis)

group_average <- group_by_iteration |>
  group_by(scenario, conf_mult, analysis) |>
  summarise(
    n_iterations = n(),
    mean_n_arm1 = mean(n_arm1),
    mean_n_arm0 = mean(n_arm0),
    mean_weight_sum_arm1 = mean(weight_sum_arm1),
    mean_weight_sum_arm0 = mean(weight_sum_arm0),
    mean_ess_arm1 = mean(ess_arm1),
    sd_ess_arm1 = sd(ess_arm1),
    mean_ess_arm0 = mean(ess_arm0),
    sd_ess_arm0 = sd(ess_arm0),
    .groups = "drop"
  )

hnum_summary <- bind_rows(all_hnum) |>
  mutate(scenario = factor(scenario, levels = scenario_names)) |>
  arrange(scenario, iteration, interval)

replicate_overall <- smd_long |>
  group_by(scenario, conf_mult, analysis, iteration) |>
  summarise(
    mean_abs_smd = mean(abs_smd),
    max_abs_smd = max(abs_smd),
    .groups = "drop"
  )

overall <- replicate_overall |>
  group_by(scenario, conf_mult, analysis) |>
  summarise(
    n_iterations = n(),
    average_mean_abs_smd = mean(mean_abs_smd),
    sd_mean_abs_smd = sd(mean_abs_smd),
    average_max_abs_smd = mean(max_abs_smd),
    q025_max_abs_smd = quantile(max_abs_smd, 0.025),
    q975_max_abs_smd = quantile(max_abs_smd, 0.975),
    .groups = "drop"
  )

write.csv(
  smd_by_iteration,
  file.path(CFG$output_dir, "balance_smd_by_iteration.csv"),
  row.names = FALSE
)
write.csv(
  smd_average,
  file.path(CFG$output_dir, "balance_smd_average.csv"),
  row.names = FALSE
)
write.csv(
  group_by_iteration,
  file.path(CFG$output_dir, "balance_group_ess_by_iteration.csv"),
  row.names = FALSE
)
write.csv(
  group_average,
  file.path(CFG$output_dir, "balance_group_ess_average.csv"),
  row.names = FALSE
)
write.csv(
  hnum_summary,
  file.path(CFG$output_dir, "common_numerator_hazard_by_iteration.csv"),
  row.names = FALSE
)
write.csv(
  overall,
  file.path(CFG$output_dir, "balance_overall_average.csv"),
  row.names = FALSE
)

plot_reps <- smd_long |>
  mutate(
    scenario_label = factor(
      as.character(scenario), levels = scenario_names,
      labels = paste(tools::toTitleCase(scenario_names), "confounding")
    ),
    covariate_label = factor(
      covariate_label,
      levels = rev(unname(covariate_labels))
    )
  )

plot_avg <- smd_average_long |>
  mutate(
    scenario_label = factor(
      as.character(scenario), levels = scenario_names,
      labels = paste(tools::toTitleCase(scenario_names), "confounding")
    ),
    covariate_label = factor(
      covariate_label,
      levels = rev(unname(covariate_labels))
    )
  )

p <- ggplot(
  plot_avg,
  aes(x = mean_abs_smd, y = covariate_label,
      color = analysis, shape = analysis)
) +
  geom_vline(
    xintercept = 0.10, linetype = "dashed",
    linewidth = 0.5, color = "grey45"
  ) +
  geom_point(
    data = plot_reps,
    aes(x = abs_smd, y = covariate_label,
        color = analysis, shape = analysis),
    position = position_jitter(height = 0.06, width = 0),
    size = 1.2, alpha = 0.16, show.legend = FALSE
  ) +
  geom_line(
    aes(group = covariate_label),
    color = "grey70", linewidth = 0.45
  ) +
  geom_point(size = 3.1) +
  facet_wrap(~ scenario_label, nrow = 1) +
  scale_color_manual(values = c(
    "Before weighting" = "#d95f02",
    "After weighting" = "#1b9e77"
  )) +
  labs(
    title = "Average Fed-CCW balance at the end of the grace period",
    subtitle = sprintf(
      "%d independent simulations/scenario; tau=%d, beta_trt=%.1f, %d sites x %d patients",
      CFG$n_reps, CFG$tau, CFG$beta_trt, CFG$K, CFG$n_per_site
    ),
    x = "Absolute standardized mean difference",
    y = NULL,
    color = NULL,
    shape = NULL,
    caption = paste(
      "Large symbols are means; faint symbols are individual iterations.",
      "Dashed line marks |SMD| = 0.10."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA)
  )

ggsave(
  file.path(CFG$output_dir, "balance_love_plot_average.png"),
  p, width = 14, height = 6.5, dpi = 180
)

message("Average balance summary across continuous covariates:")
print(overall, row.names = FALSE)
message("Wrote repeated balance diagnostics to: ",
        normalizePath(CFG$output_dir))
