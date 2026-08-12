#!/usr/bin/env Rscript
# Generate the independent numerator once, before cluster simulation.

rm(list = ls())
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
  setwd(dirname(script_path))
}

source("external_reference.R")

parse_cli <- function(x) {
  out <- list()
  for (arg in x) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) next
    z <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    out[[z[1L]]] <- paste(z[-1L], collapse = "=")
  }
  out
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
N <- if (!is.null(cli$N)) as.numeric(cli$N) else EXTERNAL_REFERENCE_N
t_star <- if (!is.null(cli[["t-star"]])) {
  as.integer(cli[["t-star"]])
} else EXTERNAL_REFERENCE_TSTAR
heterogeneity <- if (!is.null(cli$heterogeneity)) {
  cli$heterogeneity
} else EXTERNAL_REFERENCE_HETEROGENEITY
seed <- if (!is.null(cli$seed)) as.integer(cli$seed) else EXTERNAL_REFERENCE_SEED
out_file <- if (!is.null(cli$output)) cli$output else EXTERNAL_REFERENCE_FILE

message(sprintf(
  "Generating external reference: N=%g, t*=%d, heterogeneity=%s, seed=%d",
  N, t_star, heterogeneity, seed
))
ref <- simulate_external_reference_numerator(
  N = N, t_star = t_star,
  heterogeneity = heterogeneity, seed = seed
)
validate_external_reference(ref, t_star = t_star)
saveRDS(ref, out_file)

summary <- data.frame(
  interval = seq_len(ref$t_star),
  eligible = ref$eligible,
  initiated = ref$initiated,
  numerator_hazard = ref$numerator_hazard
)
write.csv(summary, "external_reference_numerator.csv", row.names = FALSE)
message("Wrote ", out_file)
message("Wrote external_reference_numerator.csv")
