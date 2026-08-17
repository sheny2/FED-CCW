# ArraySim5-1: Five Methods (fixed DGP)

Corrected finite-sample stress test comparing five estimators across five
sites. The local CCW plus survival-curve meta-analysis approach is intentionally
excluded from this simulation.

This folder starts clean: it intentionally contains no simulation results,
summary tables, plots, or completed logs from the original study. Those
outputs were generated under a DGP that allowed treatment initiation and
covariate evolution after a terminal event and must not be reused here.

## DGP correction

After a terminal event, a participant:

- is no longer eligible to initiate treatment;
- contributes no later treatment-initiation person-interval rows;
- has no later time-varying covariates generated; and
- cannot have `S` or `A_tau` changed by a post-event initiation.

The observed-data estimator and the Monte Carlo oracle use the same event-free
initiation risk-set definition. `validate_setup.R` explicitly checks that no
participant initiates after the terminal event and that local initiation
counts equal counts reconstructed directly from event-free risk sets.

## Site design

- 5 sites
- 1,000 patients per site
- total sample size: 5,000 per replicate
- site-specific treatment-initiation intercepts evenly spaced from `-4.5`
  through `-1.5`
- corresponding zero-covariate initiation probabilities range from about
  1.1% to 18.2% per interval
- every site contributes 20% of the target population

All sites share the treatment-initiation covariate slopes, covariate-transition
model, outcome model, and treatment effect. Confounding multipliers are small
(`0.5`), medium (`1.0`), and strong (`1.5`).

## Methods

1. Federated CCW
2. Pooled CCW with site fixed effects in the initiation model
3. Federated IPW without cloning
4. Federated unweighted per-protocol analysis
5. Federated landmark IPW

Fed-CCW combines weighted event/risk totals before constructing survival.

Landmark IPW targets post-landmark survival conditional on being event-free
through `tau`; it does not target the same time-zero grace-period estimand.

## Simulation grid

- `tau`: 5, 10, 15
- `beta_trt`: -1.0, -0.7, -0.5
- confounding multipliers: 0.5, 1.0, 1.5
- fixed `t_star`: 25
- 100 replicates per cell
- 3,000,000-person mixture-specific oracle per cell
- 27 array tasks
- no weight truncation (`c(0, 1)`)
- 4,000 rows per task: 100 replicates x 5 methods x 8 estimands

The oracle allocates its 3,000,000 observations equally across the five sites.

## Configuration

All simulation inputs are defined in `params.R`. Users can change site sizes,
site-specific initiation intercepts, the tau and treatment-effect grids,
confounding multipliers, replicate count, follow-up horizon, weight truncation,
oracle size, seeds, and all DGP coefficients there. The number of sites and the
SLURM array size are derived automatically from those settings.

## Validate, run, and aggregate

From this directory:

```bash
Rscript validate_setup.R
bash submit.sh
```

For a local run:

```bash
n_tasks=$(Rscript -e 'source("params.R"); cat(nrow(simulation_grid()))')
for task_id in $(seq 1 "$n_tasks"); do
  SLURM_CPUS_PER_TASK=7 Rscript run_sim_array.R "$task_id" \
    > "logs/local_${task_id}.out" \
    2> "logs/local_${task_id}.err"
done
```

`job.sh` requests 5 CPUs and 20 GB per task. `submit.sh` calculates the array
range from `params.R`, so grid changes do not require editing the SLURM script.

After all 27 tasks complete:

```bash
Rscript aggregate.R
```

Aggregation calculates RMSE from replicate-level errors and writes
`summary_all_scenarios.csv` plus bias plots for RD, RR, OR, and RMST difference.
Only federated CCW and pooled CCW report model-based confidence intervals.
