# ArraySim4-3: Ten Small Sites (fixed)

Corrected finite-sample stress test of Fed-CCW versus completing CCW locally
and averaging the ten resulting survival curves.

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

- 10 sites
- 100/300 patients per site
- total sample size: 1,000/3,000 per replicate
- site-specific treatment-initiation intercepts evenly spaced from `-4.5`
  through `-1.5`
- corresponding zero-covariate initiation probabilities range from about
  1.1% to 18.2% per interval
- every site contributes 10% of the target population

All sites share the treatment-initiation covariate slopes, covariate-transition
model, outcome model, and treatment effect. Confounding multipliers are small
(`0.5`), medium (`1.0`), and strong (`1.5`).

## Methods

1. Federated CCW
2. Pooled CCW with site fixed effects in the initiation model
3. Federated IPW without cloning
4. Federated unweighted per-protocol analysis
5. Federated landmark IPW
6. Local CCW with sample-size-weighted survival-curve aggregation

Fed-CCW combines weighted event/risk totals before constructing survival.
Local CCW constructs ten curves first and then averages them. Both use the
same common numerator and the same local denominator models.

Landmark IPW targets post-landmark survival conditional on being event-free
through `tau`; it does not target the same time-zero grace-period estimand.

## Simulation grid

- `tau`: 5, 10, 15
- `beta_trt`: -1.0, -0.7, -0.5
- confounding multipliers: 0.5, 1.0, 1.5
- fixed `t_star`: 25
- 300 replicates per cell
- 3,000,000-person mixture-specific oracle per cell
- 27 array tasks
- no weight truncation (`c(0, 1)`)
- 14,400 rows per task: 300 replicates x 6 methods x 8 estimands

The oracle allocates its 3,000,000 observations equally across the ten sites.

## Validate, run, and aggregate

From this directory:

```bash
Rscript validate_setup.R
sbatch job.sh
```

For a local run:

```bash
for task_id in {1..27}; do
  SLURM_CPUS_PER_TASK=8 Rscript run_sim_array.R "$task_id" \
    > "logs/local_${task_id}.out" \
    2> "logs/local_${task_id}.err"
done
```

`job.sh` requests 5 CPUs and 20 GB per task. No `submit.sh` wrapper is
included, so `logs/` and `results/` must exist before submission.

After all 27 tasks complete:

```bash
Rscript aggregate.R
```

Aggregation calculates RMSE from replicate-level errors and writes
`summary_all_scenarios.csv` plus bias plots for RD, RR, OR, and RMST difference.
Only federated CCW and pooled CCW report model-based confidence intervals.
