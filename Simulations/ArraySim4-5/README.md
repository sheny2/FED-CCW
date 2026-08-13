# ArraySim4-5: Sample Size and Outcome Frequency

Factorial simulation designed to determine whether the behavior of federated
landmark IPW and local CCW curve meta-analysis changes with finite site sample
size and rare versus common terminal outcomes.

This is the outcome-frequency analogue of the executed ArraySim4-4 design. It
uses the corrected event-free DGP: treatment initiation and time-varying
covariate evolution stop after a terminal event, and initiation models and
numerator counts use event-free treatment risk sets.

## Fixed settings

- 10 equally sized sites
- medium confounding (`conf_mult = 1.0`)
- treatment effect `beta_trt = -0.7`
- treatment-initiation intercept `-3.0` at every site
- follow-up horizon `t_star = 25`
- no weight truncation (`c(0, 1)`)
- 300 replicates per cell
- 5,000,000-person oracle per cell

## Factorial design

| Factor | Levels |
|---|---|
| Grace period | `tau = 4, 6, 8` |
| Patients per site | low = 300; large = 1,000 |
| Outcome frequency | rare event intercept = `-7`; common event intercept = `-3` |

This is a `3 x 2 x 2` design with 12 array tasks. Total sample size is 3,000
in the low-sample cells and 10,000 in the large-sample cells.

Pilot calibration under the fixed initiation and confounding model produced
approximately 6% cumulative mortality by `t_star=25` for intercept `-7` and
approximately 74% for intercept `-3`. Achieved mortality, initiation, landmark
sample size, and post-landmark event counts are recorded for every replicate.

## Methods

1. Federated CCW
2. Pooled CCW
3. Federated IPW without cloning
4. Federated unweighted per-protocol analysis
5. Federated landmark IPW
6. Local CCW with sample-size-weighted survival-curve aggregation

Landmark IPW targets post-landmark survival conditional on being event-free
through `tau`; it does not target the common time-zero CCW oracle. Its reported
difference from that oracle includes an estimand difference.

Each task writes 14,400 rows:

```text
300 replicates x 6 methods x 8 estimands
```

## Task mapping

| Task | tau | Patients/site | Outcome intercept |
|---:|---:|---:|---:|
| 1 | 4 | 300 | -7 |
| 2 | 6 | 300 | -7 |
| 3 | 8 | 300 | -7 |
| 4 | 4 | 1,000 | -7 |
| 5 | 6 | 1,000 | -7 |
| 6 | 8 | 1,000 | -7 |
| 7 | 4 | 300 | -3 |
| 8 | 6 | 300 | -3 |
| 9 | 8 | 300 | -3 |
| 10 | 4 | 1,000 | -3 |
| 11 | 6 | 1,000 | -3 |
| 12 | 8 | 1,000 | -3 |

Pairs differing only in finite sample size share the same oracle target and
recompute it using the same Monte Carlo size and seed.

## Validate and run

```bash
Rscript validate_setup.R
sbatch job.sh
```

`job.sh` requests 10 CPUs and 40 GB per task. For a local run:

```bash
for task_id in {1..12}; do
  SLURM_CPUS_PER_TASK=8 Rscript run_sim_array.R "$task_id" \
    > "logs/local_${task_id}.out" \
    2> "logs/local_${task_id}.err"
done
```

## Aggregate

After all 12 tasks finish:

```bash
Rscript aggregate.R
```

Outputs include:

- `summary_all_scenarios.csv` for all eight estimands;
- `design_diagnostics.csv` with achieved mortality and landmark event counts;
- `ratio_extreme_counts.csv` for extreme RR/OR estimates;
- RD and RMST-difference bias plots; and
- RR and OR log estimate/truth error plots.

RMSE is computed from a separately named replicate-level squared-error column,
avoiding the sequential `summarise()` overwrite bug identified in ArraySim4-4.
Plot facet labels are derived from `params.R`, so they cannot drift from the
executed sample sizes or event intercepts.
