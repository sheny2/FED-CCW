# ArraySim4-4: Sample Size and Treatment Initiation

Factorial simulation designed to determine whether the behavior of federated
landmark IPW and local CCW curve meta-analysis is driven by small site samples,
rare treatment initiation, or their combination.

This study uses the corrected event-free DGP from `ArraySim4-3 Ten Small Sites
fixed`: treatment initiation and time-varying covariate evolution stop after a
terminal event, and all initiation models and numerator counts use event-free
treatment risk sets.

## Factorial design

The study holds the following settings fixed:

- 10 equally sized sites
- medium confounding (`conf_mult = 1.0`)
- treatment effect `beta_trt = -0.7`
- follow-up horizon `t_star = 25`
- no weight truncation (`c(0, 1)`)
- 300 replicates per cell
- 3,000,000-person oracle per cell

It varies three factors:

| Factor | Levels |
|---|---|
| Grace period | `tau = 5, 10, 15` |
| Patients per site | low = 100; large = 1,000 |
| Initiation prevalence | low intercept = `-4.5`; high intercept = `-1.5` |

This is a `3 x 2 x 2` design with 12 array tasks. Total sample size is 1,000
in the low-sample cells and 10,000 in the large-sample cells.

Every site uses the same initiation intercept within a cell. Consequently,
the design isolates sample size and treatment prevalence rather than mixing
those factors with between-site prevalence heterogeneity. At zero covariates,
the per-interval initiation probabilities are approximately:

- low initiation: `plogis(-4.5) = 1.1%`
- high initiation: `plogis(-1.5) = 18.2%`

Actual initiation by `tau` is recorded in every replicate because cumulative
uptake also depends on covariates, survival, and grace-period length.

## Methods

1. Federated CCW
2. Pooled CCW
3. Federated IPW without cloning
4. Federated unweighted per-protocol analysis
5. Federated landmark IPW
6. Local CCW with sample-size-weighted survival-curve aggregation

The key comparisons are landmark IPW and local curve aggregation against
federated and pooled CCW across the four sample-size/initiation cells.
Landmark IPW targets post-landmark survival conditional on being event-free
through `tau`, so its difference from the common CCW oracle includes a target
difference and is not solely conventional estimator bias.

Each task writes 14,400 rows:

```text
300 replicates x 6 methods x 8 estimands
```

The estimands are `psi1`, `psi0`, RD, RR, OR, RMST1, RMST0, and RMST
difference.

## Task mapping

`run_sim_array.R` constructs the task grid from `make_study_grid()` in
`params.R`. The ordering is:

| Task | tau | Patients/site | Initiation intercept |
|---:|---:|---:|---:|
| 1 | 5 | 100 | -4.5 |
| 2 | 10 | 100 | -4.5 |
| 3 | 15 | 100 | -4.5 |
| 4 | 5 | 1,000 | -4.5 |
| 5 | 10 | 1,000 | -4.5 |
| 6 | 15 | 1,000 | -4.5 |
| 7 | 5 | 100 | -1.5 |
| 8 | 10 | 100 | -1.5 |
| 9 | 15 | 100 | -1.5 |
| 10 | 5 | 1,000 | -1.5 |
| 11 | 10 | 1,000 | -1.5 |
| 12 | 15 | 1,000 | -1.5 |

The oracle does not depend on finite sample size, so pairs of tasks that differ
only in patients per site share the same target estimand. They independently
recompute it with the same Monte Carlo size and seed.

## Validate and run

From this directory:

```bash
Rscript validate_setup.R
sbatch job.sh
```

`job.sh` requests 10 CPUs and 40 GB per task. The 1,000-patient-per-site cells
are substantially more expensive than the 100-patient-per-site cells; adjust
resources to the cluster if needed. The `logs/` and `results/` directories are
included and must exist when `sbatch` is called.

For a local run:

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

This produces:

- `summary_all_scenarios.csv`, covering all eight estimands;
- `design_diagnostics.csv`, with achieved initiation and landmark sample sizes;
- `ratio_extreme_counts.csv`, explicitly counting extreme RR/OR estimates;
- ordinary bias plots for RD and RMST difference; and
- log estimate/truth error plots for RR and OR.

RR and OR use `log(estimate / truth)` in the plots so rare landmark explosions
do not flatten every other method against zero. Raw estimates, ordinary bias,
and RMSE remain unchanged in the RDS and CSV output.
