# ArraySim4-3: Ten Small Sites

Finite-sample stress test of Fed-CCW versus completing CCW locally and
averaging the ten resulting survival curves.

## Site design

- 10 sites
- 300 patients per site
- total sample size: 3,000 per replicate
- site-specific treatment-initiation intercepts evenly spaced from `-4.5`
  through `-1.5`
- corresponding zero-covariate initiation probabilities range from about
  1.1% to 18.2% per interval
- every site contributes 10% of the target population

All sites share the treatment-initiation covariate slopes, covariate-transition
model, outcome model, and treatment effect. Confounding multipliers remain
small (`0.4`), medium (`1.0`), and strong (`2.0`).

Dividing the same total sample across ten sites forces ten nuisance-model fits
and ten completed local survival curves. Unlike the previous unbalanced study,
rare-treatment sites are not almost removed by a very small target-population
weight: each local curve receives weight `0.10` in curve aggregation.

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

## Simulation grid

- `tau`: 2, 5, 8
- `beta_trt`: -1.0, -0.7, -0.5
- confounding: small, medium, strong
- fixed `t_star`: 25
- 100 replicates per cell
- 3,000,000-person mixture-specific oracle per cell
- 27 array tasks
- no weight truncation (`c(0, 1)`)
- 4,800 rows per task

The oracle allocates its 3,000,000 observations equally across the ten sites
and uses the corresponding site initiation intercept for every observation.

## Run

From this directory:

```bash
Rscript validate_setup.R
sbatch job.sh
```

For a local eight-core run, execute array tasks sequentially:

```bash
for task_id in {1..27}; do
  SLURM_CPUS_PER_TASK=8 Rscript run_sim_array.R "$task_id" \
    > "logs/local_${task_id}.out" \
    2> "logs/local_${task_id}.err"
done
```

`job.sh` retains the previous configuration of 5 CPUs and 20 GB per task.
No `submit.sh` is included.

After all 27 tasks complete:

```bash
Rscript aggregate.R
```

The aggregation script calculates RMSE from replicate-level errors before
collapsing bias, avoiding the sequential-summary bug discovered in ArraySim4-2.
