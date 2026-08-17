# ArraySim5-2: Patient-Mix Heterogeneity

This simulation is based on ArraySim5-1 and compares the same five estimators.
The previous small/medium/strong confounding-multiplier axis is replaced by
low/moderate/high **between-site heterogeneity in baseline patient mix**.

The treatment-initiation coefficients (`DEFAULT_BETA_INIT`), outcome model,
time-varying covariate model, treatment effects, and weight settings are held
fixed across heterogeneity levels. Thus, “high” means that the sites enroll
more different patient populations; it does not mean that the conditional
confounding coefficients are larger.

## Site design

- 3 sites
- 1,000 patients per site (3,000 total per replicate)
- site-specific initiation intercepts evenly spaced from -4.5 to -1.5
- each site contributes one third of the target population and oracle

The baseline covariate distributions are:

| Heterogeneity | Site | X1 distribution | P(X2 = 1) |
|---|---:|---:|---:|
| Low | 1 | N(-0.25, 0.95^2) | 0.35 |
| Low | 2 | N(0, 1^2) | 0.40 |
| Low | 3 | N(0.25, 1.05^2) | 0.45 |
| Moderate | 1 | N(-0.75, 0.85^2) | 0.20 |
| Moderate | 2 | N(0, 1^2) | 0.40 |
| Moderate | 3 | N(0.75, 1.15^2) | 0.60 |
| High | 1 | N(-1.50, 0.70^2) | 0.10 |
| High | 2 | N(0, 1^2) | 0.40 |
| High | 3 | N(1.50, 1.30^2) | 0.70 |

With equal site sizes, all three levels preserve the pooled means
`E(X1) = 0` and `P(X2 = 1) = 0.40`. What changes is the separation between
sites (and the pooled variance of X1), making this a clean patient-mix
heterogeneity stress test.

## Methods

1. Federated CCW
2. Pooled CCW with site fixed effects in the initiation model
3. Federated IPW without cloning
4. Federated unweighted per-protocol analysis
5. Federated landmark IPW

## Simulation grid

- `tau`: 5, 10, 15
- `beta_trt`: -1.0, -0.7, -0.5
- patient-mix heterogeneity: low, moderate, high
- fixed `t_star`: 25
- fixed confounding coefficients: `DEFAULT_BETA_INIT`
- 100 replicates per cell
- 3,000,000-person scenario-specific oracle per cell
- 27 array tasks
- no weight truncation (`c(0, 1)`)
- 4,000 rows per task: 100 replicates x 5 methods x 8 estimands

Both the finite simulated data and the Monte Carlo oracle use the same
site-specific X distributions for the selected heterogeneity level.

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

After all tasks complete:

```bash
Rscript aggregate.R
```

All scenario inputs, including the full patient-mix table, are defined in
`params.R`.
