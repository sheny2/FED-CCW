# ArraySim5-4: Combined Site Heterogeneity

This simulation jointly varies three realistic sources of between-site
heterogeneity:

1. site sample size;
2. baseline patient mix through the distributions of `X1` and `X2`; and
3. residual site practice through the treatment-initiation intercept.

The conditional covariate effects in the initiation model, the time-varying
covariate process, and the outcome model are held fixed. Every level contains
five sites and 5,000 total patients.

## Joint heterogeneity settings

| Level | Site | N | X1 distribution | P(X2=1) | Initiation intercept |
|---|---:|---:|---:|---:|---:|
| Low | 1 | 900 | N(-0.250, 0.950^2) | 0.350 | -3.250 |
| Low | 2 | 950 | N(-0.125, 0.975^2) | 0.375 | -3.125 |
| Low | 3 | 1,000 | N(0, 1^2) | 0.400 | -3.000 |
| Low | 4 | 1,050 | N(0.125, 1.025^2) | 0.425 | -2.875 |
| Low | 5 | 1,100 | N(0.250, 1.050^2) | 0.450 | -2.750 |
| Moderate | 1 | 600 | N(-0.750, 0.850^2) | 0.200 | -3.750 |
| Moderate | 2 | 800 | N(-0.375, 0.925^2) | 0.300 | -3.375 |
| Moderate | 3 | 1,000 | N(0, 1^2) | 0.400 | -3.000 |
| Moderate | 4 | 1,200 | N(0.375, 1.075^2) | 0.500 | -2.625 |
| Moderate | 5 | 1,400 | N(0.750, 1.150^2) | 0.600 | -2.250 |
| High | 1 | 400 | N(-1.500, 0.700^2) | 0.100 | -4.500 |
| High | 2 | 600 | N(-0.750, 0.850^2) | 0.250 | -3.750 |
| High | 3 | 800 | N(0, 1^2) | 0.400 | -3.000 |
| High | 4 | 1,200 | N(0.750, 1.150^2) | 0.550 | -2.250 |
| High | 5 | 2,000 | N(1.500, 1.300^2) | 0.700 | -1.500 |

The gradients are intentionally aligned: later and larger sites have larger
baseline covariate values and higher residual treatment-initiation propensity.
Consequently, the high level is a joint stress test of local nuisance-model
precision, covariate overlap, and treatment-initiation positivity. The oracle
uses each level's exact site-size proportions and patient distributions, so
every cell is evaluated against its own target-population truth.

## Methods

1. **Federated CCW**: separate denominator nuisance model at each site, one
   common numerator, and central aggregation of weighted event/risk totals.
2. **Pooled CCW (fully site-stratified)**: the individual-data implementation
   of exactly the same separate site nuisance models. It should agree with
   Fed-CCW to numerical precision.
3. **Pooled CCW (site fixed effects, shared slopes)**: one conventional pooled
   person-interval logistic initiation model,

   `logit P(initiation) = interval FE + site FE + common X/L slopes`.

   This model borrows information across sites for the covariate slopes while
   retaining residual site-specific initiation propensity through site fixed
   effects. Its common-slope assumption is correct in the current DGP.
4. Federated IPW without cloning.
5. Federated unweighted per-protocol analysis.
6. Federated landmark IPW.

All three CCW implementations use the same common numerator and weighted
survival estimand. The first two differ only in data location. The third uses
a more restrictive, information-borrowing nuisance model and may have more
stable weights and lower finite-sample variance.

## Simulation grid

- `tau`: 5, 10, 15
- `beta_trt`: -1.0, -0.7, -0.5
- combined heterogeneity: low, moderate, high
- fixed `t_star`: 25
- 100 replicates per cell
- 3,000,000-person scenario-specific oracle per cell
- 27 array tasks
- no weight truncation (`c(0, 1)`)
- 4,800 result rows/task: 100 replicates x 6 methods x 8 estimands

## Validate and run

From this directory:

```bash
Rscript validate_setup.R
sbatch --array=1-27 job.sh
```

After all task files are available:

```bash
Rscript aggregate.R
```

There is no separate submission script. All settings are defined in
`params.R`.


Local run

```bash
n_tasks=$(Rscript -e 'source("params.R"); cat(nrow(simulation_grid()))')
for task_id in $(seq 1 "$n_tasks"); do
  SLURM_CPUS_PER_TASK=20 Rscript run_sim_array.R "$task_id" \
    > "logs/local_${task_id}.out" \
    2> "logs/local_${task_id}.err"
done

Rscript aggregate.R
```
