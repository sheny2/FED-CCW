# ArraySim3-3 Misspecification fixed

Four-method simulation aligned with `ArraySim3-3 Direct G-comp fixed`.

## Methods

1. Federated CCW
2. Pooled CCW
3. Direct g-computation without cloning
4. Clone-censor g-computation using stacked arm-specific clone frames

## Nuisance specifications

- `correct`: logistic models with baseline, current L, and lagged L
- `coarse_L`: logistic models using median-dichotomized current and lagged L
- `no_tv`: logistic models omitting current and lagged L
- `tree`: classification trees with the full covariate set

The specification is applied to the CCW initiation-hazard denominator and
the g-computation outcome hazard in the same scenario cell. The CCW
stabilizing numerator remains logistic, and the L transition models remain
correctly specified.

The tree implementation uses `rpart`, a recommended R package. Internal tree
cross-validation is disabled so it is not repeated inside every replicate.

## Cluster grid

- `tau`: 2, 5, 8
- treatment coefficient: -1.0, -0.7, -0.5
- confounding: small, medium, strong
- nuisance specification: correct, coarse_L, no_tv, tree
- 108 SLURM array tasks
- 100 replications per task
- 3 sites with 1,000 subjects per site
- fixed follow-up horizon `t_star = 25`

The `logs/` and `results/` directories are included, so submission can be
made directly from this directory:

```bash
sbatch job.sh
```

After all 108 tasks finish:

```bash
Rscript aggregate.R
```

Aggregation requires `dplyr` and `ggplot2`.
