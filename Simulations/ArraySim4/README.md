# ArraySim4

Cluster-ready comparison of four estimators across grace-period, treatment
effect, and time-varying confounding scenarios.

## Methods

1. Federated clone-censor-weight with time-varying IPW
2. Aligned pooled clone-censor-weight with time-varying IPW
3. Federated time-varying IPW without cloning
4. Federated naive per-protocol analysis without weighting

Fed-CCW and the no-cloning IPW comparator fit their initiation-hazard
denominator models locally at each site. They use a common stabilizing
numerator calculated from aggregated interval-level initiation and eligibility
counts. Sites then transmit weighted event and risk-set totals.

The per-protocol comparator transmits only unweighted event and risk-set counts
by observed treatment group and interval. Its central estimate is exactly the
same as the corresponding pooled crude analysis.

The two no-cloning comparators remain deliberately naive. Federation changes
the data-sharing architecture, not their susceptibility to confounding,
selection, or immortal-time bias.

## Simulation grid

- `tau`: 2, 5, 8
- `beta_trt`: -1.0, -0.7, -0.5
- confounding: small (`0.4`), medium (`1.0`), strong (`2.0`)
- fixed `t_star`: 25
- sites: 3
- patients per site: 1,000
- replicates per cell: 100
- oracle size per cell: 3,000,000
- SLURM array tasks: 27
- weight truncation: none (`c(0, 1)`)

Each task writes 3,200 rows: 100 replicates × 4 methods × 8 reported
estimands (`psi1`, `psi0`, RD, RR, OR, RMST1, RMST0, and RMST difference).

Confounding is varied by multiplying the non-intercept coefficients of the
treatment-initiation model. The intercept is held fixed so the scenarios
primarily change covariate-dependent treatment selection rather than overall
treatment prevalence.

## Validate, run, and aggregate

From this directory:

```bash
Rscript validate_setup.R
sbatch job.sh
```

No submission wrapper is required. Submit `job.sh` while the current working
directory is `ArraySim4`, because SLURM resolves the log paths before the job
script starts.

After all 27 tasks finish:

```bash
Rscript aggregate.R
```

Aggregation writes `summary_all_scenarios.csv` and bias plots for RD, RR, OR,
and RMST difference. Only Fed-CCW and pooled CCW currently report model-based
confidence intervals; the two naive comparators report point estimates.
