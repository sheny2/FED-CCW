# ArraySim4-2 Unbalance

Cluster-ready study designed to distinguish federated sufficient-statistic
aggregation from completing CCW separately at each site and averaging the
resulting survival curves.

## Deliberately unbalanced site design

| Site | Patients per replicate | Initiation intercept | Zero-covariate initiation probability per interval |
|---|---:|---:|---:|
| 1 | 200 | -4.5 | 1.1% |
| 2 | 800 | -3.0 | 4.7% |
| 3 | 2,000 | -1.5 | 18.2% |

The sites share the treatment-initiation covariate slopes, outcome model,
covariate-transition model, and treatment effect. The confounding multiplier
still changes the initiation slopes across small (`0.4`), medium (`1.0`), and
strong (`2.0`) scenarios.

This design makes the smallest site the site with the least treatment uptake.
It should have the least stable local treated-strategy curve, while the largest
site supplies most of the information. The oracle uses the same `200:800:2000`
target-population mixture and the same site-specific initiation intercepts.

## Six methods

1. Federated CCW using a common numerator and locally fitted denominators
2. Pooled CCW with a site fixed effect in its pooled denominator model
3. Federated IPW without cloning
4. Federated unweighted per-protocol analysis
5. Federated landmark IPW
6. Local CCW with sample-size-weighted survival-curve meta-analysis

All weighted federated methods use local nuisance models. All CCW methods use
the same common numerator calculated from aggregated initiation/eligibility
counts. The pooled comparator includes site in its initiation model so the
known site-specific prevalence structure is not deliberately misspecified.

The key comparison is Fed-CCW versus local CCW curve aggregation:

- Fed-CCW aggregates interval-specific weighted events and risk sets before
  constructing a survival curve.
- Local CCW constructs a complete curve at each site and then weights those
  curves by baseline site sample size.

## Simulation grid

- `tau`: 2, 5, 8
- `beta_trt`: -1.0, -0.7, -0.5
- confounding: small, medium, strong
- fixed `t_star`: 25
- 100 replicates per cell
- 3,000,000-person oracle per cell
- 27 SLURM array tasks
- no weight truncation (`c(0, 1)`)
- 4,800 output rows per task: 100 × 6 methods × 8 estimands

## Validate and run

From this directory:

```bash
Rscript validate_setup.R
sbatch job.sh
```

`job.sh` uses the previous ArraySim4-1 resource configuration: 5 CPUs and
20 GB memory per array task. No `submit.sh` is included.

After all 27 tasks complete:

```bash
Rscript aggregate.R
```

This writes `summary_all_scenarios.csv` plus bias plots for RD, RR, OR, and
RMST difference.
