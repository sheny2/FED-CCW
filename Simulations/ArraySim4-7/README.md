# ArraySim4-7: when local CCW curve meta-analysis loses support

ArraySim4-7 is a targeted stress study designed to find a setting in which
federated CCW remains close to the truth but local CCW followed by curve
meta-analysis does not. It follows the behavior seen in ArraySim4-4 through
4-6, but isolates **site fragmentation and local treatment support** rather
than adding more patient-mix heterogeneity.

## Why this setting can separate the estimators

Both estimators use the same locally fitted CCW weights and the same common
federated numerator hazard. They differ only in aggregation:

- Federated CCW sums weighted event and risk-set totals across sites first,
  then constructs one survival curve.
- Local CCW constructs a survival curve at each site and then averages those
  curves using fixed baseline sample-size weights.

With sparse initiation and very small sites, many sites have no
treated-strategy risk set late in follow-up. The current local curve estimator
sets the hazard to zero when its weighted risk set is zero, so that site's
curve remains artificially flat yet retains its full baseline weight in the
meta-average. Federated CCW instead pools the nonempty risk sets across sites.

This is a finite-sample/local-positivity stress test of this particular curve
meta-analysis rule. It is not a claim that all meta-analytic CCW estimators are
generally inconsistent.

## Design

Total sample size is fixed at 5,000 while the same homogeneous population is
split into progressively smaller sites:

| Fragmentation level | Sites | Patients/site | Total N |
|---|---:|---:|---:|
| `10_sites_x_500` | 10 | 500 | 5,000 |
| `20_sites_x_250` | 20 | 250 | 5,000 |
| `50_sites_x_100` | 50 | 100 | 5,000 |
| `100_sites_x_50` | 100 | 50 | 5,000 |

Each fragmentation level is crossed with:

| Initiation level | Intercept | Approximate initiation by tau in pilot |
|---|---:|---:|
| Moderate | -3.5 | 9.5% |
| Sparse | -4.5 | 3.6% |

Fixed settings:

- `tau = 4`, `t* = 25`
- `beta_trt = -0.7`
- 300 replicates per cell
- homogeneous baseline mix at every site: `X1 ~ N(0,1)` and `X2 ~ Bernoulli(0.4)`
- treatment initiation independent of `X` and `L` (`conf_mult = 0`)
- common event and covariate-transition models across sites
- local 1st/99th-percentile weight truncation
- 5,000,000-person Monte Carlo oracle

Removing treatment confounding is deliberate: it prevents nuisance-model
confounding bias from obscuring the aggregation/support mechanism. The CCW
denominator model remains estimated, and both aggregation methods receive the
same fitted local summaries.

## Pilot used to choose the stress cell

`pilot_search_summary.csv` records a 15-replicate exploratory screen. In the
`100 sites × 50`, sparse-initiation cell it found:

| Quantity | Federated CCW | Local curve meta |
|---|---:|---:|
| Mean RD bias | +0.028 | -0.075 |
| Mean RMST-difference bias | -0.054 | +1.097 |

An average of 56.1 of the 100 local treated-strategy risk sets were empty at
the endpoint. In the `10 × 500` sparse-initiation control, both RD biases were
near zero and almost no endpoint treated risk set was empty. These pilot
figures motivated the design; they are exploratory, not the final simulation
results.

## Task map

| Task | Fragmentation | Initiation |
|---:|---|---|
| 1 | 10 × 500 | Moderate |
| 2 | 20 × 250 | Moderate |
| 3 | 50 × 100 | Moderate |
| 4 | 100 × 50 | Moderate |
| 5 | 10 × 500 | Sparse |
| 6 | 20 × 250 | Sparse |
| 7 | 50 × 100 | Sparse |
| 8 | 100 × 50 | Sparse |

## Methods

The study focuses on three methods:

1. `fed_ccw_tvipcw`: federated aggregation of weighted event/risk totals
2. `pooled_ccw_tvipcw`: centralized benchmark
3. `local_ccw_meta`: fixed sample-size-weighted local CCW survival curves

The joint fed/local implementation fits each site's denominator model only
once. Thus their difference cannot be caused by different random fits.

## Files and workflow

- `params.R`: design constants and eight-cell grid
- `DGP_tv.R`: corrected event-free DGP and common-population oracle
- `Fed_CCW_TVIPCW.R`: estimators and shared-fit fed/local implementation
- `Simulation.R`: targeted one-replicate wrapper and support diagnostics
- `run_sim_array.R`: SLURM task runner
- `job.sh`: eight-task job array
- `validate_setup.R`: design, DGP, support, and estimator checks
- `aggregate.R`: final summaries, paired comparisons, diagnostics, and plots
- `pilot_search_summary.csv`: exploratory design-selection evidence

Run:

```bash
Rscript validate_setup.R
sbatch job.sh
```

After all eight tasks finish:

```bash
Rscript aggregate.R
```

Key outputs are `summary_all_scenarios.csv`, `support_diagnostics.csv`,
`fed_vs_local_paired.csv`, two bias boxplots, and
`bias_vs_empty_treated_support.png`.


Local run code
Rscript validate_setup.R
mkdir -p logs results

for task_id in {1..8}; do
  echo "Starting task ${task_id}"
  SLURM_CPUS_PER_TASK=4 Rscript run_sim_array.R "${task_id}" \
    2>&1 | tee "logs/local_${task_id}.log"
done