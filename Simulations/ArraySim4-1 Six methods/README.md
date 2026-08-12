# ArraySim4-1

Cluster-ready comparison of six federated or privacy-preserving estimators
across grace-period, treatment-effect, and time-varying-confounding scenarios.

## Methods

1. Federated clone-censor-weight with time-varying IPW
2. Aligned pooled clone-censor-weight with time-varying IPW
3. Federated time-varying IPW without cloning
4. Federated naive per-protocol analysis without weighting
5. Federated landmark IPW
6. Local CCW followed by sample-size-weighted survival-curve meta-analysis

The pooled CCW comparator represents the estimate that would be obtained if
patient-level records could be combined. All other methods transmit only
aggregate quantities or site-level curves.

## Federated architecture

Fed-CCW and no-cloning IPW fit treatment-initiation denominator models locally.
A common stabilizing numerator is calculated from site-level initiation and
eligibility counts. Sites then transmit weighted event and risk-set totals.

The unweighted per-protocol method transmits event and risk-set counts by
observed treatment group and interval. It exactly reproduces the corresponding
pooled crude analysis.

For landmark IPW, each site:

- restricts to patients event-free through `tau`;
- fits `A_tau` on baseline covariates, current `L_tau`, and lagged `L_(tau-1)`;
- constructs stabilized weights using a numerator calculated from aggregated
  landmark treatment counts;
- transmits weighted post-landmark event and risk-set totals.

Landmark IPW targets survival conditional on reaching `tau`, not the original
time-zero grace-period estimand. Its bias against the common oracle therefore
reflects both estimation error and this target difference.

For local CCW meta-analysis, every site completes CCW locally using the common
numerator and its own denominator model. It sends only its sample size and two
estimated survival curves. The center takes sample-size-weighted averages of
the curves and then calculates RD, RR, OR, and RMST difference. Unlike Fed-CCW,
it combines completed curves instead of event/risk totals.

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

Each task writes 4,800 rows: 100 replicates × 6 methods × 8 estimands
(`psi1`, `psi0`, RD, RR, OR, RMST1, RMST0, and RMST difference).

Confounding is varied by multiplying the non-intercept coefficients of the
treatment-initiation model. The intercept remains fixed.

## Cluster resources

`job.sh` uses the updated ArraySim4 configuration:

- 5 CPUs per task
- 20 GB memory per task
- R 4.3

No `submit.sh` wrapper is included.

## Validate, run, and aggregate

From this directory:

```bash
Rscript validate_setup.R
sbatch job.sh
```

The included `logs/` directory must exist when `sbatch` is called because
SLURM opens the output paths before starting the job.

After all 27 tasks finish:

```bash
Rscript aggregate.R
```

Aggregation writes `summary_all_scenarios.csv` and bias plots for RD, RR, OR,
and RMST difference. Fed-CCW and pooled CCW report their existing model-based
confidence intervals; the remaining four comparators report point estimates.

## Repeated covariate-balance diagnostic

Run:

```bash
Rscript balance_diagnostics.R
```

The diagnostic independently simulates 10 datasets under each small, medium,
and strong confounding setting, with `tau=5`, `beta_trt=-0.7`, three sites,
and 1,000 patients per site. These settings can be changed, for example:

```bash
Rscript balance_diagnostics.R \
  --tau=8 \
  --beta-trt=-1 \
  --n-per-site=2000 \
  --n-reps=25 \
  --seed=20260811
```

Among patients event-free through `tau`, it compares initiators by `tau` with
non-initiators using standardized mean differences for the continuous
variables: baseline `x1`, current `L1`, `L2`, and their lagged values. The
after-weighting calculation uses the same local denominator models and common
numerator as Fed-CCW. Binary `x2` remains in the weight model but is excluded
from the requested SMD summaries and plot.

Sites contribute only covariate moments (`n`, `sum(w)`, `sum(w^2)`,
`sum(w*x)`, and `sum(w*x^2)`). The center reconstructs means, variances, SMDs,
and effective sample sizes. Outputs are placed in
`balance_diagnostics_repeated/`:

- `balance_smd_by_iteration.csv`
- `balance_smd_average.csv`
- `balance_group_ess_by_iteration.csv`
- `balance_group_ess_average.csv`
- `common_numerator_hazard_by_iteration.csv`
- `balance_overall_average.csv`
- `balance_love_plot_average.png`

This is the simplified end-of-grace-period diagnostic described in the
referenced excerpt. It does not replace interval-specific balance diagnostics
for the full longitudinal weighting process.
