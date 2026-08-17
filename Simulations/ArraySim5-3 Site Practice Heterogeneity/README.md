# ArraySim5-3: Site-Practice Heterogeneity

This simulation examines treatment-initiation differences between sites that
remain after conditioning on the measured patient covariates and clinical
history. It compares the same five methods as ArraySim5-2, but replaces the
patient-mix axis with low/moderate/high heterogeneity in site-specific
treatment-initiation practices.

## What varies

All sites enroll from the same baseline distributions:

- `X1 ~ N(0, 1)`
- `P(X2 = 1) = 0.40`

The measured covariate effects in the initiation model, time-varying covariate
model, outcome model, and treatment effect are identical across sites. Only
the site-specific initiation intercept changes:

| Practice heterogeneity | Site 1 | Site 2 | Site 3 | Intercept spread |
|---|---:|---:|---:|---:|
| Low | -3.25 | -3.00 | -2.75 | 0.50 |
| Moderate | -3.75 | -3.00 | -2.25 | 1.50 |
| High | -4.50 | -3.00 | -1.50 | 3.00 |

At zero covariates, the corresponding per-interval initiation probabilities
are approximately:

| Practice heterogeneity | Site 1 | Site 2 | Site 3 |
|---|---:|---:|---:|
| Low | 3.7% | 4.7% | 6.0% |
| Moderate | 2.3% | 4.7% | 9.5% |
| High | 1.1% | 4.7% | 18.2% |

Every scenario remains centered at an initiation intercept of -3. This makes
the axis represent increasing between-site practice heterogeneity instead of
a uniform change in treatment uptake or confounding strength.

## CCW alignment

Both CCW estimators use:

- one common numerator constructed from aggregated event-free initiation and
  eligibility counts across sites;
- a separate denominator nuisance model fitted within each site;
- the same nuisance-model formula, artificial censoring, weight construction,
  truncation setting, weighted survival estimator, and influence-function
  formulas.

Federated CCW transmits only site-level quantities needed by the central
aggregator. The fully site-stratified pooled CCW benchmark has individual-level
data centrally but deliberately performs the same separate nuisance fits and
then pools the weighted outcome records. Its pooled weighted Kaplan-Meier
calculation is algebraically identical to summing site weighted event/risk
totals. Therefore their estimates and model-based standard errors should agree
up to numerical precision; this is the intended implementation benchmark.

## Methods

1. Federated CCW
2. Pooled CCW with fully site-stratified nuisance models
3. Federated IPW without cloning
4. Federated unweighted per-protocol analysis
5. Federated landmark IPW

## Simulation design

- 3 sites, 1,000 patients per site
- `tau`: 5, 10, 15
- `beta_trt`: -1.0, -0.7, -0.5
- site-practice heterogeneity: low, moderate, high
- fixed `t_star`: 25
- fixed conditional confounding coefficients: `DEFAULT_BETA_INIT`
- 100 replicates per cell
- 3,000,000-person scenario-specific oracle per cell
- 27 array tasks
- no weight truncation (`c(0, 1)`)
- 4,000 result rows/task: 100 replicates x 5 methods x 8 estimands

The oracle uses the scenario-specific site initiation intercepts and the same
equal-weight three-site target population as the finite simulations.

## Validate and run

From this directory:

```bash
Rscript validate_setup.R
sbatch --array=1-27 job.sh
```

No separate submission script is needed. After all tasks complete:

```bash
Rscript aggregate.R
```

All design inputs are defined in `params.R`.
