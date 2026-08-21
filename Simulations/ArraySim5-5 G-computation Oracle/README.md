# ArraySim5-5: G-computation oracle

ArraySim5-5 retains the complete five-site combined-heterogeneity design from
ArraySim5-4 but replaces the true-weight CCW oracle with an intervention-based
parametric Monte Carlo g-computation oracle. The redundant fully
site-stratified pooled CCW method is omitted because it is algebraically
identical to Federated CCW.

## G-computation truth

For every `(tau, beta_trt, heterogeneity)` cell, the oracle generates
`SIM_TRUTH_N` patients from the scenario-specific target population:

- site probabilities equal the configured site-size proportions;
- `X1` and `X2` follow the configured site-specific patient mix;
- `L1` and `L2` follow their true longitudinal transition equations;
- the true event hazard is evaluated under each strategy.

The two interventions are explicit:

1. **Initiate within the grace period (`g=1`)**: during intervals before
   `tau`, initiation follows the true site-specific initiation hazard given
   baseline and current time-varying covariates. Anyone still untreated is
   deterministically initiated in interval `tau`.
2. **Do not initiate (`g=0`)**: treatment is withheld throughout follow-up.

Because treatment does not affect the `L` process in this DGP, both strategies
can use the same simulated longitudinal covariate trajectories. Conditional
survival probabilities are multiplied over time and averaged over the Monte
Carlo target population. Binary event times are not simulated, which reduces
Monte Carlo noise. No estimated weights, cloning, or artificial censoring are
used to define the oracle.

This truth therefore evaluates all estimators against a clinically explicit
grace-period intervention rather than defining truth by a very large weighted
CCW analysis.

## Joint heterogeneity settings

The low, moderate, and high levels jointly vary:

1. site sample size;
2. the site-specific distributions of `X1` and `X2`; and
3. the residual site-specific treatment-initiation intercept.

Every level contains five sites and 5,000 total observed patients. All
settings remain centralized in `params.R` and are unchanged from ArraySim5-4.

## Methods

Five methods are compared:

1. Federated CCW with separate local denominator nuisance models and a common
   numerator;
2. pooled CCW with site fixed effects and shared nuisance-model slopes;
3. Federated landmark IPW;
4. Federated IPW without cloning;
5. Federated unweighted per-protocol analysis.

The aligned pooled CCW with fully site-stratified nuisance models is removed
because it returns exactly the same estimates as Federated CCW.

## Simulation grid

- grace period `tau`: 5, 10, 15;
- treatment log-effect `beta_trt`: -1.0, -0.7, -0.5;
- combined heterogeneity: low, moderate, high;
- follow-up horizon: 25 intervals;
- replicates per cell: defined by `SIM_N_REPS` in `params.R` (currently 500);
- oracle population per cell: 3,000,000;
- no weight truncation (`c(0, 1)`);
- 27 array tasks.

## Validate and run

From this directory:

```bash
Rscript validate_setup.R
mkdir -p logs results
sbatch --array=1-27 job.sh
```

After all task files are available:

```bash
Rscript aggregate.R
```

This produces `summary_all_scenarios.csv` and bias plots for RD, RR, OR, and
the RMST difference.

