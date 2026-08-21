# ArraySim5-4 weight diagnostics

This folder examines the stabilized IPCW used by CCW under the **low**,
**moderate**, and **high** combined site-heterogeneity settings in the parent
ArraySim5-4 simulation.

The default run uses 10 replicates for every combination of heterogeneity and
grace period (`tau = 5, 10, 15`) at `beta_trt = -0.7`. It uses the parent
simulation's five site sizes, patient-mix distributions, site-specific
initiation intercepts, common numerator, and no-truncation setting.

Two distinct denominator models are examined:

1. **Federated: site-specific slopes**: the actual Fed-CCW weights. Each site
   estimates its own initiation model. These are also exactly the weights used
   by the aligned site-stratified pooled CCW.
2. **Pooled site FE: shared slopes**: one pooled initiation model with site
   fixed effects and common covariate slopes.

## Run

From this diagnostics folder:

```bash
Rscript run_weight_diagnostics.R
```

The full default run fits many person-period logistic models and can take a
few minutes. A one-cell smoke run is:

```bash
Rscript run_weight_diagnostics.R \
  --n-reps=1 --heterogeneity=low --taus=5 --include-pooled=false
```

Available arguments are:

- `--n-reps=10`
- `--heterogeneity=low,moderate,high`
- `--taus=5,10,15`
- `--beta-trt=-0.7`
- `--base-seed=54001`
- `--include-pooled=true`
- `--output-dir=outputs`

## Outputs

The `outputs/` folder contains:

- `weight_by_rep_site_interval.csv`: replicate-, site-, strategy-, and
  interval-specific weight quantiles, tail proportions, risk-set size, ESS,
  ESS fraction, and design effect;
- `weight_summary_across_reps.csv`: the corresponding averages over the 10
  replicates;
- `weight_summary_at_tau.csv`: end-of-grace-period subset;
- `propensity_by_rep_site_interval.csv` and
  `propensity_summary_across_reps.csv`: denominator-hazard and
  observed-decision probability support;
- `initiation_prevalence_by_rep_site.csv`: initiation prevalence at each site;
- `model_fit_log.csv`: logistic-model warnings and errors;
- five PNG figures for ESS trajectories, weight-tail trajectories, site-level
  ESS at tau, observed-decision support, and initiation prevalence.

Weights are summarized only while each strategy clone remains in its weighted
risk set. Thus, the diagnostics correspond to the weights that actually enter
the CCW survival estimator, not to weights of already artificially censored
clones.

