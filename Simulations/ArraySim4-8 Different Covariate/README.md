# ArraySim4-8: different covariates available across sites

This small study tests whether federated estimators can use different local
covariate sets without sharing patient-level covariates. It distinguishes the
case where every site's available set is sufficient from the fundamentally
different case where active confounders are unavailable.

## Controlled comparison

All nine cells use the **same three-site DGP**, homogeneous patient mix, event
model, and treatment mechanism. Only sample size and the covariates made
available to the nuisance estimators change. Therefore, within a sample-size
level, all three availability scenarios have the same causal truth.

The treatment-initiation DGP differs by site:

| Site | Active predictors of initiation | Variables affecting the outcome |
|---:|---|---|
| 1 | `X1`, current `L1` | Full `X1`, `X2`, `L1`, `L2` history |
| 2 | `X2`, current `L2` | Full `X1`, `X2`, `L1`, `L2` history |
| 3 | `X1`, `X2`, current `L1`, current `L2` | Full history |

A variable is a treatment confounder at a site only when it predicts both
initiation and the outcome. Outcome predictors that do not enter that site's
initiation mechanism do not have to be included in its weight model for
exchangeability.

## Covariate-availability scenarios

| Scenario | Site 1 model | Site 2 model | Site 3 model | Interpretation |
|---|---|---|---|---|
| `common_full` | Full set | Full set | Full set | Harmonized control |
| `different_sufficient` | `X1`, `L1`, lagged `L1` | `X2`, `L2`, lagged `L2` | Full set | Different local sets, all sufficient |
| `different_insufficient` | `X1` only | `X2` only | `X1`, `X2` only | Active time-varying confounders unavailable |

The complete simulated variables remain in the DGP so outcomes can be
generated correctly. An estimator only receives the variables listed in its
analysis set; unavailable columns are never included in its local formula.

For pooled CCW, the analysis is restricted to the intersection of variables
available at all sites. This is the full set in `common_full`, but is empty in
the two different-variable scenarios, leaving only time and site indicators.
When common covariates exist, the pooled model includes site-by-covariate
interactions so each site's treatment mechanism can have different slopes.
This represents a conventional harmonized pooled analysis, not an oracle that
quietly accesses unavailable variables.

## Sample-size experiment

| Level | Patients/site | Sites | Total N |
|---|---:|---:|---:|
| Small | 100 | 3 | 300 |
| Medium | 500 | 3 | 1,500 |
| Large | 1,000 | 3 | 3,000 |

The full grid contains `3 sample sizes × 3 availability scenarios = 9` tasks,
with 100 replicates per task.

Other settings are held at their normal values:

- `tau = 4`, `t* = 25`
- medium confounding coefficients
- treatment-effect coefficient `beta = -0.7`
- initiation intercept `-3` at every site
- homogeneous `X1 ~ N(0,1)` and `X2 ~ Bernoulli(0.4)`
- no weight truncation
- 3,000,000-person Monte Carlo truth

## Methods

1. `fed_ccw_site_specific`: Fed-CCW with a separate valid formula per site
2. `pooled_ccw_common_set`: pooled CCW restricted to the cross-site common set
3. `fed_ipw_no_clone`: no-cloning comparator using the site-specific sets
4. `fed_perprotocol_naive`: unweighted per-protocol comparator
5. `fed_landmark_site_specific`: landmark IPW using site-specific sets
6. `local_ccw_meta_site_specific`: local CCW curves using site-specific sets,
   followed by sample-size-weighted curve meta-analysis

The experiment should answer two separate questions:

- Does Fed-CCW remain accurate when formulas differ but each site measures a
  sufficient adjustment set?
- Does increasing sample size repair the bias caused by a genuinely missing
  confounder? It should reduce variance, but it cannot restore identification.

## Task map

| Task | Sample size/site | Availability scenario |
|---:|---:|---|
| 1 | 100 | Common full |
| 2 | 500 | Common full |
| 3 | 1,000 | Common full |
| 4 | 100 | Different, sufficient |
| 5 | 500 | Different, sufficient |
| 6 | 1,000 | Different, sufficient |
| 7 | 100 | Different, insufficient |
| 8 | 500 | Different, insufficient |
| 9 | 1,000 | Different, insufficient |

## Files and workflow

- `params.R`: DGP coefficients, availability configurations, and task grid
- `DGP_tv.R`: corrected event-free DGP with site-specific initiation models
- `Fed_CCW_TVIPCW.R`: six methods with site-specific formula support
- `Simulation.R`: one-replicate wrapper
- `run_sim_array.R`: local/SLURM task runner
- `job.sh`: nine-task SLURM job
- `validate_setup.R`: structural, DGP, truth, and method checks
- `aggregate.R`: summaries, diagnostics, boundary table, and bias plots

Run on SLURM:

```bash
Rscript validate_setup.R
sbatch job.sh
```

Run locally:

```bash
for task_id in {1..9}; do
  SLURM_CPUS_PER_TASK=8 Rscript run_sim_array.R "$task_id"
done
Rscript aggregate.R
```

Important outputs are `summary_all_scenarios.csv`,
`fed_ccw_boundary_summary.csv`, `design_diagnostics.csv`, and RD/RMST bias
plots.
