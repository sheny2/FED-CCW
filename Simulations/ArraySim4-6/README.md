# ArraySim4-6: site-size balance × patient-mix heterogeneity

This study asks how the six estimators behave when sites differ in both their
sample sizes and baseline patient populations. It extends ArraySim4-4 while
retaining the corrected event-free treatment-initiation DGP.

## Factorial design

The primary design is a 2 × 2 factorial, crossed with three grace periods:

| Factor | Level | Site 1 | Site 2 | Site 3 | Total N |
|---|---|---:|---:|---:|---:|
| Site size | Balanced | 1,000 | 1,000 | 1,000 | 3,000 |
| Site size | Unbalanced | 200 | 800 | 2,000 | 3,000 |

| Patient mix | Site | `X1` distribution | `P(X2=1)` |
|---|---:|---|---:|
| Homogeneous | 1 | N(0, 1²) | 0.40 |
| Homogeneous | 2 | N(0, 1²) | 0.40 |
| Homogeneous | 3 | N(0, 1²) | 0.40 |
| Heterogeneous | 1 | N(-1.50, 0.70²) | 0.10 |
| Heterogeneous | 2 | N(0, 1²) | 0.40 |
| Heterogeneous | 3 | N(1.50, 1.30²) | 0.70 |

The heterogeneous setting is the **high** covariate-shift level in the supplied
design image. The full grid has `3 tau values × 2 size levels × 2 mix levels =
12` SLURM tasks.

Fixed settings:

- grace period `tau = 4, 6, 8`
- three sites and total `N = 3,000` in every replicate
- 300 replicates per cell
- event-model treatment coefficient `beta_trt = -0.7`
- medium treatment-confounding multiplier `1.0`
- treatment-initiation intercept `-3` at every site
- follow-up `t* = 25`
- no weight truncation (`0, 1`)
- oracle Monte Carlo size `5,000,000`

Only site size and the marginal distributions of baseline `X1` and `X2` vary.
The covariate-transition, treatment-initiation, event, and treatment-effect
coefficients remain common across sites. The pooled nuisance model includes a
site fixed effect; federated methods fit the denominator model locally.

## Important estimand detail

Each cell targets the population represented by its actual site sizes. The
oracle therefore samples sites with probabilities `n_k / sum(n_k)`, and local
CCW curve meta-analysis uses the same sample-size weighting.

Consequently, the heterogeneous/balanced target has pooled `E(X1)=0` and
`P(X2=1)=0.40`, while the heterogeneous/unbalanced target has `E(X1)=0.90` and
`P(X2=1)=0.58`. A truth difference between those cells is intentional: putting
most patients in Site 3 changes the pooled target population. The homogeneous
cells retain pooled `E(X1)=0` and `P(X2=1)=0.40` under both size allocations.

## Array task map

| Task | tau | Site sizes | Patient mix |
|---:|---:|---|---|
| 1 | 4 | Balanced | Homogeneous |
| 2 | 6 | Balanced | Homogeneous |
| 3 | 8 | Balanced | Homogeneous |
| 4 | 4 | Unbalanced | Homogeneous |
| 5 | 6 | Unbalanced | Homogeneous |
| 6 | 8 | Unbalanced | Homogeneous |
| 7 | 4 | Balanced | Heterogeneous |
| 8 | 6 | Balanced | Heterogeneous |
| 9 | 8 | Balanced | Heterogeneous |
| 10 | 4 | Unbalanced | Heterogeneous |
| 11 | 6 | Unbalanced | Heterogeneous |
| 12 | 8 | Unbalanced | Heterogeneous |

## Methods and estimands

Each replicate evaluates:

1. `fed_ccw_tvipcw`: federated CCW with a common numerator hazard
2. `pooled_ccw_tvipcw`: aligned pooled CCW
3. `fed_ipw_no_clone`: federated IPW without cloning
4. `fed_perprotocol_naive`: unweighted federated per-protocol analysis
5. `fed_landmark_ipw`: federated landmark IPW
6. `local_ccw_meta`: site-local CCW curves combined by sample-size weights

The recorded estimands are `psi1`, `psi0`, `RD`, `RR`, `OR`, `RMST1`, `RMST0`,
and `RMST_diff`.

## Files

- `params.R`: all study constants, the 2 × 2 design, and task grid
- `DGP_tv.R`: corrected event-free DGP and size-weighted oracle
- `Fed_CCW_TVIPCW.R`: all six estimators
- `Simulation.R`: one-replicate wrapper
- `run_sim_array.R`: one array task runner and design diagnostics
- `job.sh`: 12-task SLURM submission script
- `validate_setup.R`: structural, DGP, oracle, and end-to-end preflight checks
- `aggregate.R`: summaries, diagnostics, extreme-ratio counts, and plots
- `results/`: one `res_task_XXX.rds` file per task
- `logs/`: SLURM output and error logs

## Run

From this directory:

```bash
Rscript validate_setup.R
sbatch job.sh
```

After all 12 tasks finish:

```bash
Rscript aggregate.R
```

Aggregation creates:

- `summary_all_scenarios.csv`
- `design_diagnostics.csv`
- `design_targets.csv`
- `ratio_extreme_counts.csv`
- bias boxplots for `RD` and `RMST_diff`
- log-ratio-error boxplots for `RR` and `OR`

`design_diagnostics.csv` retains observed pooled and site-specific `X1`/`X2`
summaries for every replicate so the intended covariate shift can be checked
directly. RMSE is calculated from replicate-level squared errors, independently
of mean bias.
