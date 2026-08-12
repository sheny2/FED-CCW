# ArraySim3-4 Heterogeneity 2

Corrected cluster simulation for increasingly heterogeneous patient
populations across three sites. This version resolves the original
federated-versus-pooled target mismatch by using one common stabilizing
numerator hazard at every site.

## Methods

Exactly four methods are reported:

1. Updated federated clone-censor-weight with time-varying IPCW
2. Aligned pooled clone-censor-weight with time-varying IPCW
3. Time-varying IPCW without cloning
4. Naive per-protocol

G-computation is intentionally excluded from this study.

## Updated federated protocol

For each interval `m`, site `k` shares only:

- `D[k,m]`: number initiating treatment in interval `m`
- `N[k,m]`: number still eligible to initiate at the start of interval `m`

The center calculates:

```text
h_num[m] = sum_k D[k,m] / sum_k N[k,m]
```

and broadcasts this common numerator vector to all sites. Each site retains
its own conditional initiation-denominator model and returns only weighted
death and risk-set totals. The pooled comparator uses the identical common
numerator.

This ensures that federated and pooled CCW target the same distribution of
initiation times within the grace period.

## Aligned oracle

The oracle uses:

- The true site-specific conditional initiation probabilities
- The common numerator hazard calculated from aggregated site counts
- Central aggregation of weighted death and risk-set totals
- Survival-curve construction only after aggregation

It no longer averages separately stabilized site survival curves.

## Heterogeneity design

| Level | Site 1 | Site 2 | Site 3 |
|---|---|---|---|
| Low | `x1 ~ N(-0.25, 0.95^2)`, `P(x2=1)=0.35` | `x1 ~ N(0, 1^2)`, `P(x2=1)=0.40` | `x1 ~ N(0.25, 1.05^2)`, `P(x2=1)=0.45` |
| Moderate | `x1 ~ N(-0.75, 0.85^2)`, `P(x2=1)=0.20` | `x1 ~ N(0, 1^2)`, `P(x2=1)=0.40` | `x1 ~ N(0.75, 1.15^2)`, `P(x2=1)=0.60` |
| High | `x1 ~ N(-1.50, 0.70^2)`, `P(x2=1)=0.10` | `x1 ~ N(0, 1^2)`, `P(x2=1)=0.40` | `x1 ~ N(1.50, 1.30^2)`, `P(x2=1)=0.70` |

Sites are equally sized, with pooled `E[x1]=0` and `E[x2]=0.4`.
Conditional treatment, covariate-transition, outcome, and treatment-effect
models remain common across sites.

## Simulation grid

- `tau`: 2, 5, 8
- `beta_trt`: -1.0, -0.7, -0.5
- Heterogeneity: low, moderate, high
- `t_star`: 25
- Patients per site: 1,000
- Replicates per cell: 100
- Oracle size per cell: 3,000,000
- SLURM array tasks: 27

## Cluster submission

From this directory:

```bash
chmod +x submit.sh job.sh run_sim_array.R aggregate.R validate_setup.R
./submit.sh
```

If necessary, change `module load R/4.3` in `job.sh` to match the cluster.
The submission wrapper creates `logs/` and `results/` before calling SLURM.

After all tasks complete:

```bash
Rscript aggregate.R
```

This writes `summary_all_scenarios.csv` and bias figures for RD, RR, OR, and
the RMST difference.

## Preflight validation

Before cluster submission:

```bash
Rscript validate_setup.R
```

The preflight uses a small simulation to verify:

- The common numerator is identical across sites
- Federated and pooled CCW are identical for `K=1`
- Their estimates align in a heterogeneous three-site smoke test
- Exactly four requested methods are returned
- No G-computation method is included
