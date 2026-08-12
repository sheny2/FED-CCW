# ArraySim3-4: Site heterogeneity

Cluster-ready simulation derived from `ArraySim3-3 Direct G-comp fixed`.
It studies increasingly heterogeneous patient populations across three sites
while holding the conditional data-generating models and pooled baseline
means fixed.

## Methods

Exactly five methods are run:

1. Federated clone-censor-weight with time-varying IPCW
2. Pooled clone-censor-weight with time-varying IPCW
3. Direct g-computation without cloning
4. Time-varying IPCW without cloning (naive comparator)
5. Per-protocol without weighting (naive comparator)

The clone-censor/indirect g-computation variant is intentionally excluded.

## Heterogeneity design

All sites have 1,000 patients per replicate. The grid varies the baseline
distribution of continuous `x1` and binary `x2`:

| Level | Site 1 | Site 2 | Site 3 |
|---|---|---|---|
| Low | `x1 ~ N(-0.25, 0.95^2)`, `P(x2=1)=0.35` | `x1 ~ N(0, 1^2)`, `P(x2=1)=0.40` | `x1 ~ N(0.25, 1.05^2)`, `P(x2=1)=0.45` |
| Moderate | `x1 ~ N(-0.75, 0.85^2)`, `P(x2=1)=0.20` | `x1 ~ N(0, 1^2)`, `P(x2=1)=0.40` | `x1 ~ N(0.75, 1.15^2)`, `P(x2=1)=0.60` |
| High | `x1 ~ N(-1.50, 0.70^2)`, `P(x2=1)=0.10` | `x1 ~ N(0, 1^2)`, `P(x2=1)=0.40` | `x1 ~ N(1.50, 1.30^2)`, `P(x2=1)=0.70` |

The three equally sized sites retain pooled `E[x1]=0` and `E[x2]=0.4`.
Thus, differences across levels are driven by greater between-site patient-mix
separation rather than a shifted pooled mean.

## Simulation grid

- Grace period `tau`: 2, 5, 8
- Treatment log-effect `beta_trt`: -1.0, -0.7, -0.5
- Site heterogeneity: low, moderate, high
- Follow-up horizon `t_star`: 25
- Replicates per cell: 100
- Total SLURM array tasks: 27

Each task computes its own three-site oracle using 3,000,000 Monte Carlo
patients, then shares that oracle across the task's 100 replicates.

## Cluster use

From this directory:

```bash
chmod +x submit.sh job.sh run_sim_array.R aggregate.R
./submit.sh
```

`submit.sh` creates `logs/` and `results/` before calling `sbatch`. If the
cluster uses a different R module, update `module load R/4.3` in `job.sh`.

After all 27 tasks finish:

```bash
Rscript aggregate.R
```

This produces `summary_all_scenarios.csv` and four bias plots.

## Fed-versus-pooled diagnostics

Run the quick diagnostic suite for a high-heterogeneity cell:

```bash
Rscript diagnostics.R
```

For sample-size convergence and a larger oracle-alignment check:

```bash
Rscript diagnostics.R --heterogeneity=high --tau=5 --beta-trt=-0.7 \
  --run-convergence=true --convergence-reps=20 \
  --run-oracle=true --oracle-N=300000
```

Outputs are written to `diagnostics/`. The script checks one-site identity,
a fully homogeneous control, local versus common numerator hazards,
propensity-score support, weights and effective sample sizes, fitted versus
true denominator hazards, numerator/denominator decomposition, truncation
sensitivity, and optional sample-size/oracle convergence.

For a local smoke test without running the production-sized oracle, source
the four R files and call `run_once_tv()` with small `n_per_site`, `t_star`,
`mc_gcomp`, and `truth_N` values.
