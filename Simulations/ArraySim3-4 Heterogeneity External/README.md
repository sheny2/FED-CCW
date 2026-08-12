# ArraySim3-4 Heterogeneity External

One-shot federated CCW simulation under increasing site heterogeneity. Every
site uses a fixed stabilizing numerator generated from a separate external
reference cohort before the analysis begins.

## Communication

For point estimation, only one site-to-center communication is required:
sites send weighted event and risk-set sufficient statistics.

For standard errors and confidence intervals, a second round is required:
the center broadcasts global hazard/survival quantities and sites return
influence-function summaries.

The analysis cohort is never used to estimate or communicate a numerator.

## External reference numerator

The included reference is generated independently using:

- 3,000,000 patients
- Three equally sized sites
- Moderate heterogeneity patient mix
- Follow-up through `t_star=25`
- Seed `20260806`
- The same treatment-initiation and covariate-transition DGP
- No outcome data

The resulting vector is stored in:

- `external_reference_numerator.rds`: full reference object and metadata
- `external_reference_numerator.csv`: human-readable interval counts/hazards

The same numerator vector is reused across every target heterogeneity level,
grace period, treatment effect, simulation replicate, and site. This fixes
the target intervention while local denominator models adapt to each site's
patient mix.

To regenerate it exactly:

```bash
Rscript prepare_external_reference.R
```

To build a different reference:

```bash
Rscript prepare_external_reference.R \
  --N=3000000 \
  --t-star=25 \
  --heterogeneity=moderate \
  --seed=20260806
```

## Methods

Exactly four methods are reported:

1. Federated CCW using local denominators and the external numerator
2. Pooled CCW using a pooled denominator and the same external numerator
3. TV-IPCW without cloning using the external numerator
4. Naive per-protocol

G-computation is excluded.

## Simulation grid

- `tau`: 2, 5, 8
- `beta_trt`: -1.0, -0.7, -0.5
- Target heterogeneity: low, moderate, high
- `t_star`: 25
- Patients per target site: 1,000
- Replicates per cell: 100
- Oracle size per cell: 3,000,000
- SLURM array tasks: 27

The oracle uses true local denominator probabilities, the fixed external
numerator, and centrally aggregated weighted event/risk totals.

## Validate and submit

From this directory:

```bash
chmod +x prepare_external_reference.R validate_setup.R \
  submit.sh job.sh run_sim_array.R aggregate.R
Rscript validate_setup.R
./submit.sh
```

The submission script refuses to run if the external reference RDS is
missing. If necessary, update `module load R/4.3` in `job.sh`.

After all tasks finish:

```bash
Rscript aggregate.R
```

This produces `summary_all_scenarios.csv` and bias plots for RD, RR, OR, and
the RMST difference.
