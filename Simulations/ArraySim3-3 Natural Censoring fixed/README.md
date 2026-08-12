# ArraySim3-3 Natural Censoring fixed

Cluster-ready natural-censoring simulation with the method definitions aligned
to `ArraySim3-3 Direct G-comp fixed`.

## Methods

1. Federated clone-censor-weight with artificial and natural-censoring IPCW
2. Pooled clone-censor-weight
3. Clone-censor g-computation using stacked, arm-specific clone frames
4. Plain g-computation without cloning
5. IPCW without cloning
6. Naive per-protocol analysis

## Scenario grid

- `tau`: 2, 5, 8
- treatment log-odds effect: -1.0, -0.7, -0.5
- confounding: small, medium, strong
- 27 SLURM array tasks
- 100 replications per task
- 3 sites and 1,000 subjects per site
- fixed follow-up horizon `t_star = 25`

Each task writes 4,800 rows to `results/res_task_XXX.rds`.

## Cluster use

From this directory:

```bash
./submit.sh
```

`submit.sh` creates `logs/` and `results/` before submitting `job.sh`. Adjust
the R module, memory, or CPU request in `job.sh` if required by the cluster.

After all 27 tasks finish:

```bash
Rscript aggregate.R
```

The aggregation step requires the R packages `dplyr` and `ggplot2`.
