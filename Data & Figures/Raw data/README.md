# Raw data

The raw data are **not** part of this repository. To reproduce the analysis,
place the file

```
control_data_Rolke_expo.rds
```

in this folder. The scripts read it from here via `raw_data_file` in
[`config.R`](../../R%20scripts/risk_simulations_and_plotting/config.R).

Everything in `Data & Figures/` is excluded from version control (see
[`.gitignore`](../../.gitignore)), so data and generated results are never
committed.
