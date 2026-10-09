# Code for: Wintermantel et al. - Equivalence testing in pesticide risk assessment: Evaluation and practical guidance for design, analysis and interpretation

This repository contains the scripts used for the analyses in:

> Wintermantel et al. - Equivalence testing in pesticide risk assessment: Evaluation and practical guidance for design, analysis and interpretation. Environment International

## Contents

| Path | Description |
|---|---|
| `R scripts/anticlustering_randomisation/experimental_allocation.R` | Balanced allocation of test subjects (anti-clustering randomisation) |
| `R scripts/risk_simulations_and_plotting/Read original data.R` | Loads the raw data |
| `R scripts/risk_simulations_and_plotting/Functions to simulate risk assessments.R` | Functions used by simulate risk.R and Create_manuscript_figures.R scripts |
| `R scripts/risk_simulations_and_plotting/Simulate risk.R` | Runs the simulation scenarios and saves the results |
| `R scripts/risk_simulations_and_plotting/Create_manuscript_figures.R` | Builds all figures from the results |
| `R scripts/risk_simulations_and_plotting/config.R` | Central configuration: paths, `read_results`, `recompute_missing`, `n_simulations`, `seed_number` |
| `R scripts/risk_simulations_and_plotting/run_all.R` | Convenience entry point; runs the two scripts below in order |

## Data

Raw data are **not** included in this repository. They may be requested through
Bayer's transparency initiative
(<https://www.bayer.com/en/agriculture/transparency-crop-science>).

To reproduce the analysis, place the file

```
control_data_Rolke_expo.rds
```

in `Data & Figures/Raw data/`. The scripts read it from there via `raw_data_file` in
[`config.R`](R%20scripts/risk_simulations_and_plotting/config.R). See
[`Data & Figures/Raw data/README.md`](Data%20&%20Figures/Raw%20data/README.md) for details.

## Reproducibility note

Everything in `Data & Figures/` is excluded from version control (see
[`.gitignore`](.gitignore)), so neither the raw data nor any generated results or
figures are committed.

The simulation scenarios are computationally expensive: a single scenario runs 5000
simulations and can take days to weeks. The scripts therefore default to **reusing
existing results** rather than recomputing them:

- `read_results <- TRUE` — load results from `Data & Figures/Simulation results/` instead
  of recomputing them.
- `recompute_missing <- FALSE` — stop with an error when a result file is missing,
  instead of silently starting a computation that can take days to weeks.

With these defaults, `run_all.R` only reads existing results. To actually recompute
something, set `recompute_missing <- TRUE` in `config.R` deliberately.

## Running the analysis

From the repository root:

```r
source("R scripts/risk_simulations_and_plotting/run_all.R")
```

Or run the scripts individually:

```r
source("R scripts/risk_simulations_and_plotting/Simulate risk.R")
source("R scripts/risk_simulations_and_plotting/Create_manuscript_figures.R")
```

Both scripts source `config.R` themselves, so they can be run in any order and from
any working directory. Figures are written to `Data & Figures/Plots/`.

## Session information

The R session used for the analysis is recorded in
[`sessionInfo.txt`](sessionInfo.txt).
