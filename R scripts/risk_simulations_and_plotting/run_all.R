# Run the full analysis pipeline ------------------------------------------------
#
# Usage (from the repository root):
#   Rscript "R scripts/risk_simulations_and_plotting/run_all.R"
#
# The scripts are run in order:
#   1. Simulate risk.R              - runs the simulations (or loads cached results)
#   2. Create_manuscript_figures.R  - builds all figures and tables
#
# Settings such as `read_results`, `recompute_missing` and `n_simulations` are in
# config.R. Raw data are not part of this repository; see README.md.
#
# NOTE: with the default settings (`read_results = TRUE`,
# `recompute_missing = FALSE`) this script only *reads* existing results. It stops
# with an error if a result file is missing, rather than starting a computation
# that can take days to weeks. Set `recompute_missing <- TRUE` in config.R only
# when you really intend to recompute.

library(here)

script_root <- here::here("R scripts", "risk_simulations_and_plotting")

source(file.path(script_root, "config.R"))

scripts <- c(
  file.path(script_root, "Simulate risk.R"),
  file.path(script_root, "Create_manuscript_figures.R")
)

for (script in scripts) {
  message("\n=== Running ", script, " ===")
  source(script)
}

message("\nDone. Figures were written to ", plots_root)
