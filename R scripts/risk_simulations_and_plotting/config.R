# Configuration for the analysis scripts ---------------------------------------
#
# All paths are resolved relative to the repository root, so the scripts can be
# run from anywhere (via run_all.R, or by sourcing an individual script).
#
# Raw data are not part of this repository. Place the file
#   control_data_Rolke_expo.rds
# in `Data & Figures/Raw data/` before running the scripts (see README.md).

library(here)

# Guard: config.R must be sourced from within this repository. `here::here()`
# resolves to the nearest project root upwards from the working directory; if the
# R session was started in a different repository (e.g. the parent working
# folder), paths would silently point at the wrong `Data & Figures/` folder.
expected_rproj <- "equivalence-testing-and-balanced-allocation-paper.Rproj"
if (!file.exists(here::here(expected_rproj))) {
  stop(
    "config.R was sourced from the wrong project root.\n",
    "here::here() resolved to: ", here::here(), "\n",
    "Set the working directory to the repository root first, e.g.:\n",
    "  setwd(\"<path to>/equivalence-testing-and-balanced-allocation-paper\")\n",
    "or open ", expected_rproj, " in RStudio, then restart R.",
    call. = FALSE
  )
}

# Root folder containing the R scripts
scripts_root <- here::here("R scripts")

# Folders for data and generated output (not tracked by git)
output_root  <- here::here("Data & Figures")
data_root    <- file.path(output_root, "Raw data")
results_root <- file.path(output_root, "Simulation results")
plots_root   <- file.path(output_root, "Plots")

raw_data_file <- file.path(data_root, "control_data_Rolke_expo.rds")

# Analysis settings ------------------------------------------------------------
# read_results = TRUE  -> reuse existing results from `results_root`
# read_results = FALSE -> recompute them (slow: 5000 simulations per scenario)
read_results <- TRUE

# save_results = TRUE -> write newly computed results to `results_root`
save_results <- TRUE

# recompute_missing = FALSE -> stop with an error when a result file is missing
# recompute_missing = TRUE  -> recompute it instead (very slow: a single scenario
#                              can take days to weeks at n_simulations = 5000)
recompute_missing <- FALSE

n_simulations <- 5000
seed_number   <- 1

# Helper functions to build paths inside the output folders ---------------------
results_path <- function(...) file.path(results_root, ...)
plots_path   <- function(...) file.path(plots_root, ...)

# Load a result from disk, or compute it when it is not available ---------------
# `load_or_compute(name, expr)` returns the object stored in
# `<results_root>/<name>.rds` when `read_results` is TRUE and that file exists.
# Otherwise `expr` is evaluated (only then), saved when `save_results` is TRUE and
# returned. This makes `read_results = TRUE` skip the simulations for every
# scenario, not just for a subset of them.
#
# `derived = TRUE` marks objects that are cheaply derived from other results (for
# example a `bind_rows()` of already-loaded scenarios). Those are always recomputed
# when the file is absent, because doing so costs nothing.
#
# For real simulations the file is required: when it is missing and
# `recompute_missing` is FALSE the function stops instead of silently starting a
# computation that can take days to weeks.
load_or_compute <- function(name, expr, derived = FALSE) {
  path <- results_path(paste0(name, ".rds"))
  if (read_results && file.exists(path)) {
    message("Reading ", basename(path))
    return(readRDS(path))
  }
  if (!derived && !recompute_missing) {
    stop("Result file not found: ", path, "\n",
         "Set `recompute_missing <- TRUE` in config.R to recompute it ",
         "(this can take days to weeks at n_simulations = ", n_simulations, ").",
         call. = FALSE)
  }
  message("Computing ", name)
  value <- eval(substitute(expr), parent.frame())
  if (save_results) saveRDS(value, path)
  value
}

# Create the output folders if they do not exist yet
for (path in c(data_root, results_root, plots_root)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
