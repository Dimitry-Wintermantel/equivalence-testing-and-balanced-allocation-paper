# Load and prepare the original control data ------------------------------------
#
# The raw data file is not part of this repository. Place
# `control_data_Rolke_expo.rds` in `Data & Figures/Raw data/` before running this
# script (see README.md). The path is defined in config.R.

library(here)
source(here::here("R scripts", "risk_simulations_and_plotting", "config.R"))

library(tidyverse)

if (!file.exists(raw_data_file)) {
  stop("Raw data file not found: ", raw_data_file,
       "\nPlace `control_data_Rolke_expo.rds` in `Data & Figures/Raw data/` (see README.md).",
       call. = FALSE)
}

original_control_data <- readRDS(raw_data_file)
original_control_data$n_bees <- round(original_control_data$n_bees)
original_control_data$n_bees_initial <- round(original_control_data$n_bees_initial)
original_control_data$log_n_bees_initial <- log(original_control_data$n_bees_initial)
original_control_data <- rename(original_control_data, Site = Location)

original_control_data_234 <- subset(original_control_data, Assessment %in% c("2", "3", "4"))
