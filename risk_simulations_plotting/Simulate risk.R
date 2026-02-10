library(here)
library(cowplot)
library(scales)
library(tidyverse)
library(tidyr)
library(openxlsx)

source(here::here("R scripts", "Read original data.R"))
source(here::here("R scripts", "Functions to simulate risk assessments.R"))

# Set variables --------
save_results <- TRUE

seed_number <- 1
n_simulations <- 5000

alpha_d <- 0.05
alpha_e <- 0.2
alpha_Hotopp_control <- 0.1
delta_e <- 0.1

tests = c("diff", "equi", "Hotopp")

# effect_sizes_large_range <- seq(0, 0.4, 0.02)
effect_sizes_small_range <- seq(0, 0.2, 0.02)
effect_sizes_pareto <- c(0, 0.05, 0.10, 0.11, 0.12, 0.15)
effect_sizes_allocation_vs_inclusion <- c(0, 0.05, 0.07, 0.1, 0.12, 0.15) 
effect_sizes_high_risk_plot <- c(0.1, 0.12, 0.15, 0.2, 0.25, 0.3, 0.35); length(effect_sizes_high_risk_plot)
effect_sizes_high_risk <- c(0.12, 0.15, 0.2, 0.35); length(effect_sizes_high_risk_plot)
# effect_sizes_high_risk_plot <- c(seq(0.1, 0.2, 0.02), 0.35); length(effect_sizes_high_risk_plot)

# effect_sizes_allo_risk_combined <- unique(c(effect_sizes_allocation_vs_inclusion, 
#                                             effect_sizes_high_risk_plot))

# effect_sizes_tested <- sort(c(seq(0, 0.2, 0.02))); length(effect_sizes_tested)

nn_sites_high_risk_plot <- seq(4, 20, 2)
nn_sites_allocation_vs_inclusion <- c(seq(4, 20, 2), seq(25, 45, 5))

var_comps_glmm_main <- obtain_var_comps(original_data = original_control_data_234)

## Compare results using dataset simulated in parametric bootstrap to original data (permuted) ----------
### Original data
sim_5000_results_real <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = FALSE,
  n_sites = 3,
  n_colonies = 8,
  reallocate_colonies = FALSE, 
  permute = T, 
  formula = n_bees ~ Treatment*Assessment + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = alpha_e,
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = delta_e,
  effect_size = effect_sizes_small_range,
  effect_timing = "abrupt",
  tests = tests,
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number
) %>% mutate(Study = "Rolke real")

if(save_results) {
  saveRDS(sim_5000_results_real, "R output/Simulation results/sim_5000_results_real.rds")
}

### Simulated data 
sim_5000_results_simulated <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = TRUE,
  var_comps_list = var_comps_glmm_main, 
  n_sites = 3,
  n_colonies = 8,
  reallocate_colonies = FALSE, 
  permute = T, 
  formula = n_bees ~ Treatment*Assessment + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = alpha_e,
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = delta_e,
  effect_size = effect_sizes_small_range,
  effect_timing = "abrupt",
  tests = tests,
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number
) %>% mutate(Study = "Rolke simulated: Assessment 2, 3, 4")

if(save_results) {
  saveRDS(sim_5000_results_simulated, "R output/Simulation results/sim_5000_results_simulated.rds")
}

# Compare false trust rate of tests when control SD is inflated / deflated ----------

results_5000_sims_varying_control_sd_10_perc <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = T,
  var_comps_list = var_comps_glmm_main, 
  n_sites = 10,
  n_colonies = 6,
  reallocate_colonies = F, 
  formula = n_bees ~ Treatment*Assessment + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = alpha_e,
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = delta_e,
  control_sd_factor = seq(0.5, 1.5, 0.1),
  effect_size = 0.1,
  effect_timing = "abrupt",
  tests = tests,
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number) %>% 
  mutate(Study = "Rolke simulated: Assessment 2, 3, 4",
         Predictors = "Treatment * Assessment + (1|Site/Colony)") 

if(save_results) {
  saveRDS(results_5000_sims_varying_control_sd_10_perc, 
          "R output/Simulation results/results_5000_sims_varying_control_sd_10_perc.rds")
}


## Compare effectiveness of different tests and determine influence of sample size --------------
sim_5000_results_pareto_alphas <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = T,
  var_comps_list = var_comps_glmm_main, 
  n_sites = 10,
  n_colonies = 6,
  reallocate_colonies = F, 
  formula = n_bees ~ Treatment*Assessment + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = c(seq(0.01, 0.04, 0.01), seq(0.05, 0.5, 0.05)),
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = 0.1,
  effect_size = effect_sizes_pareto,
  effect_timing = "abrupt",
  tests = tests,
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number)

saveRDS(sim_5000_results_pareto_alphas, "R output/Simulation results/sim_5000_results_pareto_alphas.rds")

sim_5000_results_pareto_alphas_small_sample_size <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = T,
  var_comps_list = var_comps_glmm_main, 
  n_sites = 4,
  n_colonies = 4,
  reallocate_colonies = F, 
  formula = n_bees ~ Treatment*Assessment + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = c(seq(0.01, 0.04, 0.01), seq(0.05, 0.5, 0.05)),
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = 0.1,
  effect_size = effect_sizes_pareto,
  effect_timing = "abrupt",
  tests = tests,
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number)

saveRDS(sim_5000_results_pareto_alphas_small_sample_size, "R output/Simulation results/sim_5000_results_pareto_alphas_small_sample_size.rds")

### High risk 
sim_5000_results_high_risk <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = T,
  var_comps_list = var_comps_glmm_main, 
  n_sites = nn_sites_high_risk_plot,
  n_colonies = 6,
  reallocate_colonies = F, 
  formula = n_bees ~ Treatment*Assessment + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = alpha_e,
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = delta_e,
  effect_size = effect_sizes_high_risk,
  effect_timing = "abrupt",
  tests = tests,
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number) %>% 
  mutate(Study = "Rolke simulated: Assessment 2, 3, 4",
         Predictors = "Treatment * Assessment + (1|Site/Colony)") 

if(save_results) {
  saveRDS(sim_5000_results_high_risk, 
          "R output/Simulation results/sim_5000_results_high_risk.rds")
}

## How are risk classifications affected by accounting for the initial number of bees ----------------
## either through inclusion as a covariate and / or through a colony allocation process 
## that balances the initial number of bees per colony across sites

### Not accounting for initial number of bees per colony
sim_5000_results_n_bees_initial_not_accounted <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = T,
  var_comps_list = var_comps_glmm_main, 
  n_sites = nn_sites_allocation_vs_inclusion,
  n_colonies = 6,
  reallocate_colonies = F, 
  formula = n_bees ~ Treatment*Assessment + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = alpha_e,
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = delta_e,
  effect_size = effect_sizes_allocation_vs_inclusion,
  effect_timing = "abrupt",
  tests = tests,
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number) %>% 
  mutate(Study = "Rolke simulated: Assessment 2, 3, 4",
         Predictors = "Treatment * Assessment + (1|Site/Colony)") 

if(save_results) {
  saveRDS(sim_5000_results_n_bees_initial_not_accounted, 
          "R output/Simulation results/sim_5000_results_n_bees_initial_not_accounted.rds")
}

if(save_results) {
  
  sim_5000_results_all_not_accounted <- bind_rows(sim_5000_results_high_risk, sim_5000_results_n_bees_initial_not_accounted)
  
  saveRDS(sim_5000_results_all_not_accounted, 
          "R output/Simulation results/sim_5000_results_all_not_accounted.rds")
}

### Accounting for initial number of bees per colony 
### through balanced colony allocation process
sim_5000_results_n_bees_initial_balanced <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = T,
  var_comps_list = var_comps_glmm_main, 
  n_sites = nn_sites_allocation_vs_inclusion,
  n_colonies = 6,
  reallocate_colonies = TRUE,
  covariates = "n_bees_initial", 
  formula = n_bees ~ Treatment*Assessment + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = alpha_e,
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = delta_e,
  effect_size = effect_sizes_allocation_vs_inclusion,
  effect_timing = "abrupt",
  tests = c("equi", "diff"),
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number) %>% 
  mutate(Study = "Rolke simulated: Assessment 2, 3, 4",
         Predictors = "Treatment * Assessment + (1|Site/Colony)")

if(save_results) {
  saveRDS(sim_5000_results_n_bees_initial_balanced, 
          "R output/Simulation results/sim_5000_results_n_bees_initial_balanced.rds")
}

### Accounting for initial number of bees per colony 
### through inclusion of covariate
sim_5000_results_n_bees_initial_included <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = T,
  var_comps_list = var_comps_glmm_main, 
  n_sites = nn_sites_allocation_vs_inclusion,
  n_colonies = 6,
  reallocate_colonies = F, 
  formula = n_bees ~ Treatment*Assessment + n_bees_initial + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = alpha_e,
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = delta_e,
  effect_size = effect_sizes_allocation_vs_inclusion,
  effect_timing = "abrupt",
  tests = c("equi", "diff"),
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number)  %>% 
  mutate(Study = "Rolke simulated: Assessment 2, 3, 4",
         Predictors = "Treatment * Assessment + n_bees_initial + (1|Site/Colony)")

if(save_results) {
  saveRDS(sim_5000_results_n_bees_initial_included, 
          "R output/Simulation results/sim_5000_results_n_bees_initial_included.rds")
}

### Accounting for initial number of bees per colony 
### through balanced allocation and inclusion of covariate 
sim_5000_results_n_bees_initial_included_and_balanced <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = T,
  var_comps_list = var_comps_glmm_main, 
  n_sites = nn_sites_allocation_vs_inclusion,
  n_colonies = 6,
  reallocate_colonies = TRUE,
  covariates = "n_bees_initial", 
  formula = n_bees ~ Treatment*Assessment + n_bees_initial + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = alpha_e,
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = delta_e,
  effect_size = effect_sizes_allocation_vs_inclusion,
  effect_timing = "abrupt",
  tests = c("equi", "diff"),
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number) %>% 
  mutate(Study = "Rolke simulated: Assessment 2, 3, 4",
         Predictors = "Treatment * Assessment + n_bees_initial + (1|Site/Colony)")

if(save_results) {
  saveRDS(sim_5000_results_n_bees_initial_included_and_balanced, 
          "R output/Simulation results/sim_5000_results_n_bees_initial_included_and_balanced.rds")
}

# random reallocation
sim_5000_results_n_bees_initial_not_accounted_2 <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = T,
  var_comps_list = var_comps_glmm_main, 
  n_sites = nn_sites_allocation_vs_inclusion,
  n_colonies = 6,
  reallocate_colonies = T, 
  formula = n_bees ~ Treatment*Assessment + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = alpha_e,
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = delta_e,
  effect_size = effect_sizes_allocation_vs_inclusion,
  effect_timing = "abrupt",
  tests = c("equi", "diff"),
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number) %>% 
  mutate(Study = "Rolke simulated: Assessment 2, 3, 4",
         Predictors = "Treatment * Assessment + (1|Site/Colony)") 

if(save_results) {
  saveRDS(sim_5000_results_n_bees_initial_not_accounted_2, 
          "R output/Simulation results/sim_5000_results_n_bees_initial_not_accounted_2.rds")
}

# Refine to determine number of sites to meet 2013 guidance ------------
sim_5000_results_no_adjust_7_perc <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = T,
  var_comps_list = var_comps_glmm_main, 
  n_sites = c(40, 41, 42),
  n_colonies = 6,
  reallocate_colonies = F, 
  formula = n_bees ~ Treatment*Assessment + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = alpha_e,
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = delta_e,
  effect_size = 0.07,
  effect_timing = "abrupt",
  tests = tests,
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number) %>% 
  mutate(Study = "Rolke simulated: Assessment 2, 3, 4",
         Predictors = "Treatment * Assessment + (1|Site/Colony)") 

if(save_results) {
  saveRDS(sim_5000_results_no_adjust_7_perc, 
          "R output/Simulation results/sim_5000_results_no_adjust_7_perc.rds")
}

sim_5000_results_design_adjust_7_perc <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = T,
  var_comps_list = var_comps_glmm_main, 
  n_sites = c(29, 30, 31),
  n_colonies = 6,
  reallocate_colonies = TRUE,
  covariates = "n_bees_initial", 
  formula = n_bees ~ Treatment*Assessment + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = alpha_e,
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = delta_e,
  effect_size = effect_sizes_allocation_vs_inclusion,
  effect_timing = "abrupt",
  tests = c("equi", "diff"),
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number) %>% 
  mutate(Study = "Rolke simulated: Assessment 2, 3, 4",
         Predictors = "Treatment * Assessment + (1|Site/Colony)")

if(save_results) {
  saveRDS(sim_5000_results_design_adjust_7_perc, 
          "R output/Simulation results/sim_5000_results_design_adjust_7_perc.rds")
}

sim_5000_results_covariate_adjust <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = T,
  var_comps_list = var_comps_glmm_main, 
  n_sites = c(29,30),
  n_colonies = 6,
  reallocate_colonies = F, 
  formula = n_bees ~ Treatment*Assessment + n_bees_initial + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = alpha_e,
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = delta_e,
  effect_size = 0.07,
  effect_timing = "abrupt",
  tests = c("equi", "diff"),
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number) %>% 
  mutate(Study = "Rolke simulated: Assessment 2, 3, 4",
         Predictors = "Treatment * Assessment + n_bees_initial + (1|Site/Colony)")

if(save_results) {
  saveRDS(sim_5000_results_covariate_adjust_7_perc, 
          "R output/Simulation results/sim_5000_results_covariate_adjust_7_perc.rds")
}

sim_5000_results_combined_adjust_7_perc <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = T,
  var_comps_list = var_comps_glmm_main, 
  n_sites = c(29,30),
  n_colonies = 6,
  reallocate_colonies = TRUE,
  covariates = "n_bees_initial", 
  formula = n_bees ~ Treatment*Assessment + n_bees_initial + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = alpha_e,
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = delta_e,
  effect_size = 0.07,
  effect_timing = "abrupt",
  tests = c("equi", "diff"),
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number) %>% 
  mutate(Study = "Rolke simulated: Assessment 2, 3, 4",
         Predictors = "Treatment * Assessment + n_bees_initial + (1|Site/Colony)")

if(save_results) {
  saveRDS(sim_5000_results_combined_adjust_7_perc, 
          "R output/Simulation results/sim_5000_results_combined_adjust_7_perc.rds")
}

# Refine to determine number of sites to meet 2013 guidance ------------
sim_5000_results_no_adjust_7_perc_v2 <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = T,
  var_comps_list = var_comps_glmm_main, 
  n_sites = c(37, 38, 39),
  n_colonies = 6,
  reallocate_colonies = F, 
  formula = n_bees ~ Treatment*Assessment + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = alpha_e,
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = delta_e,
  effect_size = 0.07,
  effect_timing = "abrupt",
  tests = tests,
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number) %>% 
  mutate(Study = "Rolke simulated: Assessment 2, 3, 4",
         Predictors = "Treatment * Assessment + (1|Site/Colony)") 

if(save_results) {
  saveRDS(sim_5000_results_no_adjust_7_perc_v2, 
          "R output/Simulation results/sim_5000_results_no_adjust_7_perc_v2.rds")
}

sim_5000_results_design_adjust_7_perc_v2 <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = T,
  var_comps_list = var_comps_glmm_main, 
  n_sites = c(26, 27, 28),
  n_colonies = 6,
  reallocate_colonies = TRUE,
  covariates = "n_bees_initial", 
  formula = n_bees ~ Treatment*Assessment + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = alpha_e,
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = delta_e,
  effect_size = effect_sizes_allocation_vs_inclusion,
  effect_timing = "abrupt",
  tests = c("equi", "diff"),
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number) %>% 
  mutate(Study = "Rolke simulated: Assessment 2, 3, 4",
         Predictors = "Treatment * Assessment + (1|Site/Colony)")

if(save_results) {
  saveRDS(sim_5000_results_design_adjust_7_perc_v2, 
          "R output/Simulation results/sim_5000_results_design_adjust_7_perc_v2.rds")
}

sim_5000_results_covariate_adjust_7_perc_v2 <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = T,
  var_comps_list = var_comps_glmm_main, 
  n_sites = c(26, 27, 28),
  n_colonies = 6,
  reallocate_colonies = F, 
  formula = n_bees ~ Treatment*Assessment + n_bees_initial + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = alpha_e,
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = delta_e,
  effect_size = 0.07,
  effect_timing = "abrupt",
  tests = c("equi", "diff"),
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number) %>% 
  mutate(Study = "Rolke simulated: Assessment 2, 3, 4",
         Predictors = "Treatment * Assessment + n_bees_initial + (1|Site/Colony)")

if(save_results) {
  saveRDS(sim_5000_results_covariate_adjust_7_perc_v2, 
          "R output/Simulation results/sim_5000_results_covariate_adjust_7_perc_v2.rds")
}

sim_5000_results_combined_adjust_7_perc_v2 <- run_simulations(
  n_simulations = n_simulations,
  data = original_control_data_234,
  mimic_data = T,
  var_comps_list = var_comps_glmm_main, 
  n_sites = c(26, 27, 28),
  n_colonies = 6,
  reallocate_colonies = TRUE,
  covariates = "n_bees_initial", 
  formula = n_bees ~ Treatment*Assessment + n_bees_initial + (1|Site/Colony),
  family = "nbinom2",
  alpha_d = alpha_d,
  alpha_e = alpha_e,
  alpha_Hotopp_control = alpha_Hotopp_control,
  delta_e = delta_e,
  effect_size = 0.07,
  effect_timing = "abrupt",
  tests = c("equi", "diff"),
  assessment_scope = "both",
  parallel = T,
  seed_number = seed_number) %>% 
  mutate(Study = "Rolke simulated: Assessment 2, 3, 4",
         Predictors = "Treatment * Assessment + n_bees_initial + (1|Site/Colony)")

if(save_results) {
  saveRDS(sim_5000_results_combined_adjust_7_perc_v2, 
          "R output/Simulation results/sim_5000_results_combined_adjust_7_perc_v2.rds")
}

if (read_results) {
  sim_5000_results_no_adjust_7_perc <- readRDS("R output/Simulation results/sim_5000_results_no_adjust_7_perc_2026.rds")
  sim_5000_results_covariate_adjust_7_perc <- readRDS("R output/Simulation results/sim_5000_results_covariate_adjust_7_perc.rds")
  sim_5000_results_design_adjust_7_perc <- readRDS("R output/Simulation results/sim_5000_results_design_adjust_7_perc.rds")
  sim_5000_results_combined_adjust_7_perc  <- readRDS("R output/Simulation results/sim_5000_results_combined_adjust_7_perc.rds")
  
  sim_5000_results_no_adjust_7_perc_v2 <- readRDS("R output/Simulation results/sim_5000_results_no_adjust_7_perc_v2.rds")
  sim_5000_results_covariate_adjust_7_perc_v2 <- readRDS("R output/Simulation results/sim_5000_results_covariate_adjust_7_perc_v2.rds")
  sim_5000_results_design_adjust_7_perc_v2 <- readRDS("R output/Simulation results/sim_5000_results_design_adjust_7_perc_v2.rds")
  sim_5000_results_combined_adjust_7_perc_v2  <- readRDS("R output/Simulation results/sim_5000_results_combined_adjust_7_perc_v2.rds")
}

sim_5000_results_EFSA_GD_refinement <- bind_rows(
  sim_5000_results_no_adjust_7_perc,
  sim_5000_results_covariate_adjust_7_perc,
  sim_5000_results_design_adjust_7_perc,
  sim_5000_results_combined_adjust_7_perc,
  
  sim_5000_results_no_adjust_7_perc_v2,
  sim_5000_results_covariate_adjust_7_perc_v2,
  sim_5000_results_design_adjust_7_perc_v2,
  sim_5000_results_combined_adjust_7_perc_v2
)

if(save_results) {
  saveRDS(sim_5000_results_EFSA_GD_refinement, 
          "R output/Simulation results/sim_5000_results_EFSA_GD_refinement.rds")
}

# simu_results_n_bees_initial_included_2 <- run_simulations(
#   n_simulations = n_simulations,
#   data = original_control_data_234,
#   mimic_data = T,
#   var_comps_list = var_comps_glmm_main, 
#   n_sites = nn_sites_allocation_vs_inclusion,
#   n_colonies = 6,
#   reallocate_colonies = T, 
#   formula = n_bees ~ Treatment*Assessment + n_bees_initial + (1|Site/Colony),
#   family = "nbinom2",
#   alpha_d = alpha_d,
#   alpha_e = alpha_e,
#   alpha_Hotopp_control = alpha_Hotopp_control,
#   delta_e = delta_e,
#   effect_size = effect_sizes_allocation_vs_inclusion,
#   effect_timing = "abrupt",
#   tests = c("equi", "diff"),
#   assessment_scope = "both",
#   parallel = T,
#   seed_number = seed_number)  %>% 
#   mutate(Study = "Rolke simulated: Assessment 2, 3, 4",
#          Predictors = "Treatment * Assessment + n_bees_initial + (1|Site/Colony)")
# 
# if(save_results) {
#   saveRDS(simu_results_n_bees_initial_included_2, 
#           "R output/Simulation results/simu_results_n_bees_initial_included_2.rds")
# }


# ### vary sample size -------
# for(n_colonies in nn_colonies) {
#   simu_results_expo_one_colony <- run_simulations(
#     n_simulations = n_simulations,
#     data = original_control_data_234,
#     mimic_data = T,
#     n_sites = nn_sites,
#     n_colonies = n_colonies,
#     reallocate_colonies = F, 
#     formula = n_bees ~ Treatment*Assessment + (1|Site/Colony),
#     family = "nbinom2",
#     alpha_d = alpha_d,
#     alpha_e = alpha_e,
#     alpha_Hotopp_control = alpha_Hotopp_control,
#     delta_e = delta_e,
#     effect_size = effect_sizes_tested,
#     effect_timing = "abrupt",
#     tests = tests,
#     assessment_scope = "both",
#     parallel = T,
#     seed_number = seed_number
#   ) %>% mutate(Study = "Rolke simulated: Assessment 2, 3, 4")
#   
#   if(save_results) {
#     saveRDS(simu_results_expo_one_colony, paste0("R output/Simulation results/simu_results_expo_", n_colonies, "_colonies.rds"))
#   }
# } 
# 
# simu_results_expo <- data.frame()
# for(n_colonies in nn_colonies) {
#   simu_results_expo_one_colony <- readRDS(paste0("R output/Simulation results/simu_results_expo_", n_colonies, "_colonies.rds"))
#   simu_results_expo <- bind_rows(simu_results_expo_one_colony, simu_results_expo)
# }
# 
# if(save_results) {
#   saveRDS(simu_results_expo, "R output/Simulation results/simu_results_expo.rds")
# }
# 
# ### check why alpha error is not as wanted
# 
# simu_results_no_realloction <- run_simulations(
#   n_simulations = 1000,
#   data = original_control_data_234,
#   mimic_data = T,
#   n_sites = 10,
#   n_colonies = 6,
#   reallocate_colonies = FALSE, 
#   formula = n_bees ~ Treatment*Assessment + (1|Site/Colony),
#   family = "nbinom2",
#   alpha_d = alpha_d,
#   alpha_e = alpha_e,
#   alpha_Hotopp_control = alpha_Hotopp_control,
#   delta_e = delta_e,
#   effect_size = 0.1,
#   effect_timing = "abrupt",
#   tests = tests,
#   assessment_scope = "both",
#   parallel = T,
#   seed_number = seed_number) %>% 
#   mutate(Study = "Rolke simulated: Assessment 2, 3, 4",
#          Predictors = "Treatment * Assessment + (1|Site/Colony)") 
# 
# saveRDS(simu_results_no_realloction, "R output/Simulation results/simu_results_no_realloction.rds")
# 
# simu_results_random_realloction <- run_simulations(
#   n_simulations = 1000,
#   data = original_control_data_234,
#   mimic_data = T,
#   n_sites = 10,
#   n_colonies = 6,
#   reallocate_colonies = T, 
#   formula = n_bees ~ Treatment*Assessment + (1|Site/Colony),
#   family = "nbinom2",
#   alpha_d = alpha_d,
#   alpha_e = alpha_e,
#   alpha_Hotopp_control = alpha_Hotopp_control,
#   delta_e = delta_e,
#   effect_size = 0.1,
#   effect_timing = "abrupt",
#   tests = tests,
#   assessment_scope = "both",
#   parallel = T,
#   seed_number = seed_number) %>% 
#   mutate(Study = "Rolke simulated: Assessment 2, 3, 4",
#          Predictors = "Treatment * Assessment + (1|Site/Colony)") 
# 
# saveRDS(simu_results_random_realloction, "R output/Simulation results/simu_results_random_realloction.rds")
