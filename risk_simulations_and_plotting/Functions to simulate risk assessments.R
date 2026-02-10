library(tidyverse)
library(glmmTMB)
library(emmeans)
library(data.table)
library(doParallel)
library(foreach)
library(tibble)
library(lme4)
library(lmtest)
library(anticlust)

source(here::here("R scripts", "allocate_treatments.R"))
source(here::here("R scripts", "allocate_within_treatments.R"))

obtain_var_comps <- function(original_data, 
                             use_n_bees_initial = c("main", "no", "interaction"), 
                             use_glmm = TRUE) {
  
  # Check argument values
  use_n_bees_initial <- match.arg(use_n_bees_initial, c("main", "no",  "interaction"))
  
  if(use_glmm) {
    formula <- switch(use_n_bees_initial,
                      "no" = n_bees ~ Assessment + (1|Site/Colony),
                      "main" = n_bees ~ log(n_bees_initial) + Assessment + (1|Site/Colony),
                      "interaction" = n_bees ~ log(n_bees_initial) * Assessment + (1|Site/Colony))
    model <- glmmTMB(formula, family = nbinom2(), data = original_data)
  } else {
    formula <- switch(use_n_bees_initial,
                      "no" = log(n_bees) ~ Assessment + (1|Site/Colony),
                      "main" = log(n_bees) ~ log(n_bees_initial) + Assessment + (1|Site/Colony),
                      "interaction" = log(n_bees) ~ log(n_bees_initial) * Assessment + (1|Site/Colony))
    model <- lmer(formula, data = original_data)
  }
  
  # Extract random effect SDs
  VarCorrs <- VarCorr(model)
  if (use_glmm) {
    colony_sd <- sqrt(as.data.frame(VarCorrs[[1]]$`Colony:Site`)[[1]])
    site_sd <- sqrt(as.data.frame(VarCorrs[[1]]$Site)[[1]])
    residuals_log <- log(original_data$n_bees) - predict(model, type = "link", re.form = NULL)
    residual_sd <- sd(residuals_log)
  } else {
    VarCorrs_df <- as.data.frame(VarCorrs)
    colony_sd <- VarCorrs_df[VarCorrs_df$grp == "Colony:Site", ]$sdcor
    site_sd <- VarCorrs_df[VarCorrs_df$grp == "Site", ]$sdcor
    residual_sd <- VarCorrs_df[VarCorrs_df$grp == "Residual", ]$sdcor
  }
  
  first_assessment <- sort(unique(original_data$Assessment))[1]
  n_bees_initial <- original_data$n_bees_initial[original_data$Assessment == first_assessment]
  initial_mean <- mean(n_bees_initial, na.rm = TRUE)
  initial_sd <- sd(n_bees_initial, na.rm = TRUE)
  
  output <- list(
    model = model,
    use_n_bees_initial = use_n_bees_initial,
    use_glmm = use_glmm, 
    colony_sd = colony_sd,
    site_sd = site_sd,
    residual_sd = residual_sd,
    initial_mean = initial_mean,
    initial_sd = initial_sd)
  
  return(output)
}

calculate_sd_values <- function(n_sites_total, n_colonies) {
  # Step 1: Simulate control data
  sim_data <- mimic_control_data(
    var_comps_list = var_comps_glmm_main,
    original_data = original_control_data_234,
    n_sites_total = n_sites_total,
    n_colonies = n_colonies,
    backtransform = TRUE,
    reallocate_colonies = FALSE # First simulate unallocated structure
  )
  
  # Step 2: Add covariate (e.g., initial colony size)
  sim_data$n_bees_initial <- round(rnorm(n = nrow(sim_data), mean = 1000, sd = 100))
  sim_data$n_bees_initial[sim_data$n_bees_initial < 0] <- 0
  
  # Step 3: Allocate treatments to colonies based on covariates (colony-level)
  sim_data <- allocate_treatments(
    data = sim_data,
    treatments = c("Control", "Pesticide"),
    treatment_var = "Treatment",
    covariates = "n_bees_initial"
  ) %>%
    dplyr::select(-Site)  # temporarily remove Site to reassign after
  
  # Step 4: Assign treatments to sites (balanced reallocation)
  site_ids <- sample(unique(sim_data$Colony))  # assuming one colony per site initially
  control_sites <- assign_control_sites(site_ids)
  site_treatment_combis <- data.frame(
    Site = site_ids,
    Treatment = ifelse(site_ids %in% control_sites, "Control", "Pesticide")
  )
  
  # Step 5: Re-allocate colonies to sites (covariate-balanced)
  sim_data <- allocate_within_treatments(
    data_subjects = sim_data,
    data_groups = site_treatment_combis,
    treatment_var = "Treatment",
    group_var = "Site",
    covariates = "n_bees_initial"
  )
  
  # Step 6: Fit model and extract SDs
  model_balanced <- glmmTMB(
    n_bees ~ log(n_bees_initial) + Assessment + (1 | Site/Colony),
    data = sim_data,
    family = nbinom2()
  )
  
  output <- data.frame(
    site_sd = sqrt(VarCorr(model_balanced)$cond$Site[1]),
    colony_sd = sqrt(VarCorr(model_balanced)$cond[["Colony:Site"]][1])
  )
  
  return(output)
}

assign_control_sites <- function(site_ids, n_sites_control = NULL) {
  
  n_sites_total <- length(site_ids)
  n_sites_control <- ifelse(is.null(n_sites_control), floor(n_sites_total / 2), n_sites_control)
  control_sites <- sample(site_ids, size = n_sites_control)
  control_sites
  
}

assign_treatments <- function(data, n_sites_control = NULL) {
  site_ids <- unique(data$Site)
  # n_sites_total <- length(site_ids)
  # 
  # n_sites_control <- ifelse(is.null(n_sites_control), floor(n_sites_total / 2), n_sites_control)
  # control_sites <- sample(site_ids, size = n_sites_control)
  
  control_sites <- assign_control_sites(site_ids = site_ids, n_sites_control = n_sites_control)
  
  output <- data %>%
    mutate(Treatment = factor(ifelse(Site %in% control_sites, "Control", "Pesticide")))
  
  output
}

apply_effects <- function(data, n_sites_control = NULL, 
                          effect_size, effect_timing = "abrupt") {
  
  # Check for necessary columns
  required_columns <- c("Treatment", "Assessment", "n_bees", "Site", "Period")
  missing_columns <- setdiff(required_columns, names(data))
  
  if (length(missing_columns) > 0) {
    stop(paste("The following required columns are missing in the data:", 
               paste(missing_columns, collapse = ", ")))
  }
  
  # Ensure there's at least one row where Period == "After"
  if (sum(data$Period == "After") == 0) {
    stop("The dataset must contain at least one 'After' period in the 'Period' column.")
  }
  
  output <- data %>%
    mutate(Assessment = factor(Assessment),
           n_bees_before_effect = n_bees)
  
  
  # Apply treatment effect (abrupt or otherwise)
  if (effect_timing == "abrupt") {
    output <- output %>%
      mutate(n_bees = round(ifelse(Treatment == "Pesticide" & Period == "After", 
                                   (1 - effect_size) * n_bees, n_bees)))
  }
  
  # Additional effect_timing types can be handled here (e.g., linear effects)
  
  return(output)
}

mimic_control_data <- function(var_comps_list = NULL,
                               original_data, 
                               n_sites_total, 
                               n_colonies, 
                               use_n_bees_initial = c("main", "no", "interaction"), 
                               use_glmm = TRUE,
                               backtransform = TRUE,
                               reallocate_colonies = FALSE,
                               covariates = NULL, 
                               repetitions = 10,
                               standardize = TRUE) {
  
  use_n_bees_initial <- match.arg(use_n_bees_initial, c("main", "no", "interaction"))
  
  if (!reallocate_colonies & !is.null(covariates)) {
    warning("Covariates were specified but reallocate_colonies = FALSE. Covariates are disregarded and no reallocation is performed.")
  }
  
  # Either extract from list or compute from scratch
  if (is.null(var_comps_list)) {
    # Fit the model
    if (use_glmm) {
      formula <- switch(use_n_bees_initial,
                        "no" = n_bees ~ Assessment + (1|Site/Colony),
                        "main" = n_bees ~ log(n_bees_initial) + Assessment + (1|Site/Colony),
                        "interaction" = n_bees ~ log(n_bees_initial) * Assessment + (1|Site/Colony))
      model <- glmmTMB(formula, family = nbinom2(), data = original_data)
      VarCorrs <- VarCorr(model)
      # colony_sd <- sqrt(as.data.frame(VarCorrs[[1]]$Colony:Site)[[1]])
      colony_sd <- sqrt(as.data.frame(VarCorrs[[1]]$`Colony:Site`)[[1]])
      site_sd <- sqrt(as.data.frame(VarCorrs[[1]]$Site)[[1]])
      residuals_log <- log(original_data$n_bees) - predict(model, type = "link", re.form = NULL)
      residual_sd <- sd(residuals_log)
    } else {
      formula <- switch(use_n_bees_initial,
                        "no" = log(n_bees) ~ Assessment + (1|Site/Colony),
                        "main" = log(n_bees) ~ log(n_bees_initial) + Assessment + (1|Site/Colony),
                        "interaction" = log(n_bees) ~ log(n_bees_initial) * Assessment + (1|Site/Colony))
      model <- lmer(formula, data = original_data)
      VarCorrs_df <- as.data.frame(VarCorr(model))
      colony_sd <- VarCorrs_df[VarCorrs_df$grp == "Colony:Site", ]$sdcor
      site_sd <- VarCorrs_df[VarCorrs_df$grp == "Site", ]$sdcor
      residual_sd <- VarCorrs_df[VarCorrs_df$grp == "Residual", ]$sdcor
    }
    
    first_assessment <- sort(unique(original_data$Assessment))[1]
    n_bees_initial <- original_data$n_bees_initial[original_data$Assessment == first_assessment]
    initial_mean <- mean(n_bees_initial, na.rm = TRUE)
    initial_sd <- sd(n_bees_initial, na.rm = TRUE)
    
  } else {
    model <- var_comps_list$model
    colony_sd <- var_comps_list$colony_sd
    site_sd <- var_comps_list$site_sd
    residual_sd <- var_comps_list$residual_sd
    initial_mean <- var_comps_list$initial_mean
    initial_sd <- var_comps_list$initial_sd
  }
  
  # Create site and colony IDs
  new_sites <- paste0("F", 1:n_sites_total)
  new_colonies <- paste0(1:n_colonies)
  assessments <- unique(original_data$Assessment)
  
  IDs <- expand.grid(Site = new_sites, Colony = new_colonies) %>%
    mutate(Colony = paste(Site, Colony, sep = "_"))
  
  # Simulate initial bee counts
  if (use_n_bees_initial != "no") {
    IDs$n_bees_initial <- round(rnorm(n = nrow(IDs), mean = initial_mean, sd = initial_sd))
    IDs$n_bees_initial[IDs$n_bees_initial < 0] <- 0
  }
  
  # Reallocation
  if (reallocate_colonies) {
    if (!is.null(covariates)) {
      
      site_ids <- sample(unique(IDs$Site))
      
      # Assign treatments in balanced way
      IDs <- allocate_treatments(IDs,
                                 treatments = c("Control", "Pesticide"),
                                 treatment_var = "Treatment",
                                 covariates = covariates,
                                 repetitions = repetitions,
                                 standardize = standardize) %>% 
        dplyr::select(-Site)
      
      control_sites <- assign_control_sites(site_ids = site_ids)
      site_treatment_combis <- data.frame(
        Site = site_ids,
        Treatment = ifelse(site_ids %in% control_sites, "Control", "Pesticide")
      )
      
      IDs <- allocate_within_treatments(data_subjects = IDs,
                                        data_groups = site_treatment_combis,
                                        treatment_var = "Treatment",
                                        group_var = "Site",
                                        covariates = covariates,
                                        repetitions = repetitions, 
                                        standardize = standardize)
      
    } else {
      IDs$Site <- sample(IDs$Site, replace = FALSE)
    }
  }
  
  # Assign treatments randomly
  if(is.null(covariates)){
    IDs <- assign_treatments(IDs)
  }
  
  # Prediction grid
  new_data <- expand.grid(Assessment = assessments, Colony = IDs$Colony) %>%
    left_join(IDs, by = "Colony")
  
  # Predict fixed effects
  predicted_fixed <- predict(model, newdata = new_data, re.form = NA, allow.new.levels = TRUE)
  
  # Add random effects
  new_data <- new_data %>%
    mutate(random_site = rnorm(n_sites_total)[as.numeric(factor(Site))] * site_sd,
           random_colony = rnorm(n_sites_total * n_colonies)[as.numeric(factor(Colony))] * colony_sd,
           residual = rnorm(nrow(.), mean = 0, sd = residual_sd))
  
  new_data$log_n_bees <- predicted_fixed + new_data$random_site + new_data$random_colony + new_data$residual
  new_data$log_n_bees_no_site_effect <- predicted_fixed + new_data$random_colony + new_data$residual
  new_data$log_n_bees_fixed_effect_only <- predicted_fixed
  
  if (backtransform) {
    new_data$n_bees <- round(exp(new_data$log_n_bees)) %>% pmax(0)
    new_data$n_bees_no_site_effect <- round(exp(new_data$log_n_bees_no_site_effect)) %>% pmax(0)
    new_data$n_bees_fixed_effect_only <- round(exp(new_data$log_n_bees_fixed_effect_only)) %>% pmax(0)
  }
  
  metadata_vars <- c("Phase", "Period", "Study")
  available_metadata <- intersect(metadata_vars, names(original_data))
  if (length(available_metadata) > 0) {
    metadata_df <- original_data %>%
      dplyr::select(Assessment, all_of(available_metadata)) %>%
      distinct()
    new_data <- left_join(new_data, metadata_df, by = "Assessment")
  }
  
  # new_data <- as_tibble(new_data) %>% select(any_of(c("Treatment", "Site", "Colony", "Assessment", "n_bees", "n_bees_initial")), everything())
  
  return(new_data)
}

# mimic_control_data <- function(original_data,
#                                n_sites_total,
#                                n_colonies,
#                                use_n_bees_initial = c("main", "no", "interaction"),
#                                use_glmm = TRUE,
#                                backtransform = TRUE,
#                                reallocate_colonies = FALSE,
#                                covariates = NULL) {
# 
#   # Check argument values
#   use_n_bees_initial <- match.arg(use_n_bees_initial, c("main", "no",  "interaction"))
# 
#   if (!reallocate_colonies & !is.null(covariates)) {
#     warning("Covariates were specified but reallocate_colonies = FALSE. Covariates are disregarded and no reallocation is performed.")
#   }
# 
#   # Fit the model (either lmer or glmmTMB)
#   if (use_glmm) {
#     formula <- switch(use_n_bees_initial,
#                       "no" = n_bees ~ Assessment + (1|Site/Colony),
#                       "main" = n_bees ~ log(n_bees_initial) + Assessment + (1|Site/Colony),
#                       "interaction" = n_bees ~ log(n_bees_initial) * Assessment + (1|Site/Colony))
#     model <- glmmTMB(formula, family = nbinom2(), data = original_data)
#   } else {
#     formula <- switch(use_n_bees_initial,
#                       "no" = log(n_bees) ~ Assessment + (1|Site/Colony),
#                       "main" = log(n_bees) ~ log(n_bees_initial) + Assessment + (1|Site/Colony),
#                       "interaction" = log(n_bees) ~ log(n_bees_initial) * Assessment + (1|Site/Colony))
#     model <- lmer(formula, data = original_data)
#   }
# 
#   # Create new site and colony IDs
#   new_sites <- paste0("F", 1:n_sites_total)
#   new_colonies <- paste0(1:n_colonies)
#   assessments <- unique(original_data$Assessment)
# 
#   IDs <- expand.grid(Site = new_sites, Colony = new_colonies) %>%
#     mutate(Colony = paste(Site, Colony, sep = "_"))
# 
#   # Simulate initial bee counts if using n_bees_initial
#   if (use_n_bees_initial != "no") {
#     first_assessment <- sort(unique(original_data$Assessment))[1]
#     n_bees_initial <- subset(original_data, Assessment == first_assessment)$n_bees_initial
#     initial_mean <- mean(n_bees_initial, na.rm = TRUE)
#     initial_sd <- sd(n_bees_initial, na.rm = TRUE)
#     IDs$n_bees_initial <- round(rnorm(n = nrow(IDs), mean = initial_mean, sd = initial_sd))
#     IDs$n_bees_initial[IDs$n_bees_initial < 0] <- 0
#   }
# 
#   # Apply reallocation if requested
#   if (reallocate_colonies) {
#     if (!is.null(covariates)) {
# 
#       # Balanced reallocation
#       site_ids <- unique(IDs$Site)
#       IDs <- allocate_treatments(IDs,
#                                  treatments = site_ids,
#                                  treatment_var = "Site",
#                                  covariates = covariates)
# 
#     } else {
#       # Random reallocation
#       IDs$Site <- sample(IDs$Site, replace = FALSE)
#     }
#   }
# 
#   # Prepare prediction data
#   new_data <- expand.grid(Assessment = assessments, Colony = IDs$Colony) %>%
#     left_join(IDs, by = "Colony")
# 
#   # Predict fixed effects
#   predicted_fixed <- predict(model, newdata = new_data, re.form = NA, allow.new.levels = TRUE)
# 
#   # Extract random effect SDs
#   VarCorrs <- VarCorr(model)
# 
#   if (use_glmm) {
#     colony_sd   <- sqrt(as.data.frame(VarCorrs[[1]]$Colony:Site))[[1]]
#     site_sd <- sqrt(as.data.frame(VarCorrs[[1]]$Site))[[1]]
# 
#     predictions_including_RE <- predict(model, newdata = original_data, type = "link", re.form = NULL)
#     residuals_log_scale <- log(original_data$n_bees) - predictions_including_RE
#     residual_sd <- sd(residuals_log_scale)
# 
#   } else {
#     VarCorrs_df <- as.data.frame(VarCorrs)
#     colony_sd <- VarCorrs_df[VarCorrs_df$grp == "Colony:Site", ]$sdcor
#     site_sd   <- VarCorrs_df[VarCorrs_df$grp == "Site", ]$sdcor
#     residual_sd <- VarCorrs_df[VarCorrs_df$grp == "Residual", ]$sdcor
#   }
# 
#   # Assign random effects based on final site allocation
#   new_data <- new_data %>%
#     mutate(random_site   = rnorm(n_sites_total)[as.numeric(factor(Site))] * site_sd,
#            random_colony = rnorm(n_sites_total * n_colonies)[as.numeric(factor(Colony))] * colony_sd,
#            residual      = rnorm(nrow(.), mean = 0, sd = residual_sd))
# 
#   # Combine components
#   new_data$log_n_bees <- predicted_fixed + new_data$random_site + new_data$random_colony + new_data$residual
#   new_data$log_n_bees_no_site_effect <- predicted_fixed + new_data$random_colony + new_data$residual
#   new_data$log_n_bees_fixed_effect_only <- predicted_fixed
# 
#   if (backtransform) {
#     new_data$n_bees <- round(exp(new_data$log_n_bees))
#     new_data$n_bees[new_data$n_bees < 0] <- 0
# 
#     new_data$n_bees_no_site_effect <- round(exp(new_data$log_n_bees_no_site_effect))
#     new_data$n_bees_no_site_effect[new_data$n_bees_no_site_effect < 0] <- 0
# 
#     new_data$n_bees_fixed_effect_only <- round(exp(new_data$log_n_bees_fixed_effect_only))
#     new_data$n_bees_fixed_effect_only[new_data$n_bees_fixed_effect_only < 0] <- 0
#   }
# 
#   # Add metadata
#   metadata_vars <- c("Phase", "Period", "Study")
#   available_metadata <- intersect(metadata_vars, names(original_data))
#   if (length(available_metadata) > 0) {
#     metadata_df <- original_data %>%
#       dplyr::select(Assessment, all_of(available_metadata)) %>%
#       distinct()
# 
#     new_data <- left_join(new_data, metadata_df, by = "Assessment")
#   }
# 
#   return(new_data)
# }

# reallocate_colonies_randomly <- function(data) {
# 
#   first_assessment <- sort(unique(data$Assessment))[1]
#   data_first_assessment <- subset(data, Assessment == first_assessment)
#   data_first_assessment$Site <- sample(data_first_assessment$Site, replace = FALSE)
# 
#   assigned_sites <- data_first_assessment[, c("Site", "Colony")]
# 
#   data_other_assessments <- subset(data, Assessment != first_assessment) %>%
#     mutate(Site = NULL) %>% full_join(assigned_sites)
# 
#   output <- bind_rows(data_first_assessment, data_other_assessments)
# 
#   return(output)
# }

permute_colonies <- function(data) {
  
  first_assessment <- sort(unique(data$Assessment))[1]
  data_first_assessment <- subset(data, Assessment == first_assessment)
  treatment_site_combis <- unique(data_first_assessment[, c("Treatment", "Site")])
  data_first_assessment$Treatment <- NULL
  
  data_first_assessment$Site <- sample(data_first_assessment$Site, replace = FALSE)
  data_first_assessment <- full_join(treatment_site_combis, data_first_assessment)
  
  assigned_sites <- data_first_assessment[, c("Treatment", "Site", "Colony")]
  
  data_other_assessments <- subset(data, Assessment != first_assessment) %>% 
    dplyr::select(-Site, -Treatment) %>% full_join(assigned_sites)
  
  output <- bind_rows(data_first_assessment, data_other_assessments)
  
  return(output)
}

# reallocate_colonies_balanced <- function(data, covariates){
#   
#   first_assessment <- sort(unique(data$Assessment))[1]
#   data_first_assessment <- subset(data, Assessment == first_assessment)
#   site_ids <- unique(data_first_assessment$Site)
#   
#   data_first_assessment <- allocate_treatments(data_first_assessment, 
#                                                treatments = site_ids, 
#                                                treatment_var = "Site",
#                                                covariates = covariates) 
#   
#   assigned_sites <- data_first_assessment[, c("Site", "Colony")]
#   
#   data_other_assessments <- subset(data, Assessment != first_assessment) %>%
#     mutate(Site = NULL) %>% full_join(assigned_sites)
#   
#   output <- bind_rows(data_first_assessment, data_other_assessments)
#   
#   return(output)
# }

apply_treatment_effects <- function(data, n_sites_control = NULL,
                                    effect_size, effect_timing = "abrupt") {

  # Check for necessary columns
  required_columns <- c("Assessment", "n_bees", "Site", "Period")
  missing_columns <- setdiff(required_columns, names(data))

  if (length(missing_columns) > 0) {
    stop(paste("The following required columns are missing in the data:",
               paste(missing_columns, collapse = ", ")))
  }

  # Ensure there's at least one row where Period == "After"
  if (sum(data$Period == "After") == 0) {
    stop("The dataset must contain at least one 'After' period in the 'Period' column.")
  }

  # Assign treatments (control vs pesticide)
  unique_sites <- unique(data$Site)
  n_sites_total <- length(unique_sites)

  n_sites_control <- ifelse(is.null(n_sites_control), floor(n_sites_total / 2), n_sites_control)
  control_sites <- sample(unique_sites, size = n_sites_control)

  data_with_treatment <- data %>%
    mutate(Assessment = factor(Assessment),
           Treatment = factor(ifelse(Site %in% control_sites, "Control", "Pesticide")))

  data_with_treatment$n_bees_before_effect <- data_with_treatment$n_bees

  # Apply treatment effect (abrupt or otherwise)
  if (effect_timing == "abrupt") {
    data_with_treatment <- data_with_treatment %>%
      mutate(n_bees = round(ifelse(Treatment == "Pesticide" & Period == "After",
                                   (1 - effect_size) * n_bees, n_bees)))
  }

  # Additional effect_timing types can be handled here (e.g., linear effects)

  return(data_with_treatment)
}

# Run models -------
run_model <- function(data, formula, family = "nbinom2", ...) {
  
  # If the family is Gaussian or normal, use lmer
  if (family == "Gaussian" || family == "normal") {
    model <- lmer(formula, data = data, ...)
  } 
  # Otherwise, use glmmTMB for other families (e.g., nbinom2, poisson, etc.)
  else {
    model <- glmmTMB(formula, data = data, family = family, ...)
  }
  
  return(model)
}

# Run tests ----------
test_difference <- function(model, alpha = 0.05, 
                            delta = 0.07,
                            one_sided = TRUE, 
                            per_assessment = TRUE) {
  if (is.null(model)) return(NULL)
  
  if (one_sided == TRUE) {
    side <- "<"
  } else {
    side <- FALSE
  }
  
  by_arg <- if (per_assessment) "Assessment" else NULL
  
  
  EMM <- suppressMessages(emmeans(model, revpairwise ~ Treatment, 
                                  infer = c(TRUE, FALSE), 
                                  level = 1 - alpha, type = "response", 
                                  by = by_arg))
  
  
  stats_contrasts <- as.data.table(test(EMM$contrasts, alpha = alpha, side = side))
  setnames(stats_contrasts, old = tail(names(stats_contrasts), 1), 
           new = "p_value_diff")
  
  stats_contrasts$effect_size_estimate <- stats_contrasts$ratio - 1
  
  stats_contrasts <- stats_contrasts %>%
    mutate(risk_diff = ifelse(p_value_diff < alpha & 
                                effect_size_estimate < -delta, 1, 0))  
  
  if (!per_assessment) {
    stats_contrasts <- stats_contrasts %>%
      mutate(Assessment = "Across assessments") %>%
      dplyr::select(contrast, Assessment, everything())
  }
  
  stats_contrasts$Assessment <- as.factor(as.character(stats_contrasts$Assessment))
  
  return(stats_contrasts)
}

test_equivalence <- function(model, alpha = 0.2, delta = 0.1, 
                             one_sided = TRUE, per_assessment) {
  if (is.null(model)) return(NULL)
  level <- ifelse(one_sided, 1 - 2 * alpha, 1 - alpha)
  
  by_arg <- if (per_assessment) "Assessment" else NULL
  
  if (one_sided == TRUE) {
    side <- ">"
  } else {
    side <- FALSE
  }
  
  EMM <- suppressMessages(emmeans(model, revpairwise ~ Treatment, 
                                  infer = c(TRUE, FALSE), 
                                  level = level, type = "response", 
                                  by = by_arg))
  
  CI_contrasts <- as.data.table(EMM$contrasts)
  setnames(CI_contrasts, old = tail(names(CI_contrasts), 2), 
           new = c("CI_lower_equi", "CI_upper_equi"))
  
  # Check if the emmeans are on the log scale
  is_log_scale <- attr(EMM$emmeans, "misc")$tran == "log"
  
  # Set null based on the scale
  null_value <- if (is_log_scale) log(1 - delta) else (1 - delta)
  
  # Test contrasts with the appropriate null
  p_value_tab <- as.data.table(
    test(EMM$contrasts, null = null_value, alpha = alpha, side = side)
  )
  
  setnames(p_value_tab, old = tail(names(p_value_tab), 1), 
           new = "p_value_equi")
  
  CI_contrasts <- merge(CI_contrasts, p_value_tab)
  
  CI_contrasts <- CI_contrasts %>%
    mutate(threshold = 1 - delta,
           risk_equi = ifelse(p_value_equi >= alpha, 1, 0)) %>% 
    dplyr::select(-null, -ends_with(".ratio"))
  
  if (!per_assessment) {
    CI_contrasts <- CI_contrasts %>%
      mutate(Assessment = "Across assessments")
  }
  
  CI_contrasts$Assessment <- as.factor(as.character(CI_contrasts$Assessment))
  
  return(CI_contrasts)
}

test_Hotopp_equivalence <- function(model, alpha_Hotopp_control = 0.10,
                                    alpha = 0.2, delta = 0.1,
                                    one_sided_equi = TRUE,
                                    per_assessment) {
  if (is.null(model)) return(NULL)
  
  level_equi <- ifelse(one_sided_equi, 1 - 2 * alpha, 1 - alpha)
  by_arg <- if (per_assessment) "Assessment" else NULL
  
  EMM <- suppressMessages(emmeans(model, revpairwise ~ Treatment,
                                  infer = c(TRUE, FALSE),
                                  level = 1 - alpha_Hotopp_control,
                                  by = by_arg))
  
  EMM_tab <- as.data.table(EMM$emmeans)
  setnames(EMM_tab, old = tail(names(EMM_tab), 2),
           new = c("CI_lower", "CI_upper"))
  
  critical_score <- qnorm((1 + level_equi)/2)
  
  if (!per_assessment) {
    EMM_tab <- EMM_tab %>%
      mutate(Assessment = "Across assessments")
  }
  
  EMM_tab$Assessment <- as.factor(as.character(EMM_tab$Assessment))
  
  # Corrected pivot_wider call
  EMM_tab <- EMM_tab %>%
    pivot_wider(
      id_cols = "Assessment",
      names_from = Treatment,
      values_from = c(emmean, SE, CI_lower, CI_upper)
    )
  
  EMM_tab <- EMM_tab %>%
    mutate(
      log_ratio_Hotopp = emmean_Pesticide - CI_lower_Control,
      SE_log_ratio = sqrt(SE_Control^2 + SE_Pesticide^2),
      CI_log_ratio_lower_Hotopp = log_ratio_Hotopp - critical_score * SE_log_ratio,
      CI_log_ratio_upper_Hotopp = log_ratio_Hotopp + critical_score * SE_log_ratio,
      ratio_Hotopp = exp(log_ratio_Hotopp),
      CI_lower_Hotopp = exp(CI_log_ratio_lower_Hotopp),
      CI_upper_Hotopp = exp(CI_log_ratio_upper_Hotopp),
      threshold = 1 - delta,
      risk_equi_Hotopp = ifelse(CI_lower_Hotopp < threshold, 1, 0)
    ) %>%
    dplyr::select(Assessment, ratio_Hotopp, CI_lower_Hotopp,
                  CI_upper_Hotopp, risk_equi_Hotopp)
  
  return(EMM_tab)
  
}

run_tests <- function(model, alpha_d, alpha_e, alpha_Hotopp_control, 
                      delta_e, tests, per_assessment) {
  test_results <- list()
  
  # Adjust the suffix based on the per_assessment argument
  suffix <- if (per_assessment) "per_assessment" else "across_assessments"
  
  if ("diff" %in% tests) {
    test_results[[paste0("diff_", suffix)]] <- 
      test_difference(model, 
                      alpha = alpha_d, 
                      per_assessment = per_assessment)
  }
  
  if ("equi" %in% tests) {
    test_results[[paste0("equi_", suffix)]] <- 
      test_equivalence(model, 
                       alpha = alpha_e, 
                       delta = delta_e, 
                       per_assessment = per_assessment)
  }
  
  if ("Hotopp" %in% tests) {
    test_results[[paste0("Hotopp_", suffix)]] <- 
      test_Hotopp_equivalence(model, 
                              alpha_Hotopp_control = alpha_Hotopp_control, 
                              alpha = alpha_e, delta = delta_e, 
                              per_assessment = per_assessment)
  }
  
  return(test_results)
}

extrapolate_days_for_all_simulations <- function(minutes_for_test, n_simulations_test, n_simulations_total = 10000,
                                                 n_site_values = 1, n_colony_values = 1){
  
  days <- minutes_for_test * n_simulations_total /n_simulations_test /60 /24 * 
    n_site_values * n_colony_values
  
  days
  
}

adjust_control_sd <- function(data, control_sd_factor) {
  
  stopifnot(is.numeric(data$n_bees))  
  
  control_sd_fct <- control_sd_factor
  
  # Create new bee counts with adapted sd rounded to the next integer (Control only)
  new_data_rounded_preclamp <- data %>% 
    dplyr::group_by(Assessment) %>%
    dplyr::mutate(
      `..mu..`   = ifelse(Treatment == "Control",
                          mean(n_bees[Treatment == "Control"], na.rm = TRUE),
                          NA_real_),
      `..resid..` = ifelse(Treatment == "Control", n_bees - `..mu..`, NA_real_),
      n_bees      = round(ifelse(Treatment == "Control",
                                 `..mu..` + .env$control_sd_fct * `..resid..`,
                                 n_bees))
    ) %>%
    dplyr::ungroup()
  
  # Diagnostics on the rounded, pre-clamp values (Control only)
  n_negatives     <- sum(new_data_rounded_preclamp$Treatment == "Control" &
                           new_data_rounded_preclamp$n_bees < 0, na.rm = TRUE)
  lowest_n_bees <- suppressWarnings(
    min(new_data_rounded_preclamp$n_bees[new_data_rounded_preclamp$Treatment == "Control"], na.rm = TRUE)
  )
  
  # Now clamp to >= 0 and drop temps
  new_data <- new_data_rounded_preclamp %>%
    dplyr::mutate(n_bees = pmax(n_bees, 0)) %>%
    dplyr::select(-`..mu..`, -`..resid..`)
  
  list(
    data = new_data,
    n_negatives = n_negatives,
    lowest_n_bees = lowest_n_bees
  )
}

# Run simulations -----
run_simulations <- function(n_simulations, 
                            data = data, 
                            mimic_data = TRUE,
                            n_sites, n_colonies, 
                            var_comps_list = NULL, 
                            mimic_args = list(use_n_bees_initial = "main", 
                                              use_glmm = TRUE,
                                              backtransform = TRUE),
                            reallocate_colonies = FALSE,
                            permute = FALSE, 
                            covariates = NULL,
                            control_sd_factor = 1,
                            effect_size, 
                            effect_timing = "abrupt",
                            formula, 
                            family = "Gaussian",
                            delta_e = 0.1,
                            alpha_d = 0.05, alpha_e = 0.2, 
                            alpha_Hotopp_control = 0.1,
                            tests, assessment_scope, 
                            parallel = TRUE, 
                            cores = NULL,
                            seed_number = NULL) {
  
  # Helper function to run a single simulation
  simulate_once <- function(data,
                            mimic_data,
                            n_sites, 
                            n_colonies, 
                            var_comps_list,
                            mimic_args,
                            reallocate_colonies,
                            permute,
                            covariates,
                            control_sd_factor,
                            effect_size, 
                            effect_timing,
                            formula, 
                            family,
                            delta_e, 
                            alpha_d, 
                            alpha_e, 
                            alpha_Hotopp_control, 
                            tests, 
                            assessment_scope,
                            run_number) {
    
    # Generate the data with treatments
    if(mimic_data) {
      
      if (!is.null(var_comps_list)) {
        mimic_args <- modifyList(mimic_args, list(
          var_comps_list = var_comps_list,
          original_data = data,
          n_sites_total = 2 * n_sites,
          n_colonies = n_colonies,
          reallocate_colonies = reallocate_colonies,
          covariates = covariates
        ))
      } else {
        mimic_args <- modifyList(mimic_args, list(
          original_data = data,
          n_sites_total = 2 * n_sites,
          n_colonies = n_colonies,
          reallocate_colonies = reallocate_colonies,
          covariates = covariates
        ))
      }
      
      data_without_effects <- do.call(mimic_control_data, mimic_args) 
      
    } else {
      data_without_effects <- assign_treatments(data)  
    }
    
    if(permute){
      data_without_effects <- 
        permute_colonies(data_without_effects)
    }
    
    # # Reallocate colonies if specified
    # if (reallocate_colonies) {
    #   
    #   # Balanced reallocation if reallocate_colonies == T and a covariate(s) is/are specified
    #   if (!is.null(covariates)) {
    #     data_without_effects <- 
    #       reallocate_colonies_balanced(data_without_effects, covariates = covariates)
    #   } else {
    #     # Reallocate randomly if reallocate_colonies == T but no covariates were specified
    #     data_without_effects <- 
    #       reallocate_colonies_randomly(data_without_effects)
    #   }
    # } else {
    #   # No reallocation if reallocate_colonies == F
    #   
    #   if (!is.null(covariates)) {
    #     warning("Covariates were specified but reallocate_colonies = FALSE. Covariates are disregarded and no reallocation is performed.")
    #   }
    # }
    
    # before_effect_stats <- data_without_effects %>%
    #   group_by(Treatment, Assessment) %>%
    #   summarise(
    #     mean_n_bees_before_effect = mean(n_bees, na.rm = TRUE),
    #     sd_n_bees_before_effect = sd(n_bees, na.rm = TRUE),
    #     .groups = "drop"
    #   ) %>%
    #   pivot_wider(names_from = Treatment,
    #               values_from = c(mean_n_bees_before_effect, sd_n_bees_before_effect),
    #               names_sep = "_") %>% 
    #   mutate(
    #     diff_mean_n_bees_before_effect = mean_n_bees_before_effect_Pesticide - mean_n_bees_before_effect_Control,
    #     diff_sd_n_bees_before_effect = sd_n_bees_before_effect_Pesticide - sd_n_bees_before_effect_Control
    #   )
    
    n_negatives <- NA_integer_
    lowest_n_bees <- NA_real_ 
    
    if(control_sd_factor != 1){
      
      adj_out <- adjust_control_sd(data_without_effects, 
                                                  control_sd_factor = control_sd_factor)
      
      data_without_effects <- adj_out$data
      n_negatives <- adj_out$n_negatives
      lowest_n_bees <- adj_out$lowest_n_bees
    }
    
    before_effect_stats <- data_without_effects %>%
      group_by(Assessment) %>%
      summarise(
        mean_diff_before = mean(n_bees[Treatment == "Pesticide"], na.rm = TRUE) - 
          mean(n_bees[Treatment == "Control"], na.rm = TRUE),
        sd_diff_before = sd(n_bees[Treatment == "Pesticide"], na.rm = TRUE) - 
          sd(n_bees[Treatment == "Control"], na.rm = TRUE),
        .groups = "drop"
      )
    
    data_with_effects <- apply_effects(data = data_without_effects,
                                       effect_size = effect_size,
                                       effect_timing = effect_timing)
    
    # Run the model
    model <- run_model(data_with_effects, formula = formula, family = family)
    
    # Run tests based on assessment scope
    test_results <- list()
    
    if (assessment_scope == "per" || assessment_scope == "both") {
      test_results_per_assessment <- run_tests(
        model, 
        alpha_d = alpha_d, 
        alpha_e = alpha_e, 
        alpha_Hotopp_control = alpha_Hotopp_control, 
        delta_e = delta_e, 
        tests = tests, 
        per_assessment = TRUE
      )
      test_results_per_assessment <- suppressMessages(reduce(test_results_per_assessment, 
                                                             full_join))
      test_results <- append(test_results, list(test_results_per_assessment))
    }
    
    if (assessment_scope == "across" || assessment_scope == "both") {
      test_results_across_assessments <- run_tests(
        model, 
        alpha_d = alpha_d, 
        alpha_e = alpha_e, 
        alpha_Hotopp_control = alpha_Hotopp_control, 
        delta_e = delta_e, 
        tests = tests, 
        per_assessment = FALSE
      )
      test_results_across_assessments <- suppressMessages(reduce(test_results_across_assessments, 
                                                                 full_join))
      test_results <- append(test_results, list(test_results_across_assessments))
    }
    
    # Combine results
    combined_results <- suppressMessages(reduce(test_results, full_join))
    combined_results <- combined_results %>%
      mutate(effect_size = effect_size, effect_timing = effect_timing,
             n_sites = n_sites, n_colonies = n_colonies,
             delta_e = delta_e, 
             alpha_e = alpha_e,
             control_sd_factor = control_sd_factor,
             n_negatives = n_negatives,
             lowest_n_bees = lowest_n_bees,
             run_number = run_number)
    
    combined_results$Dataset <- if(mimic_data) "Simulated data" else "Original data"
    
    combined_results$Reallocation <- case_when(
      !reallocate_colonies ~ "none",
      reallocate_colonies & is.null(covariates) ~ "random",
      reallocate_colonies & !is.null(covariates) ~ "balanced"
    )
    
    # combined_results <- bind_cols(combined_results, before_effect_stats[rep(1, nrow(combined_results)), ])
    
    combined_results <- combined_results %>%
      left_join(before_effect_stats, by = "Assessment")
    
    
    return(combined_results)
  }
  
  start_time <- Sys.time()
  print(paste("Start time:", start_time))
  
  # Generate all combinations of several arguments to loop through
  param_combinations <- expand.grid(n_sites = n_sites, 
                                    n_colonies = n_colonies, 
                                    effect_size = effect_size,
                                    delta_e = delta_e,
                                    alpha_e = alpha_e,
                                    control_sd_factor = control_sd_factor)
  
  simulation_results <- bind_rows(lapply(1:nrow(param_combinations), function(i) {
    params <- param_combinations[i, ]
    
    if (parallel) {
      if (is.null(cores)) cores <- detectCores() - 1
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      
      library(doRNG)
      if (!is.null(seed_number)) {
        registerDoRNG(seed_number)
      }
      
      # Export necessary functions and variables 
      clusterExport(cl, list("mimic_control_data", "allocate_treatments", 
                             "allocate_within_treatments",
                             "permute_colonies", "adjust_control_sd",
                             "run_model", 
                             "run_tests", "test_difference", "test_equivalence", 
                             "test_Hotopp_equivalence", "simulate_once", "reduce", 
                             "full_join", "data", "formula", "family", "delta_e", 
                             "mimic_data", "var_comps_list", "mimic_args", "reallocate_colonies", 
                             "permute", "covariates", "effect_timing",
                             "alpha_d", "alpha_e", "alpha_Hotopp_control", "tests", 
                             "assessment_scope", "apply_effects", "assign_treatments", "assign_control_sites"), envir = environment())
      
      results <- foreach(run_number = 1:n_simulations, 
                         .packages = c("tidyverse", "glmmTMB", "lme4", "emmeans", "data.table", "anticlust")) %dopar% {
                           tryCatch({ # catch errors
                             simulate_once(
                               data = data, 
                               mimic_data = mimic_data,
                               n_sites = params$n_sites, 
                               n_colonies = params$n_colonies, 
                               var_comps_list = var_comps_list,  
                               mimic_args = mimic_args,
                               reallocate_colonies = reallocate_colonies,
                               permute = permute,
                               covariates = covariates,
                               control_sd_factor = params$control_sd_factor,
                               effect_size = params$effect_size, 
                               effect_timing = effect_timing,
                               formula = formula, 
                               family = family,
                               delta_e = params$delta_e, 
                               alpha_d = alpha_d, 
                               alpha_e = params$alpha_e, 
                               alpha_Hotopp_control = alpha_Hotopp_control, 
                               tests = tests, 
                               assessment_scope = assessment_scope,
                               run_number = run_number
                             )
                           }, error = function(e) {
                             cat("Simulation error in run", run_number, ":", conditionMessage(e), "\n")
                             return(NULL)
                           })
                         }
      
      
      stopCluster(cl)
      results <- bind_rows(results)
      
    } else {
      # Sequential execution for run_numbers
      results <- bind_rows(lapply(1:n_simulations, function(run_number) {
        simulate_once(
          data = data, 
          mimic_data = mimic_data,
          n_sites = params$n_sites, 
          n_colonies = params$n_colonies, 
          var_comps_list = var_comps_list,  
          mimic_args = mimic_args,
          reallocate_colonies = reallocate_colonies,
          permute = permute,
          covariates = covariates,
          control_sd_factor = params$control_sd_factor,
          effect_size = params$effect_size, 
          effect_timing = effect_timing,
          formula = formula, 
          family = family,
          delta_e = params$delta_e, 
          alpha_d = alpha_d, 
          alpha_e = params$alpha_e, 
          alpha_Hotopp_control = alpha_Hotopp_control, 
          tests = tests, assessment_scope = assessment_scope,
          run_number = run_number)
      }))
    }
    
    return(results)
  }))
  
  end_time <- Sys.time()
  
  run_time <- difftime(end_time, start_time, units = "auto")
  run_time_mins <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  days_for_simulation <- extrapolate_days_for_all_simulations(run_time_mins, n_simulations)
  
  print(run_time)
  print(paste("Rough estimate for 10000 runs: ", round(days_for_simulation, 2), "days"))
  
  if ("equi" %in% tests) {
    simulation_results <- as_tibble(simulation_results) %>% 
      dplyr::select(
        Assessment, effect_size, n_sites, n_colonies,
        ratio, starts_with("risk"),
        starts_with("p_value"),
        starts_with("CI_"),
        Reallocation, run_number, 
        alpha_e, delta_e, 
        control_sd_factor,
        n_negatives, lowest_n_bees,
        everything()
      )
  }
  
  return(simulation_results)
}

# Function to calculate trust rates
calculate_trust_rates <- function(simulation_results,
                                  group_by = c("effect_size", "effect_timing", 
                                               "n_sites", "n_colonies", 
                                               "Assessment")) {
  rates <- list()
  
  # Use across(all_of(group_by)) for grouping
  if ("risk_diff" %in% names(simulation_results)) {
    rates$trust_rate_diff <- simulation_results %>%
      group_by(across(all_of(group_by))) %>%
      summarise(trust_rate_diff = 1 - mean(risk_diff, na.rm = TRUE), 
                .groups = "drop")
  }
  
  if ("risk_equi" %in% names(simulation_results)) {
    rates$trust_rate_equi <- simulation_results %>%
      group_by(across(all_of(group_by))) %>%
      summarise(trust_rate_equi = 1 - mean(risk_equi, na.rm = TRUE), 
                .groups = "drop")
  }
  
  if ("risk_equi_Hotopp" %in% names(simulation_results)) {
    rates$trust_rate_equi_Hotopp <- simulation_results %>%
      group_by(across(all_of(group_by))) %>%
      summarise(trust_rate_equi_Hotopp = 1 - mean(risk_equi_Hotopp, na.rm = TRUE), 
                .groups = "drop")
  }
  
  rates_df <- suppressMessages(reduce(rates, full_join))
  
  return(rates_df)
}

prepare_trust_rate_plot_data <- function(data) {
  plot_data <- data %>% 
    pivot_longer(cols = names(data)[grepl("trust_rate_", names(data))], 
                 names_to = "Test", 
                 values_to = "trust_rate", 
                 names_prefix = "trust_rate_") %>% 
    mutate(N_sites = factor(as.character(n_sites), levels = sort(unique(n_sites))),
           Test = gsub("equi", "EFSA equivalence", Test),
           Test = gsub("diff", "Difference", Test),
           Test = gsub("EFSA equivalence_Hotopp", "Hotopp equivalence", Test),
           Effect_size = factor(as.character(-100*effect_size),
                                levels = -100*unique(effect_size)),
           True_risk = factor(ifelse(effect_size <= 0.1, "Low risk", "High risk"),
                              levels = c("Low risk", "High risk"))
    )
  return(plot_data)
}

filter_and_recode_assessments <- function(data) {
  
  across_assessments <- data %>% 
    filter(Assessment == "Across assessments") 
  
  final_assessment <- data %>% 
    filter(Assessment == 4) %>% 
    mutate(Assessment = "Final assessment")
  
  recoded_subset <- bind_rows(across_assessments, final_assessment)
  
  recoded_subset$Assessment <- factor(recoded_subset$Assessment,
                                      levels = c("Across assessments", 
                                                 "Final assessment"))
  
  recoded_subset
}

slice_simulation_results <- function(simulation_results, slices = 10, 
                                     group_by = c("effect_size", "effect_timing", 
                                                  "n_sites", "n_colonies", 
                                                  "Assessment", "Dataset", 
                                                  "Study")) {
  
  simulation_results_sliced <- simulation_results %>% 
    group_by(across(all_of(group_by))) %>% 
    mutate(Slice = (row_number() - 1) %% slices + 1) %>% 
    ungroup()
  
  simulation_results_sliced
}

calculate_and_reformat_trust_rates <- function(simulation_results,
                                               group_by = c("effect_size", "effect_timing", 
                                                            "n_sites", "n_colonies", 
                                                            "Assessment", "Dataset", 
                                                            "Study", "Reallocation"),
                                               prepare_plot_data = TRUE,
                                               refactor_assessments = TRUE,
                                               slices = NULL
) {
  
  if(is.null(slices) == FALSE){
    simulation_results <- slice_simulation_results(simulation_results, 
                                                   slices = slices, 
                                                   group_by = group_by)
    
    group_by <- c(group_by, "Slice")
  }
  
  trust_rates <- calculate_trust_rates(simulation_results, group_by = group_by)
  
  if(prepare_plot_data){
    trust_rates <- prepare_trust_rate_plot_data(trust_rates)
  }
  
  if(refactor_assessments){
    trust_rates <- filter_and_recode_assessments(trust_rates)
  }
  
  trust_rates
}

# Function to calculate false alarm and false trust rates
obtain_misjudgement_rates_by_effect_size <- function(data, delta_e = 0.1, 
                                                     IDs = c("n_sites", "n_colonies", "Assessment", "effect_timing")) {
  
  data <- data %>% mutate(effect_size_perc = as.integer(round(100*effect_size)))
  
  # Filter data into two categories: no real risk and real risk
  low_risk_data   <- data %>% filter(effect_size <= delta_e)
  high_risk_data  <- data %>% filter(effect_size > delta_e)
  
  rates_list <- list()
  
  if(is.null(low_risk_data) == FALSE) {
    
    if ("trust_rate_diff" %in% names(low_risk_data)) {
      rates_list$false_alarms_diff <-   low_risk_data  %>% mutate(false_alarm_rate_diff = 1 - trust_rate_diff) %>% 
        pivot_wider(id_cols = all_of(IDs),
                    names_from = effect_size_perc, 
                    values_from = false_alarm_rate_diff,
                    names_prefix = "false_alarm_rate_diff_")
    }
    
    if ("trust_rate_equi" %in% names(low_risk_data)) {
      rates_list$false_alarms_equi <-   low_risk_data  %>% mutate(false_alarm_rate_equi = 1 - trust_rate_equi) %>% 
        pivot_wider(id_cols = all_of(IDs),
                    names_from = effect_size_perc, 
                    values_from = false_alarm_rate_equi,
                    names_prefix = "false_alarm_rate_equi_")
    }
    
    if ("trust_rate_equi_Hotopp" %in% names(low_risk_data)) {
      rates_list$false_alarms_equi_Hotopp <-   low_risk_data  %>% mutate(false_alarm_rate_equi_Hotopp = 1 - trust_rate_equi_Hotopp) %>% 
        pivot_wider(id_cols = all_of(IDs),
                    names_from = effect_size_perc, 
                    values_from = false_alarm_rate_equi_Hotopp,
                    names_prefix = "false_alarm_rate_equi_Hotopp_")
    }
  }
  
  if(is.null(high_risk_data) == FALSE) {
    
    if ("trust_rate_diff" %in% names(high_risk_data)) {
      rates_list$false_trusts_diff <- high_risk_data %>% 
        pivot_wider(id_cols = all_of(IDs),
                    names_from = effect_size_perc, 
                    values_from = trust_rate_diff,
                    names_prefix = "false_trust_rate_diff_")
    }
    
    if ("trust_rate_equi" %in% names(high_risk_data)) {
      rates_list$false_trusts_equi <- high_risk_data %>% 
        pivot_wider(id_cols = all_of(IDs),
                    names_from = effect_size_perc, 
                    values_from = trust_rate_equi,
                    names_prefix = "false_trust_rate_equi_")
    }
    
    if ("trust_rate_equi_Hotopp" %in% names(high_risk_data)) {
      rates_list$false_trusts_equi_Hotopp <- high_risk_data %>% 
        pivot_wider(id_cols = all_of(IDs),
                    names_from = effect_size_perc, 
                    values_from = trust_rate_equi_Hotopp,
                    names_prefix = "false_trust_rate_equi_Hotopp_")
    }
  }
  
  misjudgement_rates <- suppressMessages(reduce(rates_list, full_join))
  
  return(misjudgement_rates)
}

calculate_misjudgement_rates_by_effect_size <- function(simulation_results, delta_e = 0.1, 
                                                        IDs = c("n_sites", 
                                                                "n_colonies", 
                                                                "Assessment", 
                                                                "effect_timing")) {
  
  group_by <- c("effect_size", IDs)
  
  trust_rates <- calculate_trust_rates(simulation_results, group_by = group_by)
  
  output <- obtain_misjudgement_rates_by_effect_size(trust_rates, IDs = IDs)
  
  output
}
