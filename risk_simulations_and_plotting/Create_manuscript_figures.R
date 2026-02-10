library(here)
library(ggplot2)
library(cowplot)
library(scales)
library(tidyr)
library(openxlsx)
library(rcartocolor)
library(paletteer)
library(ggrepel)
library(grid)  
library(patchwork)
library(DHARMa)

source(here::here("R scripts", "Read original data.R"))
source(here::here("R scripts", "Functions to simulate risk assessments.R"))

read_results <- TRUE

colors_tests <- c("#009E73", "#CC79A7","#000000")
colors_real_vs_simulated <- c("#56B4E9", "#E69F00")
colors_adjustment_strategy <- c("#999999", "#CC3300", "#0072B2", "#6A3D9A")
colors_effect_sizes <- c("#ffba08", "#e85d04", "#d00000", "#370617")

## Plotting functions ------------
plot_trust_vs_sites_by_effect_size <- function(data,
                                               labels = "no",
                                               hline_y = NULL,
                                               color_values = NULL,
                                               facets = . ~ Test,
                                               show_legend = FALSE,
                                               base_size = 13,
                                               ylab = "Probability pesticide is \n falsely classified as 'low risk'") {
  library(dplyr)
  library(ggplot2)
  
  # Ensure factor levels are consistent
  data <- data %>%
    mutate(Test = factor(Test, levels = unique(Test)))
  
  # Set margin space based on label setting
  left_margin <- if (labels %in% c("both", "both left")) 0.15 else 0.05
  right_margin <- if (labels %in% c("right", "both")) 0.15 else 0.05
  
  # Prepare label data if requested
  label_data <- NULL
  if (labels %in% c("right", "both")) {
    right_panel <- levels(data$Test)[length(levels(data$Test))]
    label_data <- data %>%
      filter(Test == right_panel) %>%
      group_by(Effect_size, True_risk, Test) %>%
      filter(n_sites == max(n_sites)) %>%
      ungroup() %>%
      mutate(hjust = -0.1)
  }
  if (labels == "both") {
    left_panel <- levels(data$Test)[1]
    label_data <- bind_rows(
      label_data,
      data %>%
        filter(Test == left_panel) %>%
        group_by(Effect_size, True_risk, Test) %>%
        filter(n_sites == min(n_sites)) %>%
        ungroup() %>%
        mutate(hjust = 1.1)
    )
  }
  if (labels == "both left") {
    label_data <- data %>%
      group_by(Effect_size, True_risk, Test) %>%
      filter(n_sites == min(n_sites)) %>%
      ungroup() %>%
      mutate(hjust = 1.1)
  }
  
  # Base plot
  p <- ggplot(data, aes(
    x = n_sites,
    y = trust_rate * 100,
    group = Effect_size,
    color = Effect_size,
    linetype = Effect_size
  )) +
    geom_line(linewidth = 1.2)
  
  if (!is.null(label_data)) {
    p <- p +
      geom_text(
        data = label_data,
        aes(label = paste0(Effect_size, "%"), hjust = hjust),
        size = 3.5,
        fontface = "bold",
        vjust = 0.5,
        show.legend = FALSE
      )
  }
  
  if (!is.null(hline_y)) {
    p <- p + geom_hline(yintercept = hline_y, linetype = "dashed")
  }
  
  # Facet layout
  p <- p + facet_grid(facets)
  
  # Axis + theme
  p <- p +
    theme_test(base_size = base_size) +
    labs(
      y = ylab,
      x = "Number of sites per treatment"
    ) +
    scale_y_continuous(
      labels = scales::percent_format(scale = 1),
      breaks = seq(0, 100, 10)
    ) +
    scale_x_continuous(
      breaks = unique(data$n_sites),
      expand = expansion(mult = c(left_margin, right_margin))
    ) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      panel.grid.major = element_line(color = "lightgrey", linewidth = 0.05, linetype = 3),
      strip.text.y = element_text(angle = 0, hjust = 0.5),
      strip.background = element_rect(fill = "white"),
      axis.text = element_text(size = base_size),
      axis.title = element_text(size = base_size),
      strip.text = element_text(size = base_size - 1)
    )
  
  if (!is.null(color_values)) {
    p <- p + scale_color_manual(values = color_values)
  }
  
  return(p)
}

generate_pareto_plot <- function(data,
                                 assessment_filter = "Final assessment",
                                 false_high_risk_effect_size = 0,
                                 false_low_risk_effect_size = 0.11,
                                 n_sites_filter = 10,
                                 color_values = colors_tests) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(stringr)
  
  # Filter assessment
  data_filtered <- data %>%
    filter(Assessment == assessment_filter)
  
  # Format column names as integers after multiplying by 100
  data_wide <- pivot_wider(data_filtered,
                           id_cols = c("effect_timing", "alpha_e", "n_sites", "n_colonies", "Assessment", "Test"),
                           names_from = Effect_size,
                           names_prefix = "Trust_rate_",
                           values_from = trust_rate)
  
  # Clean column names
  colnames(data_wide) <- str_replace(colnames(data_wide), "Trust_rate_-", "Trust_rate_")
  
  effect_high <- paste0("Trust_rate_", round(false_high_risk_effect_size * 100))
  effect_low  <- paste0("Trust_rate_", round(false_low_risk_effect_size * 100))
  
  # Build pareto data
  pareto_data <- data_wide %>%
    filter(Test %in% c("EFSA equivalence", "Hotopp equivalence")) %>%
    mutate(
      False_high_risk = 100 * (1 - .data[[effect_high]]),
      False_low_risk  = 100 * .data[[effect_low]],
      Correct_high_risk = 100 - False_high_risk,
      Correct_low_risk  = 100 - False_low_risk,
      Alpha_percent = paste0(100 * alpha_e, "%")
    )
  
  # Define label sets
  label_alphas <- paste0(c(1, 2, 3, 4, 5, 10, 15), "%")
  efsa_lower <- paste0(c(20, 30, 40, 50), "%")
  efsa_further_lower <- paste0(c(25, 35, 45), "%")
  hotopp_upper <- paste0(c(20, 25, 35, 45), "%")
  hotopp_further_upper <- paste0(c(30, 40, 50), "%")
  
  # Create label positions
  manual_labels_shared <- pareto_data %>%
    filter(n_sites == n_sites_filter, Alpha_percent %in% label_alphas) %>%
    mutate(label_x = False_low_risk + 2, label_y = False_high_risk + 5)
  
  manual_labels_efsa_lower <- pareto_data %>%
    filter(n_sites == n_sites_filter, Test == "EFSA equivalence", Alpha_percent %in% efsa_lower) %>%
    mutate(label_x = False_low_risk + 0.5, label_y = False_high_risk - 8)
  
  manual_labels_efsa_further_lower <- pareto_data %>%
    filter(n_sites == n_sites_filter, Test == "EFSA equivalence", Alpha_percent %in% efsa_further_lower) %>%
    mutate(label_x = False_low_risk + 0.5, label_y = False_high_risk - 14)
  
  manual_labels_hotopp_upper <- pareto_data %>%
    filter(n_sites == n_sites_filter, Test == "Hotopp equivalence", Alpha_percent %in% hotopp_upper) %>%
    mutate(label_x = False_low_risk + 0.5, label_y = False_high_risk + 5)
  
  manual_labels_hotopp_further_upper <- pareto_data %>%
    filter(n_sites == n_sites_filter, Test == "Hotopp equivalence", Alpha_percent %in% hotopp_further_upper) %>%
    mutate(label_x = False_low_risk + 0.5, label_y = False_high_risk + 10)
  
  manual_labels_all <- bind_rows(
    manual_labels_shared,
    manual_labels_efsa_lower,
    manual_labels_efsa_further_lower,
    manual_labels_hotopp_upper,
    manual_labels_hotopp_further_upper
  )
  
  # Subset for plotting
  plot_data <- pareto_data %>% filter(n_sites == n_sites_filter)
  
  if(false_low_risk_effect_size > 0){
    false_low_risk_effect_size <- -1*false_low_risk_effect_size
  }
  
  if(false_high_risk_effect_size > 0){
    false_high_risk_effect_size <- -1*false_high_risk_effect_size
  }
  
  # Plot
  ggplot(plot_data, aes(x = False_low_risk, y = False_high_risk, color = Test)) +
    geom_point(size = 4, alpha = 0.6) +
    geom_line() +
    geom_segment(data = manual_labels_all,
                 aes(x = False_low_risk, y = False_high_risk,
                     xend = label_x, yend = label_y, color = Test),
                 inherit.aes = FALSE,
                 linewidth = 0.3,
                 linetype = "dotted",
                 show.legend = FALSE) +
    geom_text(data = manual_labels_all,
              aes(x = label_x, y = label_y, label = Alpha_percent, color = Test),
              inherit.aes = FALSE,
              hjust = -0.1,
              vjust = 0.2,
              size = 3.5,
              fontface = "bold",
              show.legend = FALSE) +
    xlab(paste0("Falsely classified ", false_low_risk_effect_size * 100, "% true effect size as 'low risk' (%)")) +
    ylab(paste0("Falsely classified ", false_high_risk_effect_size * 100, "% true \neffect size as 'high risk' (%)")) +
    scale_x_continuous(labels = percent_format(scale = 1), breaks = seq(0, 100, 10)) +
    scale_y_continuous(labels = percent_format(scale = 1), breaks = seq(0, 100, 10)) +
    scale_color_manual(values = color_values) +
    coord_cartesian(clip = "off") +
    theme_test() +
    theme(
      legend.position = c(0.92, 0.92),
      legend.justification = c("right", "top"),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 13),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12)
    )
}

plot_power_vs_sites <- function(data, hline_data, x_breaks, colors, linetypes, break_y_title = FALSE) {
  
  if (break_y_title) {
    y_title <- "Probability pesticide is \nclassified as 'low risk'"
  } else {
    y_title <- "Probability pesticide is classified as 'low risk'"
  }
  
  ggplot(
    data,
    aes(
      x = n_sites,
      y = trust_rate * 100,
      color = Adjustment_strategy,
      linetype = Adjustment_strategy,
      group = Adjustment_strategy
    )
  ) +
    geom_line(linewidth = 1.2) +
    facet_grid(. ~ Effect_size_perc) +
    geom_hline(data = hline_data, aes(yintercept = hline_y), linetype = "dashed") +
    labs(y = y_title, x = "Number of sites per treatment") +
    scale_y_continuous(labels = scales::percent_format(scale = 1),
                       breaks = seq(0, 100, 20), limits = c(0, 100)) +
    scale_x_continuous(breaks = x_breaks) +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = linetypes) +
    theme_test() +
    theme(
      panel.grid.major.y = element_line(color = "lightgrey", linewidth = 0.05, linetype = 3),
      panel.grid.major.x = element_blank(),
      strip.text = element_text(angle = 0, hjust = 0.5, size = 13),
      strip.background = element_rect(fill = "white"),
      legend.position = "none",
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )
}

calculate_efsa_2013_crossing_sites <- function(data, threshold = 0.2, effect_size = 0.07, 
                                               assessment_filter = NULL) {
  data_filtered <- data %>%
    { if (!is.null(assessment_filter)) dplyr::filter(., Assessment == assessment_filter) else . } %>%
    filter(Test == "Difference" & effect_size == !!effect_size) %>%
    mutate(risk_rate = trust_rate) %>%
    arrange(Adjustment_strategy, n_sites)
  
  efsa_crossings <- data_filtered %>%
    group_by(Adjustment_strategy) %>%
    group_modify(~{
      .x <- arrange(.x, n_sites)
      ns <- .x$n_sites
      vec <- .x$risk_rate
      
      crossing_found <- FALSE
      x_cross <- NA
      
      for (i in seq_len(length(vec) - 1)) {
        y1 <- vec[i]
        y2 <- vec[i + 1]
        
        if ((y1 < threshold && y2 >= threshold) || (y1 > threshold && y2 <= threshold)) {
          x1 <- ns[i]
          x2 <- ns[i + 1]
          x_cross <- x1 + (threshold - y1) / (y2 - y1) * (x2 - x1)
          crossing_found <- TRUE
          break
        }
      }
      
      if (crossing_found) {
        tibble(n_sites_min = ceiling(x_cross))
      } else {
        tibble(n_sites_min = NA)
      }
    }) %>%
    ungroup() %>%
    filter(!is.na(n_sites_min))
  
  return(efsa_crossings)
}

add_downward_arrows_for_power_threshold <- function(plot, crossing_points, target_y) {
  plot +
    geom_segment(data = crossing_points,
                 aes(x = crossing_x, xend = crossing_x, y = target_y, yend = -Inf,
                     color = Adjustment_strategy,
                     linetype = Adjustment_strategy),
                 arrow = arrow(type = "closed", length = unit(0.1, "inches")),
                 inherit.aes = FALSE,
                 linewidth = 0.5)
}

add_accounting_legend <- function(plot, legend_position = "right", legend_ncol = 1) {
  plot +
    theme(
      legend.position = legend_position,
      legend.background = element_rect(fill = "white", color = "white", linewidth = 0.5),
      legend.box = "vertical",
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12)
    ) +
    guides(
      color = guide_legend(order = 1, title = "Adjustment for \ninitial number of bees", ncol = legend_ncol),
      linetype = guide_legend(order = 1, title = "Adjustment for \ninitial number of bees", ncol = legend_ncol)
    )
}

add_required_sites_shading <- function(plot, shaded_areas, fill_color = "#F0E442") {
  plot +
    geom_rect(data = shaded_areas,
              aes(xmin = x_min, xmax = x_max, ymin = -Inf, ymax = Inf), fill = fill_color,
              inherit.aes = FALSE,
              alpha = 0.2)
}

add_EFSA_annotations_rotated <- function(plot, shaded_areas, y_pos = 0, 
                                         y_arrow_1 = NULL, y_arrow_2 = NULL,
                                         x_arrow_1_from_right = 2, x_arrow_2_from_right = 4, 
                                         label = NULL) {
  require(ggplot2)
  require(grid)
  
  # Determine x positions from full range
  x_accounted <- min(shaded_areas$x_min, na.rm = TRUE)
  x_not_accounted <- max(shaded_areas$x_max, na.rm = TRUE)
  x_text <- mean(c(x_accounted, x_not_accounted), na.rm = TRUE)
  
  if(is.null(y_arrow_1)){
    main_label <- "Sites required under EFSA GD 2013 \nwhen initial bee number was accounted for or not"
    
    if(!is.null(label)) {
      main_label <- label
    }
    
    plot <- plot +
      annotate("text", x = x_text, y = y_pos,
               label = main_label,
               angle = 90, hjust = 0, vjust = 0.5, size = 3.5) 
    
    
    
  } else{
    main_label <- "Sites required under EFSA GD 2013 \nwhen initial bee number was accounted for  or not"
    
    if(!is.null(label)) {
      main_label <- label
    }
    
    plot <- plot +
      annotate("text", x = x_text, y = y_pos,
               label = main_label,
               angle = 90, hjust = 0, vjust = 0.5, size = 3.5) +
      
      # ← arrow pointing to the left (accounted)
      annotate("segment",
               x = x_not_accounted - x_arrow_1_from_right, xend = x_accounted,
               y = y_arrow_1, yend = y_arrow_1,
               arrow = arrow(type = "closed", length = unit(0.05, "inches")),
               color = "black", linewidth = 0.4) +
      
      # → arrow pointing to the right (not accounted)
      annotate("segment",
               x = x_not_accounted - x_arrow_2_from_right, xend = x_not_accounted,
               y = y_arrow_2, yend = y_arrow_2,
               arrow = arrow(type = "closed", length = unit(0.05, "inches")),
               color = "black", linewidth = 0.4)
  }
  plot  
}

generate_multi_panel_power_plot <- function(data,
                                            assessment_filter = "Final assessment",
                                            effect_sizes = c(0, 0.05),
                                            layout = "vertical",  # or "horizontal"
                                            show_legend = TRUE,
                                            add_EFSA_label = TRUE,
                                            EFSA_label_y_pos = 0,
                                            target_y = 80,
                                            colors = colors_adjustment_strategy,
                                            linetypes = c("solid", "dotted", "dashed", "dotdash"),
                                            legend_position = c(0.35, 0.4),
                                            x_breaks = NULL, 
                                            test = c("EFSA equivalence", "Difference", "Hotopp equivalence")) {
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  
  test <- match.arg(test)
  
  # Filter data for selected assessment
  plot_data <- data %>%
    filter(Assessment == assessment_filter)
  
  # Calculate shaded area and EFSA site thresholds
  efsa_crossings <- calculate_efsa_2013_crossing_sites(
    data,
    assessment_filter = assessment_filter
  )
  shaded_area <- efsa_crossings %>%
    summarise(x_min = min(n_sites_min), x_max = max(n_sites_min))
  
  # Filter only for EFSA equivalence tests and selected effect sizes
  plot_data_filtered <- plot_data %>%
    filter(effect_size %in% effect_sizes, Test == test)
  
  # Add hline information
  hline_data <- plot_data_filtered %>%
    mutate(hline_y = ifelse(effect_size < 0.1, target_y, 20))
  
  # Determine crossing points
  # crossing_points <- plot_data_filtered %>%
  #   mutate(trust_rate_percent = trust_rate * 100) %>%
  #   group_by(Effect_size_perc, Adjustment_strategy) %>%
  #   arrange(n_sites) %>%
  #   summarise(
  #     crossing_x = {
  #       lower_idx <- max(which(trust_rate_percent <= target_y))
  #       upper_idx <- min(which(trust_rate_percent >= target_y))
  #       if (length(lower_idx) == 0 | length(upper_idx) == 0 | lower_idx == upper_idx) NA else {
  #         x1 <- n_sites[lower_idx]; x2 <- n_sites[upper_idx]
  #         y1 <- trust_rate_percent[lower_idx]; y2 <- trust_rate_percent[upper_idx]
  #         x1 + (target_y - y1) / (y2 - y1) * (x2 - x1)
  #       }
  #     },
  #     .groups = "drop"
  #   ) %>%
  #   filter(!is.na(crossing_x))
  
  crossing_points <- plot_data_filtered %>%
    mutate(trust_rate_percent = trust_rate * 100) %>%
    group_by(Effect_size_perc, Adjustment_strategy) %>%
    arrange(n_sites, .by_group = TRUE) %>%
    summarise(
      crossing_x = {
        x <- n_sites
        y <- trust_rate_percent
        
        # Need at least two points
        if (length(y) < 2 || all(is.na(y))) {
          NA_real_
        } else {
          # Find first adjacent segment that crosses target_y in either direction
          idx <- which((y[-length(y)] < target_y & y[-1] >= target_y) |
                         (y[-length(y)] > target_y & y[-1] <= target_y))[1]
          
          if (is.na(idx)) {
            # No crossing found. Optional: if already above at first point, use first x.
            if (!is.na(y[1]) && y[1] >= target_y) x[1] else NA_real_
          } else {
            x1 <- x[idx];   x2 <- x[idx + 1]
            y1 <- y[idx];   y2 <- y[idx + 1]
            if (is.na(y1) || is.na(y2) || y2 == y1) x2
            else x1 + (target_y - y1) / (y2 - y1) * (x2 - x1)
          }
        }
      },
      .groups = "drop"
    ) %>%
    filter(!is.na(crossing_x))
  
  
  # Create a panel plot for each effect size
  plot_list <- lapply(effect_sizes, function(eff) {
    plot_df <- plot_data_filtered %>% filter(effect_size == eff)
    hline_df <- hline_data %>% filter(effect_size == eff)
    
    if(is.null(x_breaks)){
      x_breaks <- sort(unique(plot_df$n_sites))
    }
    
    crossing_df <- crossing_points %>%
      filter(grepl(paste0(eff * 100, "%"), Effect_size_perc))
    
    p <- plot_power_vs_sites(
      data = plot_df,
      hline_data = hline_df,
      x_breaks = x_breaks,
      colors = colors,
      linetypes = linetypes,
      break_y_title = TRUE
    ) %>%
      add_downward_arrows_for_power_threshold(crossing_df, target_y) %>%
      add_required_sites_shading(shaded_area)
    
    if (add_EFSA_label) {
      p <- add_EFSA_annotations_rotated(p, shaded_area,
                                        label = "Sites required \nunder EFSA GD 2013",
                                        y_pos = EFSA_label_y_pos)
    }
    
    if (eff == effect_sizes[1] && show_legend) {
      p <- add_accounting_legend(p, legend_position = legend_position)
    }
    
    return(p)
  })
  
  # Combine into a panel
  final_plot <- if (layout == "horizontal") {
    plot_grid(plotlist = plot_list, nrow = 1)
  } else {
    plot_grid(plotlist = plot_list, ncol = 1)
  }
  
  return(final_plot)
}

get_accounting_legend <- function(data,
                                  assessment_filter,
                                  eff,
                                  test,
                                  colors,
                                  linetypes,
                                  legend_position = "right",
                                  legend_ncol = 1) {
  
  plot_df <- data %>%
    dplyr::filter(
      Assessment == assessment_filter,
      effect_size == eff,
      Test == test
    )
  
  if (nrow(plot_df) == 0) {
    stop("No rows after filtering. Check Assessment / effect_size / Test.")
  }
  
  # dummy hline data (required by plot_power_vs_sites)
  hline_dummy <- tibble::tibble(hline_y = NA_real_)
  
  p <- plot_power_vs_sites(
    data = plot_df,
    hline_data = hline_dummy,
    x_breaks = sort(unique(plot_df$n_sites)),
    colors = colors,
    linetypes = linetypes,
    break_y_title = FALSE
  ) %>%
    add_accounting_legend(
      legend_position = legend_position,
      legend_ncol = legend_ncol
    )
  
  cowplot::get_legend(p)
}

plot_trust_vs_effect_size <- function(data, 
                                                color = Dataset, 
                                                color_lab = "Dataset",
                                                color_values = colors_real_vs_simulated,
                                                facet_params = Test ~ Assessment,
                                                slices = NULL) {
  
  plot_data <- calculate_and_reformat_trust_rates(data, slices = slices)
  
  plot_data <- plot_data %>% mutate(
    # Test = recode(Test, "Difference" = "Point-null"),
    Test = factor(Test, levels = c("Difference", "EFSA equivalence", "Hotopp equivalence"))
  )
  
  
  # Check if the first facetting variable is N_sites or n_sites
  first_facet <- as.character(as.formula(facet_params)[[2]])
  add_sec_axis <- first_facet %in% c("N_sites", "n_sites")
  
  # Base plot
  
  if(is.null(slices)){
    
    p0 <- ggplot(plot_data, aes(x = -100 * effect_size, y = trust_rate * 100, 
                                color = !!enquo(color))) +
      geom_line(linewidth = 1.2)
    
  } else {
    
    p0 <- ggplot(plot_data, aes(x = -100 * effect_size, y = trust_rate * 100,
                                color = !!enquo(color),
                                group = interaction(Slice, !!enquo(color), Test)
    )) +
      geom_line(alpha = 0.4,  
                linewidth = 0.6)
  }
  
  p <- p0 +
    annotate("segment", x = -7, xend = -10, y = 20, yend = 20, linetype = "dashed", color = "grey") +
    annotate("segment", x = -10, xend = -Inf, y = 20, yend = 20, linetype = "dashed", color = "black") +
    geom_vline(xintercept = -10, linetype = "dashed") +
    geom_vline(xintercept = -7, linetype = "dashed", color = "grey") +
    theme_test() +
    scale_x_reverse(labels = percent_format(scale = 1), breaks = c(-7, seq(0, -100, -5)), 
                    expand = expansion(mult = c(0, 0.1))) +
    scale_color_manual(values = color_values) +
    facet_grid(facet_params) +
    labs(
      color = color_lab,
      y = "Probability pesticide is classified as 'low risk'",
      x = "True effect size"
    ) +
    theme(
      panel.grid.major = element_line(color = "lightgrey", linewidth = 0.05, linetype = 3),
      strip.text.x = element_text(size = 14),
      strip.text.y = element_text(size = 14, angle = 0, hjust = 0.5),
      strip.background = element_rect(fill = "white"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    )
  
  # Add secondary axis conditionally
  if (add_sec_axis) {
    p <- p + scale_y_continuous(
      labels = percent_format(scale = 1),
      breaks = seq(0, 100, 20),
      sec.axis = sec_axis(~ ., name = "Number of sites per treatment", breaks = NULL)
    )
  } else {
    p <- p + scale_y_continuous(
      labels = percent_format(scale = 1),
      breaks = seq(0, 100, 20)
    )
  }
  
  return(p)
}

# False classifications of high-risk pesticides by test ------------
if(read_results){
  path_1 <- "R output/Simulation results/sim_results_2026_high_risk.rds"
  sim_results_high_risk <- readRDS(path_1)
}

## Final assessment
trust_rates_high_risk <- calculate_and_reformat_trust_rates(sim_results_high_risk, 
                                                        group_by = c("effect_size", "effect_timing", 
                                                                     "n_sites", "n_colonies", 
                                                                     "Assessment", "Study")) 

effect_size_data_final <- subset(trust_rates_high_risk, 
                                 effect_size > 0.09 &
                                   Test %in% c("EFSA equivalence", "Hotopp equivalence") &
                                   n_colonies == 6 &
                                   Assessment == "Final assessment") 

false_low_risk_by_effect_size_plot_final <- plot_trust_vs_sites_by_effect_size(
  data = effect_size_data_final,
  labels = "both left",
  hline_y = 20,
  color_values = colors_effect_sizes
)

false_low_risk_by_effect_size_plot_final

## Across assessments
trust_rates_high_risk <- calculate_and_reformat_trust_rates(sim_results_high_risk, 
                                                            group_by = c("effect_size", "effect_timing", 
                                                                         "n_sites", "n_colonies", 
                                                                         "Assessment", "Study")) 

effect_size_data_across <- subset(trust_rates_high_risk, 
                                 effect_size > 0.09 &
                                   Test %in% c("EFSA equivalence", "Hotopp equivalence") &
                                   n_colonies == 6 &
                                   Assessment == "Across assessments") 

false_low_risk_by_effect_size_plot_across <- plot_trust_vs_sites_by_effect_size(
  data = effect_size_data_across,
  labels = "both left",
  hline_y = 20,
  color_values = colors_effect_sizes
)

false_low_risk_by_effect_size_plot_across

# Pareto plots -------
if(read_results){
  sim_results_pareto_alphas <- readRDS("R output/Simulation results/sim_results_2026_pareto_alphas.rds")
}

trust_rates_pareto_alphas <- 
  calculate_and_reformat_trust_rates(sim_results_pareto_alphas, 
                                     group_by = c("effect_size", "effect_timing", 
                                                  "n_sites", "n_colonies", 
                                                  "Assessment", 
                                                  "alpha_e")) 

pareto_plot_final <- generate_pareto_plot(trust_rates_pareto_alphas,
                     assessment_filter = "Final assessment")



# Combining false_low_risk_by_effect_size_plot with pareto plot ----------

Fig_2_false_trusts_EFSA_vs_Hotopp_final <- plot_grid(
  pareto_plot_final, 
  false_low_risk_by_effect_size_plot_final, 
  align = "v", 
  axis = "lr",     
  nrow = 2,
  labels = "AUTO"
)

Fig_2_false_trusts_EFSA_vs_Hotopp_final

ggsave("R output/Plots/Figures tiff/Fig_2_false_trusts_EFSA_vs_Hotopp_final.tiff", 
       width = 7.5, height = 7.5)

# Fig_S_1_false_trusts_EFSA_vs_Hotopp_across <- plot_grid(
#   false_low_risk_by_effect_size_plot_across, 
#   pareto_plot_across, 
#   align = "v", 
#   axis = "lr",     
#   nrow = 2,
#   labels = "AUTO"
# )
# 
# Fig_S_1_false_trusts_EFSA_vs_Hotopp_across
# ggsave("R output/Plots/Figures tiff/Fig_S_1_false_trusts_EFSA_vs_Hotopp_across.tiff", 
#        width = 7.5, height = 7.5)

# Power in relation to n_sites and accounting strategy -----------
if (read_results) {
  sim_results_n_bees_initial_not_accounted <- readRDS("R output/Simulation results/sim_results_2026_n_bees_initial_not_accounted.rds")
  sim_results_n_bees_initial_balanced <- readRDS("R output/Simulation results/sim_results_2026_n_bees_initial_balanced.rds")
  sim_results_n_bees_initial_included <- readRDS("R output/Simulation results/sim_results_2026_n_bees_initial_included.rds")
  sim_results_n_bees_initial_included_and_balanced <- readRDS("R output/Simulation results/sim_results_2026_n_bees_initial_included_and_balanced.rds")
  
  sim_results_EFSA_GD_refinement <- readRDS("R output/Simulation results/sim_results_2026_EFSA_GD_refinement.rds")
  
}

sim_results_adjust_strategy <- bind_rows(
  sim_results_n_bees_initial_not_accounted,
  sim_results_n_bees_initial_balanced,
  sim_results_n_bees_initial_included,
  sim_results_n_bees_initial_included_and_balanced,
  
  sim_results_EFSA_GD_refinement
)

sim_results_adjust_strategy$risk_diff <- with(sim_results_adjust_strategy, 
                                              ifelse(p_value_diff < 0.05, 1, 0)) 

trust_rates_adjust_strategy <- calculate_and_reformat_trust_rates(
  sim_results_adjust_strategy,
  group_by = c("effect_size", "effect_timing", "n_sites", "n_colonies", "Assessment", "Reallocation", "Predictors", "Study")
) %>%
  mutate(Accounted_for_n_bees_initial = case_when(
    Reallocation %in% c("none", "random") & Predictors == "Treatment * Assessment + (1|Site/Colony)" ~ "No",
    Reallocation == "balanced" & Predictors == "Treatment * Assessment + (1|Site/Colony)" ~ "Through allocation process",
    Reallocation %in% c("none", "random") & Predictors == "Treatment * Assessment + n_bees_initial + (1|Site/Colony)" ~ "Through covariate",
    Reallocation == "balanced" & Predictors == "Treatment * Assessment + n_bees_initial + (1|Site/Colony)" ~ "Through both",
    TRUE ~ NA_character_
  )) %>%
  mutate(
    Accounted_for_n_bees_initial = factor(
      Accounted_for_n_bees_initial,
      levels = c("No", "Through allocation process", "Through covariate", "Through both")),
    Adjustment_strategy = recode(Accounted_for_n_bees_initial, 
                                 "No" = "No adjustment",
                                 "Through allocation process" = "Design adjustment",
                                 "Through covariate" = "Statistical adjustment",
                                 "Through both" = "Combined adjustment"
    ),
    Effect_size_perc = paste0("True effect size = ", Effect_size, "%"),
    Effect_size_perc = factor(Effect_size_perc, levels = unique(Effect_size_perc))
  ) 

View(trust_rates_adjust_strategy %>% filter(trust_rate > 0.8))

## Plot 
generate_multi_panel_power_plot(
  data = trust_rates_adjust_strategy,
  assessment_filter = "Across assessments",
  layout = "horizontal",
  show_legend = FALSE
)

Fig_3_power_vs_sites_by_accounting_final <- generate_multi_panel_power_plot(
  data = trust_rates_adjust_strategy,
  assessment_filter = "Final assessment", 
  x_breaks = c(seq(4, 20, 2), seq(25, 45, 5))
)

Fig_3_power_vs_sites_by_accounting_final

ggsave("R output/Plots/Figures tiff/Fig_3_power_vs_sites_by_accounting_final.tiff", 
       width = 8, height = 6)

# Fig_S_2_power_vs_sites_by_accounting_across <- generate_multi_panel_power_plot(
#   data = trust_rates_n_bees_initial,
#   assessment_filter = "Across assessments"
# )
# 
# Fig_S_2_power_vs_sites_by_accounting_across
# ggsave("R output/Plots/Figures tiff/Fig_S_2_power_vs_sites_by_accounting_across.tiff",
#        width = 8, height = 6)

# Conceptual figure on different tests ---------
# Define effect sizes and 60% CI half-widths that generate the desired outcomes
data_scenarios <- tribble(
  ~Pesticide, ~effect_size, ~ci60_halfwidth,
  # effect_size, ci60_halfwidth chosen to generate desired test outcomes
  "A",         -0.26,        0.06,   
  "B",         -0.11,        0.03,   
  "C",         -0.08,        0.06,   
  "D",         -0.07,        0.02,   
  "E",         -0.02,        0.03    
) %>%
  mutate(
    lower_60 = round(effect_size - ci60_halfwidth, 2),
    upper_60 = round(effect_size + ci60_halfwidth, 2),
    ci90_halfwidth = ci60_halfwidth * 1.9,
    lower_90 = round(effect_size - ci90_halfwidth, 2),
    upper_90 = round(effect_size + ci90_halfwidth, 2),
    
    # Significance logic
    Equivalence = ifelse(lower_60 >= -0.10, "*", "n.s."),
    Inferiority = ifelse(upper_60 < -0.10, "*", "n.s."),
    Difference = ifelse(upper_90 < 0, "*", "n.s."),
    
    # EFSA & our interpretation (optional)
    EFSA_2023 = ifelse(lower_60 < -0.10, "high risk", "low risk"),
    our_interpretation = ifelse(lower_60 < -0.10, "high risk\nnot excluded", "low risk"),
    
    # Final combined test logic with line breaks
    Combined_test = case_when(
      Inferiority == "*" & Difference == "*" ~ "high risk,\nmajor effect (>SPG)",
      Inferiority == "*" ~ "high risk",
      Equivalence == "*" & Difference == "*" ~ "low risk,\nminor effect (<SPG)",
      Equivalence == "*" ~ "low risk",
      Difference == "*" ~ "high risk \nnot excluded,\nstatistically \nsignificant effect",
      TRUE ~ "high risk\nnot excluded"
    ),
    
    Pesticide = factor(Pesticide, levels = rev(c("A", "B", "C", "D", "E")))
  )

# Plot
Test_label_x_position <- 5.4

Fig_4_conceptual_combined_tests <- ggplot(data_scenarios, 
                                         aes(y = effect_size * 100, x = Pesticide)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_90 * 100, ymax = upper_90 * 100), 
                width = 0.15) +
  geom_errorbar(aes(ymin = lower_60 * 100, ymax = upper_60 * 100), 
                width = 0.08, linewidth = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -10, linetype = "dotted") +
  coord_flip() +
  theme_test() +
  theme(axis.text = element_text(size = 11)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), 
                     breaks = seq(-40, 0, 10)) +
  expand_limits(x = c(NA, Test_label_x_position + 1.5), y = c(-40, 155)) +
  
  geom_text(x = Test_label_x_position, y = 10, label = "Equivalence", 
            hjust = 0, size = 4, fontface = "bold", angle = 90) +
  geom_text(aes(x = Pesticide, y = 10, label = Equivalence), 
            hjust = 0.5, size = 4) +
  
  geom_text(x = Test_label_x_position, y = 25, label = "Magnitude > SPG", 
            hjust = 0, size = 4, fontface = "bold", angle = 90) +
  geom_text(aes(x = Pesticide, y = 25, label = Inferiority), 
            hjust = 0.5, size = 4) +
  
  geom_text(x = Test_label_x_position, y = 40, 
            label = "Different from 0%", 
            hjust = 0, size = 4, fontface = "bold", angle = 90) +
  geom_text(aes(x = Pesticide, y = 40, label = Difference), 
            hjust = 0.5, size = 4) +
  
  geom_text(x = Test_label_x_position + 0.9, y = 85, 
            label = "Equivalence test \n(non-inferiority)", hjust = 0.5, vjust = 1,
            fontface = "bold", size = 4) +
  
  geom_text(x = Test_label_x_position + 0.2, y = 65, 
            label = "EFSA 2023\nclassification", 
            hjust = 0.5, fontface = "bold", size = 4) +
  geom_text(aes(x = Pesticide, y = 65, 
                label = EFSA_2023), size = 4, hjust = 0.5) +
  
  geom_text(x = Test_label_x_position + 0.2, y = 105, 
            label = "Our proposed\nwording", hjust = 0.5, 
            fontface = "bold", size = 4) +
  geom_text(aes(x = Pesticide, y = 105, 
                label = our_interpretation), size = 4, hjust = 0.5) +
  
  geom_text(x = Test_label_x_position + 0.9, y = 140, 
            label = "Combined \nthree-pronged \ntest", hjust = 0.5, vjust = 1, fontface = "bold", size = 4) +
  geom_text(aes(x = Pesticide, y = 140, 
                label = Combined_test), size = 4, hjust = 0.5) +
  
  labs(x = NULL, y = NULL)

Fig_4_conceptual_combined_tests 
ggsave("R output/Plots/Figures tiff/Fig_4_conceptual_combined_tests.tiff", width = 9, height = 7)

# Appendix figures ----------- 
## Fig. A1 - Results real vs simulated    --------------- 
## Compare results using dataset simulated in parametric bootstrap to non-parametric bootstrap on original data 
if(read_results){
  sim_results_real <- readRDS("R output/Simulation results/sim_results_2026_real.rds")
  sim_results_simulated <- readRDS("R output/Simulation results/sim_results_2026_simulated.rds")
}

sim_results_real_vs_simulated <- bind_rows(sim_results_real, sim_results_simulated) %>% 
  mutate(CI_span = CI_upper_equi - CI_lower_equi)

# Have a look at variation
sim_results_real_vs_simulated %>% group_by(Dataset) %>% 
  summarize(
    CI_lower_min = min(CI_lower_equi),
    CI_lower_median = median(CI_lower_equi),
    CI_lower_mean = mean(CI_lower_equi),
    CI_lower_max = max(CI_lower_equi)) # less variation in original

sim_results_real_vs_simulated %>% group_by(Dataset) %>% 
  summarize(
    CI_upper_min = min(CI_upper_equi),
    CI_upper_median = mean(CI_upper_equi),
    CI_upper_mean = mean(CI_upper_equi),
    CI_upper_max = max(CI_upper_equi)) # less variation in original

sim_results_real_vs_simulated %>% group_by(Dataset) %>% 
  summarize(
    CI_span_min = min(CI_span),
    CI_span_median = median(CI_span),
    CI_span_mean = mean(CI_span),
    CI_span_max = max(CI_span)) # clearly less variation in original

plot_trust_vs_effect_size(
  data = sim_results_real_vs_simulated
) 

Fig_A1_Trust_by_real_vs_simulated_sliced <- plot_trust_vs_effect_size(
  data = sim_results_real_vs_simulated,
  slices = 10
) 

Fig_A1_Trust_by_real_vs_simulated_sliced 
ggsave("R output/Plots/Figures tiff/Fig_A1_Trust_by_real_vs_simulated_sliced.tiff", 
       plot = Fig_A1_Trust_by_real_vs_simulated_sliced,
       width = 11, height = 7)

## Fig. A2 - Normality of initial number of bees ----------
n_bees_initial <- round(subset(original_control_data, Assessment == 1)$n_bees_initial)

# Shapiro-Wilk test
shapiro_result <- shapiro.test(n_bees_initial)
shapiro_label <- sprintf("Shapiro-Wilk\nW = %.3f\np = %.3f", 
                         shapiro_result$statistic, shapiro_result$p.value)

# Shared theme for alignment
base_theme <- theme_minimal(base_size = 12) +
  theme(plot.margin = margin(10, 10, 10, 10))

# Histogram
hist_n_bees_initial_plot <- ggplot(data.frame(x = n_bees_initial), aes(x)) +
  geom_histogram(aes(y = after_stat(density)), bins = 7, fill = "gray80", color = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(n_bees_initial), sd = sd(n_bees_initial)),
                color = "red", linewidth = 1) +
  labs(title = "Histogram with Normal Fit", x = "Values", y = "Density") +
  base_theme

# Q-Q Plot
qq_n_bees_initial_plot <- ggplot(data.frame(sample = n_bees_initial), aes(sample = sample)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, size = 4,
           label = shapiro_label, fontface = "italic") +
  labs(title = "Q-Q Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  base_theme

Fig_A2_Normality_n_bees_initial <- hist_n_bees_initial_plot + qq_n_bees_initial_plot + plot_layout(guides = "collect")

Fig_A2_Normality_n_bees_initial
ggsave("R output/Plots/Figures tiff/Fig_A2_Normality_n_bees_initial.tiff", 
       width = 8, height = 6)

## Fig. A3 - DHARMa plot of simple model --------------
model_1 <- glmmTMB(n_bees ~ Assessment + (1|Site/Colony), data = original_control_data_234, family = "nbinom2")

tiff("R output/Plots/Figures tiff/Fig_A3_DHARMa_model_residuals.tiff", width = 8, height = 7, units = "in", res = 300)
plot(simulateResiduals(fittedModel = model_1))
dev.off()

## Fig. A4 - Simularity between simulated and original data 
# done in Script Check similarity between simulated and original data 

## Fig. A5 - Trust vs effect size -------------
effect_size_data_final$Assessment <- "Final timepoint"
effect_size_data_across$Assessment <- "Across timepoints"

effect_size_data <- bind_rows(effect_size_data_final, effect_size_data_across)
effect_size_data$Assessment <- factor(effect_size_data$Assessment, levels = c("Final timepoint", "Across timepoints"))

Fig_A5_trust_vs_sites_by_effect_size <- plot_trust_vs_sites_by_effect_size(
  data = effect_size_data,
  facets = Test ~ Assessment,
  labels = "both left",
  hline_y = 20,
  color_values = colors_effect_sizes
)

Fig_A5_trust_vs_sites_by_effect_size

ggsave("R output/Plots/Figures tiff/Fig_A5_trust_vs_sites_by_effect_size.tiff", 
       width = 8.5, height = 6)

## Fig. A6 - Pareto figure ----------------
pareto_plot_across <- generate_pareto_plot(trust_rates_pareto_alphas,
                                           assessment_filter = "Across assessments")

if(read_results){
  sim_results_pareto_alphas_small_sample_size <- 
    readRDS("R output/Simulation results/sim_results_2026_pareto_alphas_small_sample_size.rds")
}

trust_rates_pareto_alphas_small_sample_size <- 
  calculate_and_reformat_trust_rates(sim_results_pareto_alphas_small_sample_size, 
                                     group_by = c("effect_size", "effect_timing", 
                                                  "n_sites", "n_colonies", 
                                                  "Assessment", 
                                                  "alpha_e")) 

pareto_plot_final_small_sample_size <- generate_pareto_plot(trust_rates_pareto_alphas_small_sample_size,
                                                            n_sites_filter = 4, 
                                                            assessment_filter = "Final assessment")

pareto_plot_across_small_sample_size <- generate_pareto_plot(trust_rates_pareto_alphas_small_sample_size,
                                                             n_sites_filter = 4, 
                                                             assessment_filter = "Across assessments")

pareto_plot_final_w_title <- pareto_plot_final + 
  ggtitle("Final assessment: \n10 sites per treatment, 6 colonies per site")
pareto_plot_across_w_title <- pareto_plot_across + 
  ggtitle("Across assessments: \n10 sites per treatment, 6 colonies per site")
pareto_plot_final_small_sample_size_w_title <- pareto_plot_final_small_sample_size + 
  ggtitle("Final assessment: \n4 sites per treatment, 4 colonies per site")
pareto_plot_across_small_sample_size_w_title <- pareto_plot_across_small_sample_size + 
  ggtitle("Across assessments: \n4 sites per treatment, 4 colonies per site")

Fig_A6_pareto <- plot_grid(pareto_plot_final_w_title, pareto_plot_across_w_title,
                           pareto_plot_final_small_sample_size_w_title, pareto_plot_across_small_sample_size_w_title,
                           align = "v", 
                           axis = "lr",     
                           nrow = 2,
                           labels = "AUTO")

Fig_A6_pareto

ggsave("R output/Plots/Figures tiff/Fig_A6_pareto.tiff", 
       width = 10.5, height = 10.5)

## Fig. A7 - Power vs sites EFSA equivalence ----------------
title_grob_final <- ggdraw() +
  draw_label("Final timepoint",
             x = 0.55, hjust = 0.5,
             fontface = "bold",
             size = 14)

title_grob_across <- ggdraw() +
  draw_label("Across timepoints",
             x = 0.55, hjust = 0.5,
             fontface = "bold",
             size = 14)

power_plot_equi_0_5_perc_final <- generate_multi_panel_power_plot(
  data = trust_rates_adjust_strategy,
  effect_sizes = c(0, 0.05),
  target_y = 80,
  assessment_filter = "Final assessment",
  show_legend = T,
  x_breaks = c(seq(4, 20, 2), seq(25, 45, 5)),
  test = "EFSA equivalence",layout = "horizontal"
) 

power_plot_equi_0_5_perc_across <- generate_multi_panel_power_plot(
  data = trust_rates_adjust_strategy,
  effect_sizes = c(0, 0.05),
  target_y = 80,
  assessment_filter = "Across assessments",
  show_legend = F,
  x_breaks = c(seq(4, 20, 2), seq(25, 45, 5)),
  test = "EFSA equivalence",layout = "horizontal"
) 

equi_7_perc_final <- plot_grid(
  title_grob_final, power_plot_equi_0_5_perc_final,
  ncol = 1,
  rel_heights = c(0.08, 1)  
)

equi_7_perc_across <- plot_grid(
  title_grob_across, power_plot_equi_0_5_perc_across,
  ncol = 1,
  rel_heights = c(0.08, 1)  
)

Fig_A7_power_vs_sites_equi_0_5_perc <- plot_grid(
  equi_7_perc_final, 
  equi_7_perc_across, 
  ncol = 1)

Fig_A7_power_vs_sites_equi_0_5_perc

ggsave("R output/Plots/Figures tiff/Fig_A7_power_vs_sites_equi_0_5_perc.tiff", 
       width = 10, height = 8)

## Fig. A8 - Power vs sites difference equivalence ----------------

legend_adjustments <- get_accounting_legend(
  data = trust_rates_adjust_strategy,
  assessment_filter = "Across assessments",
  eff = 0.07,
  test = "Difference",
  colors = colors_adjustment_strategy,
  linetypes = c("solid", "dotted", "dashed", "dotdash")
)

power_plot_diff_7_perc_final <- generate_multi_panel_power_plot(
  data = trust_rates_adjust_strategy,
  effect_sizes = c(0.07),
  target_y = 20,
  assessment_filter = "Final assessment",
  EFSA_label_y_pos = 35,
  show_legend = F, 
  x_breaks = c(seq(4, 20, 2), seq(25, 45, 5)),
  test = "Difference"
) 

power_plot_diff_7_perc_across <- generate_multi_panel_power_plot(
  data = trust_rates_adjust_strategy,
  effect_sizes = c(0.07),
  target_y = 20,
  assessment_filter = "Across assessments",
  EFSA_label_y_pos = 35,
  show_legend = F,
  x_breaks = c(seq(4, 20, 2), seq(25, 45, 5)),
  test = "Difference"
)

diff_7_perc_final <- plot_grid(
  title_grob_final, power_plot_diff_7_perc_final,
  ncol = 1,
  rel_heights = c(0.08, 1)  
)

diff_7_perc_across <- plot_grid(
  title_grob_across, power_plot_diff_7_perc_across,
  ncol = 1,
  rel_heights = c(0.08, 1)  
)

diff_7_perc_plot_wo_legend <- plot_grid(
  diff_7_perc_final, 
  diff_7_perc_across, 
  ncol = 1)

Fig_A8_power_vs_sites_diff_7_perc <- plot_grid(diff_7_perc_plot_wo_legend, 
                                               legend_adjustments, 
                                               ncol = 2, rel_widths = c(1, 0.35))

Fig_A8_power_vs_sites_diff_7_perc

ggsave("R output/Plots/Figures tiff/Fig_A8_power_vs_sites_diff_7_perc.tiff", 
       width = 8, height = 8)

# Graphics abstract --------------
# Left side
data_scenarios_graphics <- tribble(
  ~Pesticide, ~effect_size, ~ci60_halfwidth,
  # effect_size, ci60_halfwidth chosen to generate desired test outcomes
  "A",         -0.07,        0.02,   
  "B",         -0.16,        0.09
) %>%
  mutate(
    lower_60 = round(effect_size - ci60_halfwidth, 2),
    upper_60 = round(effect_size + ci60_halfwidth, 2),
    ci90_halfwidth = ci60_halfwidth * 1.9,
    lower_90 = round(effect_size - ci90_halfwidth, 2),
    upper_90 = round(effect_size + ci90_halfwidth, 2),
    Pesticide = factor(Pesticide, levels = rev(c("A", "B"))),
    Difference = c("High risk", "Low risk"),
    Equivalence = c("Low risk", "High risk \nnot \nexcluded")
  )

# Plot

windows(width = 8/2.54, height = 4/2.54)
txt9 <- 9/.pt

ggplot(data_scenarios_graphics, 
       aes(y = effect_size * 100, x = Pesticide)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower_90 * 100, ymax = upper_90 * 100), 
                width = 0.15) +
  geom_errorbar(aes(ymin = lower_60 * 100, ymax = upper_60 * 100), 
                width = 0.1, linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -10, linetype = "dotted", color = "red") +
  coord_flip() +
  theme_minimal() +
  labs(x = NULL, y = "Pesticide effect size") +
  scale_y_continuous(labels = scales::percent_format(scale = 1), 
                     breaks = seq(-10, 0, 10)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(hjust = 0.1, size = 9)) +
  expand_limits(y = c(-30, 60)) +
  geom_text(aes(x = Pesticide, y = 15, label = Difference), 
            hjust = 0.5, size = txt9, lineheight = 1) +
  geom_text(aes(x = Pesticide, y = 45, label = Equivalence), 
            hjust = 0.5, size = txt9, lineheight = 1)

ggsave("R output/Plots/Figures tiff/Graphics_abstract_left.tiff", 
       width = 8, height = 3.5, units = "cm")

