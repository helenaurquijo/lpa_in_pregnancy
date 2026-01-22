################################################################################
#                      Functions used in UKB-RAP analyses                      #
################################################################################


################################################################################
# Function to summarise variables for events counts table

summarise_variable <- function(data, varname, var_type = NULL) {
  if(is.null(var_type)) {
    # Auto-detect type
    if(is.numeric(data[[varname]])) {
      if(all(data[[varname]] %in% c(0, 1, NA), na.rm = TRUE)) {
        var_type <- "bin"
      } else {
        var_type <- "cont"
      }
    } else {
      var_type <- "cat"
    }
  }
  
  if(var_type == "cont") {
    # For continuous variables, get mean, SD, etc.
    vals <- data[[varname]]
    return(sprintf("Mean: %.2f, SD: %.2f, Range: %.2f-%.2f, N: %d", 
                  mean(vals, na.rm = TRUE), 
                  sd(vals, na.rm = TRUE),
                  min(vals, na.rm = TRUE),
                  max(vals, na.rm = TRUE),
                  sum(!is.na(vals))))
  } else if(var_type == "bin") {
    # For binary variables, get count and percentage
    count <- sum(data[[varname]] == 1, na.rm = TRUE)
    total <- sum(!is.na(data[[varname]]))
    return(sprintf("%d (%.2f%%)", count, count/total*100))
  } else {
    # For categorical, get top categories
    tbl <- table(data[[varname]], useNA = "ifany")
    top_n <- head(sort(tbl, decreasing = TRUE), 3)
    return(paste(names(top_n), top_n, sep = ": ", collapse = ", "))
  }
}

################################################################################
# Function to perform association analysis

assoc_test <- function(outcome, cat, df, exp) {
  # Determine if exposure is binary
  is_binary_exposure <- all(df[[exp]] %in% c(0, 1, NA), na.rm = TRUE)
  
  cat(sprintf("\n=== Starting association test ===\n"))
  cat(sprintf("Outcome: %s | Type: %s | Exposure type: %s\n", 
             outcome, cat, ifelse(is_binary_exposure, "binary", "continuous")))
  
  # Check if outcome exists in dataframe
  if (!outcome %in% colnames(df)) {
    warning(paste0("Outcome '", outcome, "' not found in dataframe."))
    return(NULL)
  }
  
  # Specify covariate strings
  pcs_string <- paste0("PC_", 1:40)
  geno_batch <- "geno_batch"
  
  cat("Building formula...\n")
  # Overall analysis (not stratified)
  formula <- reformulate(c(exp, geno_batch, pcs_string), response = outcome)
  
  if (cat == "bin") {
    # Check for complete separation risk
    n_cases <- sum(df[[outcome]] == 1, na.rm = TRUE)
    cat(sprintf("  Overall: %d cases out of %d (%.2f%%)\n", 
               n_cases, nrow(df), 100*n_cases/nrow(df)))
    
    if (n_cases < 5) {
      warning("Too few cases (<5) for reliable logistic regression")
      return(NULL)
    }
    
    cat("Running logistic regression model...\n")
    # Run model
    model <- glm(formula, 
                data = df, 
                family = binomial)
    
    cat("Extracting model coefficients...\n")
    # Extract model summary 
    model_summary <- summary(model)
    coef_index <- which(rownames(model_summary$coefficients) == exp)
    
    if (length(coef_index) == 0) {
      warning("Exposure variable not found in model coefficients")
      return(NULL)
    }
    
    cat("Calculating effect estimates and confidence intervals...\n")
    # Calculate statistics manually
    estimate <- model_summary$coefficients[coef_index, "Estimate"]
    std_error <- model_summary$coefficients[coef_index, "Std. Error"]
    p_value <- model_summary$coefficients[coef_index, "Pr(>|z|)"]
    
    # Calculate CI for log odds
    ci_lower <- estimate - 1.96 * std_error
    ci_upper <- estimate + 1.96 * std_error
    
    # Calculate OR and CI for OR
    or <- exp(estimate)
    or_ci_lower <- exp(ci_lower)
    or_ci_upper <- exp(ci_upper)
    
    cat("Creating results dataframe...\n")
    # Create result data frame
    result <- data.frame(
      outcome = outcome,
      b = estimate,
      se = std_error,
      OR = or,
      p_value = p_value,
      CI_lower = ci_lower,
      CI_upper = ci_upper
    )
    
    # Print result summary
    effect_desc <- if(is_binary_exposure) "comparing groups" else "per unit increase"
    cat("Results summary:\n")
    cat(sprintf("  OR = %.3f (%.3f-%.3f), p = %.4f %s\n", 
               or, or_ci_lower, or_ci_upper, p_value, effect_desc))
    
  } else if (cat == "cont") {
    cat("Running linear regression model...\n")
    model <- lm(formula, data = df)
    
    cat("Extracting model coefficients...\n")
    # Extract model summary 
    model_summary <- summary(model)
    coef_index <- which(rownames(model_summary$coefficients) == exp)
    
    if (length(coef_index) == 0) {
      warning("Exposure variable not found in model coefficients")
      return(NULL)
    }
    
    cat("Calculating effect estimates and confidence intervals...\n")
    # Calculate statistics manually
    estimate <- model_summary$coefficients[coef_index, "Estimate"]
    std_error <- model_summary$coefficients[coef_index, "Std. Error"]
    p_value <- model_summary$coefficients[coef_index, "Pr(>|t|)"]
    
    # Calculate confidence intervals
    ci_lower <- estimate - 1.96 * std_error
    ci_upper <- estimate + 1.96 * std_error
    
    cat("Creating results dataframe...\n")
    # Create result data frame 
    result <- data.frame(
      outcome = outcome,
      b = estimate,
      se = std_error,
      OR = NA,
      p_value = p_value,
      CI_lower = ci_lower,
      CI_upper = ci_upper
    )
    
    # Print result summary
    effect_desc <- if(is_binary_exposure) "comparing groups" else "per unit increase"
    cat("Results summary:\n")
    cat(sprintf("  Beta = %.3f (%.3f-%.3f), p = %.4f %s\n", 
               estimate, ci_lower, ci_upper, p_value, effect_desc))
  } else {
    return(NULL)
  }
  
  cat("=== Association test complete ===\n\n")
  return(result)
}

################################################################################
# Function for association testing with categorical exposure

cat_assoc_test <- function(outcome, cat, data, exp, ref_var) {
  
  cat(sprintf("\n=== Starting association test ===\n"))
  cat(sprintf("Outcome: %s | Type: %s | Reference: Q%s\n", outcome, cat, ref_var))
  
  # Check if outcome exists in dataframe
  if (!outcome %in% colnames(data)) {
    warning(paste0("Outcome '", outcome, "' not found in dataframe."))
    return(NULL)
  }
  
  # Specify covariate strings
  pcs_string <- paste0("PC_", 1:40)
  geno_batch <- "geno_batch"
  
  cat("Building formula...\n")
  # Set formula
  formula <- reformulate(
    c(paste0("relevel(gs_quint, ref = '", ref_var, "')"),
    paste0("PC_", 1:40), "geno_batch"),
    response = outcome)
  

  if (cat == "bin") {
    # Check for complete separation risk
    n_cases <- sum(data[[outcome]] == 1, na.rm = TRUE)
    cat(sprintf("  Overall: %d cases out of %d (%.2f%%)\n", 
               n_cases, nrow(data), 100*n_cases/nrow(data)))
    
    if (n_cases < 5) {
      warning("Too few cases (<5) for reliable logistic regression")
      return(NULL)
    }
    
    cat("Running logistic regression model...\n")
    # Run model
    model <- glm(formula, 
                data = data, 
                family = binomial)
    
    cat("Extracting model coefficients...\n")
    # Extract model summary 
    model_summary <- summary(model)
    coef_index <- 2:5
    quintile <- paste0(
      "Q",str_sub(row.names(model_summary$coefficients)[coef_index], -1)
      )
    
    if (length(coef_index) == 0) {
      warning("Exposure variable not found in model coefficients")
      return(NULL)
    }
    
    cat("Calculating effect estimates and confidence intervals...\n")
    # Calculate statistics manually
    estimate <- model_summary$coefficients[coef_index, "Estimate"]
    std_error <- model_summary$coefficients[coef_index, "Std. Error"]
    p_value <- model_summary$coefficients[coef_index, "Pr(>|z|)"]
    
    # Calculate CI for log odds
    ci_lower <- estimate - 1.96 * std_error
    ci_upper <- estimate + 1.96 * std_error
    
    # Calculate OR and CI for OR
    or <- exp(estimate)
    or_ci_lower <- exp(ci_lower)
    or_ci_upper <- exp(ci_upper)
    
    cat("Creating results dataframe...\n")
    # Create result data frame
    result <- data.frame(
      quintile = quintile,
      outcome = outcome,
      b = estimate,
      se = std_error,
      OR = or,
      p_value = p_value,
      CI_lower = ci_lower,
      CI_upper = ci_upper
    )
    
    cat("Results summary:\n")
    for(i in 1:length(quintile)) {
      cat(sprintf("  %s: OR = %.3f (%.3f-%.3f), p = %.4f\n", 
                 quintile[i], or[i], or_ci_lower[i], or_ci_upper[i], p_value[i]))
    }
    
  } else if (cat == "cont") {
    cat("Running linear regression model...\n")
    model <- lm(formula, data = data)
    
    cat("Extracting model coefficients...\n")
    # Extract model summary 
    model_summary <- summary(model)
    coef_index <- 2:5
    quintile <- paste0(
      "Q",str_sub(row.names(model_summary$coefficients)[coef_index], -1)
      )
    
    if (length(coef_index) == 0) {
      warning("Exposure variable not found in model coefficients")
      return(NULL)
    }
    
    cat("Calculating effect estimates and confidence intervals...\n")
    # Calculate statistics manually
    estimate <- model_summary$coefficients[coef_index, "Estimate"]
    std_error <- model_summary$coefficients[coef_index, "Std. Error"]
    p_value <- model_summary$coefficients[coef_index, "Pr(>|t|)"]
    
    # Calculate confidence intervals
    ci_lower <- estimate - 1.96 * std_error
    ci_upper <- estimate + 1.96 * std_error
    
    cat("Creating results dataframe...\n")
    # Create result data frame 
    result <- data.frame(
      quintile = quintile,
      outcome = outcome,
      b = estimate,
      se = std_error,
      OR = NA,
      p_value = p_value,
      CI_lower = ci_lower,
      CI_upper = ci_upper
    )
    
    cat("Results summary:\n")
    for(i in 1:length(quintile)) {
      cat(sprintf("  %s: Beta = %.3f (%.3f-%.3f), p = %.4f\n", 
                 quintile[i], estimate[i], ci_lower[i], ci_upper[i], p_value[i]))
    }
    
  } else {
    return(NULL)
  }
  
  cat("=== Association test complete ===\n\n")
  return(result)
}


################################################################################
# Plotting function for GS results

plot_strat <- function(dat, adjustment){

  cont <- dat %>%
    filter(cat == "cont") %>%
    forestplot(
      df = .,
      name = outcome_full,
      colour = .data[[adjustment]],
      estimate = b,
      logodds = FALSE,
      xlab = "Mean difference in trait (95% CI) relative to Q5",
      pvalue = FDR
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 10, "mm")
    ) + ggtitle("Continuous outcomes") +
    scale_color_manual(values = colours)
  
  bin <- dat %>%
    filter(cat == "bin") %>%
    forestplot(
      df = .,
      name = outcome_full,
      estimate = b,
      colour = .data[[adjustment]],
      logodds = TRUE,
      xlab = "Odds ratio for APPO (95% CI) relative to Q5",
      pvalue = FDR
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      plot.margin = margin(5, 5, 5, 20, "mm")
    ) + ggtitle("Binary outcomes") +
    scale_color_manual(values = colours)
  
  cowplot::plot_grid(
    bin, cont, 
    ncol = 1, 
    rel_heights = c(6.5, 1),
    align = "v",
    labels = NULL
  ) %>%
    cowplot::plot_grid(., ncol = 1, rel_heights = c(0.1, 1))
}

################################################################################
# Function to create a violin plot for lp(a) levels by genotype

create_violin_plot <- function(data, variant) {
  # Filter out NA values and convert genotype to factor with proper labels
  plot_data <- data %>%
    filter(!is.na(.data[[variant]]), !is.na(lpa_mes)) %>%
    mutate(
      genotype = factor(
        .data[[variant]],
        levels = c("0", "1", "2"),
        labels = c("WT", "HET", "HOM")
      )
    )
  
  # Check sample sizes for each genotype
  sample_sizes <- plot_data %>%
    group_by(genotype) %>%
    summarise(n = n(), .groups = "drop")
  
  print(paste("Sample sizes:", paste(sample_sizes$genotype, "=", sample_sizes$n, collapse = ", ")))
  
  # Calculate median values for each genotype
  median_values <- plot_data %>%
    group_by(genotype) %>%
    summarise(
      median_lpa = median(lpa_mes, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create the violin plot with boxplot overlay
  p <- ggplot(plot_data, aes(x = genotype, y = lpa_mes, fill = genotype)) +
    geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +  # Increased alpha and added scale parameter
    geom_boxplot(width = 0.2, fill = "white", alpha = 0.8, outlier.size = 0.3) +  # Reduced boxplot width
    # Add median values as text labels on top
    geom_text(
      data = median_values,
      aes(x = genotype, y = median_lpa, label = sprintf("%.1f", median_lpa)),
      vjust = -0.5,
      size = 3,
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    labs(
      title = paste("Lp(a) levels by", variant, "genotype"),
      x = "Genotype",
      y = "Lp(a) measurement (nmol/L)",
      caption = paste0("n(WT)=", sum(plot_data$genotype == "WT", na.rm = TRUE),
                      ", n(HET)=", sum(plot_data$genotype == "HET", na.rm = TRUE), 
                      ", n(HOM)=", sum(plot_data$genotype == "HOM", na.rm = TRUE))
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "none"
    ) +
    scale_fill_brewer(palette = "Set1")
  
  return(p)
}