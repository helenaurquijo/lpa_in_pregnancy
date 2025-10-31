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
# Function to perform association analysis - with updated output format

assoc_test <- function(outcome, cat, df, score_var, quint_var = NULL) {
  # Determine if score_var is binary
  is_binary_exposure <- all(df[[score_var]] %in% c(0, 1, NA), na.rm = TRUE)
  
  cat(paste("Analysing:", outcome, "| Type:", cat, 
            "| Exposure type:", ifelse(is_binary_exposure, "binary", "continuous"), "\n"))
  
  # Check if outcome exists in dataframe
  if (!outcome %in% colnames(df)) {
    warning(paste0("Outcome '", outcome, "' not found in dataframe."))
    return(NULL)
  }
  
  # Specify covariate strings
  pcs_string <- paste0("PC_", 1:40)
  geno_batch <- "geno_batch"
  
  if (is.null(quint_var) || is_binary_exposure) {
    # Overall analysis (not stratified by quintiles)
    # If exposure is binary, this is automatically overall only
    formula <- reformulate(c(score_var, geno_batch, pcs_string), response = outcome)
    
    if (cat == "bin") {
      # Check for complete separation risk
      n_cases <- sum(df[[outcome]] == 1, na.rm = TRUE)
      cat(sprintf("  Overall: %d cases out of %d (%.2f%%)\n", 
                 n_cases, nrow(df), 100*n_cases/nrow(df)))
      
      if (n_cases < 5) {
        warning("Too few cases (<5) for reliable logistic regression")
        return(NULL)
      }
      
      # Run model
      model <- glm(formula, 
                  data = df, 
                  family = binomial)
      
      # Extract model summary 
      model_summary <- summary(model)
      coef_index <- which(rownames(model_summary$coefficients) == score_var)
      
      if (length(coef_index) == 0) {
        warning("Score variable not found in model coefficients")
        return(NULL)
      }
      
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
      
      # Create result data frame with binary exposure flag
      result <- data.frame(
        outcome = outcome,
        quintile = "overall",
        b = estimate,
        se = std_error,
        OR = or,
        p_value = p_value,
        CI_lower = ci_lower,
        CI_upper = ci_upper
      )
      
    } else if (cat == "cont") {
      model <- lm(formula, data = df)
      
      # Extract model summary 
      model_summary <- summary(model)
      coef_index <- which(rownames(model_summary$coefficients) == score_var)
      
      if (length(coef_index) == 0) {
        warning("Score variable not found in model coefficients")
        return(NULL)
      }
      
      # Calculate statistics manually
      estimate <- model_summary$coefficients[coef_index, "Estimate"]
      std_error <- model_summary$coefficients[coef_index, "Std. Error"]
      p_value <- model_summary$coefficients[coef_index, "Pr(>|t|)"]
      
      # Calculate confidence intervals
      ci_lower <- estimate - 1.96 * std_error
      ci_upper <- estimate + 1.96 * std_error
      
      # Create result data frame 
      result <- data.frame(
        outcome = outcome,
        quintile = "overall",
        b = estimate,
        se = std_error,
        OR = NA,
        p_value = p_value,
        CI_lower = ci_lower,
        CI_upper = ci_upper
      )
    } else {
      return(NULL)
    }
    
    # Print appropriate effect description based on exposure type
    effect_type <- if(cat == "bin") "OR" else "Beta"
    effect_desc <- if(is_binary_exposure) "comparing groups" else "per SD increase"
    effect_value <- if(cat == "bin") result$OR else result$b
    
    cat(sprintf("  Overall %s = %.3f (%.3f-%.3f), p = %.4f %s\n", 
               effect_type, effect_value, result$CI_lower, result$CI_upper, 
               result$p_value, effect_desc))
    
    return(result)
  } else if (!is_binary_exposure) {
    # Only run quintile analysis for continuous exposures
    # Quintile-stratified analysis
    # Print diagnostic information about the quintiles
    quintile_counts <- table(df[[quint_var]])
    cat(sprintf("  Quintile distribution: %s\n", 
               paste(quintile_counts, collapse = ", ")))
    
    # Check if we have binary outcome and show event counts per quintile
    if (cat == "bin") {
      event_counts <- tapply(df[[outcome]], df[[quint_var]], sum, na.rm = TRUE)
      cat(sprintf("  Event counts per quintile: %s\n", 
                 paste(event_counts, collapse = ", ")))
    }
    
    # Initialize results storage
    results_list <- list()
    
    # Process each quintile using a simple for loop to ensure each is processed
    for (q in 1:5) {
      cat(sprintf("\n  Processing quintile %d...\n", q))
      
      # Subset data for this quintile
      subset_df <- df[which(df[[quint_var]] == q), ]
      
      cat(sprintf("    Subset size: %d records\n", nrow(subset_df)))
      
      # Skip if sample size is too small
      if (nrow(subset_df) < 10) {
        cat("    Skipping - sample size too small\n")
        next
      }
      
      # Specify formula
      formula <- reformulate(c(score_var, geno_batch, pcs_string), response = outcome)
      
      if (cat == "bin") {
        # Check for complete separation risk
        n_cases <- sum(subset_df[[outcome]] == 1, na.rm = TRUE)
        cat(sprintf("    Q%d: %d cases out of %d (%.2f%%)\n", 
                   q, n_cases, nrow(subset_df), 100*n_cases/nrow(subset_df)))
        
        if (n_cases < 5) {
          warning(sprintf("Too few cases (%d) for reliable logistic regression in quintile %d", n_cases, q))
          next
        }
        
        # Run model
        model <- suppressWarnings(
          glm(formula, 
              data = subset_df, 
              family = binomial)
        )
        
        # Check for convergence
        if (!model$converged) {
          warning(sprintf("Model for quintile %d did not converge - results may be unreliable", q))
        }
        
        # Extract model summary 
        model_summary <- summary(model)
        coef_index <- which(rownames(model_summary$coefficients) == score_var)
        
        if (length(coef_index) == 0) {
          warning(sprintf("Score variable not found in model coefficients for quintile %d", q))
          next
        }
        
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
        
        # Create result data frame
        this_result <- data.frame(
          outcome = outcome,
          quintile = q,
          b = estimate,  # Log odds for binary outcomes
          se = std_error, # Standard error
          OR = or,       # Odds ratio
          p_value = p_value,
          CI_lower = ci_lower,
          CI_upper = ci_upper,
          stringsAsFactors = FALSE
        )
        
        # Add to results list
        results_list[[length(results_list) + 1]] <- this_result
        
        # Print concise summary for this quintile
        cat(sprintf("    Q%d OR = %.3f (%.3f-%.3f), p = %.4f\n", 
                   q, or, or_ci_lower, or_ci_upper, p_value))
        
      } else if (cat == "cont") {
        model <- lm(formula, data = subset_df)
        
        # Extract model summary 
        model_summary <- summary(model)
        coef_index <- which(rownames(model_summary$coefficients) == score_var)
        
        if (length(coef_index) == 0) {
          warning(sprintf("Score variable not found in model coefficients for quintile %d", q))
          next
        }
        
        # Calculate statistics manually
        estimate <- model_summary$coefficients[coef_index, "Estimate"]
        std_error <- model_summary$coefficients[coef_index, "Std. Error"]
        p_value <- model_summary$coefficients[coef_index, "Pr(>|t|)"]
        
        # Calculate confidence intervals
        ci_lower <- estimate - 1.96 * std_error
        ci_upper <- estimate + 1.96 * std_error
        
        # Create result data frame
        this_result <- data.frame(
          outcome = outcome,
          quintile = q,
          b = estimate,  # Beta for continuous outcomes
          se = std_error, # Standard error
          OR = NA,       # NA for continuous outcomes
          p_value = p_value,
          CI_lower = ci_lower,
          CI_upper = ci_upper,
          stringsAsFactors = FALSE
        )
        
        # Add to results list
        results_list[[length(results_list) + 1]] <- this_result
        
        # Print concise summary for this quintile
        cat(sprintf("    Q%d Beta = %.3f (%.3f-%.3f), p = %.4f\n", 
                   q, estimate, ci_lower, ci_upper, p_value))
      }
    }
    
    # Check if any results were collected
    if (length(results_list) == 0) {
      warning(paste("No valid results for", outcome))
      return(NULL)
    }
    
    # Combine all results into a single data frame
    combined_results <- do.call(rbind, results_list)
    
    # Verify and report the number of quintiles analysed
    quintiles_found <- unique(combined_results$quintile)
    cat(sprintf("\n  Final results for %s: Analysed %d quintiles: %s\n", 
               outcome, length(quintiles_found), 
               paste(sort(quintiles_found), collapse=", ")))
    
    return(combined_results)
  } else {
    # Add a default return for any case not caught above
    warning(paste("Analysis not performed for", outcome, "- no suitable analysis branch selected"))
    return(NULL)
  }
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
      xlab = "Mean difference in trait (95% CI) per SD decrease in Lp(a) genetic score",
      pvalue = FDR
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      legend.position = "none"
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
      xlab = "Odds ratio for APPO (95% CI) per SD decrease in Lp(a) genetic score",
      pvalue = FDR
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank()
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