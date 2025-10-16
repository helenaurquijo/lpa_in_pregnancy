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
  # Minimize printing for speed
  cat(paste("Analyzing:", outcome, "|", cat, "\n"))
  
  # Check if outcome exists in dataframe
  if (!outcome %in% colnames(df)) {
    warning(paste0("Outcome '", outcome, "' not found in dataframe."))
    return(NULL)
  }
  
  # Pre-calculate PCs string once
  pcs_string <- paste0("PC_", 1:40)
  
  if (is.null(quint_var)) {
    # Overall analysis (not stratified by quintiles)
    formula <- reformulate(c(score_var, pcs_string), response = outcome)
    
    # Try-catch to prevent crashes
    tryCatch({
      if (cat == "bin") {
        # Check for complete separation risk
        n_cases <- sum(df[[outcome]] == 1, na.rm = TRUE)
        cat(sprintf("  Overall: %d cases out of %d (%.2f%%)\n", 
                   n_cases, nrow(df), 100*n_cases/nrow(df)))
        
        if (n_cases < 5) {
          warning("Too few cases (<5) for reliable logistic regression")
          return(NULL)
        }
        
        # Use control parameters to help with convergence issues
        model <- glm(formula, 
                    data = df, 
                    family = binomial,
                    control = list(maxit = 50, trace = FALSE))
        
        # Extract model summary without using tidy()
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
        
        # Create result data frame
        result <- data.frame(
          outcome = outcome,
          quintile = "overall",
          b = estimate,  # Log odds for binary outcomes
          se = std_error, # Standard error
          OR = or,       # Odds ratio
          p_value = p_value,
          CI_lower = ci_lower,
          CI_upper = ci_upper
        )
        
      } else if (cat == "cont") {
        model <- lm(formula, data = df)
        
        # Extract model summary without using tidy()
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
          b = estimate,  # Beta for continuous outcomes
          se = std_error, # Standard error
          OR = NA,       # NA for continuous outcomes
          p_value = p_value,
          CI_lower = ci_lower,
          CI_upper = ci_upper
        )
      } else {
        return(NULL)
      }
      
      # Print a single concise summary
      effect_type <- if(cat == "bin") "OR" else "Beta"
      effect_value <- if(cat == "bin") result$OR else result$b
      cat(sprintf("  Overall %s = %.3f (%.3f-%.3f), p = %.4f\n", 
                 effect_type, effect_value, result$CI_lower, result$CI_upper, result$p_value))
      
      return(result)
    }, error = function(e) {
      warning(paste("Error in overall analysis for", outcome, ":", e$message))
      return(NULL)
    })
  } else {
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
      
      # Create formula with outcome and score variable
      formula <- as.formula(paste(outcome, "~", score_var))
      
      # Try-catch to prevent crashes during model fitting
      tryCatch({
        if (cat == "bin") {
          # Check for complete separation risk
          n_cases <- sum(subset_df[[outcome]] == 1, na.rm = TRUE)
          cat(sprintf("    Q%d: %d cases out of %d (%.2f%%)\n", 
                     q, n_cases, nrow(subset_df), 100*n_cases/nrow(subset_df)))
          
          if (n_cases < 5) {
            cat("    Skipping - too few cases (<5) for reliable logistic regression\n")
            next
          }
          
          # Use more robust fitting parameters and handle warnings
          model <- suppressWarnings(
            glm(formula, 
                data = subset_df, 
                family = binomial,
                control = list(maxit = 100, trace = FALSE))
          )
          
          # Check for convergence
          if (!model$converged) {
            cat("    Warning: Model did not converge - results may be unreliable\n")
          }
          
          # Extract model summary without using tidy()
          model_summary <- summary(model)
          coef_index <- which(rownames(model_summary$coefficients) == score_var)
          
          if (length(coef_index) == 0) {
            cat("    Warning: Score variable not found in model coefficients\n")
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
          
          # Extract model summary without using tidy()
          model_summary <- summary(model)
          coef_index <- which(rownames(model_summary$coefficients) == score_var)
          
          if (length(coef_index) == 0) {
            cat("    Warning: Score variable not found in model coefficients\n")
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
      }, error = function(e) {
        cat(sprintf("    Error in quintile %d: %s\n", q, e$message))
      })
    }
    
    # Check if any results were collected
    if (length(results_list) == 0) {
      warning(paste("No valid results for", outcome))
      return(NULL)
    }
    
    # Combine all results into a single data frame
    combined_results <- do.call(rbind, results_list)
    
    # Verify and report the number of quintiles analyzed
    quintiles_found <- unique(combined_results$quintile)
    cat(sprintf("\n  Final results for %s: Analyzed %d quintiles: %s\n", 
               outcome, length(quintiles_found), 
               paste(sort(quintiles_found), collapse=", ")))
    
    return(combined_results)
  }
}