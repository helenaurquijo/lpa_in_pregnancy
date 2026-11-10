################################################################################
#                         Functions for MR analyses                            #
################################################################################

################################################################################
# Function to find and replace missing SNPs with tags

##
# This script is based off of code by Lizzy Aiton which is on Github: 
# eaiton/autoimmune-pregnancy/blob/main/functions/generate_outdat_with_proxies.R
##

replace_missing_with_proxies <- function(out, exposure_dat, outcome_dat, proxies, outcome_col = "outcome") {
  
  print(paste("Processing outcome:", out))

  # Detect all beta columns in outcome_dat
  beta_cols <- grep("^beta_", colnames(outcome_dat), value = TRUE)
  
  # Determine which columns to use for allele checking
  if (length(beta_cols) == 0) {
    # No beta columns found - use default structure
    check_cols <- "effect_allele.outcome"
    effect_col <- "effect_allele.outcome"
    other_col <- "other_allele.outcome"
  } else {
    # Beta columns exist - use beta columns for missing data check
    check_cols <- beta_cols
    
    # Check if standard allele columns exist (no suffix)
    if ("effect_allele" %in% colnames(outcome_dat)) {
      effect_col <- "effect_allele"
      other_col <- "other_allele"
    } else {
      # Fall back to .outcome suffix
      effect_col <- "effect_allele.outcome"
      other_col <- "other_allele.outcome"
    }
  }

  # Filter specific outcome dataset to available original exposure SNPs
  outcome_specific_data <- outcome_dat %>%
    filter(SNP %in% exposure_dat$SNP, .data[[outcome_col]] == out)
  
  print(paste("Number of original SNPs available for outcome", out, ":", nrow(outcome_specific_data)))

  # Remove any duplicated SNPs
  outcome_specific_data <- outcome_specific_data[!duplicated(outcome_specific_data$SNP), ]

  # Check whether outcome is available
  if (nrow(outcome_specific_data) == 0) {
    print(paste("Outcome", out, "not available"))
    return(outcome_specific_data)
  }

  # Check whether proxies are needed - check if all SNPs present for ALL beta columns
  all_snps_present <- TRUE
  for (check_col in check_cols) {
    if (sum(!is.na(outcome_specific_data[[check_col]])) < nrow(exposure_dat)) {
      all_snps_present <- FALSE
      break
    }
  }
  
  if (all_snps_present && nrow(outcome_specific_data) == nrow(exposure_dat)) {
    print(paste("No proxies needed for", out))
    return(outcome_specific_data)
  } else {
    
    # Identify proxies
    print(paste("Identifying proxies for outcome:", out))

    # Specify all available SNPs
    all_snps <- expand.grid(
      SNP = exposure_dat$SNP, temp_outcome = out
    )
    names(all_snps)[2] <- outcome_col

    # Left join SNP-out data for that outcome, NAs if missing
    snps_with_availability <- all_snps %>%
      left_join(outcome_specific_data, by = c("SNP", outcome_col))
    
    # Mark as missing if ANY of the beta columns is NA
    snps_with_availability$missing <- FALSE
    for (check_col in check_cols) {
      snps_with_availability$missing <- snps_with_availability$missing | 
        is.na(snps_with_availability[[check_col]])
    }

    # Deduplicated list of all SNPs which need a proxy
    snps_needing_proxies <- snps_with_availability %>%
      filter(missing == TRUE) %>%
      distinct(SNP)

    # Count of proxies needed 
    print("SNPs missing:")
    print(table(snps_with_availability$missing))

    # Proxies without query SNPs as the proxy
    proxies_without_query <- filter(proxies, query_rsid != rsid) %>%
      mutate(SNP = rsid)

    # List of all proxies
    proxy_snps <- unique(proxies_without_query$rsid)

    # Unfiltered SNP-out data
    unfiltered_out_data <- outcome_dat %>%
      filter(.data[[outcome_col]] == out)
    
    print(dim(unfiltered_out_data))

    # Stop if no proxies are not present in the outcome dataset at all
    if (all(!unfiltered_out_data$SNP %in% proxy_snps)) {
      print("No proxies among SNP-out data")
      return(outcome_specific_data)
    } else {
      print(paste0("Searching for proxies for ", nrow(snps_needing_proxies), " SNPs."))

      # Merge proxies with SNPs that need proxies
      proxies_with_info <- merge(x = snps_needing_proxies, y = proxies_without_query, by.x = "SNP", by.y = "query_rsid", all.x = TRUE) %>% rename(query_rsid = SNP.y)

      # Stop if none of the proxies that are needed are in the outcome dataset
      if (all(!proxies_with_info$rsid %in% unfiltered_out_data$SNP)) {
        print("No SNP-outcome associations available for any the proxies that are needed.")
        return(outcome_specific_data)
      }

      # Append outcome data for each proxy SNP
      proxy_outcome_data <- merge(x = proxies_with_info, y = unfiltered_out_data, by.x = "rsid", by.y = "SNP")

      # Exclude palindromic proxies - use the single effect/other allele columns
      filtered_proxies <- proxy_outcome_data %>%
        filter(!(.data[[effect_col]] == "A" & .data[[other_col]] == "T" |
                .data[[effect_col]] == "T" & .data[[other_col]] == "A")) %>%
        filter(!(.data[[effect_col]] == "C" & .data[[other_col]] == "G" |
                .data[[effect_col]] == "G" & .data[[other_col]] == "C"))

      # Select the proxy SNP in highest LD with the query SNP
      set.seed(9)
      selected_proxies <- filtered_proxies %>%
        group_by(SNP) %>%
        drop_na(R2) %>%
        slice_max(order_by = R2, with_ties = TRUE) %>%
        slice_sample(n = 1) %>%
        ungroup()
    }

    # Stop searching for proxies if no SNP-outcome data available
    if ((nrow(snps_needing_proxies) > 0) & (nrow(selected_proxies) == 0)) {
      print("No proxy-outcome associations available for any SNPs.")
    } else {
      
      # Harmonise proxies - use single effect/other allele columns
      harmonised_proxies <- selected_proxies
      
      # Add column specifying proxy effect allele
      harmonised_proxies$effect_allele.proxy <- apply(harmonised_proxies, 1, function(row) {
        if (row[effect_col] == substr(row["Correlated_Alleles"], 3, 3)) {
          return(substr(row["Correlated_Alleles"], 1, 1))
        } else if (row[effect_col] == substr(row["Correlated_Alleles"], 7, 7)) {
          return(substr(row["Correlated_Alleles"], 5, 5))
        } else {
          return("default_value")
        }
      })

      # Same again but for other allele
      harmonised_proxies$other_allele.proxy <- apply(harmonised_proxies, 1, function(row) {
        if (row[other_col] == substr(row["Correlated_Alleles"], 3, 3)) {
          return(substr(row["Correlated_Alleles"], 1, 1))
        } else if (row[other_col] == substr(row["Correlated_Alleles"], 7, 7)) {
          return(substr(row["Correlated_Alleles"], 5, 5))
        } else {
          return("default_value")
        }
      })

      # Treat proxy alleles as the outcome alleles
      final_harmonised_proxies <- harmonised_proxies
      final_harmonised_proxies[[effect_col]] <- final_harmonised_proxies$effect_allele.proxy
      final_harmonised_proxies[[other_col]] <- final_harmonised_proxies$other_allele.proxy
      final_harmonised_proxies$data_source.outcome <- NA
      final_harmonised_proxies$proxied <- TRUE
      final_harmonised_proxies$proxy <- final_harmonised_proxies$rsid

      # Create final dataframe which has replaced missing SNPs with their proxies
      outcome_dat_cols <- colnames(outcome_dat)
      final_harmonised_proxies <- final_harmonised_proxies %>% select(all_of(outcome_dat_cols))

      # Return message with number of SNPs which could be proxied
      print(paste0(nrow(final_harmonised_proxies), " proxies identified."))

      outcome_specific_data <- rbind(outcome_specific_data, final_harmonised_proxies)
    }

    # Return full outcome dataset with added proxies
    return(outcome_specific_data)
  }
}

################################################################################ 
# Function to harmonise LD matrix against summary data 

## Adapted from harmonise_ld_dat function in TwoSampleMR R package

harmonise_ld_dat_me<-function (x, ld) {
  
  snpnames <- do.call(rbind, strsplit(rownames(ld), split = "_")) #snpnames comes from LD matrix

  common_snps <- intersect(snpnames[, 1], x$SNP) # this creates a vector of SNPs that are in both the ld mat and the summary data, in the order they appear in the ld matrix

  ld <- ld[snpnames[, 1] %in% common_snps, snpnames[, 1] %in% common_snps, drop = FALSE] # include only SNPs in the ld matrix that are also in the summary data 

  snpnames <- snpnames[snpnames[, 1] %in% common_snps, , drop = FALSE] # include only SNPs in the snpnames df that are also in the summary data

  x <- x[x$SNP %in% common_snps, , drop = FALSE] # include only SNPs in the summary data that are also in the ld matrix

  x <- x[match(snpnames[, 1], x$SNP), , drop = FALSE] # this reorders the summary data df to match the order of the SNPs in the ld matrix, which is important for the next step

  stopifnot(all(snpnames[, 1] == x$SNP)) # this is the original check in the official function which ensures ld and x have the same SNPs in the SAME ORDER

  x$effect_allele.exposure <- as.character(x$effect_allele.exposure)
  x$other_allele.exposure <- as.character(x$other_allele.exposure)

  snpnames <- data.frame(snpnames, stringsAsFactors = FALSE)
  snpnames <- merge(subset(x, select = c(SNP, effect_allele.exposure, 
                                         other_allele.exposure)), snpnames, by.x = "SNP", by.y = "X1") # this merges the alleles from the summary data with the snpnames from the ld matrix - SNPs are reordered by alphabetical order of the SNP names

  snpnames <- snpnames[match(x$SNP, snpnames$SNP), ] #this reorders the snpnames df to match the order of the SNPs in x (summary data), which should match the order of ld as per the check above

  # at this point snpnames follows the order of x but not ld

  snpnames$keep <- ( #identifying cases where the alleles in the summary data and the ld matrix match, regardless of order
    snpnames$X2 == snpnames$effect_allele.exposure & 
                    snpnames$X3 == snpnames$other_allele.exposure) | 
    (snpnames$X3 ==  snpnames$effect_allele.exposure & 
                    snpnames$X2 == snpnames$other_allele.exposure)

  if (nrow(x) == 0) {
    message(" - none of the SNPs could be aligned to the LD reference panel")
    return(NULL)
  }
  if (any(!snpnames$keep)) {
    message(" - the following SNPs could not be aligned to the LD reference panel: \n- ", 
            paste(subset(snpnames, !keep)$SNP, collapse = "\n - "))
  }

  snpnames$flip1 <- snpnames$X2 != snpnames$effect_allele.exposure
  x <- subset(x, SNP %in% snpnames$SNP)
  temp1 <- x$effect_allele.exposure[snpnames$flip1]
  temp2 <- x$other_allele.exposure[snpnames$flip1]
  x$beta.exposure[snpnames$flip1] <- x$beta.exposure[snpnames$flip1] * -1
  x$beta.outcome[snpnames$flip1] <- x$beta.outcome[snpnames$flip1] * -1
  x$effect_allele.exposure[snpnames$flip1] <- temp2
  x$other_allele.exposure[snpnames$flip1] <- temp1

  rownames(ld) <- snpnames$SNP #SNPs in snpnames are used to RENAME the rows and columns of the ld matrix, but the order of the SNPs in the ld matrix is not changed, so they are still in their original order
  colnames(ld) <- snpnames$SNP

  if (any(!snpnames$keep)) {
    message("Removing ", sum(!snpnames$keep), " variants due to harmonisation issues")
    ld <- ld[snpnames$keep, snpnames$keep]
    x <- x[snpnames$keep, ]
  }
  if(length(ld)>1)
  {
    return(list(x = x, ld = ld))
  }
  else{
    return(NULL)
  }
}

################################################################################
# Function to use local LD matrix data as opposed to the default 1KG3

## Adapted from dat_to_MRInput function in TwoSampleMR R package

dat_to_MRInput_me <- function(dat, get_correlations = FALSE, pop = "EUR") {
  out <- plyr::dlply(
    dat, c("exposure", "outcome"), function(x) {
      x <- plyr::mutate(x)
      message("Converting:")
      message(" - exposure: ", x$exposure[1])
      message(" - outcome: ", x$outcome[1])
      if (get_correlations) {
        message(" - obtaining LD matrix")
        print(unique(x$SNP))
        ld <- NULL

        try(
          ld <- ld_matrix_local(
            unique(x$SNP),
            bfile = paste0(Sys.getenv("mrpregdir"), "/data/input/1000G/1kg.v3/EUR"),
            plink_bin = plinkbinr::get_plink_exe()
          )
        )
        out <- NULL

        if (sum(is.na(ld) == 0) && length(ld) > 1) {
          try(
            out <- harmonise_ld_dat_me(x, ld)
          )
        }

        if (is.null(out)) {
          return(NULL)
        } else {
          x <- out$x
          ld <- out$ld
          MendelianRandomization::mr_input(
            bx = x$beta.exposure,
            bxse = x$se.exposure,
            by = x$beta.outcome,
            byse = x$se.outcome,
            exposure = x$exposure[1],
            outcome = x$outcome[1],
            snps = x$SNP,
            effect_allele = x$effect_allele.exposure,
            other_allele = x$other_allele.exposure,
            eaf = x$eaf.exposure,
            correlation = ld
          )
        }
      } else {
        # run standard IVW
        MendelianRandomization::mr_input(
          bx = x$beta.exposure,
          bxse = x$se.exposure,
          by = x$beta.outcome,
          byse = x$se.outcome,
          exposure = x$exposure[1],
          outcome = x$outcome[1],
          snps = x$SNP,
          effect_allele = x$effect_allele.exposure,
          other_allele = x$other_allele.exposure,
          eaf = x$eaf.exposure
        )
      }
    }
  )

  return(out)
}

################################################################################
# Function to subset data and run MR

run_mr <- function(stu, out, gene, dat) {

  print(paste(stu, "-", gene, "-", out))

  ## Subset data
  df <- dplyr::filter(
    dat,
    .data$outcome == out &
    .data$exposure == gene
  )

  # For single study analyses
  if(length(unique(df$study))>1) {
    df <- dplyr::filter(df, .data$study == stu)
  }

  head(df)
  print(dim(df))

  # Calculate mean F-stat per analysis
  df$F <- (df$beta.exposure^2) / (df$se.exposure^2)
  F_stat <- mean(df$F)
  print(paste("Mean F-stat:", round(F_stat, 2)))

  # Calculate total variance explained (r2) per analysis
  r2 <- sum(df$info.exposure)

  ## Return null object if no SNPs selected
  if (nrow(df) == 0) {	
    res2 <- NULL
    print("No SNPs selected")

    ## Run MR Wald ratio if N SNP = 1
  } else if (nrow(df) == 1) {
    print("Running Wald ratio")
    res <- mr(df)
    res2 <- 
      mutate(
        res,
        outcome = out,
        gene.exposure = gene,
        study = stu
      )

    ## Run MR IVW accounting for LD if N SNPs > 1
  } else if (nrow(df)>1) {
    # Run generalised IVW if LD matrix could be obtained & harmonised
    print("Running IVW")
    df2 <- dat_to_MRInput_me(df, get_correlation = TRUE)
    if(!is.null(df2[[1]])) {
      res <- MendelianRandomization::mr_ivw(df2[[1]])
      res2 <- data.frame(
        study = stu,
        exposure = res$Exposure,
        outcome =  res$Outcome,
        b =     	 res$Estimate,
        se =       res$StdError,
        pval =     res$Pvalue,
        nsnp =     res$SNPs,
        het_Q =    res$Heter.Stat[1],
        het_pval = res$Heter.Stat[2],
        method =   "Generalised IVW",
        F =        F_stat,
        r2 =       r2
      )
      print("gIVW successful")
      ### Run Wald ratio (1 SNP) or standard IVW (> 1 SNP) if LD matrix could not be obtained or harmonised
    } else {
      print("WARNING: Could not be harmonised - run IVW or Wald ratio")
      #print(df)
      res <- df %>%
        mutate(rsid = SNP, pval = pval.exposure) %>%
        mr(method_list = c("mr_ivw", "mr_wald_ratio"))

      res2 <- res %>%
        mutate(
          outcome = out,
          exposure = gene,
          study = stu
        ) %>% select(-c("id.exposure", "id.outcome"))
    }
  }
  return(res2)
}

################################################################################
# Function to run LOO of study-specific MR estimates

run_loo <- function(stu_res, out) {

  # Subset data (one gene, one exposure, one outcome)
  dt <- stu_res %>% filter(outcome == out)
  idexp <- unique(dt$id.exposure)
  print(dt)

  if(nrow(dt) <= 1) {
    return(NULL)
  } else {

  # Format data for LOO
  study.loo_dat <- metafor::rma(yi = b, sei = se, data=dt, slab = study, method="FE")
  print(study.loo_dat)

  # Run LOO removing one study at a time
  study.loo_res <- metafor::leave1out(study.loo_dat)

  # Format LOO output
  study.loo_res2 <- study.loo_res %>% 
    as.data.frame %>%
	  tibble::rownames_to_column(., var = "study") %>%
    mutate(exposure = gene, outcome = out, se = estimate/zval )

  return(study.loo_res2)

  }
}

################################################################################
# Function to run LOO of study-specific MR estimates

run_loo <- function(gene, out) {

  # Subset data (one gene, one exposure, one outcome)
  dt <- stu_res %>% filter(exposure == gene & outcome == out)
  idexp <- unique(dt$id.exposure)
  print(dt)

  if(nrow(dt) <= 1) {
    return(NULL)
  } else {

  # Format data for LOO
  study.loo_dat <- metafor::rma(yi = b, sei = se, data=dt, slab = study, method="FE")
  print(study.loo_dat)

  # Run LOO removing one study at a time
  study.loo_res <- metafor::leave1out(study.loo_dat)

  # Format LOO output
  study.loo_res2 <- study.loo_res %>% 
    as.data.frame %>%
	  tibble::rownames_to_column(., var = "study") %>%
    mutate(exposure = gene, outcome = out, se = estimate/zval )

  return(study.loo_res2)

  }
}

################################################################################
# Function to calculate OR and confidence intervals for binary outcomes
calculate_or_ci <- function(dat) {
  dat %>%
    mutate(
      # Calculate OR, UCI, and LCI for binary outcomes (1 unit change)
      OR = if_else(cat == "bin", exp(b), NA_real_),
      UCI = if_else(cat == "bin", exp(b + 1.96 * se), NA_real_),
      LCI = if_else(cat == "bin", exp(b - 1.96 * se), NA_real_),
      # Calculate OR, UCI, and LCI for binary outcomes (scaled change)
      OR_scaled = if_else(cat == "bin", exp(b_scaled), NA_real_),
      UCI_scaled = if_else(cat == "bin", exp(b_scaled + 1.96 * se_scaled), NA_real_),
      LCI_scaled = if_else(cat == "bin", exp(b_scaled - 1.96 * se_scaled), NA_real_)
    )
}
