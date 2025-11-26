################################################################################
#                         Functions for coloc analyses                         #
################################################################################

################################################################################
# Function to format data for coloc
#  with conditional handling of bin/cont phenotypes

format_coloc <- function(lipid, appo) {

  lipid_df <- filter(exp_dat, exposure == lipid) %>%
    as.data.frame()

  appo_df <- outcome_dat[[appo]]
  
  df <- harmonise_data(lipid_df, appo_df, action = 1) %>% # keeping palindromic snps  
    # remove missing data
    filter(
      !is.na(beta.exposure) & 
      !is.na(beta.outcome) &
      !is.na(se.exposure) &
      !is.na(se.outcome)
    ) 
  
  if(nrow(df)<1) {return(NULL)}
    
  # incompatible snps still removed
    
  # Calculate MAF.exposure only if it doesn't already exist
  if (!"MAF.exposure" %in% colnames(df)) {
    df$MAF.exposure <- NULL
    for (j in 1:nrow(df)){
      if(df$eaf.exposure[j]>0.5)
      {df$MAF.exposure[j]=1-df$eaf.exposure[j]}
      else 
      {df$MAF.exposure[j]=df$eaf.exposure[j]}
    }
  }

  # Calculate MAF.outcome only if it doesn't already exist
  if (!"MAF.outcome" %in% colnames(df)) {
    df$MAF.outcome <- NULL
    for (j in 1:nrow(df)){
      if(df$eaf.outcome[j]>0.5)
      {df$MAF.outcome[j]=1-df$eaf.outcome[j]}
      else 
      {df$MAF.outcome[j]=df$eaf.outcome[j]}
    }
  }

  ## Conditionally add more info to dataset2 if the outcome is binary
  if (out_info$cat[out_info$outcome == appo] == "bin") {
    df$s <- df$ncase.outcome / (df$ncase.outcome + df$ncontrol.outcome)
  }
  
  return(df)
}

################################################################################
# Function to run coloc

run_coloc <- function(beta1,beta2,se1,se2,MAF1,MAF2,N1,N2,snps,s) {
  
  #For testing function
  # dataset1 <- list(beta=df$beta.exposure, varbeta=df$se.exposure^2, MAF=df$MAF.exposure,type="quant", N=N1, snp = df$SNP)
  # dataset2 <- list(beta=df$beta.outcome, varbeta=df$se.outcome^2,MAF=df$MAF.outcome, type="quant",N=N2, snp = df$SNP)
  
  # Convert the inputs in order to run in coloc function.
  dataset1 <- list(
    beta = beta1,
    varbeta = se1^2,
    MAF = MAF1,
    type = "quant",
    N = N1,
    snp = snps
  )
  
  dataset2 <- list(
    beta = beta2, 
    varbeta = se2^2,
    MAF = MAF2, 
    type = ifelse(is.null(s), "quant", "cc"),
    N = N2,
    snp = snps
  )

  ## Conditionally add more to dataset2 if s was provided (binary outcome)
  if (!is.null(s)) {
    dataset2$s <- s
  }
  
  result <- coloc.abf(dataset1, dataset2, p1=1e-4, p2=1e-4, p12=1e-5)  
  
  # List into data frame
  df <- data.frame(matrix(unlist(result$summary), nrow=1, byrow=T))
  
  # Label the columns in the data frame
  names(df) <- c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
  
  snps <- result$results %>% arrange(desc(SNP.PP.H4))
  
  res <- list(df, snps)
  
  return(res)
}

################################################################################
# Final function with output ready to save

prep_and_run_coloc <- function(lipid, appo) {
  
  # Format data
  f <- format_coloc(
    lipid,
    appo
  )
  
  if(is.null(f)) {return(NULL)}

  # Determine outcome type
  outcome_cat <- out_info$cat[out_info$outcome == appo]
  
  # Run coloc

  res <- tryCatch(
    {run_coloc(
      beta1 = f$beta.exposure,
      beta2 = f$beta.outcome,
      se1 = f$se.exposure,
      se2 = f$se.outcome,
      MAF1 = f$MAF.exposure,
      MAF2 = f$MAF.outcome,
      N1 = median(f$samplesize.exposure),
      N2 = median(f$samplesize.outcome),
      snps = f$SNP,
      s = if (outcome_cat == "bin") median(f$s) else NULL
    )},
    error = function(e) {
      message("Error in run_coloc ", e$message)
      return(NULL)
    })

  if(is.null(res)) {return(NULL)}

  res[[1]]$trait1 <- lipid
  res[[1]]$trait2 <- appo
  res[[1]]$lead_snp <- res[[2]]$snp[1]
  
  return(res[[1]])
}

################################################################################
# Function to create assoc plots

run_assocplot <- function(df) {

  ## Gene
  g <- "LPA"

  ## Biomarker name
  exp <- unique(df$exposure)

  ## Outcome name
  out <- unique(df$outcome)
  out_name <- out_info$outcome_full[out_info$outcome == out]
      
   print(paste("Doing:", g, "-", exp, "-", out))
	
  ## LD matrix (r)
	
  ### With alleles
  ld <- ld_matrix_local(
	  variants = df$SNP,
	  bfile = paste0(Sys.getenv("mrpregdir"), "/data/output/UKB_10k/ukb_eur_10k"),
	  plink_bin = plinkbinr::get_plink_exe(),
    with_alleles = T
    )

  ### Without alleles
  ld_snps <- do.call(rbind, strsplit(rownames(ld), split = "_"))[, 1]
		
  ## Format LD matrix and summary data

  ### Remove SNPs from summary data if not present in LD matrix
  df <- filter(df, SNP %in% ld_snps)

  ### Arrange SNPs in the same order between LD matrix and summary data
  ord <- match(ld_snps, df$SNP)
  ld <- ld[order(ord), order(ord)]
  ld_snps <- ld_snps[order(ord)]
  print(paste(" The order of SNPs is the same:", identical(ld_snps, df$SNP)))

  ## Format LD matrix to have same effect allele as the summary dataset
  ld_harmo <- harmonise_ld_dat(df, ld)

  ## Extract SNP info
  markers <- select(ld_harmo$x, SNP, chr.exposure, pos.exposure) %>%
    mutate(pos.exposure = as.integer(pos.exposure)) %>%
    rename(marker = SNP, chr = chr.exposure, pos = pos.exposure) 

  ## Extract z-statistics for SNP-exp and SNP-out associations
  assoc <- ld_harmo$x %>%
    rename(marker = SNP, chr = chr.exposure, pos = pos.exposure,
	  pvalue_1 = pval.exposure, pvalue_2 = pval.outcome) %>%
    select(marker, chr, pos, pvalue_1, pvalue_2) %>%
    distinct	 

  ## Select top SNP for Exposure
  topsnp <- mutate(ld_harmo$x, z_abs = abs(beta.exposure / se.exposure)) %>%
    slice_max(order_by = z_abs, n = 1, with_ties = F) %>%
    pull(SNP)

  ## Generate plot
  p <- geni.plots::fig_region_stack(
    data = assoc,
    traits = c(exp, out_name),
    corr = ld_harmo$ld,
    build = 37,
    top_marker = topsnp,
    highlights = topsnp,
    title_center = TRUE
    )

  fname <- paste0(projectdir, "/results/sensitivity_analyses/assoc_plots/stack_assoc_plot_", g, "_X_", exp, "_X_", out, ".png")

  ggsave(filename = fname, plot = p, width = 18, height = 15, units = "cm")
}