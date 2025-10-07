################################################################################
#               Extract SNP-outcome summary data for selected SNPs             #
################################################################################

##
# This script extracts SNP-outcome summary data for a list of SNPs from
# the main IV SNPs as well as its proxies in case they are required.
##


## Set up
readRenviron("scripts/config.env")
mrpregdir <- Sys.getenv("mrpregdir")

# Load modules
library(dplyr)
library(data.table)
library(purrr)

# Import list of rsids
snps_ls <- read.table("./scripts/rsids.txt") %>% pull
str(snps_ls)

# Outcome information
out_info <- read.table("./scripts/vars_metanalysis_v5_mum.csv", header=T, sep = ",") %>%
            select(Outcome:var_GWAS)

# Input data directory
dat_dir <- paste0(mrpregdir, "/results/GWAS/")


##
## Function to extract & format unadjusted outcome data

extr_dat <- function(study, prefix=NULL, suffix=NULL) {
  
 # List relevant outcomes
  if (study == "Metanalysis") {
    vars <- select(out_info, Variable_name)  
  } else if (study == "BIB-SA") {
    vars <- select(out_info, paste0("var_", "BIB.SA")) 
  } else if (study == "BIB-WE") {
    vars <- select(out_info, paste0("var_", "BIB.WE")) 
  } else if (study == "Public") {
    vars <- select(out_info, paste0("var_", "GWAS")) 
  } else {
    vars <- select(out_info, paste0("var_", toupper(study)))  
  }   
  vars <- vars %>% 
    na.omit %>% 
	pull 
  print(vars)
  
 # List pathway to files
  if (study %in% c("Metanalysis", "ALSPAC", "UKB","MOBA")) {
    path <- paste0(dat_dir, study, "/mothers/final/")
  } else {
    path <- paste0(dat_dir, study, "/mothers/")
  }
  print(path)
  
 # List files
  files <- map(vars, ~list.files(path = path, pattern = paste0(prefix, ., suffix), full.names = F)) %>% unlist
  print(files)

 for(i in 1:length(files)) {

   # File path
     f_path <- list.files(path = path, pattern = paste0("^", files[i], "$"), full.names = T)
	 print(f_path)
		
   # New file name
	new_file <- paste0("./data/", "tmp_", study, "_", files[i])
		
   # Extract header
	system(paste0("zcat ", f_path, " | awk NR==1 > ", new_file))
    
   # Extract summary data for list of rsIDs
	system(paste0("zcat ", f_path, " | grep -w -F -f ", "./scripts/rsids.txt", " >> ", new_file))
	# -w tells grep to match whole words only
	# -F search for fixed strings (plain text) rather than regular expressions
	# -f "data/rsids.txt" read search patterns from the file
 }
	
  # Format data
	sub_files <- map(files, ~paste0("tmp_", study, "_", .))
    dat <- map(sub_files, ~fread(paste0("./data/", .))) %>%
      map2(., vars, ~mutate(.x, Phenotype = .y, study = study, chr = as.numeric(chr), effect_allele = as.character(effect_allele), other_allele = as.character(other_allele), pval = as.numeric(pval))) %>%
      map(., ~filter(., SNP %in% snps_ls)) %>%
	  map(., ~mutate(., SNP = as.character(SNP))) %>%
      bind_rows 
	
# Harmonise outcome names of each study with metanalysis
	if(study %in% c("ALSPAC", "FinnGen", "MOBA", "UKB")) {
	  var_name <- paste0("var_", toupper(study))
	} else if (study == "BIB-SA") {
      var_name <- paste0("var_", "BIB.SA") 
	} else if (study == "BIB-WE") {
	  var_name <- paste0("var_", "BIB.WE")    
    } else if (study=="Public") {
	  var_name <- "var_GWAS"
	}
	if (! study == "Metanalysis") {
	dat <- select(out_info, Variable_name, all_of(var_name)) %>%
	  na.omit %>%
	  merge(., dat, by.x = var_name, by.y = "Phenotype") %>%
	  rename(Phenotype = Variable_name) %>%
	  select(-var_name)	
    }
	print(head(dat))
	return(dat)
}
 
   
##
## Extract & format outcome data   
    
## Metanalysis
ma_out_dat <- extr_dat(study = "Metanalysis", prefix = "metanalyses-R3.mum.", suffix = ".txt.gz") 

## Individual studies (remove overlapping phenotypes)

# ALSPAC
out_dat1 <- extr_dat(study = "ALSPAC", prefix = "alspac_mums_") 

# BIB-SA
out_dat2 <- extr_dat(study = "BIB-SA", prefix = "BIB_mother.southasian_") 

# BIB-WE
out_dat3 <- extr_dat(study = "BIB-WE", prefix = "BIB_mother.european_") 

# FINNGEN
out_dat4 <- extr_dat(study = "FinnGen", prefix = "FINNGEN-R12.", suffix = ".20250203.txt.gz") %>%
			rename(., eaf_case = eaf_cases, eaf_control = eaf_controls)

# MOBA
out_dat5 <- extr_dat(study = "MOBA", prefix = "moba_100k_mum_") 

# UKB
out_dat6 <- extr_dat(study = "UKB", prefix = "ukb_mums_", suffix = "_20241210.txt.gz")

# Public GWAS
out_dat7 <- extr_dat(study = "Public") 

## Combined study-specific data
stu_out_dat <- bind_rows(out_dat1, out_dat2, out_dat3, out_dat4, out_dat5, out_dat6, out_dat7)


##
## Save outcome data
write.table(ma_out_dat, paste0("./data/", "ma_out_dat.txt"), quote = F, row.names = F)
write.table(stu_out_dat, paste0("./data/", "stu_out_dat.txt"), quote = F, row.names = F)


##
## Function to extract & format adjusted outcome data (WLM estimates)

extr_wlm_dat <- function(type) {
  
  # List relevant outcomes   
    vars <- filter(out_info, !Variable_name %in% c("misc_subsamp", "s_misc_subsamp", "r_misc_subsamp", "sb_subsamp", "vpretb_all", "nvp_sev_all"))
	vars <- select(vars, Variable_name) %>% na.omit %>% pull 
	#vars <- vars[c(5:37)] # preg loss excluded from WLM
    print(vars)
  
  # List pathway to files
    path <- paste0(mrpregdir, "/results/adjusted-wlm")
	print(path)
  
  # List files 
    files <- map(vars, ~list.files(path = path, pattern = paste0("^", ., "_donuts_", type, ".txt.gz"), full.names = F)) %>% unlist
   	print(files)
	
	for(i in 1:length(files)) {

		# File path
		f_path <- list.files(path = path, pattern = paste0("^", files[i], "$"), full.names = T)
		print(f_path)
		
		# New file name
		new_file <- paste0("./data/", "tmp_", files[i])
	
		# Extract header
		system(paste0("zcat ", f_path, " | awk NR==1 > ", new_file))
    
		# Extract summary data for list of rsIDs
		system(paste0("zcat ", f_path, " | grep -w -F -f ", "./scripts/rsids.txt", " >> ", new_file))
		# -w tells grep to match whole words only
		# -F search for fixed strings (plain text) rather than regular expressions
		# -f "data/rsids.txt" read search patterns from the file
	}
	
  # Format data
   	sub_files <- map(files, ~paste0("./data/tmp_", .))
    dat <- map(sub_files, ~fread(.)) 
      
    dat <- dat %>%
         map2(., vars, ~mutate(.x, Phenotype = .y, type = type)) %>%
		 map(., ~filter(., SNP %in% snps_ls)) %>%
		 bind_rows 
	
	print(head(dat))
    return(dat)
    
}

# WLM estimates (duos)
duos_out_dat <- extr_wlm_dat(type = "duos_mum-child") %>% 
				mutate(analyses = "duos_mum-child", study = "Metanalysis")

write.table(duos_out_dat, paste0("./data/", "duos_out_dat.txt"), quote = F, row.names = F)
			
# WLM estimates (trios)
trios_out_dat <- extr_wlm_dat(type = "trios") %>% 
				 mutate(analyses = "trios", study = "Metanalysis")

write.table(trios_out_dat, paste0("./data/", "trios_out_dat.txt"), quote = F, row.names = F)

q("no")
