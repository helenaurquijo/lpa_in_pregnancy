################################################################################
#                       Functions for result plotting                          #
################################################################################

################################################################################
# Function to compare two sets of SNP effects
scatter_plot <- function(dat){
  lim <- c(-20, 170)
  
  # Calculate correlation (Pearson)
  cor_val <- cor(dat$beta.exposure, dat$beta.outcome, use = "complete.obs")
  cor_label <- paste0("r = ", sprintf("%.2f", cor_val))
  
  # Calculate confidence intervals
  dat$lSE.exposure <- dat$beta.exposure - 1.96 * dat$se.exposure
  dat$uSE.exposure <- dat$beta.exposure + 1.96 * dat$se.exposure
  dat$lSE.outcome <- dat$beta.outcome - 1.96 * dat$se.outcome
  dat$uSE.outcome <- dat$beta.outcome + 1.96 * dat$se.outcome
  
  p <- ggplot(dat, aes(x = beta.exposure, y = beta.outcome)) +
    geom_point(color = "#E69F00", alpha = 0.75) +
    geom_errorbar(aes(xmin = lSE.exposure, xmax = uSE.exposure), alpha = 0.25) +
    geom_errorbar(aes(ymin = lSE.outcome, ymax = uSE.outcome), alpha = 0.25) +
    geom_abline(color = "gray60") +
    geom_hline(yintercept = 0, color = "#FF000033") +
    geom_vline(xintercept = 0, color = "#FF000033")
    
  # Only add regression if more than one point
  if(nrow(dat) > 1) {
    p <- p + geom_smooth(method = "lm", color = "cornflowerblue", alpha = 0.25, fullrange = TRUE)
  }
  
  # Add correlation annotation to lower right
  p + coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = FALSE) +
    labs(title = "SNP association with Lp(a)",
         x = "SNP effect in Sinnott-Armstrong et al. 2021 (GWAS estimates)",
         y = "SNP effect in Shi et al. 2024 (IV estimates)") +
    theme_bw() +
    annotate("text", x = lim[2] - 0.01, y = lim[1] + 0.01, 
             label = cor_label, hjust = 1, vjust = 0, size = 4)
}


################################################################################
# Function to create and save forest plots for MR results

plot_mr <- function(dat, scaled = FALSE){

  # Set estimate and se column names and x-axis labels based on scaled parameter
  estimate_col <- if (scaled) "b_scaled" else "b"
  se_col <- if (scaled) "se_scaled" else "se"
  
  if (scaled) {
    xlab_cont <- "Mean different in trait (95% CI) per 210 nmol/L decrease in lipoprotein(a)"
    xlab_bin <- "Odds ratio for APPO risk (95% CI) per 210 nmol/L decrease in lipoprotein(a)"
  } else {
    xlab_cont <- "Mean different in trait (95% CI) per 1 nmol/L decrease in lipoprotein(a)"
    xlab_bin <- "Odds ratio for APPO risk (95% CI) per 1 nmol/L decrease in lipoprotein(a)"
  }

  cont <- dat %>%
    filter(cat == "cont") %>%
    forestplot(
      df = .,
      name = outcome_full,
      estimate = !!sym(estimate_col),
      se = !!sym(se_col),
      logodds = FALSE,
      xlab = xlab_cont,
      pvalue = FDR
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      legend.position = "none"
    ) + ggtitle("Continuous outcomes")
  
  bin <- dat %>%
    filter(cat == "bin") %>%
    forestplot(
      df = .,
      name = outcome_full,
      estimate = !!sym(estimate_col),
      se = !!sym(se_col),
      logodds =  TRUE,
      xlab = xlab_bin,
      pvalue = FDR
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank()
    ) + ggtitle("Binary outcomes")
  
  patchwork::wrap_plots(bin, cont, ncol = 1, heights = c(8, 1)) &
    theme(plot.title = element_text(hjust = 0)
  )
}

################################################################################
# Function to create and save forest plots for MR results with different labels

plot_mr_special <- function(dat, adjustment, scaled = FALSE){

  # Set estimate and se column names and x-axis labels based on scaled parameter
  estimate_col <- if (scaled) "b_scaled" else "b"
  se_col <- if (scaled) "se_scaled" else "se"
  
  if (scaled) {
    xlab_cont <- "Mean different in trait (95% CI) per 210 nmol/L decrease in lipoprotein(a)"
    xlab_bin <- "Odds ratio for APPO risk (95% CI) per 210 nmol/L decrease in lipoprotein(a)"
  } else {
    xlab_cont <- "Mean different in trait (95% CI) per 1 nmol/L decrease in lipoprotein(a)"
    xlab_bin <- "Odds ratio for APPO risk (95% CI) per 1 nmol/L decrease in lipoprotein(a)"
  }

  my_scale <- colours[1:length(unique(dat[[adjustment]]))]
  names(my_scale) <- unique(dat[[adjustment]])

  cont <- dat %>%
    filter(cat == "cont") %>%
    forestplot(
      df = .,
      name = outcome_full,
      colour = .data[[adjustment]],
      estimate = !!sym(estimate_col),
      se = !!sym(se_col),
      logodds = FALSE,
      xlab = xlab_cont,
      pvalue = FDR
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      legend.position = "none"
    ) + ggtitle("Continuous outcomes") +
    scale_color_manual(values = my_scale)
  
  bin <- dat %>%
    filter(cat == "bin") %>%
    forestplot(
      df = .,
      name = outcome_full,
      estimate = !!sym(estimate_col),
      se = !!sym(se_col),
      colour = .data[[adjustment]],
      logodds =  TRUE,
      xlab = xlab_bin,
      pvalue = FDR
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank()
    ) + ggtitle("Binary outcomes") +
    scale_color_manual(values = my_scale)
  
  patchwork::wrap_plots(bin, cont, ncol = 1, heights = c(8, 1)) &
    theme(plot.title = element_text(hjust = 0)
  )
}