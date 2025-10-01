################################################################################
#                       Functions for result plotting                          #
################################################################################

################################################################################
# Function to create and save forest plots for MR results

plot_mr <- function(dat){

  cont <- dat %>%
    filter(cat == "cont") %>%
    forestplot(
      df = .,
      name = outcome_full,
      estimate = b,
      logodds = FALSE,
      xlab = "Mean different in trait (95% CI) per 10 nmol/L decrease in lipoprotein(a)",
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
      estimate = b,
      logodds =  TRUE,
      xlab = "Odds ratio for APPO risk (95% CI) per 10 nmol/L decrease in lipoprotein(a)",
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

plot_mr_special <- function(dat, adjustment){

  my_scale <- colours[1:length(unique(dat[[adjustment]]))]
  names(my_scale) <- unique(dat[[adjustment]])

  cont <- dat %>%
    filter(cat == "cont") %>%
    forestplot(
      df = .,
      name = outcome_full,
      colour = .data[[adjustment]],
      estimate = b,
      logodds = FALSE,
      xlab = "Mean different in trait (95% CI) from decrease in lipoprotein(a)",
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
      estimate = b,
      colour = .data[[adjustment]],
      logodds =  TRUE,
      xlab = "Odds ratio for APPO risk (95% CI) from decrease in lipoprotein(a)",
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
