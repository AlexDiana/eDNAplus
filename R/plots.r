
diagnosticsPlot <- function(chain_output){
  
  chain_output_all <-
    matrix(chain_output, nrow = nchain, ncol = niter)
  
  chain_output_long <- reshape2::melt(chain_output_all)
  
  ggplot2::ggplot(data = chain_output_long, ggplot2::aes(
    x = Var2,
    y = value,
    group = Var1,
    color = factor(Var1)
  )) + ggplot2::geom_line() +
    ggplot2::xlab("Iterations") + ggplot2::ylab("Value") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 17),
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 11, face = "bold"),
      axis.text.x = ggplot2::element_text(
        size = 11,
        face = "bold",
        hjust = 1
      ),
      axis.line = ggplot2::element_line(colour = "black", size = 0.15),
      # panel.grid.minor = element_line(colour="grey", size=0.15),
      panel.grid.major = ggplot2::element_line(colour = "grey", size = 0.15),
      panel.background = ggplot2::element_rect(fill = "white", color = "black"),
      legend.position = "none"
    )
  
}

plotCoefficients <- function(modelResults, cov_num = 1, species = NULL){
  
  nchain <- MCMCparams$nchain
  nburn <- MCMCparams$nburn
  niter <- MCMCparams$niter
  nthin <- MCMCparams$nthin
  
  OTUnames <- modelResults$data_characteristics$OTUnames
  
  beta_z_output <- modelResults$params_output$beta_z_output
  
  beta_CI <- apply(beta_z_output, c(2,3), function(x){
    quantile(x, probs = c(.05,.5,.95))
  })
  
  significantSpecies <- 
    which(beta_CI[3,cov_num,] < 0 |
            beta_CI[1,cov_num,] > 0)
  
  orderSignificantSpecies <- 
    significantSpecies[order(beta_CI[2,cov_num,significantSpecies])]
  
  if(!is.null(species))  {
    orderSignificantSpecies <- speciesToAnalyze  
  }
  
  beta_CI_subset <- beta_CI[,cov_num,orderSignificantSpecies]
  colnames(beta_CI_subset) <- OTUnames[orderSignificantSpecies]
  
  subsetSpecies <- 1:length(orderSignificantSpecies)
  
  factorSub <- factor(subsetSpecies, levels = subsetSpecies)
  
  namesSpecies <- OTUnames[orderSignificantSpecies]
  
  ggplot(data = NULL, aes(x = factorSub,
                          y = beta_CI_subset[2,],
                          ymin = beta_CI_subset[1,],
                          ymax = beta_CI_subset[3,])) + geom_errorbar() + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                   axis.title = ggplot2::element_text(size = 20, face = "bold"),
                   axis.text = ggplot2::element_text(size = 9, face = "bold", angle = 0),
                   panel.grid.major = ggplot2::element_line(colour="grey", size=0.15),
                   panel.background = ggplot2::element_rect(fill = "white", color = "black")) +
    geom_hline(aes(yintercept = 0), color = "red") + 
    scale_x_discrete(name = "Species", breaks = subsetSpecies,
                     labels = namesSpecies) +
    # scale_x_continuous(breaks = subsetSpecies, name = "Species",
    # labels = colnames(OTU)[subsetSpecies]) +
    scale_y_continuous(name = "Elevation") + coord_flip()
  
}