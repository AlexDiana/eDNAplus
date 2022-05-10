
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
  
  ggplot2::ggplot(data = NULL, ggplot2::aes(x = factorSub,
                          y = beta_CI_subset[2,],
                          ymin = beta_CI_subset[1,],
                          ymax = beta_CI_subset[3,])) + geom_errorbar() + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                   axis.title = ggplot2::element_text(size = 20, face = "bold"),
                   axis.text = ggplot2::element_text(size = 9, face = "bold", angle = 0),
                   panel.grid.major = ggplot2::element_line(colour="grey", size=0.15),
                   panel.background = ggplot2::element_rect(fill = "white", color = "black")) +
    ggplot2::geom_hline(aes(yintercept = 0), color = "red") + 
    ggplot2::scale_x_discrete(name = "Species", breaks = subsetSpecies,
                     labels = namesSpecies) +
    # scale_x_continuous(breaks = subsetSpecies, name = "Species",
    # labels = colnames(OTU)[subsetSpecies]) +
    ggplot2::scale_y_continuous(name = "Elevation") + ggplot2::coord_flip()
  
}

plotBiomasses <- function(modelResults, datapoly, idxSpecies){
  
  logz_output <- modelResults$params_output$logz_output
  logz_mean <- apply(logz_output, c(2, 3), mean)
  logz_species <- logz_mean[,idxSpecies]
  logz_species <- (logz_species - min(logz_species)) / (max(logz_species) - min(logz_species))
  
  ggplot() +
    geom_polygon(data = datapoly, aes(x = x, y = y, group = id), alpha =  rep(logz_species, each = 4),
                 fill = "grey", color = "black") +
    scale_alpha_continuous(range = c(0, 1)) +
    xlab("X") + ylab("Y") +
    # geom_point(data = NULL, aes(x = siteLocations[,1], y = siteLocations[,2])) +
    scale_x_continuous(breaks = c()) +
    scale_y_continuous(breaks = c()) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.position = "none",
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 13, face = "bold", angle = 90),
          # panel.grid.major = element_line(colour="grey", size=0.015),
          panel.background = element_rect(fill = "white", color = "black")) 
  
  
  
}


