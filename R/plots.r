
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
  
  beta_CI <- apply(beta_z_output, c(3,4), function(x){
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

FPNP_plot <- function(modelResults, idxSpecies){
  
  y_col <- modelResults$data_infos$y[,,idxSpecies]
  infos <- modelResults$data_infos$infos
  
  delta_mean <- modelResults$params_output$delta_mean[,,idxSpecies,drop=F]
  delta_col <- apply(delta_mean, 2, mean)
  gamma_mean <- modelResults$params_output$gamma_mean[,,idxSpecies,drop=F]
  gamma_col <- apply(gamma_mean, 2, mean)
  theta11_mean <- modelResults$params_output$theta11_mean[,,idxSpecies,drop=F]
  theta11_col <- c(apply(theta11_mean, 2, mean), rep(0, length(delta_mean) - length(theta11_mean)))
  c_imk1_mean <- modelResults$params_output$c_imk1_mean[,,,idxSpecies,drop=F]
  c_imk1_col <- apply(c_imk1_mean, c(2,3), mean)
  c_imk2_mean <- modelResults$params_output$c_imk2_mean[,,,idxSpecies,drop=F]
  c_imk2_col <- apply(c_imk2_mean, c(2,3), mean)
  
  colnames(y_col) <- c("PCR1","PCR2","PCR2")
  theta11_col <- data.frame("PresenceProb" = theta11_col)
  delta_col <- data.frame("MeanPresence" = delta_col)
  gamma_col <- data.frame("MeanFP" = gamma_col)
  colnames(c_imk1_col) <- c("PCR1_state","PCR2_state","PCR2_state")
  colnames(c_imk2_col) <- c("PCR1_state","PCR2_state","PCR2_state")
  
  cbind(infos[,c(1,2)],
        y_col,
        theta11_col,
        delta_col,
        gamma_col,
        c_imk1_col,
        c_imk2_col
  )
  
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


plotCorrelationMatrix <- function(modelResults, 
                                  idxSpecies = NULL,
                                  numSigSpecies = 0){
  
  jointSpecies <- modelResults$data_infos$jointSpecies
  
  if(jointSpecies){
    
    Tau_output <- modelResults$params_output$Tau_output
    Tau_output <- apply(Tau_output, c(3,4), c)
    
    niter <- dim(Tau_output)[1]
    S <- dim(Tau_output)[2]
    
    # Tau_corr_output_list <- lapply(1:niter, function(i){
    #   stats::cov2cor(Tau_output[i,,])
    # })
    # 
    # Tau_corr_output <- array(NA, dim = c(niter, S, S))
    # for (iter in 1:niter) {
    #   Tau_corr_output[iter,,] <- Tau_corr_output_list[[iter]]
    # }
    
    Tau_CI <- apply(Tau_output, c(2,3), function(x){
      quantile(x, probs = c(0.025, 0.975, .5))
    })
    
    # cor_Tau <- Tau_mean
    # 
    # for (i in 1:S) {
    #   for (j in seq_len(i)) {
    #     cor_Tau[i,j] <- Tau_mean[i,j] / (sqrt(Tau_mean[i,i]) * sqrt(Tau_mean[j,j]))
    #     cor_Tau[j,i] <- cor_Tau[i,j]
    #   }
    # }
    # 
    # dimnames(cor_Tau) <- list(colnames(OTU),
    # colnames(OTU))
    
    if(!missing(speciesSubset)){
      
      Tau_CI_long <- matrix(NA, S * (S - 1) / 2, 5)
      l <- 1
      for (i in 1:S) {
        for (j in seq_len(i - 1)) {
          Tau_CI_long[l,1] <- Tau_CI[1,i,j]
          Tau_CI_long[l,2] <- Tau_CI[2,i,j]
          Tau_CI_long[l,3] <- Tau_CI[3,i,j]
          Tau_CI_long[l,4] <- i
          Tau_CI_long[l,5] <- j
          l <- l + 1
        }
      }
      
      significantCorrelations <- 
        which((Tau_CI_long[,1] < 0 & Tau_CI_long[,2] < 0) | 
                (Tau_CI_long[,1] > 0 & Tau_CI_long[,2] > 0) )
      
      if(length(significantCorrelations) == 0){
        print("No significant correlations")
      } else {
        
        Tau_CI_long_significants <- Tau_CI_long[significantCorrelations,]
        
        orderSignificantCorrelations <- order(-abs(ifelse(Tau_CI_long_significants[,3] < 0,
                                            Tau_CI_long_significants[,2],
                                            Tau_CI_long_significants[,1])))
        Tau_CI_long_ordered <- 
          Tau_CI_long_significants[order(Tau_CI_long_significants[,3]),]
        
        biggestSpecies_idx <- which((cor_Tau_long_ordered[,3] < 0 & cor_Tau_long_ordered[,2] < -0.5) | 
                                      (cor_Tau_long_ordered[,3] > 0 & cor_Tau_long_ordered[,1] > 0.5))
        
        biggestSpecies <- cor_Tau_long_ordered[biggestSpecies_idx,4:5]
        idxSpecies <- unique(as.vector(biggestSpecies))
        
      }
      
      
    } 
    
    Tau_CI_median <- Tau_CI[3,,]
    Taucorr_CI_median <- stats::cov2cor(Tau_CI_median)
    cov_to_plot <- Tau_CI[3, idxSpecies, idxSpecies]
    
    
    corrplot::corrplot(corr_to_plot)  
    
  } else {
    
    print("Species correlation not selected when running the model")
    
  }
  
}