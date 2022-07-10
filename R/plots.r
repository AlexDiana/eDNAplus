
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


#' Biomass covariates coefficients plot
#' 
#' @description Plot the covariate coefficients of the biomasses amount
#' 
#' @param modelResults Output of the function \code{runEDNA}
#' @param cov_num Index of the covariate to plot
#' @param idxSpecies Optional indexes of the species to plot. If missing, the significant species are plotted.
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return A credible intervals plot with the species on the columns
#' 
#' @examples
#' 
#' plotBiomassCoefficients(sampleResults)
#' 
plotBiomassCoefficients <- function(modelResults, cov_num = 1, idxSpecies){
  
  nchain <- MCMCparams$nchain
  nburn <- MCMCparams$nburn
  niter <- MCMCparams$niter
  nthin <- MCMCparams$nthin
  
  OTUnames <- modelResults$data_infos$OTUnames
  covName <- modelResults$data_infos$namesCovZ[cov_num]
  
  beta_z_output <- modelResults$params_output$beta_z_output
  
  beta_CI <- apply(beta_z_output, c(3,4), function(x){
    quantile(x, probs = c(.05,.5,.95))
  })
  
  beta_CI <- beta_CI[,cov_num,]
  
  significantSpecies <- 
    which(beta_CI[3,] < 0 |
            beta_CI[1,] > 0)
  
  orderSignificantSpecies <- 
    significantSpecies[order(beta_CI[2,significantSpecies])]
  
  if(!missing(idxSpecies))  {
    orderSignificantSpecies <- idxSpecies  
  }
  
  beta_CI_subset <- beta_CI[,orderSignificantSpecies,drop=F]
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
    ggtitle(covName) +
    ggplot2::scale_y_continuous(name = expression(beta)) + ggplot2::coord_flip()
  
}

#' Probability plot of true positive
#' 
#' @description Plot the credible intervals of the true positive probability
#' 
#' @param modelResults Output of the function \code{runEDNA}
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return Credible intervals plot of the true positive probabilities
#' 
#' @examples
#' 
#' plotTruePositiveProbability(sampleResults)
#' 
plotTruePositiveProbability <- function(modelResults, idxSpecies){
  
  nchain <- MCMCparams$nchain
  nburn <- MCMCparams$nburn
  niter <- MCMCparams$niter
  nthin <- MCMCparams$nthin
  
  OTUnames <- modelResults$data_infos$OTUnames
  
  p_11_output <- modelResults$params_output$p_11_output
  
  p11_CI <- apply(p_11_output, c(3), function(x){
    quantile(x, probs = c(.05,.5,.95))
  })
  
  colnames(p11_CI) <- OTUnames
  
  if(!missing(idxSpecies))  {
    speciesSubset <- idxSpecies  
  } else {
    speciesSubset <- 1:dim(p_11_output)[3]
  }
  
  p11_CI_subset <- p11_CI[,speciesSubset,drop=F]
  colnames(p11_CI_subset) <- OTUnames[speciesSubset]
  
  factorSub <- factor(speciesSubset, levels = speciesSubset)
  
  namesSpecies <- OTUnames[speciesSubset]
  
  ggplot2::ggplot() + geom_point(data = NULL, ggplot2::aes(x = factorSub,
                                            y = p11_CI_subset[2,])) + 
    geom_errorbar(data = NULL, ggplot2::aes(x = factorSub,
                                            ymin = p11_CI_subset[1,],
                                            ymax = p11_CI_subset[3,])) + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                   axis.title = ggplot2::element_text(size = 20, face = "bold"),
                   axis.text = ggplot2::element_text(size = 9, face = "bold", angle = 0),
                   panel.grid.major = ggplot2::element_line(colour="grey", size=0.15),
                   panel.background = ggplot2::element_rect(fill = "white", color = "black")) +
    ggplot2::scale_x_discrete(name = "Species", breaks = speciesSubset,
                              labels = namesSpecies) +
    # scale_x_continuous(breaks = subsetSpecies, name = "Species",
    # labels = colnames(OTU)[subsetSpecies]) +
    ggplot2::scale_y_continuous(name = "PCR true positive rate") + ggplot2::coord_flip()
  
}

#' Probability plot of false positive
#' 
#' @description Plot the credible intervals of the false positive probabilities
#' 
#' @param modelResults Output of the function \code{runEDNA}
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return Credible intervals plot of the false positive probabilities
#' 
#' @examples
#' 
#' plotTruePositiveProbability(sampleResults)
#' 
plotFalsePositiveProbability <- function(modelResults, idxSpecies){
  
  nchain <- MCMCparams$nchain
  nburn <- MCMCparams$nburn
  niter <- MCMCparams$niter
  nthin <- MCMCparams$nthin
  
  OTUnames <- modelResults$data_infos$OTUnames
  
  p_10_output <- modelResults$params_output$p_10_output
  
  p11_CI <- apply(p_10_output, c(3), function(x){
    quantile(x, probs = c(.05,.5,.95))
  })
  
  colnames(p11_CI) <- OTUnames
  
  if(!missing(idxSpecies))  {
    speciesSubset <- idxSpecies  
  } else {
    speciesSubset <- 1:dim(p_10_output)[3]
  }
  
  p11_CI_subset <- p11_CI[,speciesSubset,drop=F]
  colnames(p11_CI_subset) <- OTUnames[speciesSubset]
  
  factorSub <- factor(speciesSubset, levels = speciesSubset)
  
  namesSpecies <- OTUnames[speciesSubset]
  
  ggplot2::ggplot() + geom_point(data = NULL, ggplot2::aes(x = factorSub,
                                            y = p11_CI_subset[2,])) + 
    geom_errorbar(data = NULL, ggplot2::aes(x = factorSub,
                                            ymin = p11_CI_subset[1,],
                                            ymax = p11_CI_subset[3,])) + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                   axis.title = ggplot2::element_text(size = 20, face = "bold"),
                   axis.text = ggplot2::element_text(size = 9, face = "bold", angle = 0),
                   panel.grid.major = ggplot2::element_line(colour="grey", size=0.15),
                   panel.background = ggplot2::element_rect(fill = "white", color = "black")) +
    ggplot2::scale_x_discrete(name = "Species", breaks = speciesSubset,
                              labels = namesSpecies) +
    # scale_x_continuous(breaks = subsetSpecies, name = "Species",
    # labels = colnames(OTU)[subsetSpecies]) +
    ggplot2::scale_y_continuous(name = "PCR false positive rate") + ggplot2::coord_flip()
  
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

#' Biomasses distribution maps
#' 
#' @description Plot the maps of the biomasses
#' 
#' @param modelResults Output of the function \code{runEDNA}
#' @param datapoly Polygon file of the map
#' @param idxSpecies Index of the species to plot
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return A plot with the map of the biomasses, where a darker color signals higher biomass
#' 
#' @examples
#' 
#' plotBiomassCoefficients(sampleResults)
#' 
plotBiomassesMap <- function(modelResults, datapoly, idxSpecies){
  
  if(length(idxSpecies) > 1){
    
    print("You can only plot one species at a time")
    
  } else {
    
    logz_output <- modelResults$params_output$logz_output
    logz_mean <- apply(logz_output, c(3, 4), mean)
    logz_species <- logz_mean[,idxSpecies]
    # logz_species <- (logz_species - min(logz_species)) / (max(logz_species) - min(logz_species))
    
    datapoly2 <- datapoly
    datapoly2$LogBiomass <- rep(logz_species, each = 4)
    
    ggplot() +
      geom_polygon(data = datapoly2, aes(x = x, y = y, group = id, 
                                        alpha = LogBiomass), 
                   color = "black") +
      scale_alpha_continuous(range = c(0, 1)) +
      xlab("X") + ylab("Y") +
      # geom_point(data = NULL, aes(x = siteLocations[,1], y = siteLocations[,2])) +
      scale_x_continuous(breaks = c()) +
      scale_y_continuous(breaks = c()) +
      theme(plot.title = element_text(hjust = 0.5, size = 20),
            # legend.position = "none",
            axis.title = element_text(size = 20, face = "bold"),
            axis.text = element_text(size = 13, face = "bold", angle = 90),
            # panel.grid.major = element_line(colour="grey", size=0.015),
            panel.background = element_rect(fill = "white", color = "black")) 
      
  }
  
}

#' Species correlation matrix
#' 
#' @description Plot the species correlation matrix
#' 
#' @param modelResults Output of the function \code{runEDNA}
#' @param idxSpecies Optional indexes of the species to plot. If missing, the significant correlations are plotted.
#' @param numSigCorrs Number of significant correlation to plot
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return A plot of the median correlation matrix for the species selected
#' 
#' @examples
#' 
#' plotCorrelationMatrix(sampleResults)
#' 
plotCorrelationMatrix <- function(modelResults, 
                                  idxSpecies,
                                  numSigCorrs = 0){
  
  jointSpecies <- modelResults$data_infos$jointSpecies
  
  if(jointSpecies){
    
    Tau_output <- modelResults$params_output$Tau_output
    Tau_output <- apply(Tau_output, c(3,4), c)
    
    niter <- dim(Tau_output)[1]
    S <- dim(Tau_output)[2]

    Tau_CI <- apply(Tau_output, c(2,3), function(x){
      quantile(x, probs = c(0.025, 0.975, .5))
    })
    
    if(missing(idxSpecies)){
      
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
        
        if(missing(idxSpecies)){
          print("No significant correlations")
          opt <- options(show.error.messages = FALSE)
          on.exit(options(opt))
          stop()  
        }
        
      } else {
        
        Tau_CI_long_significants <- Tau_CI_long[significantCorrelations,]
        
        orderSignificantCorrelations <- order(-abs(ifelse(Tau_CI_long_significants[,3] < 0,
                                            Tau_CI_long_significants[,2],
                                            Tau_CI_long_significants[,1])))
        Tau_CI_long_ordered <- 
          Tau_CI_long_significants[orderSignificantCorrelations,]
        
        if(numSigCorrs > 0){
          numSigCorrs <- min(numSigCorrs, nrow(Tau_CI_long_ordered))
          idxSpecies <- unique(as.vector(Tau_CI_long_ordered[1:numSigCorrs,4:5])  )
        } else {
          idxSpecies <- unique(as.vector(Tau_CI_long_ordered[,4:5]))
        }
        
      }
      
    } 
    
    Tau_CI_median <- Tau_CI[3,,]
    Taucorr_CI_median <- stats::cov2cor(Tau_CI_median)
    corr_to_plot <- Taucorr_CI_median[idxSpecies, idxSpecies]
    
    OTUnames <- modelResults$data_infos$OTUnames
    
    rownames(corr_to_plot) <- OTUnames[idxSpecies]
    colnames(corr_to_plot) <- OTUnames[idxSpecies]
    
    return(ggcorrplot::ggcorrplot(corr_to_plot))
    
  } else {
    
    print("Species correlation not selected when running the model")
    
  }
  
}

predictNewSites <- function(){
  
  logz_output <- modelResults$params_output$logz_output
  logz_output <- apply(logz_output, c(3,4), c)
  
  if (spatialCorr) {
    # x_min <- min(X_s[,1]) - .01
    # x_max <- max(X_s[,1]) + .01
    # y_min <- min(X_s[,2]) - .01
    # y_max <- max(X_s[,2]) + .01
    # 
    # X_s_star <- as.matrix(expand.grid(seq(x_min, x_max, length.out = 20),
    #                         seq(y_min, y_max, length.out = 20)))
    
    n_star <- nrow(X_s_star)
    
    Sigma_XstarXstar <- K2(X_s_star, X_s_star, 1, l_gp)
    Sigma_XsXstar <- K2(X_s, X_s_star, 1, l_gp)
    Sigma_XstarXs <- K2(X_s_star, X_s, 1, l_gp)
    
    Sigma_star_cond <- Sigma_XstarXstar - Sigma_XstarXs %*% solve(Sigma_n) %*% Sigma_XsXstar
    chol_Sigma_star_cond <- t(chol(Sigma_star_cond))
    
    logz_star_output <- array(NA, dim = c(nchain, niter, n_star, S))
  } else {
    logz_star_output <- NULL
  }
  
  niter <- 
  
    for (iter in 1:niter) {
    Xbeta <- matrix(beta0, n, S, byrow = F) + X_z %*% beta_z
    Xbeta_star <- matrix(beta0, n_star, S, byrow = F) + X_z_star %*% beta_z
    
    mu1 <- Xbeta_star + Sigma_XstarXs %*% solve(Sigma_n) %*% (logz - Xbeta)
    Sigma1 <- Sigma_XstarXstar - Sigma_XstarXs %*% solve(Sigma_n) %*% Sigma_XsXstar
    
    if(jointSpecies){
      chol_Tau <- t(chol(Tau_params$Sigma))
    } else {
      chol_Tau <- diag(tau, nrow = S)
    }
  }
  
  
}