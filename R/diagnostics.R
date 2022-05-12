tracePlotParameters <- function(modelResults, param, idx = NULL) {
  
  eval(parse(text = paste0("param_output <- modelResults$params_output$",param,"_output")))
  
  nchain <- dim(param_output)[1]
  niter <- dim(param_output)[2]
  
  if(!is.null(idx)){
    if(length(idx) == 1){
      param_output <- param_output[,,idx,drop=F]
      param_output <-
        matrix(param_output[, , 1], nrow = nchain, ncol = niter)
    } else {
      param_output <- param_output[,,idx[1],idx[2],drop=F]
      param_output <-
        matrix(param_output[, , 1, 1], nrow = nchain, ncol = niter)
    }
  }
  
  
  param_output_long <- reshape2::melt(param_output)
  
  diagnosticsPlot <-
    ggplot2::ggplot(data = param_output_long, ggplot2::aes(
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
  
  
  plotTitle <- createPlotTitle(param_output, nchain)
  
  diagnosticsPlot <- diagnosticsPlot + ggplot2::ggtitle(plotTitle)
  
  return(diagnosticsPlot)
  
}

createPlotTitle <- function(mcmc_output, nchain) {
  eff_samplesize <- ess(mcmc_output)
  
  plotTitle <-
    paste0("Effective sample size = ", round(eff_samplesize, 1))
  
  if (nchain > 1) {
    Rhat <- compute_rhat(mcmc_output)
    
    plotTitle <-
      paste0(plotTitle, paste0(" / R.hat = ", round(Rhat, 3)))
  }
  
  plotTitle
}

compute_rhat <- function(mcmc_output) {
  mcmc_output_list <- lapply(1:nrow(mcmc_output), function(i) {
    coda::mcmc(mcmc_output[i, ])
  })
  mcmc_output_list_2 <- coda::as.mcmc.list(mcmc_output_list)
  
  Rhat <- coda::gelman.diag(mcmc_output_list_2)
  
  Rhat$psrf[2]
}

ess <- function(mcmc_output) {
  mcmc_output_list <- lapply(1:nrow(mcmc_output), function(i) {
    coda::mcmc(mcmc_output[i, ])
  })
  mcmc_output_list_2 <- coda::as.mcmc.list(mcmc_output_list)
  
  coda::effectiveSize(mcmc_output_list_2)
  
}
