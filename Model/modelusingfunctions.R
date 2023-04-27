# CPP ---------------------------------------------------------------------

library(Rcpp); library(RcppArmadillo); library(Matrix)
library(GIGrvg); library(MASS); library(ggplot2)

# library(GeneralizedHyperbolic); library(ars); library(lamW); library(clusterGeneration)
# library(ggplot2); library(MCMCpack); library(extraDistr); 
# library(Matrix); library(matrixsampling); 
# library(abind); library(beepr)

library(here)

sourceCpp(here("src","coeff.cpp"))
sourceCpp(here("src","indic_variables.cpp"))
# sourceCpp(here("Model/Functions","variable_update.cpp"))
source(here("R","functions.r"))
source(here("R","runmodel.r"))
source(here("R","plots.r"))
# R CMD INSTALL --preclean --no-multiarch --with-keep.source "functionsForForeach"

jointSpecies <- T
spatialCorr <- T

# DATA --------

# load(here("Data","Lakes.Rdata"))
load(here("Dataset","Malaise.Rdata"))
# load(here("Data","Leeches2.rda"))
# load(here("Data","Leeches.Rdata"))
# load(here("Data","Ponds.Rdata"))

subsetSpecies <- 1:(dim(y)[3])

S <- length(subsetSpecies)
# OTU <- OTU[,subsetSpecies]
y <- y[,,subsetSpecies]

S_star <- 2

# emptyTubes <- 0
# 
# X_w <- apply(X_w, 2, function(x){
#   (x - mean(x)) / sd(x)
# })

# X_w <- as.matrix(X_w)
# n <- 50
# X_z <- X_z[1:n,]
# M_site <- M_site[1:n]
# y <- y[1:sum(M_site),,]
# r <- r[1:sum(M_site),]
# OTU <- OTU[1:sum(M_site),]
# K <- K[1:sum(M_site)]
# emptySites <- emptySites[1:n]
# data_short <- data_short[1:sum(M_site),]


# RUN MODEL -------

data <- list("y" = y,
             "M_site" = M_site,
             "K" = K,
             "emptyTubes" = emptyTubes,
             "S_star" = S_star,
             "spikedSample" = spikedSample,
             "v_spikes" = v_spikes,
             "X_w" = X_w,
             "X_s" = X_s,
             "X_z" = X_z)

if(spatialCorr){
  data <- list("y" = y,
               "M_site" = M_site,
               "K" = K,
               "emptyTubes" = emptyTubes,
               "S_star" = S_star,
               "spikedSample" = spikedSample,
               "v_spikes" = v_spikes,
               "X_w" = X_w,
               "X_s" = X_s,
               "X_z" = X_z,
               "X_s_star" = X_s_star,
               "X_z_star" = X_z_star)
}


# priors
{
  a_theta0 <- 1
  b_theta0 <- 50
  
  a_p11 <- 20
  b_p11 <- 1
  
  a_p10 <- 1
  b_p10 <- 100
  
  a_tau <- 3
  b_tau <- 2
  
  lambda_Y <- 1
  
  a_sigma <- 3
  b_sigma <- 2
  # A_sigma <- .5
  
  sigma_beta <- 1
  
  sigma_mu <- 1
  sigma_gamma <- 1
  
  sigma_u <- 1
  
  sigma_beta_theta <- 1
  
  sigma_lambda <- 3
  
  # nu0 <- 1
  # Sigma0 <- diag(1, nrow = S)
  
  priors <- list("sigma_beta" = sigma_beta,
                 "sigma_lambda" = sigma_lambda,
                 "sigma_gamma" = sigma_gamma,
                 "sigma_beta_theta" = sigma_beta_theta,
                 "sigma_u" = sigma_u,
                 "a_p11" = a_p11,
                 "b_p11" = b_p11,
                 "a_p10" = a_p10,
                 "b_p10" = b_p10,
                 "a_theta0" = a_theta0,
                 "b_theta0" = b_theta0)
}

# mcmc params
{
  MCMCparams <- list("nchain" = 1,
                     "nburn" = 0,
                     "niter" = 5000,
                     "nthin" = 1,
                     "iterToAdapt" = 100000)
}

# params to update
{
  if(simulatedData){
    trueParams <- list("logz_true" = logz_true,
                       "lambda_true" = lambda_true,
                       "mu_true " = mu_true,
                       "v_true" = v_true,
                       "u_true" = u_true,
                       "r_true" = r_true,
                       "beta_w_true" = beta_w_true,
                       "beta0_true" = beta0_true,
                       "beta_z_true" = beta_z_true,
                       "delta_true" = delta_true,
                       "gamma_true" = gamma_true,
                       "c_imk_true" = c_imk_true,
                       "mu0_true" = mu0_true,
                       "n0_true" = n0_true,
                       "pi0_true" = pi0_true,
                       "mu_tilde_true" = mu_tilde_true,
                       "sigma_true" = sigma_true,
                       "beta_theta_true" = beta_theta_true,
                       "theta11_true" = theta11_true,
                       "theta10_true" = theta10_true,
                       "p11_true" = p11_true,
                       "p10_true" = p10_true)
    if(jointSpecies){
      trueParams <- c(trueParams, list("Tau_true" = Tau_true))
    } else {
      trueParams <- c(trueParams, list("tau_true" = tau_true))
    }
  } else {
    trueParams <- NULL
  }
  
  idConstraints = list(beta0equal0 = F,
                       betathetaequal0 = T)
  
  paramsUpdate <- list(
    updateAll = F,
    params =
      list(
        "lambda" = T,
        "logz" = T,
        "mu" = T,
        "v" = T,
        "u" = T,
        "uv" = T,
        "beta_theta" = T,
        "csi" = T,
        "beta_z" = T,
        "beta_w" = T,
        "tau" = T,
        "sigma" = T,
        "deltagammac" = T,
        "r" = T,
        "eta" = T,
        "lambdatilde" = T,
        "lambda0" = T,
        "p11" = T,
        "p10" = T
      ),
    correct = NULL,
    trueParams = trueParams
  )
  
  paramsUpdate <- list(
    updateAll = T,
    params = NULL,
    correct = NULL
  )
  
  # paramsUpdate <- list(
  #   updateAll = F,
  #   params =
  #     list(
  #       "lambda" = T,
  #       "logz" = T,
  #       "mu" = T,
  #       "v" = T,
  #       "u" = T,
  #       "uv" = T,
  #       "beta_theta" = T,
  #       "csi" = F,
  #       "beta_z" = T,
  #       "beta_w" = T,
  #       "tau" = T,
  #       "sigma" = T,
  #       "deltagammac" = F,
  #       "r" = F,
  #       "eta" = T,
  #       "lambdatilde" = T,
  #       "lambda0" = T,
  #       "p11" = T,
  #       "p10" = T
  #     ),
  #   correct = list(
  #     "lambda" = T,
  #     "logz" = T,
  #     "mu" = T,
  #     "v" = T,
  #     "u" = T,
  #     "uv" = T,
  #     "beta_theta" = T,
  #     "csi" = T,
  #     "beta_z" = T,
  #     "beta_w" = T,
  #     "tau" = T,
  #     "sigma" = T,
  #     "deltagammac" = T,
  #     "r" = T,
  #     "eta" = T,
  #     "lambdatilde" = T,
  #     "lambda0" = T,
  #     "p11" = T,
  #     "p10" = T
  #   ),
  #   trueParams = trueParams
  # )
}

# modelResults <- fitModel(data,
#                          priors,
#                          jointSpecies,
#                          paramsUpdate,
#                          MCMCparams)

modelResults <- fitModel(data,
                         priors,
                         jointSpecies,
                         paramsUpdate,
                         MCMCparams,
                         idConstraints)

# DIAGNOSTICS CIMK --------

j <- 1

y_imk_col <- y[,,j]
v_current <- v[,j]
colnames(y_imk_col) <- paste("PCR",1:max(K))#c("PCR1","PCR2","PCR2")
Delta <- delta[,j]
# colnames(Presence) <- "Delta"
Gamma <- gamma[,j]
# colnames(Presence) <- "Gamma"

v_current <- v[,j]

c_imk_col <- as.data.frame(c_imk[,,j])
colnames(c_imk_col) <- paste("PCR_state",1:max(K))#c("PCR1_state","PCR2_state","PCR2_state")
View(
  cbind(data_short,
        y_imk_col,
        # r,
        theta11_true[,j],
        delta_true[,j],
        v_current,
        rep(logz[,j], M_site),
        mu[j],
        Delta,
        Gamma,
        theta11[,j]
        # lambda_imk_col,
        
        # c_imk_col
        # c_imk_col2
  )
)

y[y[,1,j] > 0,,j]
y[c_imk[,1,j] == 1,,j]
which(c_imk[,1,j] == 1)

#

View(
  cbind(
    data_short,
    y_imk_col,
    theta11[,j],
    theta11_true[,j],
    delta[,j],
    delta_true[,j],
    gamma[,j],
    gamma_true[,j],
    # v_current,
    # Presence,
    # lambda_imk_col,
    v[,j],
    v_true[,j],
    logz_true[,j] + X_z %*% beta_z[,j],
    mu[j]
    # c_imk[,,j],
    # c_imk_true[,,j]
    # c_imk_col2
  )
)

j <- 1

df_toshow <- data.frame(
  y = y[,,j],
  delta[,j],
  gamma[,j],
  c_es = c_imk[,,j],
  c_mean = apply(cimk_output[1,,1,2,], 1, function(x){
    mean(x == 1)
  }),
  c_true = c_imk_true[,,j],
  expvlu = exp(v[,j] + lambda[j] + u),
  v[,j],
  v_true[,j],
  logz[,j] + X_w %*% beta_w[,j]
  
)


View(
  df_toshow  
)

# View(cbind(1:nrow(data_short),data_short,y[,,j],
#            # delta_mean[,j],
#            c_imk[,,j],
#            v[,j],
#            exp(lambda[j] + v[,j])))

c_imk[,,1] - c_imk_true[,,1]

# SIMULATED DATA ------------------------------------------------

nchain <- MCMCparams$nchain
nburn <- MCMCparams$nburn
niter <- MCMCparams$niter
nthin <- MCMCparams$nthin

# Covariates coefficients
{
  # beta z
  beta_z_output <- modelResults$params_output$beta_z_output
  
  cov_z <- 1
  j <- 0
  j <- j+1
  diagnosticsPlot( beta_z_output[,,cov_z,j]) + 
    geom_hline(aes(yintercept = beta_z_true[cov_z,j]), color = "red")  +
    ylim(c(-3,3))
  
  # beta w
  beta_w_output <- modelResults$params_output$beta_w_output
  
  cov_w <- 1
  j <- 0
  j <- j+1
  diagnosticsPlot( beta_w_output[,cov_w,j,]) + 
    geom_hline(aes(yintercept = beta_w_true[cov_w,j]), color = "red")  +
    ylim(c(-5,5))
  
  # beta theta w
  beta_theta_output <- modelResults$params_output$beta_theta_output
  
  cov_w <- 1
  j <- 0
  j <- j + 1
  diagnosticsPlot( beta_theta_output[,j,2 + cov_w,]) + 
    geom_hline(aes(yintercept = beta_theta_true[j,2 + cov_w]), color = "red")  +
    ylim(c(-2,2))
  
  
  ggsave(betas_plot, file = paste0("beta",j,"_plot.png"))
}

# Other coefficients
{
  # lambda
  lambda_output <- modelResults$params_output$lambda_output
  
  j <- 0
  j <- j+1
  diagnosticsPlot( lambda_output[,j,]) + 
    geom_hline(aes(yintercept = lambda_true[j]), color = "red", size = 2) + 
    geom_hline(aes(yintercept = lambda_prior[j]), color = "green")
  
  # beta 0
  beta0_output <- modelResults$params_output$beta0_output
  
  j <- 0
  j <- j+1
  diagnosticsPlot( beta0_output[,j,]) + 
    geom_hline(aes(yintercept = beta0_true[j]), color = "red")  #+
  # ylim(c(-3,3))
  
  # mu
  j <- 0
  j <- j+1
  diagnosticsPlot( mu_output[,j,]) + 
    geom_hline(aes(yintercept = mu_true[j]), color = "red")  #+
  # ylim(c(-3,3))
  
  # r
  j <- 0
  j <- j+1
  diagnosticsPlot( r_output[,j,]) + 
    geom_hline(aes(yintercept = r_true[j]), color = "red")  #+
  # ylim(c(0,2))
  
  # beta theta 0
  j <- 0
  j <- j + 1
  diagnosticsPlot( beta_theta_output[,j,1,]) + 
    geom_hline(aes(yintercept = beta_theta_true[j,1]), color = "red")  +
    ylim(c(-2,2))
  
  
  # beta theta 1
  j <- 0
  j <- j + 1
  diagnosticsPlot( beta_theta_output[,j,2,]) + 
    geom_hline(aes(yintercept = beta_theta_true[j,2]), color = "red")  +
    ylim(c(0,2))
  
  # tau
  j <- 0
  j <- j + 1
  diagnosticsPlot( tau_output[,j,]) + 
    geom_hline(aes(yintercept = tau_true[j]), color = "red")  #+
  # ylim(c(0,2))
  
  # sigma
  j <- 0
  j <- j + 1
  diagnosticsPlot( sigma_output[,j,]) + 
    geom_hline(aes(yintercept = sigma_true[j]), color = "red")  #+
  # ylim(c(0,2))
  
  # lambda tilde
  j <- 0
  j <- j+1
  diagnosticsPlot( mu_[,j,]) + 
    geom_hline(aes(yintercept = lambdatilde_true[j]), color = "red", size = 2) #+ 
  # geom_hline(aes(yintercept = lambda_prior[j]), color = "green")
  
}

# Individual parameters
{
  # logz
  j <- 1
  i <- 1
  i <- i+1
  diagnosticsPlot( logz_output[,i,j,] ) + 
    geom_hline(aes(yintercept = logz_true[i,j]), color = "red", size = 1) # +
  
  # v
  j <- 1
  i <- 1
  m <- 1
  l <- m + sum(M_site[seq_len(i-1)])
  diagnosticsPlot( v_output[,l,j,]) + 
    geom_hline(aes(yintercept = v_true[l,j]), color = "red")  #+
  # ylim(c(-2,2))
  l <- l + 1
  
  # u 
  i <- 1
  m <- 1
  l <- m + sum(M_site[seq_len(i-1)])
  k <- 1
  diagnosticsPlot( u_output[,l,k,]) + 
    geom_hline(aes(yintercept = u_true[l,k]), color = "red", size = 1) + 
    ylim(c(-2,2))
  l <- l + 1
  
  # v + u 
  i <- 1
  m <- 1
  j <- 1
  k <- 1
  l <- m + sum(M_site[seq_len(i-1)])
  diagnosticsPlot( v_output[,l,j,] + u_output[,l,k,]) + 
    geom_hline(aes(yintercept = v_true[l,j] + u_true[l,k]), color = "red", size = 1) + 
    ylim(c(-2,2))
  l <- l + 1
  
  # lambda ijk
  i <- 1
  m <- 1
  l <- m + sum(M_site[seq_len(i-1)])
  k <- 1
  j <- 1
  diagnosticsPlot( lambda_ijk_output[,l,k,j,])
  l <- l + 1
  
  # logz differences
  j <- 1
  i <- 1
  i <- i+1
  diagnosticsPlot( (logz_output[,i,j,] - logz_output[,i+1,j,]) ) + 
    geom_hline(aes(yintercept = (logz_true[i,j] - logz_true[i+1,j])), color = "red", size = 1) +
    ylim(c(-5,5))
  
  # v - logz
  j <- 1
  i <- 1
  m <- 1
  i <- i + 1
  m <- m + 1
  l <- m + sum(M_site[seq_len(i-1)])
  diagnosticsPlot( v_output[,l,j,] - logz_output[,i,j,]) + 
    geom_hline(aes(yintercept = v_true[l,j] - logz_true[i,j]), color = "red")  #+
  # ylim(c(-2,2))
  
  # beta theta * v
  j <- 1
  i <- 1
  i <- i+1
  diagnosticsPlot( logz_output[1,i,j,] * beta_theta_output[1,j,2,]) + 
    geom_hline(aes(yintercept = logz_true[i,j] * beta_theta_true[j,2]), color = "red") 
  
  # theta 11  
  j <- 1
  l <- 0
  l <- l + 1
  diagnosticsPlot( theta11_output[,l,j,]) + 
    geom_hline(aes(yintercept = theta11_true[l,j]), color = "red") + ylim(c(0,1)) 
  
  j <- 0
  j <- j + 1
  diagnosticsPlot( beta_theta_output[,j,2,]) + 
    geom_hline(aes(yintercept = beta_theta_true[j,2]), color = "red")  +
    ylim(c(0,1))
  
  # ggsave(logz_plot, file = paste0("beta_theta1_",j,"_plot.png"))
  
  j <- 1
  (beta_theta1_plot <- qplot(1:niter, beta_theta_output[1,j,2,], geom = "line") + 
      geom_hline(aes(yintercept = beta_theta_true[j,2]), color = "red") ) + 
    ylim(c(-2, 2))
  
  # ggsave(beta_theta1_plot, file = paste0("beta_theta_",j,"_plot.png"))
  j <- j+1
}

# amount of DNA coefficients
{
  beta_theta2_CI <- apply(beta_theta_output[,2,], 1, function(x){
    quantile(x, probs = c(.05,.5,.95))
  })
  
  # speciesOrder <- order(beta_CI[3,])
  # 
  
  
  # subsetSpecies <- speciesOrder[1:34]
  subsetSpecies <- 21:60
  
  
  ggplot(data = NULL, aes(x = subsetSpecies,
                          y = beta_theta2_CI[2,subsetSpecies],
                          ymin = beta_theta2_CI[1,subsetSpecies],
                          ymax = beta_theta2_CI[3,subsetSpecies])) + geom_errorbar() + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                   axis.title = ggplot2::element_text(size = 20, face = "bold"),
                   axis.text = ggplot2::element_text(size = 13, face = "bold", angle = 90),
                   panel.grid.major = ggplot2::element_line(colour="grey", size=0.15),
                   panel.background = ggplot2::element_rect(fill = "white", color = "black")) +
    geom_hline(aes(yintercept = 0), color = "red") + 
    scale_x_discrete(name = "Species") 
  # scale_x_continuous(breaks = subsetSpecies, name = "Species",
  #                    labels = colnames(OTU)[subsetSpecies]) +
  # scale_y_continuous(name = "Elevation")
  
}

# Tau
{
  Tau_output <- modelResults$params_output$Tau_output
  
  mean(Tau_output[1,,1,1])
  mean(Tau_output[1,,4,4])
  
  j1 <- 1
  j2 <- 3
  diagnosticsPlot( Tau_output[,,j1,j2]) + 
    geom_hline(aes(yintercept = Tau_true[j1,j2]), color = "red") +
    ylim(c(-2, 2))
  
  qplot(1:niter, Tau_output[1,,j1,j2], geom = "line") + 
    geom_hline(aes(yintercept = Tau_true[j1,j2]), color = "red")
  
  Tau_CI <- apply(Tau_output[1,,,], c(2,3), function(x){
    quantile(x, probs = c(0.025, 0.975))
  })
  
  Tau_sign <- apply(Tau_CI, c(2,3), function(x){
    sign(x[1]) == sign(x[2])
  })
  
  # confusion matrix
  
  confMatrix <- matrix(0, 2, 2)
  
  confMatrix[1,1] <- sum(Tau_sign & abs(Tau_true) > .05)
  confMatrix[1,2] <- sum(!Tau_sign & abs(Tau_true) > .05)
  confMatrix[2,1] <- sum(Tau_sign & abs(Tau_true) < .05)
  confMatrix[2,2] <- sum(!Tau_sign & abs(Tau_true) < .05)
  
  (confMatrix <- confMatrix / sum(confMatrix))
  
  # zero coefficients
  {
    Tau_corr_output <- array(NA, dim = c(nchain, niter, S, S))
    for (chain in 1:nchain) {
      for (iter in 1:niter) {
        Tau_corr_output[chain, iter,,] <- Tau_output[chain, iter,,]
        diagCorr <- diag(Tau_corr_output[chain, iter,,])
        for (j1 in 1:S) {
          for (j2 in 1:S) {
            Tau_corr_output[chain, iter,j1,j2] <- 
              Tau_corr_output[chain, iter,j1,j2] / sqrt(diagCorr[j1] * diagCorr[j2])
          }
        }
      }
    }
    
    zeroCoeffs <- which(Tau_true == 0, arr.ind = T)
    zeroCoeffs <- zeroCoeffs[zeroCoeffs[,1] > zeroCoeffs[,2],]
    
    Tau_output_zeros <- mapply(function(i, j) Tau_corr_output[,,i,j], 
                               zeroCoeffs[,1], zeroCoeffs[,2])
    CI_taus <- apply(Tau_output_zeros, 2, function(x){
      quantile(x, probs = c(0.025,0.5,0.975))
    })
    
    ggplot(data = NULL, aes(x = 1:nrow(zeroCoeffs),
                            y = CI_taus[2,],
                            ymin = CI_taus[1,],
                            ymax = CI_taus[3,])) + geom_errorbar() + 
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                     axis.title = ggplot2::element_text(size = 20, face = "bold"),
                     axis.text = ggplot2::element_text(size = 9, face = "bold", angle = 90),
                     panel.grid.major = ggplot2::element_line(colour="grey", size=0.15),
                     panel.background = ggplot2::element_rect(fill = "white", color = "black")) +
      geom_hline(aes(yintercept = 0), color = "red") + ylim(c(-1,1)) 
    # scale_x_continuous(breaks = subsetSpecies, name = "Species",
    # labels = colnames(OTU)[subsetSpecies]) +
    # scale_y_continuous(name = "Elevation") +
    
  }
  
  # nonzero coefficients
  {
    nonzeroCoeffs <- which(Tau_true != 0, arr.ind = T)
    nonzeroCoeffs <- nonzeroCoeffs[nonzeroCoeffs[,1] > nonzeroCoeffs[,2],]
    nonzerovals <- mapply(function(i, j) Tau_true[i,j], 
                          nonzeroCoeffs[,1], nonzeroCoeffs[,2])
    
    # Tau_output_nonzeros <- mapply(function(i, j) Tau_corr_output[,,i,j], 
    Tau_output_nonzeros <- mapply(function(i, j) Tau_output[,,i,j], 
                                  nonzeroCoeffs[,1], nonzeroCoeffs[,2])
    CI_taus <- apply(Tau_output_nonzeros, 2, function(x){
      quantile(x, probs = c(0.025,0.5,0.975))
    })
    
    ggplot(data = NULL, aes(x = 1:nrow(nonzeroCoeffs),
                            y = CI_taus[2,],
                            ymin = CI_taus[1,],
                            ymax = CI_taus[3,])) + geom_errorbar() + 
      geom_point(data = NULL, aes(x = 1:nrow(nonzeroCoeffs),
                                  y = nonzerovals), color = "red",
                 size = 2) + 
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                     axis.title = ggplot2::element_text(size = 20, face = "bold"),
                     axis.text = ggplot2::element_text(size = 9, face = "bold", angle = 90),
                     panel.grid.major = ggplot2::element_line(colour="grey", size=0.15),
                     panel.background = ggplot2::element_rect(fill = "white", color = "black")) +
      ylim(c(-1,1)) 
    # scale_x_continuous(breaks = subsetSpecies, name = "Species",
    # labels = colnames(OTU)[subsetSpecies]) +
    # scale_y_continuous(name = "Elevation") +
    # ylim(c(-3,3))
  }
  
}

# diagnostics
{
  lambda_output_iter[2,(iter - nburn) - 10:1]
  mu_output_iter[2,(iter - nburn) - 10:1]
  logz_output_iter[1,2,(iter - nburn) - 10:1]
  v_output_iter[2,2,(iter - nburn) - 10:1]
  u_output_iter[2,2,(iter - nburn) - 10:1]
}

# TRUE DATA ------------------------------------------------

speciesToAnalyze <- c(1,7,9,4,5,10,11,15,30)
colnames(OTU)[speciesToAnalyze]

# biomass coefficent CI 
{
  beta_CI <- apply(beta_z_output, c(2,3), function(x){
    quantile(x, probs = c(.05,.5,.95))
  })
  
  # beta_CI <- beta_CI[1,]
  
  subsetSpecies <- 1:34
  
  ggplot(data = NULL, aes(x = subsetSpecies,
                          y = beta_CI[2,subsetSpecies],
                          ymin = beta_CI[1,subsetSpecies],
                          ymax = beta_CI[3,subsetSpecies])) + geom_errorbar() + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                   axis.title = ggplot2::element_text(size = 20, face = "bold"),
                   axis.text = ggplot2::element_text(size = 13, face = "bold", angle = 90),
                   panel.grid.major = ggplot2::element_line(colour="grey", size=0.15),
                   panel.background = ggplot2::element_rect(fill = "white", color = "black")) +
    geom_hline(aes(yintercept = 0), color = "red") + 
    scale_x_continuous(breaks = subsetSpecies, name = "Species",
                       labels = colnames(OTU)[subsetSpecies]) +
    scale_y_continuous(name = "Elevation")
  
}

# biomass coefficent CI (ordered)
{
  covariate_num <- 1
  
  beta_CI <- apply(beta_z_output, c(2,3), function(x){
    quantile(x, probs = c(.05,.5,.95))
  })
  
  significantSpecies <- 
    # 1:S
    which(beta_CI[3,covariate_num,] < 0 |
            beta_CI[1,covariate_num,] > 0)
  
  significantSpecies <- significantSpecies[-c(10,19)]
  
  orderSignificantSpecies <- 
    significantSpecies[order(beta_CI[2,covariate_num,significantSpecies])]
  
  orderSignificantSpecies <- speciesToAnalyze
  
  beta_CI_subset <- beta_CI[,covariate_num,orderSignificantSpecies]
  colnames(beta_CI_subset) <- colnames(OTU)[orderSignificantSpecies]
  # speciesOrder <- order(beta_CI[1,covariate_num,])
  
  subsetSpecies <- 1:length(orderSignificantSpecies)
  
  factorSub <- factor(subsetSpecies, levels = subsetSpecies)
  
  namesSpecies <- colnames(OTU)[orderSignificantSpecies]
  
  namesSpecies <- c("Mammals - Cow",
                    "Mammals - Sheep",
                    "Mammals - Goat",
                    "Frog - Bombina Maxima",
                    "Frog - Oreolalax",
                    "Frog - Rhacophorus",
                    "Frog - Nanorana",
                    "Rodents - Shrew Hedgehog",
                    "Rodents - Oriental Squirrel")
  
  # namesSpecies <- c("Shrew Hedgehog",
  #                   "Frog - Rhacophoridae",
  #                   "Frog - Dicroglossidae Nanorana 1",
  #                   "Frog - Bombinatoridae Bombina",
  #                   "Frog - Dicroglossidae Nanorana unculuanus", #, "Anura,Dicroglossidae,Nanorana,unculuanus"
  #                   "Frog - Bombina maxima" ,
  #                   "Frog - Oreolalax jingdongensis"     ,
  #                   "Frog - Dicroglossidae Nanorana 2"    ,
  #                   "Frog - Oreolalax" ,
  #                   "Frog - Xenophrys"   ,
  #                   "Frog - Bombina 1"  ,
  #                   "Frog - Bombina 2"    ,
  #                   "Frog - Nanorana yunnanensis",
  #                   "Frog - Nanorana 3"    ,
  #                   "Goat"          ,
  #                   "Brush tailed porcupine"    ,
  #                   "Sheep"              ,
  #                   "Oriental squirrel"  )
  
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
  
  # diagnostics
  {
    j <- 0
    j <- j+1
    # diagnosticsPlot( beta_z_output[1,1,j,]) + 
    diagnosticsPlot( beta_z_output[1,1,significantSpecies[j],]) +
      ylim(c(-2,2))
    
    }
  
  
}

# coefficient 
{
  alpha_CI <- apply(alpha_output, 1, function(x){
    quantile(x, probs = c(.025,.5,.975))
  })
  
  ggplot(data = NULL, aes(x = 1:S,
                          y = alpha_CI[2,],
                          ymin = alpha_CI[1,],
                          ymax = alpha_CI[3,])) + geom_errorbar() + 
    theme_bw()
  
  #
  
  beta_theta_CI <- apply(beta_theta_output[,3,], 1, function(x){
    quantile(x, probs = c(.025,.5,.975))
  })
  
  ggplot(data = NULL, aes(x = 1:S,
                          y = beta_theta_CI[2,],
                          ymin = beta_theta_CI[1,],
                          ymax = beta_theta_CI[3,])) + geom_errorbar() + 
    theme_bw()
}

# DNA coefficients
{
  j <- 0
  j <- j+1
  diagnosticsPlot( beta_z_output[,1,j,]) + 
    ylim(c(-2,2))
  
  j <- 1
  betas_plot <- qplot(1:(niter / nthin), beta_z_output[1,j,], geom = "line") 
  # qplot(1:niter, beta_z_output[3,j,], geom = "line") + geom_hline(aes(yintercept = beta_z_true[3,j]), color = "red") + 
  # ylim(c(-2,2)) 
  
  ggsave(betas_plot, file = paste0("beta_",j,"_plot.png"))
  j <- j+1
  
  beta_z_output_all <-
    matrix(beta_z_output[, 1, j, ], nrow = nchain, ncol = niter)
  
  beta_z_output_long <- reshape2::melt(beta_z_output_all)
  
  ggplot2::ggplot(data = beta_z_output_long, ggplot2::aes(
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
  j <- j + 1
  
}

# General coefficients
{
  # lambda
  j <- 0
  j <- j+1
  diagnosticsPlot( lambda_output[,j,]) + 
    ylim(c(0,10))
  
  # beta 0
  j <- 0
  j <- j+1
  diagnosticsPlot( beta0_output[,j,]) +
    ylim(c(-3,3))
  
  # beta z
  j <- 0
  j <- j+1
  diagnosticsPlot( beta_z_output[,1,j,]) + 
    ylim(c(-2,2))
  
  # alpha
  j <- 0
  j <- j+1
  diagnosticsPlot( alpha_output[,j,]) + 
    ylim(c(-2,2))
  
  # r
  j <- 0
  j <- j+1
  diagnosticsPlot( r_output[,j,]) + 
    ylim(c(0,20))
  
  # beta theta 0
  j <- 0
  j <- j + 1
  diagnosticsPlot( beta_theta_output[,j,1,]) + 
    ylim(c(-4,0))
  
  # beta theta 1
  j <- 0
  j <- j + 1
  diagnosticsPlot( beta_theta_output[,j,2,]) + 
    ylim(c(0,2))
  
  # tau
  j <- 0
  j <- j + 1
  diagnosticsPlot( tau_output[,j,]) + 
    ylim(c(0,2))
  
  # sigma
  j <- 0
  j <- j + 1
  diagnosticsPlot( sigma_output[,j,]) + 
    ylim(c(0,10))
  
}

# individual parameters
{
  # logz
  j <- 1
  i <- 1
  i <- i+1
  diagnosticsPlot( logz_output[,i,j,]) + 
    ylim(c(-2,2))
  
  # beta theta * v
  j <- 1
  i <- 1
  i <- i+1
  diagnosticsPlot( logz_output[1,i,j,] * beta_theta_output[1,j,2,]) + 
    geom_hline(aes(yintercept = logz_true[i,j] * beta_theta_true[j,2]), color = "red") 
  
  # theta 11  
  j <- 1
  l <- 0
  l <- l + 1
  diagnosticsPlot( theta11_output[,l,j,]) +  ylim(c(0,1)) 
  
  
  
  #
  j <- 1
  i <- 1
  m <- 1
  l <- m + sum(M_site[seq_len(i-1)])
  diagnosticsPlot( v_output[,l,j,]) 
  
  
  j <- 0
  j <- j + 1
  diagnosticsPlot( beta_theta_output[,j,2,]) + 
    geom_hline(aes(yintercept = beta_theta_true[j,2]), color = "red")  +
    ylim(c(0,1))
  
  # ggsave(logz_plot, file = paste0("beta_theta1_",j,"_plot.png"))
  
  j <- 1
  (beta_theta1_plot <- qplot(1:niter, beta_theta_output[1,j,2,], geom = "line") + 
      geom_hline(aes(yintercept = beta_theta_true[j,2]), color = "red") ) + 
    ylim(c(-2, 2))
  
  # ggsave(beta_theta1_plot, file = paste0("beta_theta_",j,"_plot.png"))
  j <- j+1
}

# correlation
{
  Tau_corr_output_list <- lapply(1:niter, function(i){
    cov2cor(Tau_output[1,i,,])
  })
  
  Tau_corr_output <- array(NA, dim = c(niter, S, S))
  for (iter in 1:niter) {
    Tau_corr_output[iter,,] <- Tau_corr_output_list[[iter]]
  }
  
  Tau_CI <- apply(Tau_corr_output, c(2,3), function(x){
    quantile(x, probs = c(0.025, 0.975, .5))
  })
  
  Tau_mean <- apply(Tau_corr_output, c(2,3), mean)
  
  cor_Tau <- Tau_mean
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
  
  cor_Tau_long <- matrix(NA, S * (S - 1) / 2, 5)
  l <- 1
  for (i in 1:S) {
    for (j in seq_len(i - 1)) {
      cor_Tau_long[l,1] <- Tau_CI[1,i,j]
      cor_Tau_long[l,2] <- Tau_CI[2,i,j]
      cor_Tau_long[l,3] <- Tau_CI[3,i,j]
      cor_Tau_long[l,4] <- i
      cor_Tau_long[l,5] <- j
      l <- l + 1
    }
  }
  
  significantCorrelations <- 
    which((cor_Tau_long[,1] < 0 & cor_Tau_long[,2] < 0) | 
            (cor_Tau_long[,1] > 0 & cor_Tau_long[,2] > 0) )
  cor_Tau_long_significants <- cor_Tau_long[significantCorrelations,]
  cor_Tau_long_ordered <- 
    cor_Tau_long_significants[order(cor_Tau_long_significants[,3]),]
  
  biggestSpecies_idx <- which((cor_Tau_long_ordered[,3] < 0 & cor_Tau_long_ordered[,2] < -0.3) | 
                                (cor_Tau_long_ordered[,3] > 0 & cor_Tau_long_ordered[,1] > 0.3))
  
  biggestSpecies <- cor_Tau_long_ordered[biggestSpecies_idx,4:5]
  uniquebiggestSpecies <- unique(as.vector(biggestSpecies))
  
  uniquebiggestSpecies <- uniquebiggestSpecies[-5]
  
  uniquebiggestSpecies <- speciesToAnalyze
  
  corr_biggest <- cor_Tau[uniquebiggestSpecies, uniquebiggestSpecies]
  
  
  
  namesSpecies <- colnames(OTU)[uniquebiggestSpecies]
  
  namesSpecies <- c("Mammals - Cow",
                    "Mammals - Sheep",
                    "Mammals - Goat",
                    "Frog - Bombina Maxima",
                    "Frog - Oreolalax",
                    "Frog - Rhacophorus",
                    "Frog - Nanorana",
                    "Rodents - Shrew Hedgehog",
                    "Rodents - Oriental Squirrel")
  
  colnames(corr_biggest) <- namesSpecies
  rownames(corr_biggest) <- namesSpecies
  
  library(corrplot)
  corrplot(corr_biggest)  
}

# true positives
{
  p_CI <- apply(p_10_output, 2, function(x){
    quantile(x, probs = c(0.025, 0.5, 0.975))
  })
  
  ggplot(data = NULL, aes(x = 1:S,
                          y = p_CI[2,],
                          ymin = p_CI[1,],
                          ymax = p_CI[3,])) + geom_errorbar() + geom_point() + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                   axis.title = ggplot2::element_text(size = 20, face = "bold"),
                   axis.text = ggplot2::element_text(size = 9, face = "bold", angle = 0),
                   panel.grid.major = ggplot2::element_line(colour="grey", size=0.15),
                   panel.background = ggplot2::element_rect(fill = "white", color = "black")) +
    # geom_hline(aes(yintercept = 0), color = "red") + 
    scale_x_discrete(name = "", breaks = subsetSpecies,
                     labels = namesSpecies) +
    # scale_x_continuous(breaks = subsetSpecies, name = "Species",
    # labels = colnames(OTU)[subsetSpecies]) +
    scale_y_continuous(name = "PCR false positive rate", 
                       limits = c(0, .03)
    ) + coord_flip()
}

# tau sigma
{
  tau_CI <- sapply(1:S, function(j){
    quantile(log(sqrt(Tau_output[,,j,j])), 
             probs = c(0.025, 0.5, 0.975))
  })
  
  
  sigma_CI <- apply(log(sigma_output), 2, function(x){
    quantile(x, probs = c(0.025, 0.5, 0.975))
  })
  
  ggplot() + 
    geom_errorbar(data = NULL, aes(x = 1:S,
                                   y = tau_CI[2,],
                                   ymin = tau_CI[1,],
                                   ymax = tau_CI[3,])) + 
    geom_errorbar(data = NULL, aes(x = 1:S,
                                   y = sigma_CI[2,],
                                   ymin = sigma_CI[1,],
                                   ymax = sigma_CI[3,]), color = "blue") + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                   axis.title = ggplot2::element_text(size = 20, face = "bold"),
                   axis.text = ggplot2::element_text(size = 9, face = "bold", angle = 0),
                   panel.grid.major = ggplot2::element_line(colour="grey", size=0.15),
                   panel.background = ggplot2::element_rect(fill = "white", color = "black")) +
    # geom_hline(aes(yintercept = 0), color = "red") + 
    scale_x_discrete(name = "Species", breaks = subsetSpecies,
                     labels = namesSpecies) +
    # scale_x_continuous(breaks = subsetSpecies, name = "Species",
    # labels = colnames(OTU)[subsetSpecies]) +
    scale_y_continuous(name = "Log Variance") + coord_flip()
}

# biomasses over time (francesco)
{
  
  years <- c(-8212, -7358, -5580, -4435, -4020, -3603, -3350, -3207, -3071, 
             -2906, -2639, -2394, -2134, -1945, -1767, -1594, -1453, -1377, 
             -1280, -1129, -1047, -930, -829, -676, -590, -495, -388, -282, 
             -176, -51, 58, 162, 297, 407, 563, 713, 863, 1027, 1124, 1235, 
             1339, 1429, 1504, 1607, 1779, 1975, 2006)
  
  biom_CI <- apply(logz_output, c(2,3), function(x){
    quantile(x, probs = c(0.025, 0.5, 0.975))
  })
  
  setwd("/cluster/home/osr/ad625/eDNA/Project/Model/Results/Biomass")
  
  for (j in 1:S) {
    
    currentplot <- ggplot(data = NULL, aes(x = years,
                                           y = biom_CI[2,,j],
                                           ymin = biom_CI[1,,j],
                                           ymax = biom_CI[3,,j])) + 
      geom_point() + geom_errorbar() + geom_line() + 
      scale_x_continuous(breaks = years,
                         labels = years,
                         name = "Years") +
      ylab("Log-Biomass") + 
      ggtitle(colnames(OTU)[j]) + 
      theme(plot.title = element_text(hjust = 0.5, size = 17),
            axis.title = element_text(size = 16, face = "bold"),
            axis.text.y = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(size = 11, face = "bold", angle = 90, hjust = 1),
            axis.line = element_line(colour="black", size=0.15),
            # panel.grid.minor = element_line(colour="grey", size=0.15),
            panel.grid.major = element_line(colour="grey", size=0.15),
            panel.background = element_rect(fill = "white", color = "black"))
    
    ggsave(plot = currentplot, 
           filename = paste0(j," - ",colnames(OTU)[j],".jpeg"))
  }
  
}

# BIOMASSES -------

biom_mean <- apply(logz_output, c(2,3), mean)




# SIMULATED DATA ------------------------------------------------

# Covariates coefficients
{
  # z
  cov_z <- 1
  j <- 0
  j <- j+1
  diagnosticsPlot( beta_z_output[,cov_z,j,]) +  
    ylim(c(-3,3))
  
  # beta w
  cov_w <- 1
  j <- 0
  j <- j+1
  diagnosticsPlot( beta_w_output[,cov_w,j,]) + 
    ylim(c(-5,5))
  
  # beta theta w
  cov_w <- 1
  j <- 0
  j <- j + 1
  diagnosticsPlot( beta_theta_output[,j,2 + cov_w,]) + 
    ylim(c(-2,2))
  
}

# Other coefficients
{
  # lambda
  j <- 0
  j <- j+1
  diagnosticsPlot( lambda_output[,j,]) 
  
  # beta 0
  j <- 0
  j <- j+1
  diagnosticsPlot( beta0_output[,j,]) 
  # ylim(c(-3,3))
  
  # r
  j <- 0
  j <- j+1
  diagnosticsPlot( r_output[,j,]) 
  # ylim(c(0,2))
  
  # beta theta 0
  j <- 0
  j <- j + 1
  diagnosticsPlot( beta_theta_output[,j,1,]) 
  # ylim(c(-2,2))
  
  # beta theta 1
  j <- 0
  j <- j + 1
  diagnosticsPlot( beta_theta_output[,j,2,]) + 
    geom_hline(aes(yintercept = beta_theta_true[j,2]), color = "red")  #+
  # ylim(c(0,2))
  
  # tau
  j <- 0
  j <- j + 1
  diagnosticsPlot( tau_output[,j,]) + 
    geom_hline(aes(yintercept = tau_true[j]), color = "red")  #+
  # ylim(c(0,2))
  
  # sigma
  j <- 0
  j <- j + 1
  diagnosticsPlot( sigma_output[,j,]) 
  # ylim(c(0,2))
  
}

# Individual parameters
{
  # logz
  j <- 1
  i <- 1
  i <- i+1
  diagnosticsPlot( logz_output[,i,j,] ) 
  
  # v
  j <- 1
  i <- 1
  m <- 1
  l <- m + sum(M_site[seq_len(i-1)])
  diagnosticsPlot( v_output[,l,j,]) 
  # ylim(c(-2,2))
  l <- l + 1
  
  # u 
  i <- 1
  m <- 1
  l <- m + sum(M_site[seq_len(i-1)])
  k <- 1
  diagnosticsPlot( u_output[,l,k,]) + 
    ylim(c(-2,2))
  l <- l + 1
  
  # lambda ijk
  i <- 1
  m <- 1
  l <- m + sum(M_site[seq_len(i-1)])
  k <- 1
  j <- 1
  diagnosticsPlot( lambda_ijk_output[,l,k,j,])
  l <- l + 1
  
  # logz differences
  j <- 1
  i <- 1
  i <- i+1
  diagnosticsPlot( (logz_output[,i,j,] - logz_output[,i+1,j,]) ) + 
    geom_hline(aes(yintercept = (logz_true[i,j] - logz_true[i+1,j])), color = "red", size = 1) +
    ylim(c(-5,5))
  
  # v - logz
  j <- 1
  i <- 1
  m <- 1
  i <- i + 1
  m <- m + 1
  l <- m + sum(M_site[seq_len(i-1)])
  diagnosticsPlot( v_output[,l,j,] - logz_output[,i,j,]) + 
    geom_hline(aes(yintercept = v_true[l,j] - logz_true[i,j]), color = "red")  #+
  # ylim(c(-2,2))
  
  # beta theta * v
  j <- 1
  i <- 1
  i <- i+1
  diagnosticsPlot( logz_output[1,i,j,] * beta_theta_output[1,j,2,]) + 
    geom_hline(aes(yintercept = logz_true[i,j] * beta_theta_true[j,2]), color = "red") 
  
  # theta 11  
  j <- 1
  l <- 0
  l <- l + 1
  diagnosticsPlot( theta11_output[,l,j,]) + 
    geom_hline(aes(yintercept = theta11_true[l,j]), color = "red") + ylim(c(0,1)) 
  
  j <- 0
  j <- j + 1
  diagnosticsPlot( beta_theta_output[,j,2,]) + 
    geom_hline(aes(yintercept = beta_theta_true[j,2]), color = "red")  +
    ylim(c(0,1))
  
  # ggsave(logz_plot, file = paste0("beta_theta1_",j,"_plot.png"))
  
  j <- 1
  (beta_theta1_plot <- qplot(1:niter, beta_theta_output[1,j,2,], geom = "line") + 
      geom_hline(aes(yintercept = beta_theta_true[j,2]), color = "red") ) + 
    ylim(c(-2, 2))
  
  # ggsave(beta_theta1_plot, file = paste0("beta_theta_",j,"_plot.png"))
  j <- j+1
}

# amount of DNA coefficients
{
  beta_theta2_CI <- apply(beta_theta_output[,2,], 1, function(x){
    quantile(x, probs = c(.05,.5,.95))
  })
  
  # speciesOrder <- order(beta_CI[3,])
  # 
  
  
  # subsetSpecies <- speciesOrder[1:34]
  subsetSpecies <- 21:60
  
  
  ggplot(data = NULL, aes(x = subsetSpecies,
                          y = beta_theta2_CI[2,subsetSpecies],
                          ymin = beta_theta2_CI[1,subsetSpecies],
                          ymax = beta_theta2_CI[3,subsetSpecies])) + geom_errorbar() + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                   axis.title = ggplot2::element_text(size = 20, face = "bold"),
                   axis.text = ggplot2::element_text(size = 13, face = "bold", angle = 90),
                   panel.grid.major = ggplot2::element_line(colour="grey", size=0.15),
                   panel.background = ggplot2::element_rect(fill = "white", color = "black")) +
    geom_hline(aes(yintercept = 0), color = "red") + 
    scale_x_discrete(name = "Species") 
  # scale_x_continuous(breaks = subsetSpecies, name = "Species",
  #                    labels = colnames(OTU)[subsetSpecies]) +
  # scale_y_continuous(name = "Elevation")
  
}

# Tau
{
  j1 <- 4
  j2 <- 5
  diagnosticsPlot( Tau_output[,,j1,j2]) + 
    geom_hline(aes(yintercept = Tau_true[j1,j2]), color = "red") +
    ylim(c(-2, 2))
  
  qplot(1:niter, Tau_output_iter[j1,j2,], geom = "line") + 
    geom_hline(aes(yintercept = Tau_true[j1,j2]), color = "red")
  
  Tau_CI <- apply(Tau_output_iter, c(2,3), function(x){
    quantile(x, probs = c(0.025, 0.975))
  })
  
  Tau_sign <- apply(Tau_CI, c(2,3), function(x){
    sign(x[1]) == sign(x[2])
  })
  
  # confusion matrix
  
  confMatrix <- matrix(0, 2, 2)
  
  confMatrix[1,1] <- sum(Tau_sign & abs(Tau_true) > .05)
  confMatrix[1,2] <- sum(!Tau_sign & abs(Tau_true) > .05)
  confMatrix[2,1] <- sum(Tau_sign & abs(Tau_true) < .05)
  confMatrix[2,2] <- sum(!Tau_sign & abs(Tau_true) < .05)
  
  (confMatrix <- confMatrix / sum(confMatrix))
  
  # zero coefficients
  {
    Tau_corr_output <- array(NA, dim = c(nchain, niter, S, S))
    for (chain in 1:nchain) {
      for (iter in 1:niter) {
        Tau_corr_output[chain, iter,,] <- Tau_output[chain, iter,,]
        diagCorr <- diag(Tau_corr_output[chain, iter,,])
        for (j1 in 1:S) {
          for (j2 in 1:S) {
            Tau_corr_output[chain, iter,j1,j2] <- 
              Tau_corr_output[chain, iter,j1,j2] / sqrt(diagCorr[j1] * diagCorr[j2])
          }
        }
      }
    }
    
    zeroCoeffs <- which(Tau_true == 0, arr.ind = T)
    zeroCoeffs <- zeroCoeffs[zeroCoeffs[,1] > zeroCoeffs[,2],]
    
    Tau_output_zeros <- mapply(function(i, j) Tau_corr_output[,,i,j], 
                               zeroCoeffs[,1], zeroCoeffs[,2])
    CI_taus <- apply(Tau_output_zeros, 2, function(x){
      quantile(x, probs = c(0.025,0.5,0.975))
    })
    
    ggplot(data = NULL, aes(x = 1:nrow(zeroCoeffs),
                            y = CI_taus[2,],
                            ymin = CI_taus[1,],
                            ymax = CI_taus[3,])) + geom_errorbar() + 
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                     axis.title = ggplot2::element_text(size = 20, face = "bold"),
                     axis.text = ggplot2::element_text(size = 9, face = "bold", angle = 90),
                     panel.grid.major = ggplot2::element_line(colour="grey", size=0.15),
                     panel.background = ggplot2::element_rect(fill = "white", color = "black")) +
      geom_hline(aes(yintercept = 0), color = "red") + ylim(c(-1,1)) 
    # scale_x_continuous(breaks = subsetSpecies, name = "Species",
    # labels = colnames(OTU)[subsetSpecies]) +
    # scale_y_continuous(name = "Elevation") +
    
  }
  
  # nonzero coefficients
  {
    nonzeroCoeffs <- which(Tau_true != 0, arr.ind = T)
    nonzeroCoeffs <- nonzeroCoeffs[nonzeroCoeffs[,1] > nonzeroCoeffs[,2],]
    nonzerovals <- mapply(function(i, j) Tau_true[i,j], 
                          nonzeroCoeffs[,1], nonzeroCoeffs[,2])
    
    # Tau_output_nonzeros <- mapply(function(i, j) Tau_corr_output[,,i,j], 
    Tau_output_nonzeros <- mapply(function(i, j) Tau_output[,,i,j], 
                                  nonzeroCoeffs[,1], nonzeroCoeffs[,2])
    CI_taus <- apply(Tau_output_nonzeros, 2, function(x){
      quantile(x, probs = c(0.025,0.5,0.975))
    })
    
    ggplot(data = NULL, aes(x = 1:nrow(nonzeroCoeffs),
                            y = CI_taus[2,],
                            ymin = CI_taus[1,],
                            ymax = CI_taus[3,])) + geom_errorbar() + 
      geom_point(data = NULL, aes(x = 1:nrow(nonzeroCoeffs),
                                  y = nonzerovals), color = "red",
                 size = 2) + 
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                     axis.title = ggplot2::element_text(size = 20, face = "bold"),
                     axis.text = ggplot2::element_text(size = 9, face = "bold", angle = 90),
                     panel.grid.major = ggplot2::element_line(colour="grey", size=0.15),
                     panel.background = ggplot2::element_rect(fill = "white", color = "black")) +
      ylim(c(-1,1)) 
    # scale_x_continuous(breaks = subsetSpecies, name = "Species",
    # labels = colnames(OTU)[subsetSpecies]) +
    # scale_y_continuous(name = "Elevation") +
    # ylim(c(-3,3))
  }
  
}

# diagnostics
{
  lambda_output_iter[2,(iter - nburn) - 10:1]
  mu_output_iter[2,(iter - nburn) - 10:1]
  logz_output_iter[1,2,(iter - nburn) - 10:1]
  v_output_iter[2,2,(iter - nburn) - 10:1]
  u_output_iter[2,2,(iter - nburn) - 10:1]
}

# DIAGNOSTICS -------------------------------------------------------------

# beta z

j <- 3
qplot(1:niter, beta_z_output[2,j,], geom = "line") + geom_hline(aes(yintercept = beta_z_true[2,j]), color = "red") + 
  ylim(c(-2,2))
qplot(1:niter, beta_z_output[3,j,], geom = "line") + geom_hline(aes(yintercept = beta_z_true[3,j]), color = "red") + 
  ylim(c(-2,2)) 
j <- j+1

# beta p11

qplot(1:niter, beta_p11_output[1,], geom = "line") + geom_hline(aes(yintercept = beta_p11_true[1]), color = "red") + 
  ylim(c(-2,2))
qplot(1:niter, beta_p11_output[2,], geom = "line") + geom_hline(aes(yintercept = beta_p11_true[2]), color = "red") + 
  ylim(c(-2,2))

# logz
i <- 47
m <- 1
m + sum(M_site[seq_len(i-1)])
j <- 2
y[m + sum(M_site[seq_len(i-1)]), ,]
round(w_true[m + sum(M_site[seq_len(i-1)]),], 3)
u[m + sum(M_site[seq_len(i-1)]),]; u_true[m + sum(M_site[seq_len(i-1)]),] 
c_imk_true[m + sum(M_site[seq_len(i-1)]), ,j]

qplot(1:niter, logz_output[i,j,], geom = "line") + geom_hline(aes(yintercept = logz_true[i,j])) + 
  ylim(c(-5,5))

# v
qplot(1:niter, v_output[m + sum(M_site[seq_len(i-1)]),j,], geom = "line") + 
  geom_hline(aes(yintercept = v_true[m + sum(M_site[seq_len(i-1)]),j]))# + ylim(c(-5,5))

# p_10
j <- 2
qplot(1:niter, p_10_output[j,], geom = "line") + 
  geom_hline(aes(yintercept = p10_true[j])) 

# beta 0




j <- 4
qplot(1:niter, beta_theta11_output[j,1,], geom = "line") + 
  geom_hline(aes(yintercept = beta_theta11_true[j,1]), color = "red") + 
  ylim(c(-3,3))
qplot(1:niter, beta_theta11_output[j,2,], geom = "line") + 
  geom_hline(aes(yintercept = beta_theta11_true[j,2]), color = "red") + 
  ylim(c(-2,2))




# lambda0 and lambdatilde
qplot(1:niter, lambda_output, geom = "line")
qplot(1:niter, lambda0_output, geom = "line")

# beta0
qplot(1:niter, beta_z_output[1,j,], geom = "line") + geom_hline(aes(yintercept = beta_z_true[1,j])) + 
  ylim(c(-2,2)) 

qplot(1:niter, beta_z_output[1,j,] + eta_output[j,], geom = "line") + 
  geom_hline(aes(yintercept = beta_z_true[1,j] + eta_true[j])) + 
  ylim(c(-2,2)) 




qplot(1:niter, beta_p11_output[1,], geom = "line") + geom_hline(aes(yintercept = beta_p11_true[1])) + 
  ylim(c(-2,2))
qplot(1:niter, beta_p11_output[2,], geom = "line") + geom_hline(aes(yintercept = beta_p11_true[2])) + 
  ylim(c(-2,2))
qplot(1:niter, beta_p11_output[3,], geom = "line") + geom_hline(aes(yintercept = beta_p11_true[3])) + 
  ylim(c(-2,2))

qplot(1:niter, sigma_output[3,], geom = "line")
qplot(1:niter, tau_output[3,], geom = "line")


qplot(1:niter, lambdatilde_output, geom = "line")



qplot(1:niter, lambda_output, geom = "line") + geom_hline(aes(yintercept = lambda_true))
qplot(1:niter, u_output[1,3,], geom = "line")
qplot(1:niter, lambda_output * u_output[1,3,], geom = "line")

# DIAGNOSTICS TRUE MODEL -------------------------------------------------------------

library(ggplot2)
lambda_output <- modelResults$params_output$lambda_output
beta_z_output <- modelResults$params_output$beta_z_output


# lambda

j <- 1
diagnosticsPlot(lambda_output[,j,]) + ylim(c(-1,10))
j <- j + 1

# beta z

j <- 1
cov_num <- 1
diagnosticsPlot(beta_z_output[,cov_num,j,]) +  ylim(c(-3,3))
j <- j+1

# alpha

j <- 1
diagnosticsPlot( alpha_output[,j,]) + 
  ylim(c(-2,2))
j <- j + 1


# logz
i <- 1
j <- 1
diagnosticsPlot( logz_output[,i,j,])  + ylim(c(-10,10))
# i <- i+1
j <- j+1

# tau
j <- 1
diagnosticsPlot( tau_output[,j,]) + ylim(c(0,10))
j <- j + 1

# sigma
j <- 1
diagnosticsPlot( sigma_output[,j,]) + ylim(c(0,10))
j <- j + 1

# r
j <- 1
diagnosticsPlot( r_output[,j,]) + ylim(c(0,10))
j <- j + 1

# v
qplot(1:niter, v_output[m + sum(M_site[seq_len(i-1)]),j,], geom = "line")#+ ylim(c(-1.5,1.5))

# p_10
j <- 2
qplot(1:niter, p_10_output[j,], geom = "line") #+ 
geom_hline(aes(yintercept = p10_true[j])) 

# beta theta 1
j <- 1
diagnosticsPlot( beta_theta_output[,j,1,]) + 
  ylim(c(-6,-1.5))
j <- j + 1

# beta theta 2
j <- 1
diagnosticsPlot( beta_theta_output[,j,2,]) #+
# ylim(c(0,3))
j <- j + 1


# MORE OUTPUT ------

j <- 1
i <- 1
m <- 2
k <- 1
mean_mu <- exp(lambda[j] + 
                 v[m + sum(M_site[seq_len(i-1)]),j] + 
                 u[m + sum(M_site[seq_len(i-1)]),k])

# sample

beta0_mean <- apply(beta0_output, 2, mean)
tau_mean <- apply(tau_output, 2, mean)
sigma_mean <- apply(sigma_output, 2, mean)
lambda_mean <- apply(lambda_output, 2, mean)

beta0_mean <- beta0_true
tau_mean <- tau_true
sigma_mean <- sigma_true
lambda_mean <- lambda_true

j <- 1
nsamples <- 100000
sampledRead <- sapply(1:nsamples, function(i){
  logz_current <- rnorm(1, beta0_mean[j], sd = tau_mean[j])
  v_current <- rnorm(1, logz_current, sd = sigma_mean[j])
  currentMean <- exp(lambda_mean[j] + v_current)
  # currentMean
  rnbinom(1, size = r_nb[j], mu = currentMean)
})

ggplot(data = NULL, aes(x = sampledRead[sampledRead < 1000],
                        y = ..density..)) + 
  geom_histogram(color = "black", fill = "cornsilk", 
                 bins = 40)

ggplot(data = NULL, aes(x = log(sampledRead[sampledRead < 5000] + .01), 
                        y = ..density..)) + 
  geom_histogram(color = "black", fill = "cornsilk", 
                 bins = 40)

#
ggplot() + geom_histogram(data = NULL, 
                          aes(x = log(as.vector(y[,,j])[as.vector(y[,,j]) > 10] + .01), 
                              y = ..density..),
                          color = "black", fill = "green")
j <- j + 1

#

ggplot() + geom_histogram(data = NULL, 
                          aes(x = log(as.vector(y[,,j])[as.vector(y[,,j]) < 5000 & 
                                                          as.vector(y[,,j]) > 0] + .01), 
                              y = ..density..),
                          color = "black", fill = "green") + 
  geom_histogram(data = NULL, 
                 aes(x = log(sampledRead[sampledRead < 5000] + .01), 
                     y = ..density..),
                 color = "black", fill = "red")

#

sort(apply(y, 3, function(x){
  sum(x[,1] > 20)
}))
