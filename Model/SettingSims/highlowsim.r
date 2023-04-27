# CPP ---------------------------------------------------------------------

# library(Rcpp); library(RcppArmadillo); library(Matrix)
# library(GIGrvg); library(MASS); library(ggplot2)
# 
# library(GeneralizedHyperbolic); library(ars); library(lamW); library(clusterGeneration)
# library(ggplot2); library(MCMCpack); library(extraDistr); 
# library(Matrix); library(matrixsampling); 
# library(abind); library(beepr)

library(here)
# library(eDNAPlus)

setwd("/cluster/home/osr/ad625/eDNA/eDNAPlus/eDNAplus")
sourceCpp("src/coeff.cpp")
sourceCpp("src/indic_variables.cpp")
# sourceCpp(here("Model/Functions","variable_update.cpp"))
source("R/functions.r")
source("R/runmodel.R")

# data settings -----

# simulation params
{
  beta0_mean <- 0
  beta_theta_0_mean <- -1.5
  sigmas <- .5
  taus <- .5
  r_0 <- 100
  p11_0 <- .95
  p10_0 <- .02
  theta_10 <- .02
  
  sd_beta_theta_0 <- .001
  sigma_beta0 <- 1
  sigma_u <- 1
  sigma_gamma <- 1
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

jointSpecies <- F
spatialCorr <- F

nsims <- 10

ldiff_accuracy <- matrix(NA, nrow(nMK_grid), nsims)


simsSettings <- list("tau" = tau_true[1],
                     "sigma" = sigma_true[1],
                     "beta_theta" = beta_theta_0_mean,
                     "M" = M_site[1],
                     "n" = n,
                     "K" = K[1])

M_grid <- c(1, 2, 3, 4, 5)
K_grid <- c(1, 2, 3, 4)

nMK_grid <- expand.grid(M_grid, K_grid)


for (sim in 1:nsims) {
  
  # DATA --------
  
  {
    simulatedData <- T
    # settings 
    {
      S <- 40
      n <- 300#500
      ncov_z <- 0
      ncov_w <- 0
      
      M_site <- rep(5, n)
      emptyTubes <- 0
      
      K <- rep(5, sum(M_site) + emptyTubes)
      
      S_star <- 0 # numOfSpikes
    }
    
    
    
    # design
    {
      
      data_short <- data.frame(Site = rep(1:n, M_site),
                               Sample = as.vector((sapply(M_site, function(x){
                                 seq_len(x)
                               }))))
      
      if(emptyTubes > 0){
        data_short <- rbind(data_short,
                            data.frame(Site = 0,
                                       Sample = 1:emptyTubes))
        
      }
      
      X_z_0 <- matrix(rnorm(n * ncov_z), n, ncov_z)
      if(ncov_z > 0) colnames(X_z_0) <- paste("Z",1:ncov_z)
      
      X_w_0 <- matrix(rnorm(sum(M_site) * ncov_w), sum(M_site), ncov_w)
      if(ncov_w > 0) colnames(X_w_0) <- paste("W",1:ncov_w)
      
      v_spikes0 <- matrix(0, nrow = sum(M_site) + emptyTubes, ncol = S_star)
      spikedSample0 <- rbind(matrix(1, nrow = sum(M_site), ncol = S_star),
                             matrix(rbinom(emptyTubes * S_star, 1, .5), nrow = emptyTubes, ncol = S_star))
      
    }
    
    beta0_true <- rep(0, S)#rnorm(S, beta0_mean, sd = sigma_beta0)#
    
    beta_z_true <- matrix(sample(c(1,-1), (ncov_z) * S, replace = T), 
                          nrow = ncov_z, ncol = S, byrow = T)
    
    sigma_true <- rep(sigmas, S)#pmin(rhcauchy(S, 2), .1)
    
    if(jointSpecies) {
      
      sizeBlocks <- 3
      
      Tau_list <- lapply(1:(S/sizeBlocks), function(j){
        cormat <- matrix(1, sizeBlocks, sizeBlocks)
        repeat {
          for (i in 2:sizeBlocks) {
            for (j in seq_len(i-1)) {
              cormat[i,j] <- sample(c(-1,1) * .8, size = 1)
              cormat[j,i] <- cormat[i,j]
            }
          }
          if(all(eigen(cormat)$values > 0)){
            break
          }
        }
        
        cormat <- cormat * taus
      })
      Tau_true <- as.matrix(bdiag(Tau_list))
      # Tau_true <- diag(taus, nrow = S)
    } else {
      tau_true <- rep(taus, S)#pmin(rhcauchy(S, .5), 1)
    }
    
    beta_w_true <- matrix(sample(c(1,-1), (ncov_w) * S, replace = T), 
                          nrow = ncov_w, ncol = S, byrow = T)
    
    beta_theta_true <- cbind(
      rnorm(S, mean = beta_theta_0_mean, sd = sd_beta_theta_0),  # baseline
      rep(1, S), # DNA coefficient
      t(matrix(sample(c(1,-1), (ncov_w) * S, replace = T), 
               nrow = ncov_w, ncol = S, byrow = T))) # covariate coefficient
    
    theta10_true <- rep(theta_10, S)
    
    lambda_true <- rnorm(S + S_star, 9, sd = 1)
    r_true <- rep(r_0, S + S_star)#pmax(rgamma(S, 1, .2), 10)
    # lambdatilde_true <- rep(.3, S + S_star)
    
    mu_true <- rnorm(S, sd = 1)
    
    mu0_true <- 5
    n0_true <- 10
    pi0_true <- .9
    
    mu_tilde_true <- 100
    n_tilde_true <- 100000
    
    p11_true <- rep(p11_0, S + S_star)#rbeta(S, a_p11, b_p11)
    p10_true <- rep(p10_0, S + S_star)
    
    # simulation
    highLowEffect <- c(rep(-.5, n / 2), rep(.5, n / 2))
    if(jointSpecies){
      logz_true <- matrix(NA, n, S)
      for (i in 1:n) {
        Xb_i <-  c(1, X_z_0[i,]) %*% rbind(beta0_true, beta_z_true)
        logz_true[i,] <- mvrnorm(1, Xb_i, Tau_true)
        
      }
    } else {
      # logz_true <- X_z %*% beta_z_true +
      logz_true <- cbind(1, X_z_0) %*% rbind(beta0_true, beta_z_true) +
        sapply(1:S, function(j) rnorm(n, 0, sd = tau_true[j])) +
        matrix(highLowEffect, n, S)
      
    }
    
    theta11_true <- matrix(NA, nrow = sum(M_site), ncol = S)
    for (j in 1:S) {
      X_wbeta_theta <- cbind(1, rep(logz_true[,j], M_site), X_w_0) %*% beta_theta_true[j,]
      # X_wbeta_theta <- cbind(1, rep(exp(logz_true[,j]), M_site), X_w) %*% beta_theta_true[j,]
      theta11_true[,j] <- logistic(X_wbeta_theta)
    }
    
    delta_true <- array(NA, dim = c(sum(M_site) + emptyTubes, S + S_star))
    for (i in 1:n) {
      for (m in 1:M_site[i]) {
        for (j in 1:S) {
          delta_true[m + sum(M_site[seq_len(i-1)]),j] <- 
            rbinom(1, 1, theta11_true[m + sum(M_site[seq_len(i-1)]),j])
        }
        for (j in seq_len(S_star)) {
          delta_true[m + sum(M_site[seq_len(i-1)]),S + j] <- spikedSample[m + sum(M_site[seq_len(i-1)]), j]  
        }
      }
    }
    # empty tubes
    for (m in seq_len(emptyTubes)) {
      delta_true[sum(M_site) + m,] <- 0
      for(j in seq_len(S_star)){
        delta_true[sum(M_site) + m, S + j] <- spikedSample[sum(M_site) + m, j]
      }
    }
    
    gamma_true <- array(NA, dim = c(sum(M_site) + emptyTubes, S + S_star))
    for (i in 1:n) {
      for (m in 1:M_site[i]) {
        for (j in 1:S) {
          if(delta_true[m + sum(M_site[seq_len(i-1)]),j] == 0){
            gamma_true[m + sum(M_site[seq_len(i-1)]),j] <- rbinom(1, 1, theta10_true[j])  
          } else if(delta_true[m + sum(M_site[seq_len(i-1)]),j] == 1){
            gamma_true[m + sum(M_site[seq_len(i-1)]),j] <- 0
          }
        }
        for (j in seq_len(S_star)) {
          gamma_true[m + sum(M_site[seq_len(i-1)]),S + j] <- 0
        }
      }
    }
    # empty tubes
    for (m in seq_len(emptyTubes)) {
      gamma_true[sum(M_site) + m,] <- 0
    }
    
    Xw_betaw <- X_w_0 %*% beta_w_true
    
    # log amount of DNA
    v_true <- matrix(NA, sum(M_site) + emptyTubes, S + S_star)
    for (i in 1:n) {
      for (m in 1:M_site[i]) {
        for (j in 1:S) {
          if(delta_true[m + sum(M_site[seq_len(i-1)]),j] == 1){
            v_true[m + sum(M_site[seq_len(i-1)]),j] <- 
              rnorm(1, logz_true[i,j] +
                      Xw_betaw[m + sum(M_site[seq_len(i-1)]),j], 
                    sigma_true[j]) 
          } else if (gamma_true[m + sum(M_site[seq_len(i-1)]),j] == 1){
            v_true[m + sum(M_site[seq_len(i-1)]),j] <- 
              rnorm(1, mu_true[j], 
                    sd = sigma_gamma)
          } 
        }
        for (j in seq_len(S_star)) {
          v_true[m + sum(M_site[seq_len(i-1)]), S + j] <- 
            v_spikes[m + sum(M_site[seq_len(i-1)]),j]
        }
      }
    }
    for (m in seq_len(emptyTubes)) {
      for (j in seq_len(S_star)) {
        v_true[sum(M_site) + m, S + j] <- 
          v_spikes[sum(M_site) + m,j]
      }
    }
    
    # DNA detections
    c_imk_true <- array(NA, dim = c(sum(M_site) + emptyTubes, max(K), S + S_star))
    for (i in 1:n) {
      for (m in 1:M_site[i]) {
        numRep <- K[m + sum(M_site[seq_len(i-1)])]
        for (k in 1:numRep) {
          for (j in 1:(S+S_star)) {
            if(delta_true[m + sum(M_site[seq_len(i-1)]),j] == 1 | 
               gamma_true[m + sum(M_site[seq_len(i-1)]),j] == 1){
              c_imk_true[m + sum(M_site[seq_len(i-1)]),k,j] <- rbinom(1, 1, p11_true[j])
            } else {
              c_imk_true[m + sum(M_site[seq_len(i-1)]),k,j] <- 2 * rbinom(1, 1, p10_true[j])
            }
          }
        }
      }
    }
    # empty tubes
    for (m in seq_len(emptyTubes)) {
      numRep <- K[sum(M_site) + m]
      for (k in 1:numRep) {
        for (j in 1:(S+S_star)) {
          if(delta_true[m + sum(M_site),j] == 1 | 
             gamma_true[m + sum(M_site),j] == 1){
            c_imk_true[sum(M_site) + m,k,j] <- rbinom(1, 1, p11_true[j])
          } else {
            c_imk_true[sum(M_site) + m,k,j] <- 2 * rbinom(1, 1, p10_true[j])
          }
        }
        # for (j in 1:(S+S_star)) {
        #   c_imk_true[sum(M_site) + m,k,j] <- 2 * rbinom(1, 1, p10_true[j])
        # }
      }
    }
    
    # PCR noise
    u_true <- matrix(NA, sum(M_site) + emptyTubes, max(K))
    for (i in 1:n) {
      for (m in 1:M_site[i]) {
        numRep <- K[m + sum(M_site[seq_len(i-1)])]
        for (k in 1:numRep) {
          u_true[m + sum(M_site[seq_len(i-1)]),k] <- rnorm(1,0, sigma_u)
        }
      }
    }
    # empty tubes
    for (m in seq_len(emptyTubes)) {
      numRep <- K[sum(M_site) + m]
      for (k in 1:numRep) {
        u_true[sum(M_site) + m,k] <- rnorm(1, 0, sigma_u)
      }
    }
    u_true <- u_true - mean(u_true)
    
    OTUnames <- paste0("OTU_",1:(S+S_star))
    
    y <- array(NA, dim = c(sum(M_site) + emptyTubes, max(K), S + S_star),
               dimnames = list(NULL, NULL, OTUnames))
    for (i in 1:n) {
      # print(i)
      for (m in 1:M_site[i]) {
        numRep <- K[m + sum(M_site[seq_len(i-1)])]
        for (k in 1:numRep) {
          u_imk <- u_true[m + sum(M_site[seq_len(i-1)]), k]
          for (j in 1:(S + S_star)) {
            if (c_imk_true[m + sum(M_site[seq_len(i - 1)]), k, j] == 0) {
              y[m + sum(M_site[seq_len(i - 1)]), k, j] <- 
                ifelse(runif(1) < pi0_true, 0, 1 + rnbinom(1, mu = mu0_true, size = n0_true))
            } else if (c_imk_true[m + sum(M_site[seq_len(i - 1)]), k, j] == 2) {
              y[m + sum(M_site[seq_len(i - 1)]), k, j] <-
                # runif(1, 0, exp(lambda_true[j]))
                rnbinom(1, mu = mu_tilde_true, size = n_tilde_true)
              # rpois(1, lambdatilde_true[j] * exp(lambda_true[j]))
            } else {
              mean_true <- exp(v_true[m + sum(M_site[seq_len(i - 1)]), j] + 
                                 lambda_true[j] + 
                                 u_imk)
              # rpois(1, mean_true)
              y[m + sum(M_site[seq_len(i - 1)]), k, j] <-
                rnbinom(1, size = r_true[j], mu = mean_true)
            }
          }
        }
      }
    }
    # empty tubes
    for (m in seq_len(emptyTubes)) {
      numRep <- K[sum(M_site) + m]
      for (k in 1:numRep) {
        u_imk <- u_true[sum(M_site) + m, k]
        for (j in 1:(S + S_star)) {
          if (c_imk_true[sum(M_site) + m, k, j] == 0) {
            y[sum(M_site) + m, k, j] <- 
              ifelse(runif(1) < pi0_true, 0, 1 + rnbinom(1, mu = mu0_true, size = n0_true))
          } else if (c_imk_true[sum(M_site) + m, k, j] == 2) {
            y[sum(M_site) + m, k, j] <-
              # rpois(1, lambdatilde_true[j] * exp(lambda_true[j]))
              rnbinom(1, mu = mu_tilde_true, size = n_tilde_true)
          } else {
            mean_true <- exp(v_true[sum(M_site) + m, j] + 
                               lambda_true[j] + 
                               u_imk)
            y[sum(M_site) + m, k, j] <-
              rnbinom(1, size = r_true[j], mu = mean_true)
          }
        }
      }
    }
    
    {
      y0 <- y
      c_imk_true_0 <- c_imk_true
      v_true0 <- v_true
      gamma_true0 <- gamma_true
      delta_true0 <- delta_true
      p11_true0 <- p11_true
      p10_true0 <- p10_true
      lambda_true0 <- lambda_true
      r_true0 <- r_true
      mutilde_true0 <- mu_tilde_true
      ntilde_true0 <- n_tilde_true
    }
    
  }
  
  
  # PARAMS -----
  
  
  
  # mcmc params
  {
    MCMCparams <- list("nchain" = 1,
                       "nburn" = 5000,
                       "niter" = 5000,
                       "nthin" = 1,
                       "iterToAdapt" = Inf)
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
      updateAll = T,
      params = NULL,
      correct = NULL
    )
    
  }
  
  
  # SIMS -----------
  
  for (setting_idx in 1:nrow(nMK_grid)) {
    # for (setting_idx in 21:nrow(nMK_grid)) {
    
    print(paste0("Sim = ",sim," - Setting idx = ",setting_idx))
    
    n_start <- nrow(X_z_0)
    n <- n_start
    M_current <- nMK_grid[setting_idx, 1]
    K_current <- nMK_grid[setting_idx, 2]
    
    # cut data
    {
      idx_sites <- c(seq_len(n / 2), n_start / 2 + seq_len(n / 2))
      X_z <- X_z_0[idx_sites,,drop = F]
      M_site <- rep(M_current, n)
      K <- rep(K_current, sum(M_site))
      idx_obs <- which(data_short$Site %in% idx_sites & data_short$Sample <= M_current)
      y <- y0[idx_obs,seq_len(K[1]),,drop = F]
      X_w <- X_w_0[idx_obs,,drop = F]
      spikedSample <- spikedSample0[idx_obs,]
      v_spikes <- v_spikes0[idx_obs,]
    }
    
    # run mcmc
    {
      data <- list("y" = y,
                   "M_site" = M_site,
                   "K" = K,
                   "emptyTubes" = emptyTubes,
                   "S_star" = S_star,
                   "spikedSample" = spikedSample,
                   "v_spikes" = v_spikes,
                   "X_w" = X_w,
                   "X_z" = X_z)
      
      
      modelResults <- fitModel(data,
                               priors,
                               jointSpecies,
                               spatialCorr,
                               paramsUpdate,
                               MCMCparams,
                               idConstraints)
      
    }
    
    # save output
    {
      # summarise results
      {
        # u_mean <- apply(u_output, c(2,3), mean)
        # (bias <- mean(abs(u_mean - u_true)))
        # u_biases[s_idx] <- bias
        # 
        # u_vars[s_idx] <- mean(apply(u_output, c(2,3), sd))
        
        #
        
        # beta_z_output <- modelResults$params_output$beta_z_output
        # 
        # beta_mean <- apply(beta_z_output, c(3,4), mean)
        # (bias <- mean(abs(beta_mean - beta_z_true)))
        # beta_biases[setting_idx] <- bias
        # 
        # beta_vars[setting_idx] <- mean(apply(beta_z_output, c(3,4), sd))
        # 
        #
        
        logz_output <- modelResults$params_output$logz_output
        
        l_sign <- array(NA, dim = c(S, (n / 2) * (n / 2)))
        
        for (j in 1:S) {
          # print(j)
          l <- 1
          for (i1 in seq_len(n / 2)) {
            
            for (i2 in (n/2) + seq_len(n / 2)) {
              
              trueSign <- as.numeric(logz_true[i1,j] > logz_true[i2,j])
              estimatedSign <- mean(logz_output[,,i1,j] > logz_output[,,i2,j])
              
              # y[i1,,j]
              # y[i2,,j]
              # ggplot() + geom_point(data = NULL, aes(x = 1:5000, 
              #                                        y = logz_output[1,,i1,j])) +
              #   geom_point(data = NULL, aes(x = 1:5000, 
              #                               y = logz_output[1,,i2,j]), color = "red")
              
              # ggplot() + geom_point(data = NULL, aes(x = 1:5000, 
              # y = logz_output[1,,i1,j] - logz_output[1,,i2,j])) #+
              # geom_point(data = NULL, aes(x = 1:5000, 
              # y = logz_output[1,,i2,j]), color = "red")
              
              l_sign[j,l] <- (trueSign - estimatedSign)^2
              
              l <- l + 1
            }
            
          }
        }
        
        ldiff_accuracy[setting_idx,sim] <- mean(l_sign)
        
      }
    }
    
  }
  
  
  
}


# RESULTS -----

ldiff_accuracy <- c(0.236428961861022, 0.225282255359467, 0.210897160256622, 0.199922942099956, 
                    0.189949491251378, 0.233306003169956, 0.219025160809378, 0.204190063294756, 
                    0.191989378721511, 0.181179756233244, 0.230729007304444, 0.215981208905778, 
                    0.201415897926667, 0.187330608709956, 0.175728344457778, 0.2286606804628, 
                    0.212785329137111, 0.196833056697022, 0.185326364830711, 0.171935724184
)

settings_df <- as.data.frame(nMK_grid)
colnames(settings_df) <- c("M","K")

plot_data <- cbind(settings_df, data.frame(vals = ldiff_accuracy))

plot_data$K <- as.factor(plot_data$K)
# plot_data$M <- as.factor(plot_data$M)

var_settings <- paste0(
  plot_data$M," - ",plot_data$K
)

plot_data <- cbind(plot_data, var_settings)

# plot_data <- plot_data[plot_data$n == max(plot_data$n),]

ggplot(data = 
         plot_data,
       aes(x = M, 
           y = vals,
           group = K,
           linetype = K)) + 
  geom_line() + geom_point() + #facet_grid(cols = vars(K)) +
  # coord_cartesian(ylim = c(0, max(plot_data$vals))) +
  scale_x_continuous(name = "M", breaks = unique(plot_data$M)) + 
  scale_y_continuous(name = "Brier score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
        axis.line = element_line(colour="black", size=0.15),
        # panel.grid.minor = element_line(colour="grey", size=0.15),
        panel.grid.major = element_line(colour="grey", size=0.15),
        panel.background = element_rect(fill = "white", color = "black")) + 
  ggtitle(title)
# PLOTS ----------------

# beta bias
ggplot(data = NULL, aes(x = S_star_grid, 
                        y = beta_biases)) + 
  geom_line() + geom_point() + 
  coord_cartesian(ylim = c(0, max(beta_biases))) + 
  scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
        axis.line = element_line(colour="black", size=0.15),
        # panel.grid.minor = element_line(colour="grey", size=0.15),
        panel.grid.major = element_line(colour="grey", size=0.15),
        panel.background = element_rect(fill = "white", color = "black"))

# beta vars
ggplot(data = NULL, aes(x = S_star_grid, 
                        y = beta_vars)) + 
  geom_line() + geom_point() + 
  coord_cartesian(ylim = c(0, max(beta_vars))) + 
  scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
        axis.line = element_line(colour="black", size=0.15),
        # panel.grid.minor = element_line(colour="grey", size=0.15),
        panel.grid.major = element_line(colour="grey", size=0.15),
        panel.background = element_rect(fill = "white", color = "black"))

# l diff biases
ggplot(data = NULL, aes(x = S_star_grid, 
                        y = ldiff_biases)) + 
  geom_line() + geom_point() + 
  coord_cartesian(ylim = c(0, max(ldiff_biases))) + 
  scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
        axis.line = element_line(colour="black", size=0.15),
        # panel.grid.minor = element_line(colour="grey", size=0.15),
        panel.grid.major = element_line(colour="grey", size=0.15),
        panel.background = element_rect(fill = "white", color = "black"))

# l diff vars
ggplot(data = NULL, aes(x = S_star_grid, 
                        y = ldiff_vars)) + 
  geom_line() + geom_point() + 
  coord_cartesian(ylim = c(0, max(ldiff_vars))) + 
  scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
        axis.line = element_line(colour="black", size=0.15),
        # panel.grid.minor = element_line(colour="grey", size=0.15),
        panel.grid.major = element_line(colour="grey", size=0.15),
        panel.background = element_rect(fill = "white", color = "black"))

# u biases
ggplot(data = NULL, aes(x = S_star_grid, 
                        y = u_biases)) + 
  geom_line() + geom_point() + 
  coord_cartesian(ylim = c(0, max(u_biases))) + 
  scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
        axis.line = element_line(colour="black", size=0.15),
        # panel.grid.minor = element_line(colour="grey", size=0.15),
        panel.grid.major = element_line(colour="grey", size=0.15),
        panel.background = element_rect(fill = "white", color = "black"))

# u vars
ggplot(data = NULL, aes(x = S_star_grid, 
                        y = u_vars)) + 
  geom_line() + geom_point() + 
  coord_cartesian(ylim = c(0, max(u_vars))) + 
  scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
        axis.line = element_line(colour="black", size=0.15),
        # panel.grid.minor = element_line(colour="grey", size=0.15),
        panel.grid.major = element_line(colour="grey", size=0.15),
        panel.background = element_rect(fill = "white", color = "black"))

setwd("/cluster/home/osr/ad625/eDNA/Project/Model/Spike-in Sims")
save(simsSettings,
     beta_biases, beta_vars, ldiff_biases, ldiff_vars, u_biases, u_vars, 
     file = paste0("Results_",
                   "tau=",simsSettings$tau,"_",
                   "sigma=",simsSettings$sigma,"_",
                   "sigma_u=",simsSettings$sigma_u,"_",
                   "M=",simsSettings$M,"_",
                   "n=",simsSettings$n,"_",
                   "K=",simsSettings$K,
                   ".rda"))


