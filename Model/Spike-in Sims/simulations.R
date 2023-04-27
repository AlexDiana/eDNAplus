# CPP ---------------------------------------------------------------------

library(Rcpp); library(RcppArmadillo); library(MASS)
library(GeneralizedHyperbolic); library(ars); library(lamW); library(clusterGeneration)
library(ggplot2); library(MCMCpack); library(extraDistr); library(GIGrvg)#library(nleqslv)
library(Matrix); library(matrixsampling); library(here); library(abind); library(beepr)

sourceCpp(here("Model/Functions","coeff.cpp"))
sourceCpp(here("Model/Functions","indic_variables.cpp"))
# sourceCpp(here("Model/Functions","variable_update.cpp"))
source(here("Model/Functions","functions.r"))
# R CMD INSTALL --preclean --no-multiarch --with-keep.source "functionsForForeach"

jSDM <- F
betaThetaEqual1 <- T

# DATA --------

# load(here("Data","Leeches.Rdata"))
load(here("Data","Ponds.Rdata"))

subsetSpecies <- 1:10

S <- length(subsetSpecies)
OTU <- OTU[,subsetSpecies]
y <- y[,,subsetSpecies]

S_star <- 0
emptyTubes <- 0

# n <- 50
# X_z <- X_z[1:n,]
# M_site <- M_site[1:n]
# y <- y[1:sum(M_site),,]
# r <- r[1:sum(M_site),]
# OTU <- OTU[1:sum(M_site),]
# K <- K[1:sum(M_site)]
# emptySites <- emptySites[1:n]
# data_short <- data_short[1:sum(M_site),]

# simulate base data  ----------

{
  
  # SIMULATE DATA -----------------------------------------------------------
  
  # design
  {
    n <- 1000
    ncov_z <- 1
    ncov_w <- 0
    
    M_site <- rep(1, n)
    emptyTubes <- 0
    
    K <- rep(3, sum(M_site) + emptyTubes)
    
    PCR_spiked <- rep(T, sum(K))
    S_star <- 10 # numOfSpikes
    
    data_short <- data.frame(Site = rep(1:n, M_site),
                             Sample = 1:sum(M_site))
    
    if(emptyTubes > 0){
      data_short <- cbind(data_short,
                          data.frame(Site = 0,
                                     Sample = 1:emptyTubes))
      
    }
    
    X_z <- matrix(runif(n * ncov_z), n, ncov_z)
    X_w <- matrix(runif(sum(M_site) * ncov_w), sum(M_site), ncov_w)
    v_spikes <- matrix(0, nrow = sum(M_site), ncol = S_star)
    
  }
  
  # prior
  {
    beta0_mean <- 0
    beta_theta_0_mean <- 5
    sigmas <- 1
    taus <- 1
    r_0 <- 1000
    p11_0 <- 1
    p10_0 <- .1
    theta_10 <- .0
    
    sd_beta_theta_0 <- .001
    sigma_beta0 <- .001
    sigma_u <- 1
  }
  
  simulatedData <- T
  
  beta0_true <- rnorm(S, beta0_mean, sd = sigma_beta0)#
  
  beta_z_true <- matrix(sample(c(1,-1), (ncov_z) * S, replace = T), 
                        nrow = ncov_z, ncol = S, byrow = T)
  
  sigma_true <- rep(sigmas, S)#pmin(rhcauchy(S, 2), .1)
  
  if(jSDM) {
    
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
  } else {
    tau_true <- rep(taus, S)#pmin(rhcauchy(S, .5), 1)
  }
  
  beta_w_true <- matrix(sample(c(1,-1), (ncov_w) * S, replace = T), 
                        nrow = ncov_w, ncol = S, byrow = T)
  
  beta_theta_true <- cbind(rnorm(S, mean = beta_theta_0_mean, sd = sd_beta_theta_0),  # baseline
                           rep(1, S), # DNA coefficient
                           t(matrix(sample(c(1,-1), (ncov_w) * S, replace = T), 
                                    nrow = ncov_w, ncol = S, byrow = T))) # covariate coefficient
  
  theta10_true <- rep(theta_10, S)
  
  lambda_true <- rnorm(S + S_star, log(1000), sd = .5)
  r_true <- rep(r_0, S + S_star)#pmax(rgamma(S, 1, .2), 10)
  lambdatilde_true <- rep(.3, S + S_star)
  
  mu_true <- rnorm(S, sd = 1)
  
  mu0_true <- 5
  n0_true <- 10
  pi0_true <- .9
  
  p11_true <- rep(p11_0, S + S_star)#rbeta(S, a_p11, b_p11)
  p10_true <- rep(p10_0, S + S_star)
  
  # simulation
  
  if(jSDM){
    logz_true <- matrix(NA, n, S)
    for (i in 1:n) {
      Xb_i <-  c(1, X_z[i,]) %*% rbind(beta0_true, beta_z_true)
      logz_true[i,] <- mvrnorm(1, Xb_i, Tau_true)
      
    }
  } else {
    # logz_true <- X_z %*% beta_z_true +
    logz_true <- cbind(1, X_z) %*% rbind(beta0_true, beta_z_true) +
      sapply(1:S, function(j) rnorm(n, 0, sd = tau_true[j]))
  }
  
  theta11_true <- matrix(NA, nrow = sum(M_site), ncol = S)
  for (j in 1:S) {
    X_wbeta_theta <- cbind(1, rep(exp(logz_true[,j]), M_site), X_w) %*% beta_theta_true[j,]
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
        if(PCR_spiked[m + sum(M_site[seq_len(i-1)])]){
          delta_true[m + sum(M_site[seq_len(i-1)]),S + j] <- 1  
        } else {
          delta_true[m + sum(M_site[seq_len(i-1)]),S + j] <- 0  
        }
      }
    }
  }
  # empty tubes
  for (m in seq_len(emptyTubes)) {
    delta_true[sum(M_site) + m,] <- 0
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
  
  Xw_betaw <- X_w %*% beta_w_true
  
  # log amount of DNA
  v_true <- matrix(NA, sum(M_site) + emptyTubes, S + S_star)
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        if(delta_true[m + sum(M_site[seq_len(i-1)]),j] == 1){
          v_true[m + sum(M_site[seq_len(i-1)]),j] <- 
            rnorm(1, logz_true[i,j], sigma_true[j]) +
            Xw_betaw[m + sum(M_site[seq_len(i-1)]),j] 
        } else if (gamma_true[m + sum(M_site[seq_len(i-1)]),j] == 1){
          v_true[m + sum(M_site[seq_len(i-1)]),j] <- 
            rnorm(1, mu_true[j], 
                  sd = 1)
        } 
        # else {
        #   v_true[m + sum(M_site[seq_len(i-1)]),j] <- 0
        # }
      }
      for (j in seq_len(S_star)) {
        v_true[m + sum(M_site[seq_len(i-1)]), S + j] <- 
          v_spikes[m + sum(M_site[seq_len(i-1)]),j]
      }
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
        c_imk_true[sum(M_site) + m,k,j] <- 2 * rbinom(1, 1, p10_true[j])
      }
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
  u_true <- u_true - mean(u_true)
  # empty tubes
  for (m in seq_len(emptyTubes)) {
    numRep <- K[sum(M_site) + m]
    for (k in 1:numRep) {
      u_true[sum(M_site) + m,k] <- rnorm(1, 0, sigma_u)
    }
  }
  
  y <- array(NA, dim = c(sum(M_site) + emptyTubes, max(K), S + S_star))
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
              rpois(1, lambdatilde_true[j] * exp(lambda_true[j]))
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
            rpois(1, lambdatilde_true[j] * exp(lambda_true[j]))
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
    lambdatilde_true0 <- lambdatilde_true
    v_spikes0 <- v_spikes
  }
  
}

# SIMS -----------

S_star_grid <- c(0, 1, 2, 5, 10)

beta_biases <- rep(NA, length(S_star_grid))
beta_vars <- rep(NA, length(S_star_grid))

ldiff_biases <- rep(NA, length(S_star_grid))
ldiff_vars <- rep(NA, length(S_star_grid))

u_biases <- rep(NA, length(S_star_grid))
u_vars <- rep(NA, length(S_star_grid))

for (s_idx in seq_along(S_star_grid)) {
  
  S_star <- S_star_grid[s_idx]
  
  # cut data
  {
    y <- y0[,,1:(S+S_star),drop = F]
    c_imk_true <- c_imk_true_0[,,1:(S+S_star),drop = F]
    v_true <- v_true0[,1:(S+S_star),drop = F]
    gamma_true <- gamma_true0[,1:(S+S_star)]
    delta_true <- delta_true0[,1:(S+S_star)]
    p11_true <- p11_true0[1:(S+S_star)]
    p10_true <- p10_true0[1:(S+S_star),drop = F]
    lambda_true <- lambda_true0[1:(S+S_star)]
    r_true <- r_true0[1:(S+S_star)]
    lambdatilde_true <- lambdatilde_true0[1:(S+S_star)]
    v_spikes <- v_spikes0[,seq_len(S_star),drop = F]
  }
  
  # run mcmc
  {
    # CLEAN DATA ---------
    
    ncov_z <- ncol(X_z)
    ncov_w <- ncol(X_w)
    
    # PRIOR -------------------------------------------------------------------
    
    a_theta0 <- 1
    b_theta0 <- 20
    
    a_p11 <- 1
    b_p11 <- 1
    
    a_p10 <- 1
    b_p10 <- 20
    
    a_tau <- 3
    b_tau <- 2
    # A_tau <- .5
    lambda_Y <- 2
    
    a_sigma <- 3
    b_sigma <- 2
    # A_sigma <- .5
    
    # ggplot() + stat_function(fun  = dinvgamma,
    # args = list(alpha = 3, beta = 2)) + xlim(c(0,3))
    
    sigma_beta <- 1
    
    sigma_mu <- 1
    sigma_gamma <- 1
    
    sigma_u <- 1
    
    sigma_beta_theta <- 1
    sigma_alpha <- 1
    
    prior_r <- 10
    r_var <- 25
    
    sigma_lambda <- 3
    
    nu0 <- 1
    Sigma0 <- diag(1, nrow = S)
    
    # MCMC --------------------------------------------------------------------
    
    nchain <- 1
    nburn <- 5
    niter <- 5000
    nthin <- 1
    
    iterToAdapt <- 100000
    
    simulatedData <- F
    mixedScenario <- T
    
    # params to update
    {
      updateAll <- F
      if(updateAll){
        updateLambda <- T; correctLambda <- !updateLambda
        updateL <- T; correctL <- !updateL
        updateMu <- T; correctMu <- !updateMu
        
        updateV <- F; correctV <- !updateV
        updateU <- F; correctU <- !updateU
        updateUV <- T; correctU <- !updateUV; correctV <- !updateUV
        
        updateLambdaijk <- T
        updateR <- T; correctR <- !updateR
        
        updateBeta_w <- T; correctBeta_w <- !updateBeta_w
        updateBeta_z <- T; correctBeta_z <- !updateBeta_z
        
        updateSigma <- T; correctSigma <- !updateSigma
        updateTau <- T; correctTau <- !updateTau
        
        updateDeltaGammaC <- T; correctDeltaGammaC <- !updateDeltaGammaC#F
        
        updateLambdaTilde <- T; correctLambdaTilde <- !updateLambdaTilde
        updateLambda0 <- T; correctLambda0 <- !updateLambda0
        
        updateP11 <- T; correctP11 <- !updateP11
        updateP10 <- T; correctP10 <- !updateP10
        
        updateBetaTheta <- T; correctBetaTheta <- !updateBetaTheta
        updateTheta10 <- T; correctTheta10 <- !updateTheta10
      } else {
        
        updateLambda <- T; correctLambda <- !updateLambda
        updateLambdaJoint <- F; correctLambda <- !updateLambdaJoint
        updateL <- T; correctL <- !updateL
        updateMu <- T; correctMu <- !updateMu
        
        updateV <- F; correctV <- !updateV
        updateU <- F; correctU <- !updateU
        updateUV <- T; correctU <- !updateUV; correctV <- !updateUV
        
        updateLambdaijk <- T
        updateR <- T; correctR <- !updateR
        
        updateBeta_w <- T; correctBeta_w <- !updateBeta_w
        updateBeta_z <- T; correctBeta_z <- !updateBeta_z
        
        updateSigma <- T; correctSigma <- !updateSigma
        updateTau <- T; correctTau <- !updateTau
        
        updateDeltaGammaC <- F; correctDeltaGammaC <- !updateDeltaGammaC#F
        
        updateLambdaTilde <- F; correctLambdaTilde <- !updateLambdaTilde
        updateLambda0 <- F; correctLambda0 <- !updateLambda0
        
        updateP11 <- F; correctP11 <- !updateP11
        updateP10 <- F; correctP10 <- !updateP10
        
        updateBetaTheta <- F; correctBetaTheta <- !updateBetaTheta
        updateTheta10 <- F; correctTheta10 <- !updateTheta10
      }
      
    }
    
    # output
    {
      # trueIters <- niter / nthin
      # lambda0_output <- rep(NA, niter)
      mu0_output <- matrix(NA, nchain, niter)
      n0_output <- matrix(NA, nchain, niter)
      lambda_output <- array(NA, c(nchain, S + S_star, niter))
      lambdatilde_output <- array(NA, c(nchain, S + S_star, niter))
      u_output <- array(NA, c(nchain, sum(M_site) + emptyTubes, max(K), niter))
      zeta_output <- array(NA, c(nchain, niter, sum(M_site)))
      logz_output <- array(NA, dim = c(nchain, n, S, niter))
      v_output <- array(NA, dim = c(nchain, sum(M_site) + emptyTubes, S + S_star, niter))
      beta_w_output <- array(NA, c(nchain, ncov_w, S, niter))
      sigma_output <- array(NA, c(nchain, S, niter))
      if(jSDM){
        Tau_output <- array(NA, dim = c(nchain, niter, S, S))  
      } else {
        tau_output <- array(NA, c(nchain, S, niter))
      }
      beta_theta_output <- array(NA, dim = c(nchain, S, 2 + ncov_w, niter))
      beta0_output <- array(NA, dim = c(nchain, S, niter))
      mu_output <- array(NA, dim = c(nchain, S, niter))
      beta_z_output <- array(NA, dim = c(nchain, ncov_z, S, niter))
      lambda_ijk_output <- array(NA, dim = c(nchain, sum(M_site) + emptyTubes, max(K), S + S_star, niter))
      theta11_output <- array(NA, dim = c(nchain, sum(M_site), S, niter))
      p_11_output <- array(NA, dim = c(nchain, S + S_star, niter))
      p_10_output <- array(NA, dim = c(nchain, S + S_star, niter))
      r_output <- array(NA, dim = c(nchain, S + S_star, niter))
      delta_output <- array(NA, dim = c(nchain, sum(M_site) + emptyTubes, S + S_star, niter))
      gamma_output <- array(NA, dim = c(nchain, sum(M_site) + emptyTubes, S + S_star, niter))
      cimk_output <- array(NA, dim = c(nchain, sum(M_site) + emptyTubes, max(K), S + S_star, niter))
    }
    
    for (chain in 1:nchain) {
      
      # chain output
      {
        mu0_output_iter <- rep(NA, niter)
        n0_output_iter <- rep(NA, niter)
        lambda_output_iter <- matrix(NA, nrow = S + S_star, ncol = niter)
        lambdatilde_output_iter <- matrix(NA, nrow = S + S_star, ncol = niter)
        u_output_iter <- array(NA, dim = c(sum(M_site) + emptyTubes, max(K), niter))
        zeta_output_iter <- matrix(NA, niter, sum(M_site))
        logz_output_iter <- array(NA, dim = c(n, S, niter))
        v_output_iter <- array(NA, dim = c(sum(M_site) + emptyTubes, S + S_star, niter))
        beta_w_output_iter <- array(NA, c(ncov_w, S, niter))
        sigma_output_iter <- matrix(NA, S, niter)
        if(jSDM){
          Tau_output_iter <- array(NA, dim = c(niter, S, S))
        } else {
          tau_output_iter <- matrix(NA, S, niter)
        }
        beta_theta_output_iter <- array(NA, dim = c(S, 2 + ncov_w, niter))
        beta0_output_iter <- array(NA, dim = c(S, niter))
        mu_output_iter <- array(NA, dim = c(S, niter))
        beta_z_output_iter <- array(NA, dim = c(ncov_z, S, niter))
        lambda_ijk_output_iter <- array(NA, dim = c(sum(M_site) + emptyTubes, max(K), S + S_star, niter))
        theta11_output_iter <- array(NA, dim = c(sum(M_site), S, niter))
        p_11_output_iter <- array(NA, dim = c(S + S_star, niter))
        p_10_output_iter <- array(NA, dim = c(S + S_star, niter))
        r_output_iter <- matrix(NA, S + S_star, niter)
        delta_output_iter <- array(NA, dim = c(sum(M_site) + emptyTubes, S + S_star, niter))
        gamma_output_iter <- array(NA, dim = c(sum(M_site) + emptyTubes, S + S_star, niter))
        cimk_output_iter <- array(NA, dim = c(sum(M_site) + emptyTubes, max(K), S  + S_star, niter))
      }
      
      # starting values
      {
        # set lambda prior
        {
          # lambda_prior <- apply(OTU, 2, function(x){
          #   log(mean(x[x > quantile(x[x > 0], probs = c(0.1))]))
          #   # mean(x[x > 50])
          #   # 0
          # })
          # print("change!!")
          lambda_prior <- apply(y, 3, function(x){
            log(mean(x[x > quantile(x[x > 0], probs = c(0.1))]))
            # mean(x[x > 50])
            # 0
          })
        }
        
        if(mixedScenario){
          
          if(correctLambda){
            
            lambda <- lambda_true
            
          } else 
          {
            
            # lambda <- apply(OTU, 2, function(x){
            #   rnorm(1, log(mean(x[x > quantile(x[x > 0], probs = c(0.1))])), sd = .05)
            #   # mean(x[x > 50])
            #   # 0
            # })
            
            lambda <- apply(y, 3, function(x){
              log(mean(x[x > quantile(x[x > 0], probs = c(0.1))]))
              # mean(x[x > 50])
              # 0
            })
            
          }
          
          if(correctL){
            
            logz <- logz_true
            
          } else 
          {
            logz <- matrix(0, n, S)
          }
          
          if(correctMu){
            mu <- mu_true
          } else 
          {
            mu <- rep(0, S)
          }
          
          if(correctV){
            v <- v_true
          } else 
          {
            v <- matrix(0, sum(M_site) + emptyTubes, S)
            if(S_star > 0){
              v <- cbind(v, rbind(v_spikes,
                                  matrix(0, nrow = emptyTubes,
                                         ncol = S_star)))
            }
          }
          
          if(correctU){
            u <- u_true
          } else 
          {
            u <- matrix(0, sum(M_site) + emptyTubes, max(K))
          }
          
          if(correctR){
            r_nb <- r_true
          } else 
          {
            r_nb <- rpois(S + S_star, 10)
          }
          
          if(correctBeta_w){
            beta_w <- beta_w_true
          } else 
          {
            beta_w <- matrix(0, ncov_w, S)
          }
          
          if(correctBeta_z){
            beta0 <- beta0_true
            beta_z <- beta_z_true
          } else
          {
            beta0 <- rep(0, S)
            beta_z <- matrix(0, nrow = ncov_z, ncol = S)
          }
          
          if(correctDeltaGammaC){
            c_imk <- c_imk_true
            
            delta <- delta_true
            
            gamma <- gamma_true
          } else 
          {
            c_imk <- array(NA, dim = c(sum(M_site) + emptyTubes, max(K), S + S_star))
            for (i in 1:n) {
              for (m in 1:M_site[i]) {
                numRep <- K[m + sum(M_site[seq_len(i-1)])]
                for (k in 1:numRep) {
                  for (j in 1:(S + S_star)) {
                    
                    # if(y[m + sum(M_site[seq_len(i-1)]),k,j] > 100){
                    if(y[m + sum(M_site[seq_len(i-1)]),k,j] > sample(c(100,500), 1)){
                      c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- 1
                      # } else if(y[m + sum(M_site[seq_len(i-1)]),k,j] > 10){
                    } else if(y[m + sum(M_site[seq_len(i-1)]),k,j] > sample(c(10,50), 1)){
                      c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- 2
                    } else {
                      c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- 0
                    }
                    
                  }
                }
              }
            }
            for (m in seq_len(emptyTubes)) {
              if(y[sum(M_site) + m,k,j] > sample(c(100,500), 1)){
                c_imk[sum(M_site) + m,k,j] <- 1
                # } else if(y[m + sum(M_site[seq_len(i-1)]),k,j] > 10){
              } else if(y[sum(M_site) + m,k,j] > sample(c(10,50), 1)){
                c_imk[sum(M_site) + m,k,j] <- 2
              } else {
                c_imk[sum(M_site) + m,k,j] <- 0
              }
            }
            
            delta <- array(NA, dim = c(sum(M_site) + emptyTubes, S + S_star))
            for (i in 1:n) {
              for (m in 1:M_site[i]) {
                for (j in 1:S) {
                  if(sum(c_imk[m + sum(M_site[seq_len(i-1)]),
                               1:K[m + sum(M_site[seq_len(i-1)])],j]) > 0){
                    delta[m + sum(M_site[seq_len(i-1)]),j] <- 1  
                  } else {
                    delta[m + sum(M_site[seq_len(i-1)]),j] <- 0
                  }
                  
                }
                for (j in seq_len(S_star)) {
                  delta[m + sum(M_site[seq_len(i-1)]),S + j] <- 1
                }
              }
            }
            for (m in seq_len(emptyTubes)) {
              delta[sum(M_site) + m,] <- 0
            }
            
            gamma <- array(NA, dim = c(sum(M_site) + emptyTubes, S + S_star))
            for (i in 1:n) {
              for (m in 1:M_site[i]) {
                for (j in 1:(S + S_star)) {
                  gamma[m + sum(M_site[seq_len(i-1)]),j] <- 0 
                }
              }
            }
            for (m in seq_len(emptyTubes)) {
              gamma[sum(M_site) + m,] <- 0 
            }
            
          }
          
          if(correctLambda0){
            mu0 <- mu0_true 
            n0 <- n0_true 
            pi0 <- pi0_true 
          } else 
          {
            mu0 <- 2
            n0 <- 1
            pi0 <- .99
          }
          
          if(correctLambdaTilde){
            lambdatilde <- lambdatilde_true
          } else 
          {
            lambdatilde <- rep(.5, S + S_star)
          }
          
          if(correctSigma){
            sigma <- sigma_true
          } else 
          {
            sigma <- rep(1, S)
          }
          
          if(correctTau){
            if(jSDM){
              invTau <- solve(Tau_true) # Omega_Y 
              Tau <- Tau_true # Sigma_y #solve(Omega)
              tau_NG <- matrix(1, nrow = S, ncol = S)
              lambda_NG <- 1
              gamma_NG <- 1
            } else {
              tau <- tau_true
            }
          } else 
          {
            if(jSDM){
              invTau <- diag(1, nrow = S) # Omega_Y 
              Tau <- diag(1, nrow = S) # Sigma_y #solve(Omega)
              tau_NG <- matrix(1, nrow = S, ncol = S)
              lambda_NG <- 1
              gamma_NG <- 1
            } else {
              tau <- rep(1, S)  
            }
          }
          
          if(correctP11){
            p_11 <- p11_true
          } else 
          {
            p_11 <- rep(0.9, S + S_star)
          }
          
          if(correctP10){
            p_10 <- p10_true
          } else 
          {
            p_10 <- rep(0.01, S + S_star)
          }
          
          if(correctBetaTheta){
            beta_theta <- beta_theta_true
            theta11 <- theta11_true
          } else 
          {
            # beta_theta <- cbind(rep(0, S), pmax(rnorm(S, 1), 0), rep(0, S))
            beta_theta <- cbind(rep(0, S), pmax(rnorm(S, 1), 0), matrix(0, nrow = S, ncol = ncov_w))
            theta11 <- matrix(NA, nrow = sum(M_site), ncol = S)
            for (i in 1:n) {
              for (m in 1:M_site[i]) {
                for (j in 1:S) {
                  theta11[m + sum(M_site[seq_len(i-1)]),j] <- 
                    logistic(c(1, logz[i,j], X_w[m + sum(M_site[seq_len(i-1)]),]) %*% 
                               beta_theta[j,])
                }
              }
            }
          }
          
          if(correctTheta10){
            theta10 <- theta10_true
          } else 
          {
            theta10 <- rep(.01, S)
          }
          
          # other
          {
            zeta <- rep(0, n)
            # zeta <- rep(0, sum(M_site))
            zeta_z <- rep(0, S)
            
            lambda_ijk <- array(NA, c(sum(M_site) + emptyTubes, max(K), S + S_star))
            for (i in 1:n) {
              # print(i)
              for (m in 1:M_site[i]) {
                numRep <- K[m + sum(M_site[seq_len(i-1)])]
                for (k in 1:numRep) {
                  u_imk <- u[m + sum(M_site[seq_len(i-1)]), k]
                  for (j in 1:(S + S_star)) {
                    lambda_ijk[m + sum(M_site[seq_len(i - 1)]), k, j] <-  exp(v[m + sum(M_site[seq_len(i - 1)]), j] + 
                                                                                lambda[j] + 
                                                                                u_imk)
                  }
                }
              }
            }
            for (m in seq_len(emptyTubes)) {
              numRep <- K[sum(M_site) + m]
              for (k in 1:numRep) {
                u_imk <- u[sum(M_site) + m, k]
                for (j in 1:(S + S_star)) {
                  lambda_ijk[sum(M_site) + m, k, j] <-  exp(v[sum(M_site) + m, j] + 
                                                              lambda[j] + 
                                                              u_imk)
                }
              }
            }
            
            A <- array(NA, dim = c(sum(M_site), S, (3 * 2^max(K))))
            
            for (i in 1:n) {
              for (m in 1:M_site[i]) {
                currentK <- K[m + sum(M_site[seq_len(i-1)])]
                A[m + sum(M_site[seq_len(i-1)]),,] <- 1 / (3 * 2^currentK)
              }
            }
            
            pi_hat <- array(0, dim = c(sum(M_site), 3 * 2^max(K), S))
          }
          
        } else {
          if(simulatedData){
            
            # intensities
            {
              
              lambda <- lambda_true
              
              v <- v_true
              
              u <- u_true
              
              lambda_ijk <- array(NA, c(sum(M_site), max(K), S))
              for (i in 1:n) {
                # print(i)
                for (m in 1:M_site[i]) {
                  numRep <- K[m + sum(M_site[seq_len(i-1)])]
                  for (k in 1:numRep) {
                    u_imk <- u_true[m + sum(M_site[seq_len(i-1)]), k]
                    for (j in 1:S) {
                      lambda_ijk[m + sum(M_site[seq_len(i - 1)]), k, j] <-  exp(v[m + sum(M_site[seq_len(i - 1)]), j] + 
                                                                                  lambda[j] + 
                                                                                  u_imk)
                    }
                  }
                }
              }
              
              logz <- logz_true
              mu <- mu_true
              
              mu0 <- mu0_true 
              n0 <- n0_true 
              pi0 <- pi0_true 
              
              r_nb <- r_true
              
              lambdatilde <- lambdatilde_true
              
              beta0 <- beta0_true
              beta_z <- beta_z_true
              alpha <- alpha_true
              
              zeta <- rep(0, sum(M_site))
              zeta_z <- rep(0, S)
            }
            
            # indicator variables
            {
              c_imk <- c_imk_true
              
              delta <- delta_true
              
              gamma <- gamma_true
              
            }
            
            # probabilities
            {
              beta_theta <- beta_theta_true
              theta11 <- theta11_true
              theta10 <- theta10_true
              
              p_11 <- p11_true
              
              p_10 <- p10_true
            }
            
            # variances
            {
              sigma <- sigma_true
              # tau <- tau_true
              Tau <- Tau_true
              
              a_sigma <- rep(1, S)
              a_tau <- rep(1, S)
              
              # graphical horseshoe parameters
              {
                invTau <- diag(1, nrow = S) # Omega_Y 
                Tau <- diag(1, nrow = S) # Sigma_y #solve(Omega)
                # Lambdasq <- matrix(1, nrow = n_s, ncol = n_s)
                # Nu <- matrix(1, nrow = n_s, ncol = n_s)
                # tausq <- 1
                # csi <- 1  
                tau_NG <- matrix(1, nrow = S, ncol = S)
                lambda_NG <- 1
                gamma_NG <- 1
                lambda_Y <- 2
              }
              
              # sigma_u <- sigma_u_true
            }
            
          } 
          else {
            
            # intensities
            {
              
              # lambda <- apply(OTU, 2, function(x){
              #   rnorm(1, log(mean(x[x > quantile(x[x > 0], probs = c(0.1))])), sd = .05)
              #   # mean(x[x > 50])
              #   # 0
              # })
              
              print("change!!")
              lambda_prior <- apply(y, 2, function(x){
                log(mean(x[x > quantile(x[x > 0], probs = c(0.1))]))
                # mean(x[x > 50])
                # 0
              })
              
              logz <- matrix(0, n, S)
              mu <- rep(0, S)
              v <- matrix(0, sum(M_site), S)
              
              u <- matrix(0, sum(M_site), max(K))
              
              mu0 <- 1
              n0 <- 1
              pi0 <- .99
              
              r_nb <- rpois(S, 10)
              
              lambdatilde <- rep(.5, S)
              
              beta0 <- rep(0, S)
              alpha <- rep(0, S)
              beta_z <- matrix(0, nrow = ncov_z, ncol = S)
              
              zeta <- rep(0, sum(M_site))
              zeta_z <- rep(0, S)
              
              lambda_ijk <- array(NA, c(sum(M_site), max(K), S))
              for (i in 1:n) {
                # print(i)
                for (m in 1:M_site[i]) {
                  numRep <- K[m + sum(M_site[seq_len(i-1)])]
                  for (k in 1:numRep) {
                    u_imk <- u[m + sum(M_site[seq_len(i-1)]), k]
                    for (j in 1:S) {
                      lambda_ijk[m + sum(M_site[seq_len(i - 1)]), k, j] <-  exp(v[m + sum(M_site[seq_len(i - 1)]), j] + 
                                                                                  lambda[j] + 
                                                                                  u_imk)
                    }
                  }
                }
              }
              
            }
            
            # indicator variables
            {
              c_imk <- array(NA, dim = c(sum(M_site), max(K), S))
              for (i in 1:n) {
                for (m in 1:M_site[i]) {
                  numRep <- K[m + sum(M_site[seq_len(i-1)])]
                  for (k in 1:numRep) {
                    for (j in 1:S) {
                      if(emptySites[i] == 1){
                        c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- 0
                      } else {
                        # if(y[m + sum(M_site[seq_len(i-1)]),k,j] > 100){
                        if(y[m + sum(M_site[seq_len(i-1)]),k,j] > sample(c(100,500), 1)){
                          c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- 1
                          # } else if(y[m + sum(M_site[seq_len(i-1)]),k,j] > 10){
                        } else if(y[m + sum(M_site[seq_len(i-1)]),k,j] > sample(c(10,50), 1)){
                          c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- 2
                        } else {
                          c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- 0
                        }
                      }
                    }
                  }
                }
              }
              
              delta <- array(NA, dim = c(sum(M_site), S))
              for (i in 1:n) {
                for (m in 1:M_site[i]) {
                  for (j in 1:S) {
                    if(emptySites[i] == 0){
                      if(sum(c_imk[m + sum(M_site[seq_len(i-1)]),
                                   1:K[m + sum(M_site[seq_len(i-1)])],j]) > 0){
                        delta[m + sum(M_site[seq_len(i-1)]),j] <- 1  
                      } else {
                        delta[m + sum(M_site[seq_len(i-1)]),j] <- 0
                      }
                    } else {
                      delta[m + sum(M_site[seq_len(i-1)]),j] <- 0
                    }
                  }
                }
              }
              
              gamma <- array(NA, dim = c(sum(M_site), S))
              for (i in 1:n) {
                for (m in 1:M_site[i]) {
                  for (j in 1:S) {
                    gamma[m + sum(M_site[seq_len(i-1)]),j] <- 0 
                  }
                }
              }
              
              # all_probs <- array(0, dim = c(sum(M_site), S, 2^(max(K)) + 3^(max(K))))
              # index_allprobs <- 0
            }
            
            # probabilities
            {
              # beta_theta <- cbind(rep(0, S), rep(1, S), rep(0, S))
              beta_theta <- cbind(rep(0, S), pmax(rnorm(S, 1), 0), rep(0, S))
              
              theta11 <- matrix(NA, nrow = sum(M_site), ncol = S)
              for (i in 1:n) {
                for (m in 1:M_site[i]) {
                  for (j in 1:S) {
                    theta11[m + sum(M_site[seq_len(i-1)]),j] <- 
                      logistic(c(1, logz[i,j], r[m + sum(M_site[seq_len(i-1)])]) %*% 
                                 beta_theta[j,])
                  }
                }
              }
              
              theta10 <- rep(.01, S)
              
              p_11 <- rep(0.9, S)
              
              p_10 <- rep(0.01, S)
            }
            
            # variances
            {
              sigma <- rep(1, S)
              tau <- rep(1, S)
              
              a_sigma <- rep(1, S)
              a_tau <- rep(1, S)
              
              # graphical horseshoe parameters
              {
                invTau <- diag(1, nrow = S) # Omega_Y 
                Tau <- diag(1, nrow = S) # Sigma_y #solve(Omega)
                # Lambdasq <- matrix(1, nrow = n_s, ncol = n_s)
                # Nu <- matrix(1, nrow = n_s, ncol = n_s)
                # tausq <- 1
                # csi <- 1  
                tau_NG <- matrix(1, nrow = S, ncol = S)
                lambda_NG <- 1
                gamma_NG <- 1
              }
              
            }
            
          }
        }
        
        
      }
      
      for (iter in 1:(nburn + nthin * niter)) {
        
        print(paste0("S_star = ", S_star))
        if(iter <= nburn){
          print(paste0("Chain = ",chain," - Burn-in Iteration = ",iter))  
        } else {
          print(paste0("Chain = ",chain," - Iteration = ",iter - nburn))
        }  
        
        # LAMBDA ----------------------------------------------
        
        if(updateLambda){
          # print("Update lambda")
          
          # lambda <- update_lambda(beta0, mu, lambda,
          #                              sigma_beta, sigma_mu,
          #                              exp(lambda_prior), sigma_lambda,
          #                              S_star, betaThetaEqual1)
          
          list_lambda <- update_lambda_CP(beta0, mu, lambda,
                                          sigma_beta, sigma_mu,
                                          lambda_prior, sigma_lambda,
                                          S_star)
          lambda <- list_lambda$lambda
          beta0 <- list_lambda$beta0
          mu <- list_lambda$mu
          v <- list_lambda$v
          logz <- list_lambda$logz
          
        }
        
        # LAMBDA JOINT -------
        
        if(updateLambdaJoint){
          print("Update lambda")
          
          # lambda <- update_lambda(beta0, mu, lambda,
          #                              sigma_beta, sigma_mu,
          #                              exp(lambda_prior), sigma_lambda,
          #                              S_star, betaThetaEqual1)
          
          list_lambda <- update_lambdajoint_cpp(lambda, v, u, lambda_ijk, c_imk, delta, gamma, sigma, 
                                                logz, mu, r_nb, sigma_gamma, M_site, X_w, beta_w, lambda_prior, 
                                                sigma_lambda, K, emptyTubes, S_star)
          lambda <- list_lambda$lambda
          v <- list_lambda$v
          
          lambda <- update_lambda_spikeins(mu, 
                                           lambda,
                                           u,
                                           r_nb,
                                           lambda_prior, 
                                           sigma_lambda,
                                           S_star)
          
        }
        
        # MU -----------------------------------
        
        if(updateMu){
          # print("Update mu")
          
          list_mu <- update_mu_cpp(mu, lambda, delta, gamma, sigma, sigma_gamma, beta0,
                                   beta_z, logz, v, beta_theta, M_site, sigma_mu, S_star)
          mu <- list_mu$mu
          # v <- list_mu$v  
        }
        
        # LOG_Z ------------------------------
        
        if(updateL){
          
          if(jSDM){
            logz <- update_logz_JSDM_cpp(logz, beta0, X_z, beta_z, mu, 
                                         v, lambda, beta_theta, X_w, beta_w, 
                                         Tau, delta, gamma, sigma, 
                                         M_site, S_star)
          } else {
            logz <- update_logz_cpp(logz, beta0, X_z, beta_z,
                                    mu, v, lambda, beta_theta,
                                    X_w, beta_w, tau, delta, gamma, sigma,
                                    M_site,
                                    S_star)
          }
          
        }
        
        # LAMBDA IJK ----
        
        if(updateLambdaijk){
          # print("Update lambda ijk")
          lambda_ijk <- update_lambdaijk(lambda, lambda_ijk, v, u, r_nb, c_imk, M_site, y, K,
                                         S_star, emptyTubes)
        }
        
        # VU ---------
        
        if(updateUV){
          list_uv <- update_uv_poisgamma_cpp(u, v, logz, lambda, X_z, beta_theta, beta_z, beta0, 
                                             r_nb, mu, lambda_ijk, c_imk, delta, gamma, sigma, sigma_gamma, 
                                             sigma_u, M_site, X_w, beta_w, K, S_star)
          u <- list_uv$u
          v <- list_uv$v  
          lambda <- list_uv$lambda
        }
        
        # V ------------------------------------------
        
        if(updateV){
          print("Update v")
          
          v[1:sum(M_site),] <- update_v_poisgamma_cpp(v, logz, 
                                                      lambda,  X_z, 
                                                      beta_theta, u, beta_z,
                                                      beta0, r_nb, mu, lambda_ijk,
                                                      c_imk, delta, gamma, sigma, 
                                                      sigma_gamma, M_site, 
                                                      X_w, beta_w,
                                                      K, S_star)
        }
        
        # U --------------------------------------
        
        if(updateU){
          print("Update u")
          
          list_u <- update_u_poisgamma_cpp(v, zeta, u, lambda, beta0, beta_z, logz,
                                           mu, lambda_ijk, r_nb, X_w, beta_w, c_imk,
                                           delta, gamma, sigma_u, beta_theta, sigma, 
                                           sigma_gamma, M_site, 
                                           K, S_star, emptyTubes)
          u <- list_u$u
          # zeta <- list_u$zeta
          # v <- list_u$v
          # logz <- list_u$logz2
          # print(zeta)
          
        }
        
        # BETA_Z --------------
        
        if(updateBeta_z){
          # print("Update beta")
          if(jSDM){
            list_beta_z <- update_betaz_CP_corr(logz, Tau, X_z, sigma_beta, S_star)  
          } else {
            list_beta_z <- update_betaz_CP(logz, lambda, tau, X_z, sigma_beta, S_star)
          }
          
          beta0 <- list_beta_z$beta0
          beta_z <- list_beta_z$beta_z
          
        }
        
        # BETA_W ---------
        
        if(updateBeta_w){
          # print("Update betaw")
          if(ncov_w > 0){
            beta_w <- update_betaw_cpp(beta_w, v, delta, logz, X_w, sigma, sigma_alpha, M_site)  
          }
        }
        
        # TAU ----------------------------------------------------------
        
        if(updateTau){
          
          if(jSDM){
            # list_Tau <- update_Tau_Wish(nu0, Sigma0, logz, beta0, betaz, X_z)
            
            # invTau <- list_Tau$invSigma
            # Tau <- list_Tau$Sigma
            list_Tau <- update_Tau(invTau, Tau, tau_NG, lambda_NG, gamma_NG, lambda_Y,
                                   logz, beta0, beta_z, X_z)
            
            invTau <- list_Tau$Omega
            Tau <- list_Tau$Sigma
            tau_NG <- list_Tau$tau
            lambda_NG <- list_Tau$lambda
            gamma_NG <- list_Tau$gamma
          } else {
            tau <- update_tau_cpp(tau, logz, X_z, beta_z, beta0,
                                  a_tau, rep(b_tau, S)
                                  # a_tau = .5, b_tau = 1 / a_tau
            )
            # list_Tau <- update_Tau_Wish(nu0, Sigma0, logz, beta0, betaz, X_z)
            # 
            # invTau <- list_Tau$invSigma
            # Tau <- list_Tau$Sigma
            # tau <- sqrt(diag(Tau))
          }
          
        }
        
        # UPDATE a TAU ----------------------------------------------
        
        # a_tau <- sapply(1:S, function(j) {
        #   rinvgamma_cpp(.5, 1 / A_tau + 1 / tau[j]^2)
        # })
        
        # SIGMA ----------------------------------------------------------
        
        if(updateSigma){
          print("Update sigma")
          # sigma <- rep(1, S)
          sigma <- update_sigma_cpp(sigma, lambda, beta_z, beta0,
                                    mu, logz, v, X_w, beta_w,
                                    delta, gamma, beta_theta,
                                    a_sigma, rep(b_sigma, S),
                                    # a_sigma = .5, b_sigma = 1 / a_sigma,
                                    M_site, S_star)
          
        }
        
        # UPDATE a SIGMA ----------------------------------------------
        
        # a_sigma <- sapply(1:S, function(j) {
        #   rinvgamma_cpp(.5, 1 / A_sigma + 1 / sigma[j]^2)
        # })
        
        # R - NB ---------------------
        
        if(updateR){
          # print("Update r")
          
          r_nb <- update_r_nb_cpp(r_nb, prior_r, sqrt(r_var), lambda, u,  
                                  v, y, delta, gamma, c_imk, M_site, K, 
                                  optimStep = 
                                    # T,
                                    F,
                                  # ((iter - 1) %% 5 == 0),
                                  sd_r_proposal = .05)
          
        }
        
        # DELTA/C/D --------------------------------------------------------
        
        if(updateDeltaGammaC){
          # print("Update delta-gamma-c")
          
          if(iter < iterToAdapt){
            
            list_deltagammac <- update_delta_c_d_rjmcmc(y, v, lambda, r_nb,
                                                        M_site, K, lambdatilde,
                                                        mu0, n0, pi0, u, logz, X_w, beta_w,
                                                        sigma, mu, sigma_gamma, v_sd = .5,
                                                        p_11, p_10, theta11, theta10, emptyTubes,
                                                        S_star)
            # list_deltagammac <- update_delta_foreach(y, v, lambda, r_nb,
            #                                          M_site, K, lambdatilde,
            #                                          mu0, n0, pi0, u, logz, X_w, beta_w,
            #                                          sigma, mu, sigma_gamma, v_sd = .5,
            #                                          p_11, p_10, theta11, theta10, emptySites)
            delta <- list_deltagammac$delta
            gamma <- list_deltagammac$gamma
            c_imk <- list_deltagammac$c_imk
            v <- list_deltagammac$v
            
          } else {
            
            list_deltagammac <- update_delta_c_d_proposals(c_imk, delta, gamma, A, y, v, lambda, r_nb,
                                                           M_site, K, lambdatilde,
                                                           mu0, n0, pi0, u, logz, r, alpha,
                                                           sigma, mu, sigma_gamma, v_sd = .5,
                                                           p_11, p_10, theta11, theta10, emptySites)
            delta <- list_deltagammac$delta
            gamma <- list_deltagammac$gamma
            c_imk <- list_deltagammac$c_imk
            v <- list_deltagammac$v
            
            
          }
          
        }
        
        if(iter %% 100 == -1){
          
          # pi_hat <- apply(deltagammac_output_iter[,,,1:(iter - 1)], c(1,2,3), mean)
          # 
          # pi_hat <- deltagammac_output_iter[,,,1:(iter - 1)]
          
          for (i in 1:n) {
            for (m in 1:M_site[i]) {
              currentK <- K[m + sum(M_site[seq_len(i-1)])]
              for (j in 1:S) {
                currentProbs <- pi_hat[m + sum(M_site[seq_len(i-1)]),1:(3 * 2^currentK),j]
                newProbs <- .95 * currentProbs + .05 * rep(1 / (3 * 2^currentK), 3 * 2^currentK)
                A[m + sum(M_site[seq_len(i-1)]), j, 1:(3 * 2^currentK)] <- newProbs
              }
            }
          }
          
        }
        
        
        # for (i in 1:n) {
        #   for (m in 1:M_site[i]) {
        #     currentK <- K[m + sum(M_site[seq_len(i-1)])]
        #     for (j in 1:S) {
        #       index_deltagammac <- convertDeltaIndexes(delta[m + sum(M_site[seq_len(i-1)]),j], 
        #                                                gamma[m + sum(M_site[seq_len(i-1)]),j], 
        #                                                c_imk[m + sum(M_site[seq_len(i-1)]),1:currentK,j], 
        #                                                currentK)
        #       deltagammac_output_iter[m + sum(M_site[seq_len(i-1)]),
        #                               index_deltagammac + 1,
        #                               j,
        #                               iter] <- 1
        #     }
        #   }
        # }
        
        
        # deltagammac_current <- updateDeltaGammaC_iter(delta, gamma, c_imk,
        #                                                            K, n, M_site, S)
        # 
        # pi_hat <- (pi_hat * (iter - 1) + deltagammac_current ) / iter
        
        # LAMBDA 0 ---------------------------------------------------------
        
        if(updateLambda0){
          print("Update lambda0")
          
          list_lambda0 <- update_lambda0_mixt(y, c_imk, mu0, n0, sd_mu0 = .05, sd_n0 = .05)
          mu0 <- list_lambda0$mu0
          n0 <- list_lambda0$n0
          pi0 <- list_lambda0$pi0
        }
        
        # LAMBDA TILDE ---------------------------------------------------------
        
        if(updateLambdaTilde){
          print("Update lambda tilde")
          
          lambdatilde <- update_lambdatilde(y, c_imk, lambda)
        }
        
        # P11 --------------------------------------------------------------
        
        if(updateP11){
          print("Update p11")
          
          p_11 <- update_p_11_cpp(p_11, delta, gamma, c_imk, M_site,
                                  K, a_p11, b_p11)
        }
        
        # P10 --------------------------------------------------------------
        
        if(updateP10){
          print("Update p10")
          
          p_10 <- update_p_10_cpp(p_10, delta, gamma, c_imk, M_site, K, a_p10, b_p10, emptyTubes)
        }
        
        # UPDATE BETA THETA11 --------------------------------------------------------------
        
        if(updateBetaTheta){
          print("Update beta theta")
          
          list_beta_theta <- update_betatheta11_cpp(logz, beta_theta, theta11, delta,
                                                    X_w, M_site,
                                                    b_theta11 = c(0,1,rep(0, ncov_w)),
                                                    B_theta11 = diag(sigma_beta_theta, nrow = 2 + ncov_w),
                                                    !betaThetaEqual1)
          beta_theta <- list_beta_theta$beta_theta
          theta11 <- list_beta_theta$theta11
          
        }
        
        # UPDATE THETA10 --------------------------------------------------------------
        
        if(updateTheta10){
          print("Update theta10")
          
          theta10 <- update_theta10_cpp(theta10, delta, gamma, M_site, K, emptySites,
                                        a_theta0, b_theta0)
        }
        
        # WRITE RESULTS -----------------------------------------------------------
        
        if(iter > nburn & (iter - nburn) %% nthin == 0){
          trueIter <- (iter - nburn) / nthin
          
          # lambda0_output[trueIter] <- lambda0
          mu0_output_iter[trueIter] <- mu0
          n0_output_iter[trueIter] <- n0
          lambdatilde_output_iter[,trueIter] <- lambdatilde
          lambda_output_iter[,trueIter] <- lambda
          beta_w_output_iter[,,trueIter] <- beta_w
          zeta_output_iter[trueIter,] <- zeta
          u_output_iter[,,trueIter] <- u
          v_output_iter[,,trueIter] <- v
          logz_output_iter[,,trueIter] <- logz
          if(jSDM){
            Tau_output_iter[trueIter,,] <- Tau
          } else {
            tau_output_iter[,trueIter] <- tau
          }
          sigma_output_iter[,trueIter] <- sigma
          beta_theta_output_iter[,,trueIter] <- beta_theta
          lambda_ijk_output_iter[,,,trueIter] <- lambda_ijk
          beta_z_output_iter[,,trueIter] <- beta_z
          beta0_output_iter[,trueIter] <- beta0
          mu_output_iter[,trueIter] <- mu
          theta11_output_iter[,,trueIter] <- theta11
          p_11_output_iter[,trueIter] <- p_11
          p_10_output_iter[,trueIter] <- p_10
          r_output_iter[,trueIter] <- r_nb
          
          delta_output_iter[,,trueIter] <- delta
          gamma_output_iter[,,trueIter] <- gamma
          cimk_output_iter[,,,trueIter] <- c_imk
          
          # deltagammac_output_iter[,,1,trueIter] <- delta
          # deltagammac_output_iter[,,2,trueIter] <- gamma
          
        }
        
      }
      
      # copy chains output
      {
        mu0_output[chain,] <- mu0_output_iter
        n0_output[chain,] <- n0_output_iter
        lambda_output[chain,,] <- lambda_output_iter
        lambdatilde_output[chain,,] <- lambdatilde_output_iter
        u_output[chain,,,] <- u_output_iter
        zeta_output[chain,,] <- zeta_output_iter
        logz_output[chain,,,] <- logz_output_iter
        v_output[chain,,,] <- v_output_iter
        beta_w_output[chain,,,] <- beta_w_output_iter
        sigma_output[chain,,] <- sigma_output_iter
        if(jSDM) {
          Tau_output[chain,,,] <- Tau_output_iter  
        } else {
          tau_output[chain,,] <- tau_output_iter
        }
        beta_theta_output[chain,,,] <- beta_theta_output_iter
        beta0_output[chain,,] <- beta0_output_iter
        mu_output[chain,,] <- mu_output_iter
        beta_z_output[chain,,,] <- beta_z_output_iter
        lambda_ijk_output[chain,,,,] <- lambda_ijk_output_iter
        theta11_output[chain,,,] <- theta11_output_iter
        p_11_output[chain,,] <- p_11_output_iter
        p_10_output[chain,,] <- p_10_output_iter
        r_output[chain,,] <- r_output_iter
        delta_output[chain,,,] <- delta_output_iter
        gamma_output[chain,,,] <- gamma_output_iter
        cimk_output[chain,,,,] <- cimk_output_iter
      }
      
    }
  }
  
  # save output
  {
    # summarise results
    {
      u_mean <- apply(u_output, c(2,3), mean)
      (bias <- mean(abs(u_mean - u_true)))
      u_biases[s_idx] <- bias
      
      u_vars[s_idx] <- mean(apply(u_output, c(1,2), var))
      
      #
      
      beta_mean <- apply(beta_z_output, c(2,3), mean)
      (bias <- mean(abs(beta_mean - beta_z_true)))
      beta_biases[s_idx] <- bias
      
      beta_vars[s_idx] <- mean(apply(beta_output, c(1,2), var))
      
      #
      
      # computeL()
      
      l_vars <- array(NA, dim = c(S, n - 1))
      l_bias <- array(NA, dim = c(S, n - 1))
      
      for (j in 1:S) {
        print(j)
        l <- 1
        for (i in seq_len(n-1)) {
          l_vars[j,i] <- var(logz_output[,i,j,] - logz_output[,i + 1,j,])
          l_bias[j,i] <- mean(abs((logz_output[,i,j,] - logz_output[,i + 1,j,]) -
                                    (logz_true[i,j] - logz_true[i + 1,j])))
        }
        # for (i1 in 1:n) {
        #   for (i2 in seq_len(i1-1)) {
        #     l_vars[j,l] <- var(l_output[i1,j,] - l_output[i2,j,])
        #     l_bias[j,l] <- mean((l_output[i1,j,] - l_output[i2,j,]) -
        #                           (l_true[i1,j] - l_true[i2,j]))
        #     l <- l + 1
        #   }
        # }
      }
      
      ldiff_biases[s_idx] <- mean(l_bias)
      
      ldiff_vars[s_idx] <- mean(l_vars)
      
    }
  }
  
}

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
save(beta_biases, beta_vars, ldiff_biases, ldiff_vars, u_biases, u_vars, file = "results_tau1_sigma1_n1000_m1_k1.rda")
