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

load(here("Data","Lakes.Rdata"))
# load(here("Data","Leeches2.rda"))
# load(here("Data","Leeches.Rdata"))
# load(here("Data","Ponds.Rdata"))

subsetSpecies <- 1:(dim(y)[3])

S <- length(subsetSpecies)
OTU <- OTU[,subsetSpecies]
y <- y[,,subsetSpecies]

S_star <- 0

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

# SIMULATE DATA -----------------------------------------------------------

# design
{
  S <- 2
  n <- 1000
  ncov_z <- 1
  ncov_w <- 0
  
  M_site <- rep(1, n)
  emptyTubes <- 0
  
  K <- rep(1, sum(M_site) + emptyTubes)
  
  PCR_spiked <- rep(F, sum(K))
  S_star <- 0 # numOfSpikes
  
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
  beta_theta_0_mean <- 0
  sigmas <- 1
  taus <- 1
  r_0 <- 100
  p11_0 <- .95
  p10_0 <- .02
  theta_10 <- .02
  
  sd_beta_theta_0 <- .000001
  sigma_beta0 <- 1
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

beta_theta_true <- cbind(
  rnorm(S, mean = beta_theta_0_mean, sd = sd_beta_theta_0),  # baseline
  rep(1, S), # DNA coefficient
  t(matrix(sample(c(1,-1), (ncov_w) * S, replace = T), 
           nrow = ncov_w, ncol = S, byrow = T))) # covariate coefficient

theta10_true <- rep(theta_10, S)

lambda_true <- rnorm(S + S_star, log(1000), sd = .5)
r_true <- rep(r_0, S + S_star)#pmax(rgamma(S, 1, .2), 10)
# lambdatilde_true <- rep(.3, S + S_star)

mu_true <- rnorm(S, sd = 1)

mu0_true <- 5
n0_true <- 10
pi0_true <- .9

mu_tilde_true <- 5
n_tilde_true <- 100000

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
  X_wbeta_theta <- cbind(1, rep(logz_true[,j], M_site), X_w) %*% beta_theta_true[j,]
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
# log amount of DNA
v_true <- matrix(NA, sum(M_site), S + S_star)
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

all((which(c_imk_true[,,1] == 1) %in% 
       which(delta_true[,1] == 1 | gamma_true[,1] == 1)))

all((which(c_imk_true[,,2] == 1) %in% 
       which(delta_true[,2] == 1 | gamma_true[,2] == 1)))

# CLEAN DATA ---------

ncov_z <- ncol(X_z)
ncov_w <- ncol(X_w)

# PRIOR -------------------------------------------------------------------

a_theta0 <- 1
b_theta0 <- 20

a_p11 <- 20
b_p11 <- 1

a_p10 <- 1
b_p10 <- 100

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
sigma_gamma <- .5

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
nburn <- 5000
niter <- 5000
nthin <- 1

iterToAdapt <- 100000

simulatedData <- F
mixedScenario <- T

# params to update
{
  updateAll <- T
  if(updateAll){
    updateLambda_CP <- T; correctLambda <- !updateLambda_CP
    updateLambda_NP <- F; correctLambda <- !updateLambda_CP
    
    updateBeta_w <- T; correctBeta_w <- !updateBeta_w
    updateBeta_z <- T; correctBeta_z <- !updateBeta_z
    updateMu <- T; correctMu <- !updateMu
    
    updateL <- T; correctL <- !updateL
    
    updateV <- T; correctV <- !updateV
    updateU <- T; correctU <- !updateU
    updateUV <- T; correctU <- !updateUV; correctV <- !updateUV
    updateLambdaijk <- T; correctLambdaijk <- !updateLambdaijk
    
    updateSigma <- T; correctSigma <- !updateSigma
    updateTau <- T; correctTau <- !updateTau
    
    updateR <- T; correctR <- !updateR
    
    updateDeltaGammaC <- T; correctDeltaGammaC <- !updateDeltaGammaC#F
    
    updateLambdaTilde <- T; correctLambdaTilde <- !updateLambdaTilde
    updateLambda0 <- T; correctLambda0 <- !updateLambda0
    
    updateP11 <- T; correctP11 <- !updateP11
    updateP10 <- T; correctP10 <- !updateP10
    
    updateBetaTheta <- T; correctBetaTheta <- !updateBetaTheta
    updateTheta10 <- T; correctTheta10 <- !updateTheta10
  } else {
    
    updateLambda_CP <- F; correctLambda <- !updateLambda_CP
    updateLambda_NP <- F; #correctLambda <- !updateLambda
    # updateLambdaJoint <- F; correctLambda <- !updateLambdaJoint
    updateL <- F; correctL <- !updateL
    updateMu <- F; correctMu <- !updateMu
    
    updateV <- T; correctV <- !updateV
    updateU <- F; correctU <- !updateU
    updateUV <- F; correctU <- !updateUV; correctV <- !updateUV
    
    updateLambdaijk <- T
    updateR <- F; correctR <- !updateR
    
    updateBeta_w <- F; correctBeta_w <- !updateBeta_w
    updateBeta_z <- T; correctBeta_z <- !updateBeta_z
    
    updateSigma <- F; correctSigma <- !updateSigma
    updateTau <- F; correctTau <- !updateTau
    
    updateDeltaGammaC <- T; correctDeltaGammaC <- !updateDeltaGammaC#F
    
    updateLambdaTilde <- T; correctLambdaTilde <- !updateLambdaTilde
    updateLambda0 <- T; correctLambda0 <- !updateLambda0
    
    updateP11 <- F; correctP11 <- !updateP11
    updateP10 <- F; correctP10 <- !updateP10
    
    updateBetaTheta <- T; correctBetaTheta <- !updateBetaTheta
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
  mutilde_output <- array(NA, c(nchain, S + S_star, niter))
  u_output <- array(NA, c(nchain, sum(M_site), max(K), niter))
  # zeta_output <- array(NA, c(nchain, niter, sum(M_site)))
  logz_output <- array(NA, dim = c(nchain, n, S, niter))
  # v_output <- array(NA, dim = c(nchain, sum(M_site), S + S_star, niter))
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
  # lambda_ijk_output <- array(NA, dim = c(nchain, sum(M_site) + emptyTubes, max(K), S + S_star, niter))
  theta11_output <- array(NA, dim = c(nchain, sum(M_site), S, niter))
  p_11_output <- array(NA, dim = c(nchain, S + S_star, niter))
  p_10_output <- array(NA, dim = c(nchain, S + S_star, niter))
  r_output <- array(NA, dim = c(nchain, S + S_star, niter))
  # delta_output <- array(NA, dim = c(nchain, sum(M_site) + emptyTubes, S + S_star, niter))
  # gamma_output <- array(NA, dim = c(nchain, sum(M_site) + emptyTubes, S + S_star, niter))
  # cimk_output <- array(NA, dim = c(nchain, sum(M_site) + emptyTubes, max(K), S + S_star, niter))
}

for (chain in 1:nchain) {
  
  # chain output
  {
    mu0_output_iter <- rep(NA, niter)
    n0_output_iter <- rep(NA, niter)
    lambda_output_iter <- matrix(NA, nrow = S + S_star, ncol = niter)
    mutilde_output_iter <- matrix(NA, nrow = S + S_star, ncol = niter)
    u_output_iter <- array(NA, dim = c(sum(M_site), max(K), niter))
    # zeta_output_iter <- matrix(NA, niter, sum(M_site))
    logz_output_iter <- array(NA, dim = c(n, S, niter))
    # v_output_iter <- array(NA, dim = c(sum(M_site), S + S_star, niter))
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
    # lambda_ijk_output_iter <- array(NA, dim = c(sum(M_site), max(K), S + S_star, niter))
    theta11_output_iter <- array(NA, dim = c(sum(M_site), S, niter))
    p_11_output_iter <- array(NA, dim = c(S + S_star, niter))
    p_10_output_iter <- array(NA, dim = c(S + S_star, niter))
    r_output_iter <- matrix(NA, S + S_star, niter)
    # delta_output_iter <- array(NA, dim = c(sum(M_site) + emptyTubes, S + S_star, niter))
    # gamma_output_iter <- array(NA, dim = c(sum(M_site) + emptyTubes, S + S_star, niter))
    # cimk_output_iter <- array(NA, dim = c(sum(M_site) + emptyTubes, max(K), S  + S_star, niter))
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
      # lambda_prior <- lambda_true
      lambda_prior <- apply(y, 3, function(x){
        log(mean(x[x > quantile(x[x > 0], probs = c(0.1), na.rm = T)], na.rm = T))
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
        
        lambda <- lambda_prior
        
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
        v <- matrix(0, sum(M_site), S)
        if(S_star > 0){
          v <- cbind(v, v_spikes)
        }
      }
      
      if(correctU){
        u <- u_true
      } else 
      {
        u <- matrix(0, sum(M_site), max(K))
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
          if(y[sum(M_site) + m,k,j] > sample(c(10,50), 1)){
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
        mu_tilde <- mu_tilde_true
        n_tilde <- n_tilde_true
      } else
      {
        mu_tilde <- 100
        n_tilde <- 100000
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
        # beta_theta <- cbind(rep(0, S), pmax(rnorm(S, 1), 0), matrix(0, nrow = S, ncol = ncov_w))
        # theta11 <- matrix(NA, nrow = sum(M_site), ncol = S)
        # for (i in 1:n) {
        #   for (m in 1:M_site[i]) {
        #     for (j in 1:S) {
        #       theta11[m + sum(M_site[seq_len(i-1)]),j] <- 
        #         logistic(c(1, logz[i,j], X_w[m + sum(M_site[seq_len(i-1)]),]) %*% 
        #                    beta_theta[j,])
        #     }
        #   }
        # }
        beta_theta <- cbind(rep(0, S),pmax(rnorm(S, 1), 0), matrix(0, nrow = S, ncol = ncov_w))
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
      
      if(correctLambdaijk){
        lambda_ijk <- y[1:sum(M_site), 1:max(K), 1:(S + S_star),drop = F]
      } else {
        lambda_ijk <- array(NA, c(sum(M_site), max(K), S + S_star))
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
      }
      
      # other
      {
        
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
          
          # lambdatilde <- lambdatilde_true
          
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
          
          # lambdatilde <- rep(.5, S)
          
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
    # print(cbind(theta10, p_10))
    if(iter <= nburn){
      print(paste0("Chain = ",chain," - Burn-in Iteration = ",iter))  
    } else {
      print(paste0("Chain = ",chain," - Iteration = ",iter - nburn))
    }  
    
    # LAMBDA ----------------------------------------------
    
    if(updateLambda_CP){
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
    
    # LAMBDA ----------------------------------------------
    
    if(updateLambda_NP){
      # print("Update lambda")
      
      # lambda <- update_lambda(beta0, mu, lambda,
      #                              sigma_beta, sigma_mu,
      #                              exp(lambda_prior), sigma_lambda,
      #                              S_star, betaThetaEqual1)
      
      lambda <- update_lambda_NP(lambda_ijk, c_imk, mu,
                                 r_nb, v, u,
                                 lambda_prior, sigma_lambda)
      
    }
    
    # LAMBDA JOINT -------
    
    # if(updateLambdaJoint){
    #   print("Update lambda")
    #   
    #   # lambda <- update_lambda(beta0, mu, lambda,
    #   #                              sigma_beta, sigma_mu,
    #   #                              exp(lambda_prior), sigma_lambda,
    #   #                              S_star, betaThetaEqual1)
    #   
    #   list_lambda <- update_lambdajoint_cpp(lambda, v, u, lambda_ijk, c_imk, delta, gamma, sigma, 
    #                                         logz, mu, r_nb, sigma_gamma, M_site, X_w, beta_w, lambda_prior, 
    #                                         sigma_lambda, K, emptyTubes, S_star)
    #   lambda <- list_lambda$lambda
    #   v <- list_lambda$v
    #   
    #   lambda <- update_lambda_spikeins(mu, 
    #                                    lambda,
    #                                    u,
    #                                    r_nb,
    #                                    lambda_prior, 
    #                                    sigma_lambda,
    #                                    S_star)
    #   
    # }
    
    
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
                                     S_star)
    }
    
    # VU ---------
    
    if(updateUV){
      list_uv <- update_uv_poisgamma_cpp(u, v, logz, lambda, X_z, beta_theta, beta_z, beta0, 
                                         r_nb, mu, lambda_ijk, c_imk, delta, gamma, sigma, sigma_gamma, 
                                         sigma_u, M_site, X_w, beta_w, K, S_star)
      u <- list_uv$u
      v <- list_uv$v  
      lambda <- list_uv$lambda
      # lambdatilde <- list_uv$lambdatilde
    }
    
    # V ------------------------------------------
    
    if(updateV){
      # print("Update v")
      
      v <- update_v_poisgamma_cpp(v, logz, 
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
      # print("Update u")
      
      list_u <- update_u_poisgamma_cpp(v, u, lambda, beta0, beta_z, logz,
                                       mu, lambda_ijk, r_nb, X_w, beta_w, c_imk,
                                       delta, gamma, sigma_u, beta_theta, sigma, 
                                       sigma_gamma, M_site, 
                                       K, S_star)
      u <- list_u$u
      # lambda <- list_u$lambda
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
        beta_w <- update_betaw_cpp(beta_w, v, delta, logz, X_w, sigma, sigma_beta, M_site)  
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
      # print("Update sigma")
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
        
        v_pres <- (delta == 1) | (gamma == 1)
        list_deltagammac <- update_delta_c_d_rjmcmc(v_pres, y, v, lambda, r_nb,
                                                    M_site, K, #lambdatilde,
                                                    mu0, n0, pi0, mu_tilde, 
                                                    n_tilde, u, logz, X_w, beta_w,
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
        
        # update_delta_c_d_rjmcmc(delta, gamma, y, v, lambda, r_nb,
        #                         M_site, K, #lambdatilde,
        #                         mu0, n0, pi0, mu_tilde,
        #                         n_tilde, u, logz, X_w, beta_w,
        #                         sigma, mu, sigma_gamma, v_sd = 1,
        #                         p_11, p_10, theta11, theta10, emptyTubes,
        #                         S_star)$c_imk[75,1,2]
        # 
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
    
    # LAMBDA IJK ----
    
    if(updateLambdaijk){
      # print("Update lambda ijk")
      lambda_ijk <- update_lambdaijk(lambda, lambda_ijk, v, u, r_nb, c_imk, M_site, y, K,
                                     S_star)
    }
    
    # LAMBDA 0 ---------------------------------------------------------
    
    if(updateLambda0){
      # print("Update lambda0")
      
      list_lambda0 <- update_lambda0_mixt(y, c_imk, mu0, n0, sd_mu0 = .05, sd_n0 = .05)
      mu0 <- list_lambda0$mu0
      n0 <- list_lambda0$n0
      pi0 <- list_lambda0$pi0
    }
    
    # LAMBDA TILDE ---------------------------------------------------------
    
    if(updateLambdaTilde){
      # print("Update lambda tilde")
      
      # lambdatilde <- update_lambdatilde(y, c_imk, lambda)
      list_lambda_tilde <- update_lambda_tilde_NB(y, c_imk, mu_tilde, n_tilde, sd_mu0 = 1, sd_n0 = 1)
      mu_tilde <- list_lambda_tilde$mu_tilde
      n_tilde <- list_lambda_tilde$n_tilde
    }
    
    # P11 --------------------------------------------------------------
    
    if(updateP11){
      # print("Update p11")
      
      p_11 <- update_p_11_cpp(p_11, delta, gamma, c_imk, M_site,
                              K, a_p11, b_p11)
    }
    
    # P10 --------------------------------------------------------------
    
    if(updateP10){
      # print("Update p10")
      
      p_10 <- update_p_10_cpp(p_10, delta, gamma, c_imk, M_site, K, a_p10, b_p10, emptyTubes)
    }
    
    # UPDATE BETA THETA11 --------------------------------------------------------------
    
    if(updateBetaTheta){
      # print("Update beta theta")
      
      list_beta_theta <- update_betatheta11_cpp(logz, 
                                                beta_theta, 
                                                theta11, 
                                                delta[1:sum(M_site),],
                                                X_w, M_site,
                                                b_theta11 = c(1,rep(0, ncov_w)),
                                                B_theta11 = diag(sigma_beta_theta, nrow = 1 + ncov_w),
                                                !betaThetaEqual1)
      # b_theta11 = c(0,1,rep(0, ncov_w)),
      # B_theta11 = diag(sigma_beta_theta, nrow = 2 + ncov_w),
      # !betaThetaEqual1)
      beta_theta <- list_beta_theta$beta_theta
      theta11 <- list_beta_theta$theta11
      
    }
    
    # UPDATE THETA10 --------------------------------------------------------------
    
    if(updateTheta10){
      # print("Update theta10")
      
      theta10 <- update_theta10_cpp(theta10, delta, gamma, M_site, 
                                    a_theta0, b_theta0)
      
    }
    
    # WRITE RESULTS -----------------------------------------------------------
    
    if(iter > nburn & (iter - nburn) %% nthin == 0){
      trueIter <- (iter - nburn) / nthin
      
      # lambda0_output[trueIter] <- lambda0
      mu0_output_iter[trueIter] <- mu0
      n0_output_iter[trueIter] <- n0
      mutilde_output_iter[,trueIter] <- mu_tilde
      lambda_output_iter[,trueIter] <- lambda
      beta_w_output_iter[,,trueIter] <- beta_w
      # zeta_output_iter[trueIter,] <- zeta
      u_output_iter[,,trueIter] <- u
      # v_output_iter[,,trueIter] <- v
      logz_output_iter[,,trueIter] <- logz
      if(jSDM){
        Tau_output_iter[trueIter,,] <- Tau
      } else {
        tau_output_iter[,trueIter] <- tau
      }
      sigma_output_iter[,trueIter] <- sigma
      beta_theta_output_iter[,,trueIter] <- beta_theta
      # lambda_ijk_output_iter[,,,trueIter] <- lambda_ijk
      beta_z_output_iter[,,trueIter] <- beta_z
      beta0_output_iter[,trueIter] <- beta0
      mu_output_iter[,trueIter] <- mu
      theta11_output_iter[,,trueIter] <- theta11
      p_11_output_iter[,trueIter] <- p_11
      p_10_output_iter[,trueIter] <- p_10
      r_output_iter[,trueIter] <- r_nb
      
      # delta_output_iter[,,trueIter] <- delta
      # gamma_output_iter[,,trueIter] <- gamma
      # cimk_output_iter[,,,trueIter] <- c_imk
      
      # deltagammac_output_iter[,,1,trueIter] <- delta
      # deltagammac_output_iter[,,2,trueIter] <- gamma
      
    }
    
  }
  
  # copy chains output
  {
    mu0_output[chain,] <- mu0_output_iter
    n0_output[chain,] <- n0_output_iter
    lambda_output[chain,,] <- lambda_output_iter
    mutilde_output[chain,,] <- mutilde_output_iter
    u_output[chain,,,] <- u_output_iter
    # zeta_output[chain,,] <- zeta_output_iter
    logz_output[chain,,,] <- logz_output_iter
    # v_output[chain,,,] <- v_output_iter
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
    # lambda_ijk_output[chain,,,,] <- lambda_ijk_output_iter
    theta11_output[chain,,,] <- theta11_output_iter
    p_11_output[chain,,] <- p_11_output_iter
    p_10_output[chain,,] <- p_10_output_iter
    r_output[chain,,] <- r_output_iter
    # delta_output[chain,,,] <- delta_output_iter
    # gamma_output[chain,,,] <- gamma_output_iter
    # cimk_output[chain,,,,] <- cimk_output_iter
  }
  
}

# PLOT FUNCTION -------------------

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

# ERROR DIAGNOSTICS -------

j <- 1
i <- 1
m <- 1

chain <- 2
startIter <- 2250
lastiter <- 2500
ggplot() + 
  # geom_point(data = NULL, aes(x = startIter:lastiter, 
  #                             y = lambda_output[chain,j,startIter:lastiter]),
  #            color = "purple") + 
  geom_point(data = NULL, aes(x = startIter:lastiter, 
                              y = logz_output[chain,i,j,startIter:lastiter]),
             color = "blue") + 
  # geom_point(data = NULL, aes(x = startIter:lastiter, 
  #                             y = v_output[chain,i,j,startIter:lastiter]),
  #            color = "blue") + 
  # geom_point(data = NULL, aes(x = startIter:lastiter, 
  #                             y = tau_output[chain,j,startIter:lastiter]),
  #            color = "red") + 
  geom_point(data = NULL, aes(x = startIter:lastiter, 
                              y = sigma_output[chain,j,startIter:lastiter]),
             color = "red") + 
  geom_point(data = NULL, aes(x = startIter:lastiter, 
                              y = beta_theta_output[chain,j,2,startIter:lastiter]),
             color = "black") 
# geom_point(data = NULL, aes(x = startIter:lastiter, 
#                             y = beta_z_output[chain,1,j,startIter:lastiter]),
#            color = "green")
# +
#   geom_point(data = NULL, aes(x = startIter:lastiter, 
#                            y = lambda_output[chain,j,startIter:lastiter]),
# color = "green") 


# SPIKE-IN ANALYSIS -----

# covariates coefficients
{
  ncov_z <- 1
  beta_vars <- apply(beta_z_output[,ncov_z,,], 1, function(x){
    sd(x)
    # quantile(x, probs = c(0.025,0.975))
  })
  
  beta_means <- apply(beta_z_output[,ncov_z,,], 1, function(x){
    mean(x)
    # quantile(x, probs = c(0.025,0.975))
  })
  
  mean(beta_means - beta_z_true[ncov_z,])
  
  mean(beta_vars) 
  
  # BIAS
  # M = 2
  # S = 50 # 0 ->  0.1541339  # 5 -> 0.01238269 # 10 -> 0.04416577
  
  # VARIANCE
  # M = 2
  # S = 50 # 0 ->  0.158871  # 5 -> 0.1373485 # 10 -> 0.1433103
  
}

# biomass differences
{
  l_vars <- array(NA, dim = c(S, n * (n-1) / 2))
  l_bias <- array(NA, dim = c(S, n * (n-1) / 2))
  
  for (j in 1:S) {
    print(j)
    l <- 1
    for (i1 in 1:n) {
      for (i2 in seq_len(i1-1)) {
        l_vars[j,l] <- sd(logz_output[,i1,j,] - logz_output[,i2,j,])
        l_bias[j,l] <- mean((logz_output[,i1,j,] - logz_output[,i2,j,]) - 
                              (logz_true[i1,j] - logz_true[i2,j]))
        l <- l + 1
      }
    }
  }
  
  mean(l_vars)
  
  mean(l_bias)
  
  # BIAS
  # M = 2
  # S = 50 # 0 ->    # 5 ->  2.47284e-05 # 10 -> -0.0009967025
  
  # VARIANCE
  # M = 2
  # S = 50 # 0 ->    # 5 -> 0.2767587 # 10 -> 0.2780609
}

# plots on u
{
  u_pcis <- apply(u_output, c(2,3), function(x){
    quantile(x, probs = c(0.025, 0.5, 0.975))
  })
  
  subsetPoints <- 1:20
  pcr <- 1
  ggplot(data = NULL, aes(x = subsetPoints,
                          y = u_true[subsetPoints, pcr],
                          ymin = u_pcis[1,subsetPoints,pcr],
                          ymax = u_pcis[3,subsetPoints,pcr])) +
    geom_point(size = 2, color = "red") + geom_errorbar() +
    theme_bw() + xlab("Samples") + ylab("u") + ylim(c(-4,4)) +
    scale_x_discrete(breaks = subsetPoints, labels = subsetPoints)
  
  (average_width <- mean(apply(u_pcis, c(2,3), function(x){
    x[3] - x[1]
  })))
  # 1.569874
}

# plots on lambda for spike-in
{
  lambda_pcis <- apply(lambda_ijk_output, c(2, 3, 4), function(x){
    quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = T)
  })
  
  pcr <- 1
  species <- S + 1
  subsetPoints <- which(delta_true[,species] == 1)#1:30
  ggplot(data = NULL, aes(x = subsetPoints,
                          y = y[subsetPoints,pcr,species],
                          ymin = lambda_pcis[1,subsetPoints,pcr,species],
                          ymax = lambda_pcis[3,subsetPoints,pcr,species])) + 
    geom_point(size = 2, color = "red") +
    geom_errorbar() +
    theme_bw() + xlab("Samples") + ylab("u") + ylim(c(0,10000)) + 
    scale_x_discrete(breaks = subsetPoints, labels = subsetPoints)
  
  (average_width <- mean(apply(lambda_pcis[,,,1:(S + S_star)], c(2,3), function(x){
    x[3] - x[1]
  })))
  (average_width <- mean(apply(lambda_pcis[,,,1:S], c(2,3), function(x){
    x[3] - x[1]
  })))
  (average_width <- mean(apply(lambda_pcis[,,,S + 1:S_star], c(2,3), function(x){
    x[3] - x[1]
  })))
}

# betaz credible intervals
{
  ncov_z <- 1
  
  beta_CI <- apply(beta_z_output, c(2,3), function(x){
    quantile(x, probs = c(.05,.5,.95))
  })
  
  (average_width <- mean(apply(beta_CI[,ncov_z,], 2, function(x){
    x[3] - x[1]
  })))
  # 1.141678
}

# lambda ijk traceplot
{
  i <- 1
  m <- 1
  l <- m + sum(M_site[seq_len(i-1)])
  k <- 1
  j <- 1
  diagnosticsPlot( lambda_ijk_output[,l,k,j,]) 
  l <- 1
  
}

# DIAGNOSTICS CIMK --------

theta11_mean <- apply(theta11_output, c(2,3), mean)
delta_mean <- apply(delta_output, c(2,3), mean)
cimk1_mean <- apply(cimk_output, c(2,3,4), function(x){
  mean(x == 1)
})
cimk2_mean <- apply(cimk_output, c(2,3,4), function(x){
  mean(x == 2)
})
# cimk1_mean <- apply(cimk_output, c(2,3,4), function(x){
# mean(x == 1)
# })
logz_mean <- apply(logz_output, c(2,3), mean)
v_mean <- apply(v_output, c(2,3), mean)

j <- 1
# estimatedProbs <- rep(NA, sum(M_site))
# for (i in 1:n) {
#   for (m in 1:M_site[i]) {
#     estimatedProbs[m + sum(M_site[seq_len(i-1)])] <- 
#       mean(logistic(beta_theta_output[,j,1,] + beta_theta_output[,j,2,] * logz_output[,i,j,] + 
#                       beta_theta_output[,j,3,] * r[m + sum(M_site[seq_len(i-1)])]))
#   }
# }
# 
# 
# lambda_ijk_mean <- matrix(NA, nrow = sum(M_site), ncol = S)
# for (i in 1:n) {
#   for (m in 1:M_site[i]) {
#     for (k in 1:K[m + sum(M_site[seq_len(i-1)])]) {
#       lambda_ijk_mean <- mean(exp(lambda_output[,j,] + ))
#     }
#   }
# }
# 
c_imk_col <- as.data.frame(c_imk[,,j])
# c_imk_col <- cimk1_mean[,,j]
# c_imk_col2 <- cimk2_mean[,,j]

colnames(c_imk_col) <- paste("PCR_state",1:max(K))#c("PCR1_state","PCR2_state","PCR2_state")
lambda_imk_col <- as.data.frame(lambda_ijk[,,j])
colnames(lambda_imk_col) <- paste("PCR_state",1:max(K))
y_imk_col <- y[,,j]
v_current <- v[,j]
gamma_col <- gamma[,j]
colnames(y_imk_col) <- paste("PCR",1:max(K))#c("PCR1","PCR2","PCR2")
Presence <- delta[,j]
# colnames(Presence) <- "DNA_state"
View(
  cbind(data_short,
        y_imk_col,
        # rep(logz_mean[,j], M_site),
        # r,
        theta11_mean[,j],
        theta11_true[,j],
        delta_mean[,j],
        delta_true[,j],
        # v_current,
        # Presence,
        # lambda_imk_col,
        
        c_imk_col,
        gamma_col
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

# Covariates coefficients
{
  # z
  cov_z <- 1
  j <- 0
  j <- j+1
  diagnosticsPlot( beta_z_output[,cov_z,j,]) + 
    geom_hline(aes(yintercept = beta_z_true[cov_z,j]), color = "red")  +
    ylim(c(-3,3))
  
  # beta w
  cov_w <- 1
  j <- 0
  j <- j+1
  diagnosticsPlot( beta_w_output[,cov_w,j,]) + 
    geom_hline(aes(yintercept = beta_w_true[cov_w,j]), color = "red")  +
    ylim(c(-5,5))
  
  # beta theta w
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
  j <- 0
  j <- j+1
  diagnosticsPlot( lambda_output[,j,]) + 
    geom_hline(aes(yintercept = lambda_true[j]), color = "red", size = 2) #+ 
  # geom_hline(aes(yintercept = lambda_prior[j]), color = "green")
  
  # beta 0
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
