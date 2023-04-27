library(Matrix)

# SIMULATION --------

jointSpecies <- T
spatialCorr <- T

# design
{
  S <- 10
  n <- 100
  ncov_z <- 0
  ncov_w <- 0
  
  M_site <- rep(2, n)
  emptyTubes <- 5
  
  K <- rep(3, sum(M_site) + emptyTubes)
  
  S_star <- 2 # numOfSpikes
  
  data_short <- data.frame(Site = rep(1:n, M_site),
                           Sample = 1:sum(M_site))
  
  if(emptyTubes > 0){
    data_short <- rbind(data_short,
                        data.frame(Site = 0,
                                   Sample = 1:emptyTubes))
    
  }
  
  if(spatialCorr){
    X_s <- cbind(runif(n, 0, 1), runif(n, 0, 1))  
  } else {
    X_s <- NULL
  }
  
  # ggplot(data = NULL, aes(x = X_s[,1],
  #                         y = X_s[,2])) + geom_point()
  
  X_z <- matrix(runif(n * ncov_z), n, ncov_z)
  if(ncov_z > 0) colnames(X_z) <- paste("Z",1:ncov_z)
  
  X_w <- matrix(runif(sum(M_site) * ncov_w), sum(M_site), ncov_w)
  if(ncov_w > 0) colnames(X_w) <- paste("W",1:ncov_w)
  
  v_spikes <- matrix(0, nrow = sum(M_site) + emptyTubes, ncol = S_star)
  spikedSample <- rbind(matrix(1, nrow = sum(M_site), ncol = S_star),
                        matrix(rbinom(emptyTubes * S_star, 1, .5), nrow = emptyTubes, ncol = S_star))
  
}

# prior
{
  beta0_mean <- 0
  beta_theta_0_mean <- -1.5
  sigmas <- .5
  taus <- .5
  l_gp <- .1
  r_0 <- 100
  p11_0 <- .95
  p10_0 <- .02
  theta_10 <- .02
  
  sd_beta_theta_0 <- .5
  sigma_beta0 <- 1
  sigma_u <- 1
  sigma_gamma <- 1
}

simulatedData <- T

beta0_true <- rep(0, S)#rnorm(S, beta0_mean, sd = sigma_beta0)#

beta_z_true <- matrix(sample(c(1,-1), (ncov_z) * S, replace = T), 
                      nrow = ncov_z, ncol = S, byrow = T)

sigma_true <- rep(sigmas, S)#pmin(rhcauchy(S, 2), .1)

if(jointSpecies) {
  
  # sizeBlocks <- 3
  # 
  # Tau_list <- lapply(1:(S/sizeBlocks), function(j){
  #   cormat <- matrix(1, sizeBlocks, sizeBlocks)
  #   repeat {
  #     for (i in 2:sizeBlocks) {
  #       for (j in seq_len(i-1)) {
  #         cormat[i,j] <- sample(c(-1,1) * .8, size = 1)
  #         cormat[j,i] <- cormat[i,j]
  #       }
  #     }
  #     if(all(eigen(cormat)$values > 0)){
  #       break
  #     }
  #   }
  #   
  #   cormat <- cormat * taus
  # })
  
  sizeBlocks <- c(3,2,3,2)
  
  Tau_list <- lapply(1:length(sizeBlocks), function(j){
    cormat <- matrix(1, sizeBlocks[j], sizeBlocks[j])
    repeat {
      for (i in 2:sizeBlocks[j]) {
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

lambda_true <- rnorm(S + S_star, 15, sd = .5)
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

XB <- cbind(1, X_z) %*% rbind(beta0_true, beta_z_true)

if(!spatialCorr){
  Sigma <- diag(1, nrow = n, n)
} else {
  Sigma <- K2(X_s, X_s, 1, l_gp)  
} 

if(!jointSpecies){
  Tau <- diag(tau_true, nrow = S)
} else {
  Tau <- Tau_true
}
logz_true <- rmtrnorm(XB, Sigma, Tau)

# # old
# if(jointSpecies){
#   logz_true <- matrix(NA, n, S)
#   for (i in 1:n) {
#     Xb_i <-  c(1, X_z[i,]) %*% rbind(beta0_true, beta_z_true)
#     logz_true[i,] <- mvrnorm(1, Xb_i, Tau_true)
#     
#   }
# } else {
#   # logz_true <- X_z %*% beta_z_true +
#   logz_true <- cbind(1, X_z) %*% rbind(beta0_true, beta_z_true) +
#     sapply(1:S, function(j) rnorm(n, 0, sd = tau_true[j]))
# }

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

Xw_betaw <- X_w %*% beta_w_true

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

# write data ----

setwd("~/eDNAPlus/Dataset/Simulated")

OTUnames <- paste0("OTU_",1:(S+S_star))

data_infos <- data.frame(Site = c(rep(1:n, M_site),
                                  rep("empty", emptyTubes)),
                         Sample = 1:(sum(M_site) + emptyTubes),
                         Replicates = K)

PCR1 <- y[,1,1:S]
colnames(PCR1) <- OTUnames[1:S]
PCR2 <- y[,2,1:S]
colnames(PCR2) <- OTUnames[1:S]
PCR3 <- y[,3,1:S]
colnames(PCR3) <- OTUnames[1:S]
PCR <- cbind(data_infos, PCR1, PCR2, PCR3)

X_z_1 <- data.frame(Site = 1:n,
                      X_z)

X_z_all <- X_z_1[match(data_infos$Site, X_z_1$Site),-1]

X_w_withempty <- rbind(X_w, matrix(NA, nrow = emptyTubes, ncol = ncol(X_w)))
X_w_all <- data.frame(Sample = 1:(sum(M_site) + emptyTubes),
                      X_w_withempty)
X_w_all <- X_w_all[,-1]

PCR_X <- cbind(PCR, X_z_all, X_w_all)

write.csv(PCR_X, file = "PCR.csv", row.names = F)

PCR1_spike <- y[,1,S + seq_len(S_star)]
colnames(PCR1_spike) <- OTUnames[S + seq_len(S_star)]
PCR2_spike <- y[,2,S + seq_len(S_star)]
colnames(PCR2_spike) <- OTUnames[S + seq_len(S_star)]
PCR3_spike <- y[,3,S + seq_len(S_star)]
colnames(PCR3_spike) <- OTUnames[S + seq_len(S_star)]
PCR_spike <- cbind(PCR1_spike, PCR2_spike, PCR3_spike)

write.csv(PCR_spike, file = "PCR_spike.csv", row.names = F)

