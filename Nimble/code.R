
library(tidyverse); library(nimble); library(here)

nimbleCode_occupancy <- nimbleCode({
  
  if(ncov_psi > 0){
    Xbetapsi[1:n, 1:S] <-
      (X_psi[1:n,1:ncov_psi] %*% beta_psi[1:ncov_psi,1:S])[1:n,1:S]
  }
  
  # if(ncov_theta > 0){
  #   Xbetatheta[1:N, 1:S] <-
  #     (X_theta[1:N,1:ncov_theta] %*% beta_theta[1:ncov_theta,1:S])[1:N,1]
  # }
  
  for(j in 1:S){
    for(i in 1:n) {
      
      logit(psi[i,j]) <- Xbetapsi[i,j] 
      z[i,j] ~ dbern(psi[i,j])
      
      for(m in 1:M[i]){
        
        logit(theta[sumM[i] + m, j]) <- beta0_theta[j]
        delta[sumM[i] + m, j] ~ dbern(theta[sumM[i] + m,j] * z[i,j])
        
        y[sumM[i] + m, j] ~ 
          dbin(delta[sumM[i] + m, j] * p[j] + 
                 (1 - delta[sumM[i] + m,j]) * q[j], K)
      }
      
    }  
  }
  
  for(j in 1:S){
    
    beta0_theta[j] ~ dnorm(0, 1) 
    
    for(i in seq_len(ncov_psi)){
      beta_psi[i,j] ~ dnorm(0, 3)
    }
    
  }
  
})

load(here("Dataset","data_malaise.rda"))

PCR_table <- data$PCR_table
y1 <- PCR_table[,1:50]
y2 <- PCR_table[,51:100]
y3 <- PCR_table[,101:150]
y <- y1 + y2 + y3
y <- y > 0
mode(y) <- "numeric"

X_z <- data$X_z %>% 
  rename(Site = SiteName)

data_all <- data$infos %>% 
  left_join(., X_z) %>% 
  cbind(., y) %>% 
  filter(Site != "empty")

X_psi <- data_all %>% dplyr::select(be30,lg_DistRoad)
X_psi <- cbind(1, X_psi)

M_df <- data_all %>% group_by(Site) %>% 
  summarise(M = n())
M <- M_df$M
n <- nrow(M_df)

K <- 3
S <- ncol(y)
ncov_psi <- ncol(X_psi)

eDNAData <- list(y = y,
                 X_psi = X_psi)

constants <- list(n = n,
                  M = M,
                  N = sum(M),
                  S = S,
                  sumM = c(0, cumsum(M)[-n]),
                  K = K,
                  ncov_psi = ncov_psi)

# starting points
{
  delta_start <- y > 0
  mode(delta_start) <- "numeric"
  
  z <- t(sapply(1:n, function(i){
    sapply(1:S, function(j){
      as.numeric(any(delta_start[sum(M[seq_len(i-1)]) + 1:M[i], j] > 0)  )  
    })
  }))
  
  startingVariables <- list(
    z = z,
    delta = delta_start,
    beta_psi = matrix(0, ncov_psi, S),
    beta0_theta = rep(0, S),
    p = rep(.9, S),
    q = rep(.1, S))
}

MCMCsamples <- nimbleMCMC(code = nimbleCode_occupancy, 
                          constants = constants,
                          data = eDNAData, 
                          inits = startingVariables,
                          nchains = 1, niter = 10000,
                          summary = TRUE, 
                          monitors = c('beta_psi',
                                       'beta0_theta',
                                       'p','q'))

# occModel <- nimbleModel(code = nimbleCode_occupancy, constants = constants,
#                     data = eDNAData, inits = startingVariables, calculate = T)
# CoccModel <- compileNimble(occModel)
# 
# occConf <- configureMCMC(occModel, print = TRUE)
# occMCMC <- buildMCMC(occConf)
# CoccMCMC <- compileNimble(occMCMC)
# MCMCsamples <- runMCMC(CoccMCMC, niter = 1000)

# OUTPUT --------

samples <- MCMCsamples$samples
str(samples)

# plot covariates -------

beta_psi_samples <- samples[,grep("beta_psi\\[", colnames(samples))]

beta_psi_CI <- apply(beta_psi_samples, 2, function(x){
  quantile(x, probs = c(0.025, 0.5, 0.975))
})

ncov_psi <- 3

beta_psi_CI_array <- array(beta_psi_CI, dim = c(3, ncov_psi, S))

beta_psi_CI_cov1 <- t(beta_psi_CI_array[,1,])

data_plot <- beta_psi_CI_cov1 %>% 
  as.data.frame() %>% 
  mutate(Species = row_number()) 
# %>% 
#   pivot_longer(cols = starts_with("V"), 
#                names_to = "CI",
               # values_to = "value")

ggplot(data_plot, aes(x = Species,
                      y = V2, 
                      ymin = V1,
                      ymax = V3)) + geom_errorbar()


# false positives ---------

p_samples <- samples[,grep("p\\[", colnames(samples))]
q_samples <- samples[,grep("q\\[", colnames(samples))]

