# simulate data
simulate_delta_v_from_l <- function(logz, sigma_gamma, sigma, 
                                    theta10, beta_theta, X_w, beta_w){
  
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
  
  delta <- array(NA, dim = c(sum(M_site), S))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        delta[m + sum(M_site[seq_len(i-1)]),j] <- 
          rbinom(1, 1, theta11[m + sum(M_site[seq_len(i-1)]),j])  
      } 
      
    }
  }
  
  
  gamma <- array(NA, dim = c(sum(M_site), S))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        if(delta[m + sum(M_site[seq_len(i-1)]),j] == 0){
          gamma[m + sum(M_site[seq_len(i-1)]),j] <- rbinom(1, 1, theta10[j])  
        } else if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1){
          gamma[m + sum(M_site[seq_len(i-1)]),j] <- 0
        }
      }
    }
  }
  
  v <- matrix(NA, sum(M_site), S)
  for (i in 1:n) {
    
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1){
          v[m + sum(M_site[seq_len(i-1)]),j] <- 
            rnorm(1, logz[i,j], sigma[j]) +
            beta_w[,j] * X_w[m + sum(M_site[seq_len(i-1)]),] 
        } else if (gamma[m + sum(M_site[seq_len(i-1)]),j] == 1){
          v[m + sum(M_site[seq_len(i-1)]),j] <- 
            rnorm(1, mu[j], 
                  sd = sigma_gamma)
        } else {
          v[m + sum(M_site[seq_len(i-1)]),j] <- 0
        }
      }
    }
    
  }
  
  list("delta" = delta,
       "gamma" = gamma,
       "v" = v)
}

# params ----------

jointSpecies <- T
spatialCorr <- T

# fixed values
{
  n <- 2
  M_site <- rep(1, n)
  
  prior_logz <- 0
  prior_var <- 1
  
  S <- 3
  
  ncov_z <- 1
  ncov_w <- 1
  
  X_z <- matrix(rnorm(n), nrow = n)
  X_w <- matrix(rnorm(sum(M_site)), nrow = sum(M_site))
  
  sizeBlocks <- 3
  
  taus <- 1
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
  Tau <- as.matrix(bdiag(Tau_list))
  
  df_t <- n + 10
  Sigma_n <- rWishart(1, df_t, diag(1, nrow = n) / df_t)[,,1]
  # Sigma_n <- diag(1, nrow = n)#rWishart(1, df_t, diag(1, nrow = n) / df_t)[,,1]
  invSigma_n <- solve(Sigma_n)
  
  lambda <- rep(10, S)
  
  mu <- rep(0, S)
  v <- matrix(0, sum(M_site), S)
  
  a_theta0 <- 1
  b_theta0 <- 20
  theta10_true <- rbeta(S, a_theta0, b_theta0)
  
  beta0 <- rep(0, S)
  beta_z <- matrix(0, nrow = ncov_z, ncol = S)
  beta_w <- matrix(0, nrow = ncov_w, ncol = S)
  
  logz <- matrix(0, n, S)
  
  beta_theta <- cbind(rep(0, S), rep(1, S), rep(0, S))
  
  delta <- array(NA, dim = c(sum(M_site), S))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
          delta[m + sum(M_site[seq_len(i-1)]),j] <- 
            logistic(c(1, logz[i,j], X_w[m + sum(M_site[seq_len(i-1)]),]) %*% 
                       beta_theta[j,])
        
          
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
  
  theta10 <- rep(.01, S)
  
  sigma <- rep(.5, S)
  
  sigma_gamma <- 1
  
  S_star <- 0
  emptyTubes <- 0
}

# starting values
{
  logz <- matrix(0, n, S)
  
  inv_chol_Sigma_n <- solve(t(chol(Sigma_n)))
  inv_chol_Sigma_S <- solve(t(chol(Tau)))
}

# mcmc -------------

niter <- 20000

# output
{
  logz_output <- array(NA, c(n, S, niter))
  logz_vec_output <- array(NA, c(n * S, niter))
}

for (iter in 1:niter) {
  
  print(iter)
  
  # simulate data
  list_logz <- simulate_delta_v_from_l(logz, sigma_gamma, sigma, theta10, beta_theta, X_w, beta_w)
  v <- list_logz$v
  delta <- list_logz$delta
  gamma <- list_logz$gamma
  
  # update l
  # list_logz_bar <- update_logz_corr(beta0, beta_z, lambda, mu, logz, v, Tau, delta, X_z)
  logz <- update_logz_joint_cpp(logz, beta0, X_z, beta_z, mu,
                               v, lambda, beta_theta, X_w, beta_w,
                               Sigma_n, invSigma_n, Tau, delta, gamma, sigma,
                               M_site, S_star, emptyTubes)
  # logz <- update_logz_joint_fast_cpp(logz, beta0, X_z, beta_z, mu,
  #                                    v, lambda, beta_theta, X_w, beta_w,
  #                                    Sigma_n, Tau, inv_chol_Sigma_n, inv_chol_Sigma_S,
  #                                    delta, gamma, sigma,
  #                                    M_site, S_star, emptyTubes)
  
  # output
  logz_output[,,iter] <- logz
  logz_vec_output[,iter] <- as.vector(logz)
  
}

#  MEAN

i <- 1
j <- 1

prior_mean <- c(1, X_z[i,]) %*% c(beta0[j], beta_z[,j]) 

library(ggplot2)
qplot(1:niter, logz_output[i,j,]) + 
  geom_hline(aes(yintercept = prior_mean), color = "red")


# POST COV

Sigma_tilde_hat <- cov(t(logz_vec_output))
Sigma_tilde <- kronecker(Tau, Sigma_n)

Sigma_tilde_hat
Sigma_tilde

max(abs(Sigma_tilde - Sigma_tilde_hat))

i <- 1
j <- 1
ggplot() + 
  geom_histogram(data = NULL, aes(x = logz_output[i,j,], y = ..density..)) +  
  stat_function(fun = dnorm, 
                args = list(mean = prior_mean, sd = sqrt(Sigma_tilde[i,i])))


c(prior_mean, mean(logz_output[i,j,]))
c(sqrt(Tau[j,j]),sd(logz_output[i,j,]))
c(Tau[j,j],var(logz_output[i,j,]))

ggplot() + 
  geom_histogram(data = NULL, aes(x = logz_output[i,j,], y = ..density..)) +  
  stat_function(fun = dnorm, args = list(mean = prior_mean, sd = sqrt(Tau[j,j])))
  # stat_function(fun = dnorm, args = list(mean = prior_mean, sd = sd(logz_output[i,j,])))

round(cov(t(logz_output[i,,])), 3) #- Tau

