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

# fixed values
{
  n <- 3
  M_site <- rep(1, n)
  
  prior_logz <- 0
  prior_var <- 1
  
  S <- 99
  
  ncov_w <- 1
  
  X_z <- matrix(rnorm(n), nrow = n)
  X_w <- matrix(rnorm(sum(M_site)), nrow = sum(M_site))
  
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
  Tau <- as.matrix(bdiag(Tau_list))
  
  lambda <- rep(10, S)
  
  mu <- rep(0, S)
  v <- matrix(0, sum(M_site), S)
  
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
  
}

# starting values
{
  logz <- matrix(0, n, S)
}

# mcmc -------------

niter <- 10000

# output
{
  logz_output <- array(NA, c(n, S, niter))
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
  logz <- update_logz_corr_cpp(logz, beta0, X_z, beta_z, mu,
                               v, lambda, beta_theta, X_w, beta_w,
                               Tau, delta, gamma, sigma,
                               M_site, S_star, emptyTubes)
  
  # output
  logz_output[,,iter] <- logz
  
}

# 

i <- 1
j <- 1

prior_mean <- c(1, X_z[i,]) %*% c(beta0[j], beta_z[,j]) 

library(ggplot2)
qplot(1:niter, logz_output[i,j,]) + 
  geom_hline(aes(yintercept = prior_mean), color = "red")

c(prior_mean, mean(logz_output[i,j,]))
c(sqrt(Tau[j,j]),sd(logz_output[i,j,]))
c(Tau[j,j],var(logz_output[i,j,]))

ggplot() + 
  geom_histogram(data = NULL, aes(x = logz_output[i,j,], y = ..density..)) +  
  stat_function(fun = dnorm, args = list(mean = prior_mean, sd = sqrt(Tau[j,j])))
  # stat_function(fun = dnorm, args = list(mean = prior_mean, sd = sd(logz_output[i,j,])))

round(cov(t(logz_output[i,,])), 3) #- Tau

