# simulate data
simulate_v_from_mu <- function(mu, sigma_gamma, sigma, 
                                    delta, gamma, v_spikes){
  
  # log amount of DNA
  v <- matrix(NA, sum(M_site), S + S_star)
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1){
          # v[m + sum(M_site[seq_len(i-1)]),j] <- 
          #   rnorm(1, logz_true[i,j], sigma_true[j]) +
          #   Xw_betaw[m + sum(M_site[seq_len(i-1)]),j] 
        } else if (gamma[m + sum(M_site[seq_len(i-1)]),j] == 1){
          v[m + sum(M_site[seq_len(i-1)]),j] <- 
            rnorm(1, mu[j], 
                  sd = sigma_gamma)
        } 
        # else {
        #   v_true[m + sum(M_site[seq_len(i-1)]),j] <- 0
        # }
      }
      for (j in seq_len(S_star)) {
        v[m + sum(M_site[seq_len(i-1)]), S + j] <- 
          v_spikes[m + sum(M_site[seq_len(i-1)]),j]
      }
    }
  }
  
  v
}

# params ----------

# fixed values
{
  prior_logz <- 0
  prior_var <- 1
  
  S <- 2
  lambda <- rep(10, S)
  
  mu <- rep(0, S)
  v <- matrix(0, sum(M_site), S)
  
  theta10_true <- rbeta(S, a_theta0, b_theta0)
  
  beta0 <- rep(0, S)
  alpha <- rep(0, S)
  beta_z <- matrix(1, nrow = ncov_z, ncol = S)
  
  logz <- matrix(0, n, S)
  
  beta_theta <- cbind(rep(0, S), rep(1, S), rep(0, S))
  
  delta <- array(NA, dim = c(sum(M_site), S))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        if(emptySites[i] == 0){
          delta[m + sum(M_site[seq_len(i-1)]),j] <- 
            logistic(c(1, logz[i,j], r[m + sum(M_site[seq_len(i-1)])]) %*% 
                       beta_theta_true[j,])
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
  
  
}

# starting values
{
  logz <- matrix(0, n, S)
}

# mcmc -------------

niter <- 3000

# output
{
  logz_output <- array(NA, c(n, S, niter))
}

for (iter in 1:niter) {
  
  print(iter)
  
  # simulate data
  list_logz <- simulate_delta_v_from_l(logz, sigma_gamma, sigma, theta10, beta_theta, r, alpha)
  v <- list_logz$v
  delta <- list_logz$delta
  
  # update l
  list_logz_bar <- update_logz_corr(beta0, beta_z, lambda, mu, logz, v, Tau, delta, X_z)
  # list_logz_bar <- update_logz_cpp(logz, beta0, X_z[,-1,drop=F], beta_z,
  #                                  mu, v, lambda, beta_theta,
  #                                  r, alpha, diag(Tau), delta, gamma, sigma, M_site,
  #                                  emptySites, sigma_prop = (.01)^2)
  logz <- list_logz_bar$logz
  v <- list_logz_bar$v
  
  # output
  logz_output[,,iter] <- logz
  
}

# 

i <- 1
j <- 1

prior_mean <- X_z[i,] %*% c(beta0[j], beta_z[,j]) 

library(ggplot2)
qplot(1:niter, logz_output[i,j,]) + 
  geom_hline(aes(yintercept = prior_mean), color = "red")

c(prior_mean, mean(logz_output[i,j,]))

ggplot() + 
  geom_histogram(data = NULL, aes(x = logz_output[i,j,], y = ..density..)) + #xlim(c(prior_mean,3)) + 
  stat_function(fun = dnorm, args = list(mean = prior_mean, sd = sqrt(Tau[j,j])))

cov(t(logz_output[i,,]))

