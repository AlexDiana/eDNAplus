# simulate data
simulate_delta_v_from_l <- function(logz, sigma_gamma, sigma, 
                                    theta10, beta_theta, r, alpha){
  
  theta11 <- matrix(NA, nrow = sum(M_site), ncol = S)
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        theta11[m + sum(M_site[seq_len(i-1)]),j] <- 
          logistic(c(1, exp(logz[i,j]), r[m + sum(M_site[seq_len(i-1)])]) %*% 
                     beta_theta[j,])
      }
    }
  }
  
  delta <- array(NA, dim = c(sum(M_site), S))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        if(emptySites[i] == 0){
          delta[m + sum(M_site[seq_len(i-1)]),j] <- 
            rbinom(1, 1, theta11[m + sum(M_site[seq_len(i-1)]),j])  
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
    if(emptySites[i] == 0){
      for (m in 1:M_site[i]) {
        for (j in 1:S) {
          if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1){
            v[m + sum(M_site[seq_len(i-1)]),j] <- 
              rnorm(1, logz[i,j], sigma[j]) +
              alpha[j] * r[m + sum(M_site[seq_len(i-1)])] 
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
  }
  
  list("delta" = delta,
       "gamma" = gamma,
       "v" = v)
}

# params ----------

# fixed values
{
  prior_logz <- 0
  prior_var <- 1
  
  S <- 2
  
  r <- rep(0, sum(M_site))
  
  Tau <- rcorrmatrix(S)
  # Tau <- diag(1, nrow = S)
  
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
  
  theta10 <- rep(.01, S)
  
  sigma <- rep(1, S)
  
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
  gamma <- list_logz$gamma
  
  # update l
  # list_logz_bar <- update_logz_corr(beta0, beta_z, lambda, mu, logz, v, Tau, delta, X_z)
  list_logz_bar <- update_logz_cpp(logz, beta0, X_z, beta_z,
                                   mu, v, lambda, beta_theta,
                                   r, alpha, diag(Tau), delta, 
                                   gamma, sigma, M_site,
                                   emptySites, sigma_prop = (.01)^2)
  logz <- list_logz_bar$logz
  v <- list_logz_bar$v
  
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

ggplot() + 
  geom_histogram(data = NULL, aes(x = logz_output[i,j,], y = ..density..)) + #xlim(c(prior_mean,3)) + 
  stat_function(fun = dnorm, args = list(mean = prior_mean, sd = sqrt(Tau[j,j])))

cov(t(logz_output[i,,]))

