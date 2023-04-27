# simulate data
simulate_v_from_sigma <- function(logz, sigma, X_w, beta_w){
  
  v <- matrix(NA, sum(M_site), S)
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1){
          v[m + sum(M_site[seq_len(i-1)]),j] <- 
            rnorm(1, logz[i,j], sigma[j]) +
            beta_w[,j] * X_w[m + sum(M_site[seq_len(i-1)]),] 
        } else if (gamma[m + sum(M_site[seq_len(i-1)]),j] == 1){
          # v[m + sum(M_site[seq_len(i-1)]),j] <- 
          #   rnorm(1, mu[j], 
          #         sd = sigma_gamma)
        } else {
          # v[m + sum(M_site[seq_len(i-1)]),j] <- 0
        }
      }
    }
    
  }
  
  v
}

# params ----------

# fixed values
{
  n <- 100
  M_site <- rep(10, n)
  
  ncov_w <- 1
  
  S <- 2
  
  logz <- matrix(0, n, S)
  
  X_w <- matrix(rnorm(sum(M_site)), nrow = sum(M_site))
  beta_w <- matrix(1, nrow = ncov_w, ncol = S)
  
  lambda <- rep(10, S)
}

# starting values
{
  a_sigma <- 10
  b_sigma <- 10
  sigma <- rep(1, S)
}

# mcmc -------------

niter <- 3000

# output
{
  sigma_output <- array(NA, c(S, niter))
}

for (iter in 1:niter) {
  
  print(iter)
  
  # simulate data
  v <- simulate_v_from_sigma(logz, sigma, X_w, beta_w)
  
  # update l
  sigma <- update_sigma_cpp(sigma, lambda, beta_z, beta0,
                            mu, logz, v, X_w, beta_w,
                            delta, gamma, beta_theta,
                            a_sigma, rep(b_sigma, S),
                            # a_sigma = .5, b_sigma = 1 / a_sigma,
                            M_site, S_star)
  
  # output
  sigma_output[,iter] <- sigma
  
}

# 

j <- 2

library(ggplot2)
ggplot() + 
  geom_histogram(data = NULL, aes(x = sigma_output[j,]^2, y = ..density..)) + #xlim(c(prior_mean,3)) + 
  stat_function(fun = dinvgamma, args = list(alpha = a_sigma, 
                                             beta = b_sigma))

cov(t(logz_output[i,,]))

