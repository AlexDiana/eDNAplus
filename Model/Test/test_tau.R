# simulate data
simulate_l_from_tau <- function(beta0, beta_z,
                                tau, X_z){
  
  logz <- cbind(1, X_z) %*% rbind(beta0, beta_z) +
    sapply(1:S, function(j) rnorm(n, 0, sd = tau[j]))
  
  logz
}

# params ----------

# fixed values
{
  n <- 100
  M_site <- rep(10, n)
  
  ncov_z <- 1
  
  S <- 2
  
  X_z <- matrix(rnorm(n), nrow = n)
  
  beta0 <- rep(0, S)
  beta_z <- matrix(1, nrow = ncov_z, ncol = S)
  
  
}

# starting values
{
  a_tau <- 10
  b_tau <- 10
  tau <- rep(1, S)
}

# mcmc -------------

niter <- 3000

# output
{
  tau_output <- array(NA, c(S, niter))
}

for (iter in 1:niter) {
  
  print(iter)
  
  # simulate data
  logz <- simulate_l_from_tau(beta0, beta_z, tau, X_z)
  
  # update l
  tau <- update_tau_cpp(tau, logz, X_z, beta_z, beta0,
                        a_tau, rep(b_tau, S)
                        # a_tau = .5, b_tau = 1 / a_tau
  )
  
  # output
  tau_output[,iter] <- tau
  
}

# 

j <- 1

library(ggplot2)
ggplot() + 
  geom_histogram(data = NULL, aes(x = tau_output[j,]^2, y = ..density..)) + #xlim(c(prior_mean,3)) + 
  stat_function(fun = dinvgamma, args = list(alpha = a_tau, 
                                             beta = b_tau))

cov(t(logz_output[i,,]))

