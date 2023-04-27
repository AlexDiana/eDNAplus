# simulate data
simulate_logz_from_tau <- function(tau, 
                                   X_z, 
                                   beta0, 
                                   beta_z){
  
  logz <- X_z %*% rbind(beta0, beta_z) +
    sapply(1:S, function(j) rnorm(n, 0, sd = tau[j]))
  
  logz
}

# starting values ----------

# fixed values
{
  n <- 100
  emptySites <- rep(0, n)
  
  beta0 <- rep(0, S)
  beta_z <- matrix(0, nrow = ncov_z, ncol = S)
  
  X_z <- cbind(1, matrix(runif(n * ncov_z), n, ncov_z))
}

# params
{
  tau <- rep(1, S)
}

# mcmc -------------

niter <- 5000

# output
{
  tau_output <- array(NA, c(S, niter))
}

for (iter in 1:niter) {
  
  print(iter)
  
  # simulate data ---------
  
  logz <- simulate_logz_from_tau(tau, 
                                 X_z, 
                                 beta0, 
                                 beta_z)
  
  # update tau --------
  
  tau <- update_tau_cpp(tau, logz, X_z[,-1,drop=F], beta_z, beta0,
                        a_tau = 10, b_tau = rep(10, S), emptySites
                        # a_tau = .5, b_tau = 1 / a_tau
  )
  
  
  # output
  tau_output[,iter] <- tau
  
}

# 

tau_samples <- sapply(1:1000, function(i){
  rinvgamma_cpp(10, 10)
}) 

j <- 1
library(ggplot2)
ggplot() + geom_histogram(data = NULL, aes(x = tau_output[j,]^2,
                                           y = ..density..)) + 
  stat_function(fun = dinvgamma, args = list(alpha = 10, 
                                             beta = 10)) + 
  xlim(c(0,3))
