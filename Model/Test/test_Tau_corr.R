# simulate data
simulate_l_from_Tau <- function(X_z, beta0, beta_z, Tau){
  
  logz <- matrix(NA, n, S)
  for (i in 1:n) {
    Xb_i <-  c(1, X_z[i,]) %*% rbind(beta0, beta_z)
    logz[i,] <- mvrnorm(1, Xb_i, Tau)
  }
  
  logz
}

# params ----------

jointSpecies <- T

# fixed values
{
  n <- 200
  
  ncov_z <- 1
  
  S <- 3
  
  X_z <- matrix(rnorm(n), nrow = n)
  
  beta0 <- rep(0, S)
  beta_z <- matrix(1, nrow = ncov_z, ncol = S)
  
  lambda_Y <- 1
  Tau_Priors <- list("lambda_Y" = lambda_Y)
  
}

# starting values
{
  Tau_params <- list("Omega" = diag(1, nrow = S),
                     "lambdasq" = matrix(1, nrow = S, ncol = S),
                     "tausq" = 1,
                     "nu" = matrix(1, nrow = S, ncol = S),
                     "csi" = 1,
                     "Sigma" = diag(1, nrow = S))
}

# mcmc -------------

niter <- 20000

# output
{
  Tau_output <- array(NA, c(S, S, niter))
}

for (iter in 1:niter) {
  
  print(iter)
  
  # simulate data
  logz <- simulate_l_from_Tau(X_z, beta0, beta_z, Tau_params$Omega)
  
  # update l
  Tau_params <- update_Tau(X_z, logz, beta0, beta_z,
                           Tau_params, Tau_Priors,
                           jointSpecies)
  
  # output
  Tau_output[,,iter] <- Tau_params$Omega
  
}

# DIAGONAL ------

j <- 3

library(ggplot2)
ggplot() + 
  geom_histogram(data = NULL, aes(x = Tau_output[j,j,], y = ..density..), bins = 100) + #xlim(c(prior_mean,3)) + 
  stat_function(fun = dexp, args = list(rate = lambda_Y / 2)) + 
  coord_cartesian(xlim = c(0, 10))

sapply(1:S, function(j){
  mean( Tau_output[j,j,])
})

(1 / (lambda_Y / 2))

qplot(1:niter, Tau_output[j,j,])

# OFF-DIAGONAL ------

j1 <- 2
j2 <- 3

library(ggplot2)
ggplot() + 
  geom_histogram(data = NULL, aes(x = Tau_output[j1,j2,], y = ..density..), bins = 100) + #xlim(c(prior_mean,3)) + 
  stat_function(fun = dexp, args = list(rate = lambda_Y)) + 
  coord_cartesian(xlim = c(-10, 10))

qplot(1:niter, Tau_output[j1,j2,])

