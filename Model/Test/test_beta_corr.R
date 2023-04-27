simulate_logz_from_beta <- function(X_z, beta0, beta_z, Sigma){
  
  logz <- matrix(NA, n, S)
  for (i in 1:n) {
    logz[i,] <- mvrnorm(1, mu = as.vector(beta0 + X_z[i,] %*% beta_z), Sigma)
  }
  
  logz
}

# DATA -------

library(MASS); library(Matrix); library(ggplot2)

n <- 100
S <- 5

ncov_z <- 1
X_z <- matrix(runif(ncov_z * n), n, ncov_z)

beta0 <- rep(0, S)
beta_z <- matrix(0, ncov_z, S)

sigma_beta <- 1

Sigma <- rWishart(1, S + 1, diag(1, nrow = S))[,,1]

# MCMC ----

niter <- 1000

beta_z_output <- array(NA, dim = c(ncov_z, S, niter))

for (iter in 1:niter) {

  print(iter)
  # logz
  logz <- simulate_logz_from_beta(X_z, beta0, beta_z, Sigma)
  
  # update 
  list_beta_z <- update_betaz_CP_corr(beta0, beta_z, logz, Sigma, X_z, sigma_beta, T)
  beta0 <- as.vector(list_beta_z$beta0)
  beta_z <- list_beta_z$beta_z
  
  beta_z_output[,,iter] <- beta_z
}

cov(t(beta_z_output[1,,]))

qplot(1:niter, beta_z_output[1,1,])
