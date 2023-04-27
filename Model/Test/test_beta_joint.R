simulate_logz_from_beta <- function(X_z, beta0, beta_z, Sigma, Tau){

  XB <- cbind(1, X_z) %*% rbind(beta0, beta_z)

  logz <- rmtrnorm(XB, Sigma, Tau)

  logz
}

# DATA -------

library(MASS); library(Matrix); library(ggplot2)
library(eDNAPlus)

n <- 100
S <- 5

ncov_z <- 1
X_z <- matrix(runif(ncov_z * n), n, ncov_z)

beta0 <- rep(0, S)
beta_z <- matrix(0, ncov_z, S)

# beta_mean <- matrix(0, ncov_z + 1, S)
# beta_mean <- matrix(rnorm((ncov_z + 1) * S), ncov_z + 1, S)
sigma_beta <- 1

# Sigma_n <- diag(1, nrow = n)
# Sigma_S <- diag(1, nrow = S)
Sigma_n <- rWishart(1, n + 500, diag(1, nrow = n) / (n + 500))[,,1]
Sigma_S <- rWishart(1, S + 50, diag(1, nrow = S) / (S + 50))[,,1]

invSigma_S <- solve(Sigma_S)
invSigma_n <- solve(Sigma_n)

inv_chol_Sigma_n <- solve(t(chol(Sigma_n)))
inv_chol_Sigma_S <- solve(t(chol(Sigma_S)))

# MCMC ----

niter <- 20000

beta_z_output <- array(NA, dim = c(ncov_z, S, niter))

for (iter in 1:niter) {
  
  print(iter)
  # logz
  logz <- simulate_logz_from_beta(X_z, beta0, beta_z, Sigma_n, Sigma_S)
  
  # update 
  # list_beta_z <- update_betaz_CP_joint(beta0, beta_z, logz, 
  #                                      inv_chol_Sigma_n, inv_chol_Sigma_S,
  #                                      X_z, sigma_beta, beta_mean,
  #                                      T)
  list_beta_z <- update_betaz_CP_joint(beta0, beta_z, logz,
                                       invSigma_S, invSigma_n,
                                       X_z, sigma_beta, #beta_mean,
                                       T)
  beta0 <- as.vector(list_beta_z$beta0)
  beta_z <- list_beta_z$beta_z
  
  beta_z_output[,,iter] <- beta_z
}

apply(beta_z_output[1,,], 1, mean)
beta_mean[2,]

cov(t(beta_z_output[1,,]))

j <- 2
qplot(1:niter, beta_z_output[1,j,], geom = "line") + geom_hline(aes(yintercept = beta_mean[2,j]))
qplot(1:niter, beta_z_output[1,2,])
