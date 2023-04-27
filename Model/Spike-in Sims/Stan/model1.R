library(here); library("rstan") 
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
setwd("~/eDNAPlus/Model/Spike-in Sims/Stan")

# DATA ------

n <- 2
S <- 10
S_star <- 5
M <- 4
K <- 3

# prior
{
  sigma <- 1
  tau <- 100
  sigma_u <- 1
  sigma_y <- 1
  sigma_lambda <- 1
}

logz_true <- sapply(1:S, function(j) rnorm(n, 0, sd = tau))

v_star <- matrix(0, n * M, S_star)

v_true <- matrix(NA, n * M, S)
for (i in 1:n) {
  for (m in 1:M) {
    for (j in 1:S) {
      v_true[m + (i - 1)*M,j] <- 
        rnorm(1, logz_true[i,j], 
              sigma) 
    }
  }
}

u_true <- matrix(NA, n * M, K)
for (i in 1:n) {
  for (m in 1:M) {
    for (k in 1:K) {
      u_true[m + (i - 1)*M,k] <- 
        rnorm(1, 0, sigma_u) 
    }
  }
}

lambda_true <- rep(NA, S + S_star)
for (j in 1:(S + S_star)) {
  lambda_true[j] <- rnorm(1, 0, sigma_lambda)
}

y <- array(NA, dim = c(n * M, K, S + S_star))
for (i in 1:n) {
  for (m in 1:M) {
    for (k in 1:K) {
      for (j in 1:S) {
        y[m + (i - 1)*M,k,j] <- 
          rnorm(1, lambda_true[j] + v_true[(i - 1)*M + m,j] + u_true[(i - 1)*M + m,k], 
                sigma_y) 
      }  
      for (j in 1:S_star) {
        y[m + (i - 1)*M,k,S + j] <- 
          rnorm(1, lambda_true[S + j] + v_star[(i - 1)*M + m,j] + u_true[(i - 1)*M + m,k], 
                sigma_y) 
      }  
    }
  }
}

# FIT ------



edna_dat <- list(n = n,
                 S = S,
                 S_star = S_star,
                 M = M,
                 K = K,
                 y = y,
                 v_star = v_star,
                 sigma = sigma,
                 sigma_y = sigma_y,
                 sigma_u = sigma_u,
                 sigma_lambda = sigma_lambda,
                 tau = tau
)

# MCMC ----

# occ_model <- stanc(file = 'model.stan', model_name = "occmodel")
# mod <- stan_model(stanc_ret = occ_model, verbose = TRUE)
fit <- sampling(mod, data = edna_dat, chains = 1, iter = 5000, verbose = T)
# fit <- stan(file = 'occupancy.stan', data = occ_dat, verbose = T)

matrix_of_draws <- as.matrix(fit)

var( matrix_of_draws[,1] - matrix_of_draws[,2])
var( matrix_of_draws[,3] - matrix_of_draws[,4])
var( matrix_of_draws[,5] - matrix_of_draws[,6])
var( matrix_of_draws[,7] - matrix_of_draws[,8])
var( matrix_of_draws[,9] - matrix_of_draws[,10])

sigma_ratio <- sigma_u^2 / sigma_y^2
2 * (1 / M) * (sigma^2 + (sigma_y^2 / K) * (1 + (sigma_ratio / (sigma_ratio * S_star + 1))))


colnames(matrix_of_draws)
niter <- nrow(matrix_of_draws)

sd(matrix_of_draws[,6])

library(ggplot2)
j <- 3
ggplot(data = NULL) + geom_point(aes(x = 1:niter, 
                                     y = matrix_of_draws[,(j - 1)*(ncov_z + 1) + 2])) + 
  geom_hline(aes(yintercept = as.vector(t(beta_z_true))[j])) 

qplot(1:niter, matrix_of_draws[,6])

quantile(matrix_of_draws[,1], probs = c(0.025, 0.975))
