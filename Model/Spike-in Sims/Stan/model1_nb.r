library(here); library("rstan") 
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
setwd("~/eDNAPlus/Model/Spike-in Sims/Stan")
occ_model <- stanc(file = 'model_nb.stan', model_name = "occmodel")
mod <- stan_model(stanc_ret = occ_model, verbose = TRUE)

# DATA ------

# study design
{
  # n <- 2
  # S <- 2
  # S_star <- 5
  # M <- 4
  # K <- 3  
}

# data generating process
{
  # sigma <- 1
  # tau <- 1
  sigma_u <- 1
  # r <- 1000
  sigma_lambda <- 1
  mean_lambda <- 7
}

nsims <- 5

S_star_grid <- c(3, 2, 1, 0)

settings <- expand.grid(M = c(1,3,5),
                        K = c(1,3),
                        tau = c(.5, 1),
                        sigma = c(.5, 1))

data <- expand.grid(M = c(1,3,5),
                    K = c(1,3),
                    tau = c(.5, 1),
                    sigma = c(.5, 1),
                    S_star = S_star_grid)

var_ratio_all <- rep(NA, nrow = nrow(data))

for (idx_setting in 1:nrow(settings)) {
  
  var_ratios <- matrix(NA, nsims, length(S_star_grid))
  mae_ratios <- matrix(NA, nsims, length(S_star_grid))
  
  for (sim in 1:nsims) {
    
    # data generating parameters
    {
      n <- 5
      S <- 3
      S_star <- max(S_star_grid)
      M <- settings$M[idx_setting]
      K <- settings$K[idx_setting]
      sigma <- settings$sigma[idx_setting]
      tau <- settings$tau[idx_setting]
      sigma_u <- 1
      r <- 1000
      sigma_lambda <- 1
    }
    
    # simulate data
    {
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
        lambda_true[j] <- rnorm(1, mean_lambda, sigma_lambda)
      }
      
      y <- array(NA, dim = c(n * M, K, S + S_star))
      for (i in 1:n) {
        for (m in 1:M) {
          for (k in 1:K) {
            for (j in 1:S) {
              y[m + (i - 1)*M,k,j] <- 
                rnbinom(1, mu = exp(lambda_true[j] + v_true[(i - 1)*M + m,j] + u_true[(i - 1)*M + m,k]), 
                        size = r) 
            }  
            for (j in 1:S_star) {
              y[m + (i - 1)*M,k,S + j] <- 
                rnbinom(1, mu = exp(lambda_true[S + j] + v_star[(i - 1)*M + m,j] + u_true[(i - 1)*M + m,k]), 
                        size = r) 
            }  
          }
        }
      }
      
      y_0 <- y
      v_star_0 <- v_star
    }
    
    for (idx_s_star in seq_along(S_star_grid)) {
      
      S_star <- S_star_grid[idx_s_star]
      
      # cut data 
      {
        y <- y_0[,,seq_len(S + S_star),drop = F]
        v_star <- v_star_0[,seq_len(S_star), drop = F]
      }
      
      edna_dat <- list(n = n,
                       S = S,
                       S_star = S_star,
                       M = M,
                       K = K,
                       y = y,
                       v_star = v_star,
                       sigma_u = sigma_u,
                       mean_lambda = mean_lambda,
                       sigma = rep(sigma,S),
                       r = rep(r, S + S_star),
                       sigma_lambda = sigma_lambda)
      
      
      # fit model
      fit <- sampling(mod, data = edna_dat, chain = 1, iter = 4000, verbose = T)
      
      matrix_of_draws <- as.matrix(fit)
      niter <- nrow(matrix_of_draws)     
      
      colnames(matrix_of_draws)
      
      i <- 4
      j <- 2
      
      l <- (j - 1)*n + i
      colnames(matrix_of_draws)[l]
      qplot(1:niter, matrix_of_draws[,l] - matrix_of_draws[,l+i]) + 
        geom_hline(aes(yintercept = logz_true[i,j] - logz_true[i+1,j]))
      
      i <- 1
      m <- 1
      k <- 3
      l <- (S * n + S * n * M) + (k - 1) * (M * n) + (i - 1) * M + m
      qplot(1:niter, matrix_of_draws[,l]) + 
        geom_hline(aes(yintercept = u_true[(i - 1) * M + m,k] ))
      # qplot(1:niter, matrix_of_draws[,l] - matrix_of_draws[,l+i]) + 
      #   geom_hline(aes(yintercept = u_true[(i - 1) * M + m,k] - logz_true[(i - 1) * M + m + 1,k]))
      # qplot(1:niter, matrix_of_draws[,2] - matrix_of_draws[,3])
      
      biomass_differences_var <- matrix(NA, n - 1, S)
      for (j in 1:S) {
        for (i in 1:(n-1)) {
          biomass_differences_var[i,j] <- 
            var( matrix_of_draws[,(j - 1)*n + i] - matrix_of_draws[,(j - 1)*n + (i+1)])
        }
      }
      
      biomass_differences_mae <- matrix(NA, n - 1, S)
      for (j in 1:S) {
        for (i in 1:(n-1)) {
          biomass_differences_mae[i,j] <- 
            abs( mean(matrix_of_draws[,(j - 1)*n + i] - matrix_of_draws[,(j - 1)*n + (i+1)]) - 
                   (logz_true[i,j] - logz_true[i+1,j]))
        }
      }
      
      
      
      var_ratios[sim, idx_s_star] <- mean(biomass_differences_var)
      mae_ratios[sim, idx_s_star] <- mean(biomass_differences_mae)
      
    }
    
  }
  
  var_ratios <- apply(var_ratios, 1, function(x){
    x / x[length(x)]
  })
  
  var_ratios_mean <- apply(var_ratios, 2, mean)
  
  # data
}

# ----------

# data generating process
{
  sigma <- 1
  tau <- 1
  sigma_u <- 1
  r <- 1000
  sigma_lambda <- 1
  
}

for (sim in 1:nsims) {
  
  # simulate data ---
  
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
            rnbinom(1, mu = exp(lambda_true[j] + v_true[(i - 1)*M + m,j] + u_true[(i - 1)*M + m,k]), 
                    size = r) 
        }  
        for (j in 1:S_star) {
          y[m + (i - 1)*M,k,S + j] <- 
            rnbinom(1, mu = exp(lambda_true[S + j] + v_star[(i - 1)*M + m,j] + u_true[(i - 1)*M + m,k]), 
                    size = r) 
        }  
      }
    }
  }
  
  edna_dat <- list(n = n,
                   S = S,
                   S_star = S_star,
                   M = M,
                   K = K,
                   y = y,
                   v_star = v_star,
                   sigma = sigma,
                   r = r,
                   sigma_u = sigma_u,
                   sigma_lambda = sigma_lambda,
                   tau = tau)
  
  # MCMC ----
  
  fit <- sampling(mod, data = edna_dat, chains = 1, iter = 5000, verbose = T)
  
  matrix_of_draws <- as.matrix(fit)
  
  
  var( matrix_of_draws[,1] - matrix_of_draws[,2])
  var( matrix_of_draws[,3] - matrix_of_draws[,4])
  var( matrix_of_draws[,5] - matrix_of_draws[,6])
  var( matrix_of_draws[,7] - matrix_of_draws[,8])
  var( matrix_of_draws[,9] - matrix_of_draws[,10])
  
  
}

# prior
{
  
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
                 r = r,
                 sigma_u = sigma_u,
                 sigma_lambda = sigma_lambda,
                 tau = tau)

# MCMC ----

# occ_model <- stanc(file = 'model_nb.stan', model_name = "occmodel")
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
