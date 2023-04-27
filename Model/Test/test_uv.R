# simulate data
simulate_y_from_uv <- function(v, 
                              u, 
                              lambda, 
                              c_imk, 
                              r_nb,
                              M_site,
                              K,
                              emptyTubes,
                              S, S_star){
  
  y <- array(NA, dim = c(sum(M_site) + emptyTubes, max(K), S + S_star))
  for (i in 1:n) {
    # print(i)
    for (m in 1:M_site[i]) {
      numRep <- K[m + sum(M_site[seq_len(i-1)])]
      for (k in 1:numRep) {
        u_imk <- u[m + sum(M_site[seq_len(i-1)]), k]
        for (j in 1:(S + S_star)) {
          if (c_imk[m + sum(M_site[seq_len(i - 1)]), k, j] == 0) {
            # y[m + sum(M_site[seq_len(i - 1)]), k, j] <- 
            #   ifelse(runif(1) < pi0_true, 0, 1 + rnbinom(1, mu = mu0_true, size = n0_true))
          } else if (c_imk[m + sum(M_site[seq_len(i - 1)]), k, j] == 2) {
            # y[m + sum(M_site[seq_len(i - 1)]), k, j] <-
            #   # runif(1, 0, exp(lambda_true[j]))
            #   rnbinom(1, mu = mu_tilde_true, size = n_tilde_true)
            # rpois(1, lambdatilde_true[j] * exp(lambda_true[j]))
          } else {
            mean_true <- exp(v[m + sum(M_site[seq_len(i - 1)]), j] + 
                               lambda[j] + 
                               u_imk)
            # rpois(1, mean_true)
            y[m + sum(M_site[seq_len(i - 1)]), k, j] <-
              rnbinom(1, size = r_nb[j], mu = mean_true)
          }
        }
      }
    }
  }
  
  y
}

# starting values ----------

# fixed values
{
  n <- 2
  M_site <- rep(1, n)
  K <- rep(5, sum(M_site))
  
  emptyTubes <- 0
  
  S <- 5
  
  ncov_z <- 1
  ncov_w <- 1
  
  lambda <- rep(10, S)
  
  logz <- matrix(0, n, S)
  mu <- rep(0, S)
  
  beta0 <- rep(0, S)
  beta_z <- matrix(0, nrow = ncov_z, ncol = S)
  beta_w <- matrix(0, nrow = ncov_w, ncol = S)
  
  X_z <- matrix(runif(n * ncov_z), n, ncov_z)
  X_w <- matrix(runif(sum(M_site) * ncov_w), sum(M_site), ncov_w)
  
  sigma <- pmin(rhcauchy(S, .5), 1)
  
  delta <- array(NA, dim = c(sum(M_site) + emptyTubes, S + S_star))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        delta[m + sum(M_site[seq_len(i-1)]),j] <- 1#rbinom(1, 1, .8)
      }
      for (j in seq_len(S_star)) {
        if(PCR_spiked[m + sum(M_site[seq_len(i-1)])]){
          delta[m + sum(M_site[seq_len(i-1)]),S + j] <- 1  
        } else {
          delta[m + sum(M_site[seq_len(i-1)]),S + j] <- 0  
        }
      }
    }
  }
  # empty tubes
  for (m in seq_len(emptyTubes)) {
    delta[sum(M_site) + m,] <- 0
  }
  
  gamma <- array(NA, dim = c(sum(M_site) + emptyTubes, S + S_star))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        if(delta[m + sum(M_site[seq_len(i-1)]),j] == 0){
          gamma[m + sum(M_site[seq_len(i-1)]),j] <- rbinom(1, 1, .2)  
        } else if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1){
          gamma[m + sum(M_site[seq_len(i-1)]),j] <- 0
        }
      }
      for (j in seq_len(S_star)) {
        gamma[m + sum(M_site[seq_len(i-1)]),S + j] <- 0
      }
    }
  }
  # empty tubes
  for (m in seq_len(emptyTubes)) {
    gamma[sum(M_site) + m,] <- 0
  }
  
  c_imk <- array(NA, dim = c(sum(M_site) + emptyTubes, max(K), S + S_star))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      numRep <- K[m + sum(M_site[seq_len(i-1)])]
      for (k in 1:numRep) {
        for (j in 1:(S+S_star)) {
          if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1 | 
             gamma[m + sum(M_site[seq_len(i-1)]),j] == 1){
            c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- 1#rbinom(1, 1, .8)
          } else {
            c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- 2 * rbinom(1, 1, .1)
          }
        }
      }
    }
  }
  
  beta_theta <- cbind(rep(0, S), rep(1, S), rep(0, S))
  
  sigma_u <- 1
  
  r_nb <- rep(10, S)
}

# params
{
  v <- matrix(NA, sum(M_site), S)
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1 | 
           gamma[m + sum(M_site[seq_len(i-1)]),j] == 1){
          v[m + sum(M_site[seq_len(i-1)]),j]  <-  0
        }
      }
    }
  }
  
  
  u <- matrix(0, sum(M_site), max(K))
  
}

# mcmc -------------

niter <- 20000

# output
{
  v_output <- array(NA, c(sum(M_site), S, niter))
  u_output <- array(NA, c(sum(M_site), max(K), niter))
}

for (iter in 1:niter) {
  
  print(iter)
  
  # simulate data ---------
  
  y <- simulate_y_from_uv(v, 
                          u, 
                          lambda, 
                          c_imk, 
                          r_nb,
                          M_site,
                          K,
                          emptyTubes,
                          S, S_star)
  
  # update lambda_ijk --------
  
  lambda_ijk <- update_lambdaijk(lambda, lambda_ijk, v, u, r_nb, 
                                 c_imk, M_site, y, K,
                                 S_star)
  
  # update v ---------------
  
  list_uv <- update_uv_poisgamma_cpp(u, v, logz, lambda, X_z, beta_theta, beta_z, beta0,
                                     r_nb, mu, lambda_ijk, c_imk, delta, gamma, sigma, sigma_gamma,
                                     sigma_u, M_site, X_w, beta_w, K, S_star)
  u <- list_uv$u
  v <- list_uv$v
  lambda <- list_uv$lambda
  
  v <- update_v_poisgamma_cpp(v, logz,
                              lambda,  X_z,
                              beta_theta, u, beta_z,
                              beta0, r_nb, mu, lambda_ijk,
                              c_imk, delta, gamma, sigma,
                              sigma_gamma, M_site,
                              X_w, beta_w,
                              K, S_star)

  list_u <- update_u_poisgamma_cpp(v, u, lambda, beta0, beta_z, logz,
                                   mu, lambda_ijk, r_nb, X_w, beta_w, c_imk,
                                   delta, gamma, sigma_u, beta_theta, sigma,
                                   sigma_gamma, M_site,
                                   K, S_star)
  u <- list_u$u
  lambda <- list_u$lambda
  
  # output
  v_output[,,iter] <- v
  u_output[,,iter] <- u
  
}

# v
library(ggplot2)

l <- 1
j <- 1

ggplot() + geom_histogram(data = NULL, aes(x = v_output[l,j,],
                                           y = ..density..)) + 
  stat_function(fun = dnorm, args = list(mean = 0, 
                                         sd = sigma[j])) #+ 
  # xlim(c(-3,3))

c(mean(v_output[l,j,]), 0)
c(var(v_output[l,j,]), sigma[j]^2)


# v mean

v_mean_output <- apply(v_output, c(1,3), mean)

l <- 1

current_var <- sum(sigma^2) / S^2

ggplot() + geom_histogram(data = NULL, aes(x = v_mean_output[l,],
                                           y = ..density..)) + 
  stat_function(fun = dnorm, args = list(mean = 0, 
                                         sd = sqrt(current_var))) 


c(var(v_mean_output[l,]), current_var)
c(mean(v_mean_output[l,]), 0)

# u
l <- 1
k <- 1

ggplot() + geom_histogram(data = NULL, aes(x = u_output[l,k,],
                                           y = ..density..)) + 
  stat_function(fun = dnorm, args = list(mean = 0, 
                                         sd = sigma_u)) 
  # xlim(c(-3,3))

c(var(u_output[l,k,]), sigma_u * sigma_u)

# qplot(1:niter, u_output[l,k,])

# u mean

u_mean_output <- apply(u_output, c(1,3), mean)

l <- 1

current_var <- sigma_u * sigma_u / K[l]

ggplot() + geom_histogram(data = NULL, aes(x = u_mean_output[l,],
                                           y = ..density..)) + 
  stat_function(fun = dnorm, args = list(mean = 0, 
                                         sd = sqrt(current_var))) 


c(var(u_mean_output[l,]), current_var)
c(mean(u_mean_output[l,]), 0)

#

qplot(1:niter, v_output[1,1,])
qplot(1:niter, u_output[1,1,])
