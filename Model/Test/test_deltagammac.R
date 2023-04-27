# simulate data
simulate_y_v_from_deltagammac <- function(delta,
                                          gamma,
                                          c_imk,
                                          lambda,
                                          mu_tilde,
                                          n_tilde,
                                          u,
                                          logz,
                                          sigma,
                                          mu,
                                          sigma_gamma,
                                          pi0,
                                          mu0, 
                                          n0, 
                                          r_nb){
  
  v <- matrix(NA, sum(M_site), S)
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1){
          v[m + sum(M_site[seq_len(i-1)]),j] <- 
            rnorm(1, logz[i,j], sigma[j]) 
        } else if (gamma[m + sum(M_site[seq_len(i-1)]),j] == 1){
          v[m + sum(M_site[seq_len(i-1)]),j] <- 
            rnorm(1, mu[j], 
                  sd = sigma_gamma)
        } else {
          v[m + sum(M_site[seq_len(i-1)]),j] <- 
            NA
        }
      }
    }
  }
  
  y <- array(NA, dim = c(sum(M_site), max(K), S))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      numRep <- K[m + sum(M_site[seq_len(i-1)])]
      for (k in 1:numRep) {
        u_imk <- u[m + sum(M_site[seq_len(i-1)]), k]
        for (j in 1:S) {
          if (c_imk[m + sum(M_site[seq_len(i - 1)]), k, j] == 0) {
            y[m + sum(M_site[seq_len(i - 1)]), k, j] <- 
              ifelse(runif(1) < pi0, 0, 1 + rnbinom(1, mu = mu0, size = n0))
          } else if (c_imk[m + sum(M_site[seq_len(i - 1)]), k, j] == 2) {
            y[m + sum(M_site[seq_len(i - 1)]), k, j] <-
              rnbinom(1, mu = mu_tilde, size = n_tilde)
          } else {
            mean_true <- exp(v[m + sum(M_site[seq_len(i - 1)]), j] + 
                               lambda[j] + 
                               u_imk)
            y[m + sum(M_site[seq_len(i - 1)]), k, j] <-
              rnbinom(1, size = r_nb[j], mu = mean_true)
          }
        }
      }
    }
  }
  
  list("y" = y,
       "v" = v)
}

# starting values ----------

# data params
{
  S <- 5
  n <- 5
  M_site <- rep(2, n)
  K <- rep(2, sum(M_site))
  
  ncov_w <- 0
  
  emptyTubes <- 0
  X_w <- matrix(runif(sum(M_site) * ncov_w), sum(M_site), ncov_w)
}

# fixed values
{
  lambda <- rep(10, S)
  
  mu0 <- 100
  n0 <- 1
  pi0 <- .5
  
  mu_tilde <- 100
  n_tilde <- 100
  
  logz <- matrix(0, n, S)
  mu <- rep(0, S)
  u <- matrix(0, sum(M_site), max(K))
  
  beta0 <- rep(0, S)
  beta_z <- matrix(0, nrow = ncov_z, ncol = S)
  beta_w <- matrix(0, nrow = ncov_w, ncol = S)
  
  v <- matrix(0, sum(M_site), S)
  
  sigma <- rep(.5, S)
  
  theta11 <- matrix(runif(sum(M_site) * S), 
                    nrow = sum(M_site), ncol = S)
  theta10 <- runif(S, 0, .05)

  p_11 <- rbeta(S, a_p11, b_p11)
  p_10 <- rep(0.02, S)
  
  r_nb <- rep(10, S)
}

# starting values
{
  delta <- matrix(NA, 
                  nrow = sum(M_site), ncol = S)
  gamma <- matrix(NA, 
                  nrow = sum(M_site), ncol = S)
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        delta[m + sum(M_site[seq_len(i-1)]),j] <- rbinom(1, 1, theta11[m + sum(M_site[seq_len(i-1)]),j])
        if(delta[m + sum(M_site[seq_len(i-1)]),j] == 0){
          gamma[m + sum(M_site[seq_len(i-1)]),j] <- rbinom(1, 1, theta10[j])
        }
      }
    }
  }
  
  c_imk <- array(0, dim = c(sum(M_site), max(K), S))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      numRep <- K[m + sum(M_site[seq_len(i-1)])]
      for (k in 1:numRep) {
        for (j in 1:S) {
          if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1 | 
             gamma[m + sum(M_site[seq_len(i-1)]),j] == 1){
            c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- rbinom(1, 1, p_11[j])
          } else {
            c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- 2 * rbinom(1, 1, p_10[j])
          }
        }
      }
    }
  }
  
}

# mcmc -------------

niter <- 10000

# output
{
  delta_output <- array(NA, dim = c(sum(M_site), S, niter))
  gamma_output <- array(NA, dim = c(sum(M_site), S, niter))
  cimk_output <- array(NA, dim = c(sum(M_site), max(K), S, niter))
}

for (iter in 1:niter) {
  
  print(iter)
  
  # simulate data ---------
  
  list_yv <- simulate_y_v_from_deltagammac(delta,
                                           gamma,
                                           c_imk,
                                           lambda,
                                           mu_tilde,
                                           n_tilde,
                                           u,
                                           logz,
                                           sigma,
                                           mu,
                                           sigma_gamma,
                                           pi0,
                                           mu0, 
                                           n0, 
                                           r_nb)
  y <- list_yv$ys
  v <- list_yv$v
  
  # update delta gamma c --------
  
  v_pres <- (delta == 1) | (gamma == 1)
  list_deltagammac <- update_delta_c_d_rjmcmc_old(v_pres, 
                                                  y, 
                                                  v, 
                                                  lambda, 
                                                  r_nb,
                                                  M_site, K, 
                                                  mu0, n0, pi0, 
                                                  mu_tilde,
                                                  n_tilde, u, logz, 
                                                  X_w, beta_w,
                                                  sigma, mu, 
                                                  sigma_gamma, v_sd = .5,
                                                  p_11, p_10, theta11, 
                                                  theta10, emptyTubes,
                                                  S_star)
  delta <- list_deltagammac$delta
  gamma <- list_deltagammac$gamma
  c_imk <- list_deltagammac$c_imk
  v <- list_deltagammac$v
  
  # output
  delta_output[,,iter] <- delta
  gamma_output[,,iter] <- gamma
  cimk_output[,,,iter] <- c_imk
  
}

# 
delta_mean <- apply(delta_output, c(1,2), mean)
c_imk_mean <- apply(cimk_output, c(1,2,3), mean)



i <- 4
m <- 1
k <- 1
j <- 1

delta_mean[m + sum(M_site[seq_len(i-1)]), j]
theta11[m + sum(M_site[seq_len(i-1)]), j]

theta11[m + sum(M_site[seq_len(i-1)]), j] * p_11[j]
c_imk_mean[m + sum(M_site[seq_len(i-1)]), k, j]
qplot(1:niter, cumsum(delta_output[m + sum(M_site[seq_len(i-1)]), j,]) / 1:niter)
