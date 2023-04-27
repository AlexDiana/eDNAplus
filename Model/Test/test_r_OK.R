# simulate data
simulate_y_from_r <- function(r, c_imk, u, v, lambda, mu0, n0, lambdatilde){
  
  y <- array(NA, dim = c(sum(M_site), max(K), S))
  for (i in 1:n) {
    # print(i)
    for (m in 1:M_site[i]) {
      numRep <- K[m + sum(M_site[seq_len(i-1)])]
      for (k in 1:numRep) {
        u_imk <- u[m + sum(M_site[seq_len(i-1)]), k]
        for (j in 1:S) {
          if (c_imk[m + sum(M_site[seq_len(i - 1)]), k, j] == 0) {
            y[m + sum(M_site[seq_len(i - 1)]), k, j] <- 
              # rpois(1, lambda0_true)
              rnbinom(1, mu = mu0, size = n0)
          } else if (c_imk[m + sum(M_site[seq_len(i - 1)]), k, j] == 2) {
            y[m + sum(M_site[seq_len(i - 1)]), k, j] <-
              rpois(1, lambdatilde[j] * exp(lambda[j]))
          } else {
            mean_true <- exp(v[m + sum(M_site[seq_len(i - 1)]), j] + 
                               lambda[j] + 
                               u_imk)
            # rpois(1, mean_true)
            y[m + sum(M_site[seq_len(i - 1)]), k, j] <-
              rnbinom(1, size = r[j], mu = mean_true)
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
  prior_r <- 10
  prior_var <- 1
  
  S <- 2
  
  lambda <- rep(10, S)
  
  logz <- matrix(0, n, S)
  mu <- rep(0, S)
  v <- matrix(0, sum(M_site), S)
  
  u <- matrix(0, sum(M_site), max(K))
  
  mu0 <- 1
  n0 <- 1
  
  lambdatilde <- rep(.5, S)
  
  beta0 <- rep(0, S)
  alpha <- rep(0, S)
  beta_z <- matrix(0, nrow = ncov_z, ncol = S)
  
  zeta <- 0
  zeta_z <- rep(0, S)
  
  c_imk <- array(NA, dim = c(sum(M_site), max(K), S))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      numRep <- K[m + sum(M_site[seq_len(i-1)])]
      for (k in 1:numRep) {
        for (j in 1:S) {
          if(emptySites[i] == 1){
            c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- 0
          } else {
            c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- rbinom(1, 1, .9)
          }
        }
      }
    }
  }
  
  delta <- array(NA, dim = c(sum(M_site), S))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        if(emptySites[i] == 0){
          if(sum(c_imk[m + sum(M_site[seq_len(i-1)]),
                       1:K[m + sum(M_site[seq_len(i-1)])],j]) > 0){
            delta[m + sum(M_site[seq_len(i-1)]),j] <- 1  
          } else {
            delta[m + sum(M_site[seq_len(i-1)]),j] <- 0
          }
        } else {
          delta[m + sum(M_site[seq_len(i-1)]),j] <- 0
        }
      }
    }
  }
  
  gamma <- array(NA, dim = c(sum(M_site), S))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        gamma[m + sum(M_site[seq_len(i-1)]),j] <- 0 
      }
    }
  }
  
  beta_theta <- cbind(rep(0, S), rep(1, S), rep(0, S))
}

# params
{
  r_nb <- rep(prior_r, S)
}

# mcmc -------------

niter <- 1000

# output
{
  r_output <- matrix(NA, S, niter)
}

for (iter in 1:niter) {

  print(iter)
  
  # simulate data
  y <- simulate_y_from_r(r_nb, c_imk, u, v, lambda, mu0, n0, lambdatilde)
      
  # update r
  r_nb <- update_r_nb_cpp(r_nb, prior_r, sqrt(r_var), lambda, X_z[,-1,drop=F], 
                          beta_theta, u, beta_z, beta0, mu, 
                          v, logz, y, delta, gamma, c_imk, M_site, K, 
                          # optimStep = ((iter - 1) %% 5 == 0),
                          optimStep = T,
                          sd_r_proposal = .1)
  
  # output
  r_output[,iter] <- r_nb
  
}

# 

ggplot() + 
  geom_histogram(data = NULL, aes(x = r_output[1,], y = ..density..)) + #xlim(c(prior_mean,3)) + 
  stat_function(fun = dnorm, args = list(mean = prior_r, sd = sqrt(r_var)))

library(ggplot2)
qplot(1:niter, r_output[1,]) + geom_hline(aes(yintercept = prior_r))
