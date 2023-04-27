# simulate data
simulate_data_from_lambda <- function(lambda, 
                                      logz_bar,
                                      beta_theta, 
                                      sigma_beta,
                                      sigma_mu){
  
  theta11 <- matrix(NA, nrow = sum(M_site), ncol = S)
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        theta11[m + sum(M_site[seq_len(i-1)]),j] <- 
          logistic(c(1, logz_bar[i,j] - lambda[j], X_w[m + sum(M_site[seq_len(i-1)]),]) %*% 
                     beta_theta[j,])
      }
    }
  }
  
  delta <- array(NA, dim = c(sum(M_site), S))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        delta[m + sum(M_site[seq_len(i-1)]),j] <- 
          rbinom(1, 1, theta11[m + sum(M_site[seq_len(i-1)]),j])  
        
      }
    }
  }
  
  beta_bar <- sapply(1:S, function(i){
    rnorm(1, lambda[j], sigma_beta)
  })
 
  mu_bar <- sapply(1:S, function(i){
    rnorm(1, lambda[j], sigma_mu)
  })
  
  list("delta" = delta,
       "mu_bar" = mu_bar,
       "beta_bar" = beta_bar)
}

# starting values ----------

# fixed values
{
  n <- 100
  M_site <- rep(5, n)
  K <- rep(1, sum(M_site))
  S <- 2
  
  ncov_z <- 1
  ncov_w <- 1
  
  lambda <- rep(10, S)
  
  u <- matrix(0, sum(M_site), max(K))
  v <- matrix(0, sum(M_site), S)
  
  c_imk <- array(NA, dim = c(sum(M_site), max(K), S))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      numRep <- K[m + sum(M_site[seq_len(i-1)])]
      for (k in 1:numRep) {
        for (j in 1:S) {
          c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- rbinom(1, 1, .9)
        }
      }
    }
  }
  
  r_nb <- rep(10, S)
  
  beta_z <- matrix(0, nrow = ncov_z, ncol = S)
  beta_w <- matrix(0, nrow = ncov_w, ncol = S)
  
  X_z <- matrix(runif(n * ncov_z), n, ncov_z)
  X_w <- matrix(runif(sum(M_site) * ncov_w), sum(M_site), ncov_w)
  
  logz_bar <- matrix(0, n, S)
  beta_theta <- cbind(rep(0, S), rep(1, S), rep(0, S))
}

# params
{
  lambda <- rep(10, S)
  
  lambda_prior <- rep(5, S)
  sigma_lambda <- 2
}

# mcmc -------------

niter <- 5000

# output
{
  lambda_output <- matrix(NA, niter, S)
}

for (iter in 1:niter) {
  
  print(iter)
  
  # simulate data ---------
  
  list_datalambda <- simulate_data_from_lambda(lambda, 
                                               logz_bar,
                                               beta_theta, 
                                               sigma_beta,
                                               sigma_mu)
  mu_bar <- list_datalambda$mu_bar
  beta_bar <- list_datalambda$beta_bar
  delta <- list_datalambda$delta
  
  # update lambda --------
  
  df_t <- 3
  
  for (j in 1:S) {
    
    # nonPCRcounts <- as.vector(y[,,j])[as.vector(c_imk[,,j]) == 2]
    # nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
    # a <- length(nonPCRcounts) * lambdatilde[j]
    
    Xwbetatheta <- beta_theta[j,1] + X_w %*% beta_theta[j,-c(1,2)]
    
    X_l <- rep(logz_bar[,j], M_site)
    
    optim_fun <- function(x){
      -logdpost_cpp(x, X_l, 
                    beta_theta[j,2], Xwbetatheta, delta[,j], 
                    beta_bar[j], sigma_beta, mu_bar[j], sigma_mu, lambda_prior[j], 
                    sigma_lambda)
    }
    
    lambda_star <- optimize(optim_fun, c(lambda[j] - 10, lambda[j] + 10))$minimum
    
    sd_star <-  1 / sqrt(-der2_logdpost_cpp(lambda_star, X_l, 
                                            beta_theta[j,2], Xwbetatheta, 
                                            delta[,j], 
                                            beta_bar[j], sigma_beta, mu_bar[j], sigma_mu, lambda_prior[j], 
                                            sigma_lambda))
    
    lambda_new <- rt2(lambda_star, sd_star, 3)
    
    logproposal_ratio <-  log(dt2(lambda[j], lambda_star, sd_star, df_t)) - 
      log(dt2(lambda_new, lambda_star, sd_star, df_t))
    
    logposterior <- logdpost_cpp(lambda[j], X_l, 
                                 beta_theta[j,2], Xwbetatheta, delta[,j], 
                                 beta_bar[j], sigma_beta, mu_bar[j], sigma_mu, lambda_prior[j], 
                                 sigma_lambda)
    logposterior_star <- logdpost_cpp(lambda_new, X_l, 
                                      beta_theta[j,2], Xwbetatheta, delta[,j], 
                                      beta_bar[j], sigma_beta, mu_bar[j], sigma_mu, lambda_prior[j], 
                                      sigma_lambda)
    
    (mh_ratio <- exp(logposterior_star - logposterior + logproposal_ratio))
    
    if(runif(1) < mh_ratio){
      lambda[j] <- lambda_new
    }
    
  }
  
  
  # output
  lambda_output[iter,] <- lambda
  
}

# 

j <- 1

library(ggplot2)
ggplot() + geom_histogram(data = NULL, aes(x = lambda_output[,j],
                                           y = ..density..)) + 
  stat_function(fun = dnorm, args = list(mean = lambda_prior[j], 
                                          sd = sigma_lambda)) + 
  xlim(c(0,15))

j <- 2
c(sd(lambda_output[,j]),sigma_lambda)
c(mean(lambda_output[,j]),lambda_prior[j])

