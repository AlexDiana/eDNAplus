# simulate data
simulate_data_from_lambda <- function(lambda, 
                                      v,
                                      u, 
                                      c_imk,
                                      r_nb){
  
  y <- array(NA, dim = c(sum(M_site), max(K), S))
  for (i in 1:n) {
    # print(i)
    for (m in 1:M_site[i]) {
      numRep <- K[m + sum(M_site[seq_len(i-1)])]
      for (k in 1:numRep) {
        u_imk <- u[m + sum(M_site[seq_len(i-1)]), k]
        for (j in 1:S) {
           if (c_imk[m + sum(M_site[seq_len(i - 1)]), k, j] == 1) {
             mean_true <- exp(v[m + sum(M_site[seq_len(i-1)]), j] + 
                                lambda[j] + 
                                u_imk)
             
             y[m + sum(M_site[seq_len(i-1)]), k, j] <-
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
  n <- 100
  S <- 2
  
  M_site <- M_site[1:n]
  K <- rep(1, length(M_site))
  
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
}

# params
{
  lambda <- rep(10, S)
  
  lambda_prior <- rep(5, S)
  sigma_lambda <- 2
}

# mcmc -------------

niter <- 1000

# output
{
  lambda_output <- matrix(NA, niter, S)
}

for (iter in 1:niter) {
  
  print(iter)
  
  # simulate data ---------
  
  y <- simulate_data_from_lambda(lambda, 
                                 v,
                                 u, 
                                 c_imk,
                                 r_nb)
  
  # lambda ijk
  
  lambda_ijk <- update_lambdaijk(lambda, lambda_ijk, v, u, r_nb, c_imk, M_site, y, K,
                                 S_star)
  
  # update lambda --------
  
  lambda <- update_lambda_NP(lambda_ijk, c_imk, mu,
                             r_nb, v, u,
                             lambda_prior, sigma_lambda)
  
  
  # output
  lambda_output[iter,] <- lambda
  
}

# 

a_lambda_prior <- lambda_prior[j]^2 / sigma_lambda^2
b_lambda_prior <- lambda_prior[j] / sigma_lambda^2

library(ggplot2)
ggplot() + geom_histogram(data = NULL, aes(x = lambda_output[,j],
                                           y = ..density..)) + 
  stat_function(fun = dgamma, args = list(shape = a_lambda_prior, 
                                         rate = b_lambda_prior)) + 
  xlim(c(0,20))
