# simulate data
simulate_delta_from_betatheta <- 
  # function(beta_theta, logz, X_w){
  function(theta11){
  
  # theta11 <- matrix(NA, nrow = sum(M_site), ncol = S)
  # for (i in 1:n) {
  #   for (m in 1:M_site[i]) {
  #     for (j in 1:S) {
  #       theta11[m + sum(M_site[seq_len(i - 1)]), j] <-
  #         logistic(c(1, logz[i, j], X_w[m + sum(M_site[seq_len(i - 1)]), ]) %*%
  #                    beta_theta[j, ])
  #     }
  #   }
  # }
  
  delta <- array(NA, dim = c(sum(M_site), S))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        delta[m + sum(M_site[seq_len(i - 1)]), j] <-
          rbinom(1, 1, theta11[m + sum(M_site[seq_len(i - 1)]), j])
        
      }
    }
  }
  
  delta
}

# params ----------

# fixed values
{
  n <- 250
  M_site <- rep(2, n)
  
  S <- 2
  
  ncov_w <- 1
  X_w <- matrix(runif(sum(M_site) * ncov_w), ncol = ncov_w)
  
  logz <- matrix(runif(n * S), n, S)
  
}

# starting values
{
  beta_theta <- cbind(rep(0, S), pmax(rnorm(S, 1), 0), matrix(0, S, ncov_w))
  theta11 <- matrix(NA, nrow = sum(M_site), ncol = S)
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        theta11[m + sum(M_site[seq_len(i-1)]),j] <- 
          logistic(c(1, logz[i,j], X_w[m + sum(M_site[seq_len(i-1)]),]) %*% 
                     beta_theta[j,])
      }
    }
  }
}

# prior
{
  b_theta11 <- c(1,rep(0, ncov_w))
  B_theta11 <- diag(sigma_beta_theta, nrow = 1 + ncov_w)
}

# mcmc -------------

niter <- 10000

# output
{
  beta_theta_output <- array(NA, dim = c(S, 2 + ncov_w, niter))
}

for (iter in 1:niter) {
  
  print(iter)
  
  # simulate data
  delta <- simulate_delta_from_betatheta(
    # beta_theta, logz, X_w
    theta11
    )
  
  # update l
  list_beta_theta <- update_betatheta11_cpp(logz,
                                            beta_theta,
                                            theta11,
                                            delta[1:sum(M_site),],
                                            X_w, M_site,
                                            b_theta11,
                                            B_theta11,
                                            F)
  beta_theta <- list_beta_theta$beta_theta
  theta11 <- list_beta_theta$theta11
  
  # output
  beta_theta_output[,,iter] <- beta_theta
  
}


library(ggplot2)
# qplot(1:niter, beta_theta_output[i,j,]) +
  # geom_hline(aes(yintercept = 0), color = "red")

j <- 1
l <- 2
ggplot() + 
  geom_histogram(data = NULL, aes(x = beta_theta_output[j,l,], y = ..density..)) + #xlim(c(prior_mean,3)) + 
  stat_function(fun = dnorm, args = list(mean = b_theta11[l - 1], sd = sqrt(B_theta11[l-1,l-1])))

c(mean(beta_theta_output[j, l, ]), b_theta11[l - 1])
c(var(beta_theta_output[j, l, ]), B_theta11[l - 1, l - 1])

