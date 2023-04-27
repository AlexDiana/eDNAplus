mean_v <- rep(0, sum(M_site))
var_v <- rep(0, sum(M_site))
vbar_im <- rep(0, sum(M_site))

Xw_beta <- X_w %*% beta_w

for (i in 1:n) {
  for (m in 1:M_site[i]) {
    numSamples = 0
    for (j in 1:S) {
      if(delta[m + sum(M_site[seq_len(i-1)]), j] == 1){
        vbar_im[m + sum(M_site[seq_len(i-1)])] <- vbar_im[m + sum(M_site[seq_len(i-1)])] + 
          v[m + sum(M_site[seq_len(i-1)]), j];
        mean_v[m + sum(M_site[seq_len(i-1)])] <- mean_v[m + sum(M_site[seq_len(i-1)])] + 
          logz[i, j] + Xw_beta[m + sum(M_site[seq_len(i-1)]), j]
        var_v[m + sum(M_site[seq_len(i-1)])] <- var_v[m + sum(M_site[seq_len(i-1)])]  + 
          sigma[j] * sigma[j]
        numSamples <- numSamples + 1
      } else if(gamma[m + sum(M_site[seq_len(i-1)]), j] == 1){
        vbar_im[m + sum(M_site[seq_len(i-1)])] <- vbar_im[m + sum(M_site[seq_len(i-1)])] + 
          v[m + sum(M_site[seq_len(i-1)]), j];
        mean_v[m + sum(M_site[seq_len(i-1)])] <- mean_v[m + sum(M_site[seq_len(i-1)])] + 
         mu[j]
        var_v[m + sum(M_site[seq_len(i-1)])] <- var_v[m + sum(M_site[seq_len(i-1)])]  + 
          sigma_gamma * sigma_gamma
        numSamples <- numSamples + 1
      }
    }
    vbar_im[m + sum(M_site[seq_len(i-1)])] <- vbar_im[m + sum(M_site[seq_len(i-1)])] / numSamples
    mean_v[m + sum(M_site[seq_len(i-1)])] <- mean_v[m + sum(M_site[seq_len(i-1)])] / numSamples
    var_v[m + sum(M_site[seq_len(i-1)])] <- var_v[m + sum(M_site[seq_len(i-1)])] / (numSamples * numSamples)
  }
}

ubar_im <- rep(NA, sum(M_site))
for (i in 1:n) {
  for (m in 1:M_site[i]) {
    ubar_im[m + sum(M_site[seq_len(i-1)])] = mean(u[m + sum(M_site[seq_len(i-1)]),])
  }
}

vtilde <- matrix(NA, sum(M_site), S)
for (i in 1:n) {
  for (m in 1:M_site[i]) {
    for(j in 1:S){
      if(delta[m + sum(M_site[seq_len(i-1)]), j] == 1 | 
         gamma[m + sum(M_site[seq_len(i-1)]), j] == 1){
        vtilde[m + sum(M_site[seq_len(i-1)]), j] = v[m + sum(M_site[seq_len(i-1)]), j] - 
          vbar_im[m + sum(M_site[seq_len(i-1)])]
      }
    }
  }
}

utilde <- matrix(NA, sum(M_site), max(K))
for (i in 1:n) {
  for (m in 1:M_site[i]) {
    for(k in 1:K[m + sum(M_site[seq_len(i-1)])]){
      utilde[m + sum(M_site[seq_len(i-1)]), k] <-  u[m + sum(M_site[seq_len(i-1)]), k] - 
        ubar_im[m + sum(M_site[seq_len(i-1)])]
    }
  }
}

for (i in 1:n) {
  for (m in 1:M_site[i]) {
    
    uv_current <- c(ubar_im[m + sum(M_site[seq_len(i-1)])], 
                    vbar_im[m + sum(M_site[seq_len(i-1)])])
    
    m_v <-  mean_v[m + sum(M_site[seq_len(i-1)])]
    v_v <- var_v[m + sum(M_site[seq_len(i-1)])]
    
    m_u <- 0
    v_u <- sigma_u^2 / K[m + sum(M_site[seq_len(i-1)])]
    
    a_1 <- 0
    a_2 <- 0
    a_3 <- 0
    a_4 <- 0
    
    l2 <- 0
    r_all <- rep(NA, S * K[m + sum(M_site[seq_len(i-1)])])
    lambdas <- rep(NA, S * K[m + sum(M_site[seq_len(i-1)])])
    x_present <- rep(NA, S * K[m + sum(M_site[seq_len(i-1)])])
    l3 <- 0
    r_all2 <- rep(NA, S * K[m + sum(M_site[seq_len(i-1)])])
    lambdas2 <- rep(NA, S * K[m + sum(M_site[seq_len(i-1)])])
    x_present2 <- rep(NA, S * K[m + sum(M_site[seq_len(i-1)])])
    for (k in 1:K[m + sum(M_site[seq_len(i-1)])]) {
      for (j in 1:S) {
        if(c_imk[m + sum(M_site[seq_len(i-1)]), k, j] == 1){
          l2 <- l2 + 1
          r_all[l2] = r_nb[j]
          x_present[l2] = lambda[j] + vtilde[m + sum(M_site[seq_len(i-1)]), j] + 
            utilde[m + sum(M_site[seq_len(i-1)]), k]
          lambdas[l2] = lambda_ijk[m + sum(M_site[seq_len(i-1)]), k, j]
          a_1 <- a_1 + ( - r_nb[j] * exp(- x_present[l2]) * lambdas[l2])
          a_2 <- a_2 + r_nb[j]
        }
      }
      
      for (j in S + seq_len(S_star)) {
        if(c_imk[m + sum(M_site[seq_len(i-1)]), k, j] == 1){
          l3 <- l3 + 1
          r_all2[l3] = r_nb[j]
          x_present2[l3] = lambda[j] + v[m + sum(M_site[seq_len(i-1)]), j] + 
            utilde[m + sum(M_site[seq_len(i-1)]), k]
          lambdas2[l3] = lambda_ijk[m + sum(M_site[seq_len(i-1)]), k, j]
          a_3 <- a_3 + ( - r_nb[j] * exp(- lambda[j] - utilde[m + sum(M_site[seq_len(i-1)]), k] - 
                                           v[m + sum(M_site[seq_len(i-1)]), j]) * 
                           lambdas2[l3])
          a_4 <- a_4 + r_nb[j]
        }
      }
    }
    r_all <- r_all[1:l2]
    lambdas <- lambdas[1:l2]
    x_present <- x_present[1:l2]
    r_all2 <- r_all2[1:l3]
    lambdas2 <- lambdas2[1:l3]
    x_present2 <- x_present2[1:l3]
    
    ubar_star = -log(-a_4/a_3)
    vbar_star = log(a_1*a_4/(a_2*a_3))
    mu_lik <- c(ubar_star, vbar_star)
    
    Sigma_lik <- solve(- H_fjoint_noprior(ubar_star, vbar_star, a_1, a_3))
    
    # sum(uv_current)
    # 
    # ubar_star + vbar_star
    
    mu_prior <- c(mean_v[m + sum(M_site[seq_len(i-1)])], 0)
    D_var <- diag(c(var_v[m + sum(M_site[seq_len(i-1)])], (sigma_u / K)^2) , 2)
    
    Sigma_star <- solve(solve(Sigma) + solve(D_var))
          
    mu <- solve(Sigma_star) %*% mu_star + solve(D_var) %*% mu_prior
    mu_star <- Sigma_star %*% mu
    
    mu_new <- mvrnorm(1, mu_star, Sigma_star)
    
    loglik_new <- logpost_gamma_uv(mu_new[1], mu_new[2], 
                                   lambdas, x_present, r_all, 
                                   lambdas2, x_present2, r_all2)
    loglik_current <- logpost_gamma_uv(uv_current[1], uv_current[2], 
                                       lambdas, x_present, r_all, 
                                       lambdas2, x_present2, r_all2)

    loglik <- loglik_new - loglik_current
    
    logprior_new <- dnorm(mu_new[1], m_u, sqrt(v_u), log = T) +
      dnorm(mu_new[1], m_v, sqrt(v_v), log = T)
    logprior_current <- dnorm(uv_current[1], m_u, sqrt(v_u), log = T) +
      dnorm(uv_current[1], m_v, sqrt(v_v), log = T)
    
    logprior <- logprior_new - logprior_current
    
    logproposal_new <- dmvnorm_cpp(mu_new, mu_star, Sigma_star, 1)   
    logproposal_current <- dmvnorm_cpp(uv_current, mu_star, Sigma_star, 1)   
    
    logproposal <- logproposal_current - logproposal_new
      
    if(runif(1) < exp(loglik + logprior + logproposal)){
      c(ubar_im[m + sum(M_site[seq_len(i-1)])], 
        vbar_im[m + sum(M_site[seq_len(i-1)])])
    }
    
  }
}


Sigma_all <- matrix(0, n * 2, n * 2)
for (i in 1:n) {
  Sigma_all[(i-1)*2 + 1:2,
            (i-1)*2 + 1:2] <- matrix(c(1,-.5,-.5,1), 2, 2)
}

A_mat <- rep(c(0,1), times = n)

Sigma_all %*% A_mat

invMat <- solve(A_mat %*% Sigma_all %*% A_mat)

uv <- rep(NA, 2*n)
for (i in 1:n) {
  uv[(i-1)*2 + 1:2] <- c(vbar_im[i], ubar_im[i])
}

Axme <- A_mat %*% uv

uv_star <- uv - (Sigma_all %*% A_mat) %*% invMat %*% Axme

A_mat %*% uv_star

uv_star[1] + uv_star[2]
uv[1] + uv[2]

uv[1:2]
uv_star[1:2]
