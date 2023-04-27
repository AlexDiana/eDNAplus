list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, 
                                gamma, beta_theta, M_site)
beta_bar <- list_CP_cpp$beta_bar
mu_bar <- list_CP_cpp$mu_bar
logz_bar <- list_CP_cpp$logz_bar
v_bar <- list_CP_cpp$v_bar

Xb <- X_z %*% beta_z

for (i in 1:n) {
  print(i)
  Xb_i <-  Xb[i,,drop=F]
  
  logz_star_mean <- rep(0, S)
  logz_star_sd <- rep(0, S)
  
  logz_bar_current <- logz_bar[i,]
  
  # invSigma <- matrix(0, S, S)
  
  for (j in 1:S) {
    
    logzbar_current_j <- logz_bar[i,j]
    
    idxPresent <- which(delta[1:M_site[i] + sum(M_site[seq_len(i-1)]),j] == 1)
    
    v_samples <- v_bar[idxPresent + sum(M_site[seq_len(i-1)]),j] - 
      r[idxPresent + sum(M_site[seq_len(i-1)])] * alpha[j]
    
    y_all <- delta[1:M_site[i] + sum(M_site[seq_len(i-1)]),j]
    v_all <- beta_theta[j, 1] +
      r[1:M_site[i] + sum(M_site[seq_len(i-1)])] * beta_theta[j, 3]
    x_all <- beta_theta[j, 2] / exp(lambda[j])
    
    prior_mean <- beta_bar[j] + Xb_i[j]
    prior_sd <- sqrt(Tau[j,j])
    
    logz_star_mean[j] <- findzero_cpp(logzbar_current_j - 15,
                                      logzbar_current_j + 15,
                                      .01,
                                      x_all,
                                      y_all, v_all, 
                                      v_samples,
                                      prior_sd,
                                      sigma[j], prior_mean)
    
    logz_star_sd[j] <- sqrt(1 / (-h_f_cpp(logz_star_mean[j],
                                                 x_all,
                                                 y_all, v_all, 
                                                 v_samples,
                                          prior_var,
                                          sigma[j], prior_mean)))
    
    l_grid <- seq(2.5, 7.5, length.out = 100)
    logf_values <- sapply(l_grid, function(l){
      # logposterior_logz_cpp(v_samples, sigma[j],
      #                       0, l,
      #                       prior_mean, prior_sd) +
      #   logposterior_logz_logistic_cpp(y_all,
      #                                  x_all,
      #                                  v_all,
      #                                  0,
      #                                  l)
      logf_cpp_correct(l, x_all, sigma[j], v_samples, v_all, y_all, prior_sd, prior_mean)
    })
     
    qplot(l_grid, logf_values)
    # qplot(l_grid, exp(logf_values - max(logf_values)))
    
    
    # invSigma[j,j] <- sqrt(1 /  logz_star_sd[j]^2)
    
  }
  invSigma <- diag(1 / logz_star_sd^2)
  print(invSigma)
  # post_cov <- solve(invTau + invSigma)
  # post_mu <- (invSigma %*% logz_star_mean) + invTau %*% t(Xb_i)
  # 
  # post_mean <- post_cov %*% post_mu
  # 
  # logz_star <- mvrnorm(1, post_mean, 2 * post_cov)
  #  
  # logposterior_ratios <- 0
  # 
  # for (j in 1:S) {
  #   
  #   logzbar_current_j <- logz_bar[i,j]
  #   
  #   idxPresent <- which(delta[1:M_site[i] + sum(M_site[seq_len(i-1)]),j] == 1)
  #   
  #   v_samples <- v_bar[idxPresent + sum(M_site[seq_len(i-1)]),j] - 
  #     r[idxPresent + sum(M_site[seq_len(i-1)])] * alpha[j]
  #   
  #   y_all <- delta[1:M_site[i] + sum(M_site[seq_len(i-1)]),j]
  #   v_all <- beta_theta[j, 1] +
  #     r[1:M_site[i] + sum(M_site[seq_len(i-1)])] * beta_theta[j, 3]
  #   x_all <- beta_theta[j, 2] / exp(lambda[j])
  #   
  #   logpost_diff <- 
  #     logposterior_logz_both_cpp(logz_star[j], v_samples, sigma[j],
  #                                y_all, x_all,
  #                                v_all) -
  #     logposterior_logz_both_cpp(logzbar_current_j, v_samples, sigma[j],
  #                                y_all, x_all,
  #                                v_all)
  #   
  #   logposterior_ratios <- logposterior_ratios + logpost_diff
  #   # print(logpost_diff)
  # }
  # 
  # logprior_ratio <- 
  #   dmvnorm_cpp(logz_star, Xb_i, Tau, 1) - 
  #   dmvnorm_cpp(logz_bar_current, Xb_i, Tau, 1);
  # 
  # logproposal_ratio <- 
  #   dmvnorm_cpp(logz_bar_current, post_mean, 2 * post_cov, 1) - 
  #   dmvnorm_cpp(logz_star, post_mean, 2 * post_cov, 1)
  # 
  # mh_ratio <- exp(logposterior_ratios +  logprior_ratio + logproposal_ratio)
  
}
