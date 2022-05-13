# UTILITY --------

logistic <- function(x){
  1 / (1 + exp(-x))
}

# LAMBDA ----------

update_lambda_CP <- function(beta0, beta_z, logz, 
                             mu, lambda, v, u, lambda_ijk, r_nb,
                             c_imk, delta, gamma, X_w, beta_theta, 
                             M_site, sigma_beta, sigma_mu,
                             lambda_prior, sigma_lambda,
                             S_star, emptyTubes){
  
  df_t  = 3
  
  S <- length(mu)
  
  list_CP_cpp <- convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, 
                                   delta, gamma, beta_theta, M_site, S_star,
                                   emptyTubes)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  beta_theta_bar <- list_CP_cpp$beta_theta_bar
  
  # update paramters
  
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
    
    # lambda_star <- optimize(optim_fun, c(lambda[j] - 10, lambda[j] + 10), tol = .001)$minimum
    # try two mins
    # optim1 <- 
    # lambda_star_1 <- optimize(optim_fun, c(lambda[j] - 1, lambda[j] + 1), tol = .001)$minimum
    # lambda_obj_1 <- optimize(optim_fun, c(lambda[j] - 1, lambda[j] + 1), tol = .001)$objective
    # lambda_star_2 <- optimize(optim_fun, c(lambda[j] - 4, lambda[j] + 4), tol = .01)$minimum
    # lambda_obj_2 <- optimize(optim_fun, c(lambda[j] - 4, lambda[j] + 4), tol = .01)$objective
    # lambda_star_3 <- optimize(optim_fun, c(lambda[j] - 10, lambda[j] + 10), tol = .01)$minimum
    # lambda_obj_3 <- optimize(optim_fun, c(lambda[j] - 10, lambda[j] + 10), tol = .01)$objective
    # if(lambda_obj_1 < lambda_obj_2 & lambda_obj_1 < lambda_obj_3){
    #   lambda_star <- lambda_star_1
    # } else if(lambda_obj_2 < lambda_obj_1 & lambda_obj_2 < lambda_obj_3){
    #   lambda_star <- lambda_star_2
    # } else {
    #   lambda_star <- lambda_star_3
    # }
    lambda_star <- optimize(optim_fun, c(lambda[j] - 10, lambda[j] + 10))$minimum
    
    sd_star <-  1 / sqrt(-der2_logdpost_cpp(lambda_star, X_l, 
                                            beta_theta[j,2], Xwbetatheta, 
                                            delta[,j], 
                                            beta_bar[j], sigma_beta, mu_bar[j], sigma_mu, lambda_prior[j], 
                                            sigma_lambda))
    
    # lambda_grid <- seq(lambda[j] - 1, lambda[j] + 1, length.out = 100)
    # # lambda_grid <- seq(lambda_star - 01, lambda_star + 01, length.out = 100)
    # y_grid <- sapply(lambda_grid, function(x){
    #   log(dt2(x, lambda_star, sd_star, df_t))
    #   # optim_fun(x)
    #   logdpost_cpp(x, X_l,
    #                beta_theta[j,2], Xwbetatheta, delta[,j],
    #                nonPCRcounts, #lambdatilde[j],
    #                beta_bar[j], sigma_beta, mu_bar[j], sigma_mu, lambda_prior[j],
    #                sigma_lambda)
    # })
    # 
    # # qplot(lambda_grid, y_grid) +
    # qplot(lambda_grid, exp(y_grid - max(y_grid)) / sum(exp(y_grid - max(y_grid)))) +
    #   geom_vline(aes(xintercept = lambda_true[j])) +
    #   geom_vline(aes(xintercept = lambda[j]), color = "red") +
    #   geom_vline(aes(xintercept = lambda_star), color = "green")
    
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
  
  for (j in seq_len(S_star)) {
    
    a_lambda_prior <- lambda_prior[S + j]^2 / sigma_lambda^2
    b_lambda_prior <- lambda_prior[S + j] / sigma_lambda^2
    
    # nonPCRcounts <- as.vector(y[,,j])[as.vector(c_imk[,,j]) == 2]
    # nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
    
    # psi_gig <- length(nonPCRcounts) * lambdatilde[S + j] + b_lambda_prior
    psi_gig <- b_lambda_prior
    
    PCRcounts <- as.vector(lambda_ijk[,,S + j])[as.vector(c_imk[,,S + j]) == 1]
    PCRcounts <- PCRcounts[!is.na(PCRcounts)]
    
    vpu <- u 
    vpu_PCR <- as.vector(vpu)[as.vector(c_imk[,,S + j]) == 1]
    vpu_PCR <- vpu_PCR[!is.na(vpu_PCR)]
    
    chi_gig <- sum(PCRcounts * r_nb[S + j] / exp(vpu_PCR)) 
    
    lambda_gig <- a_lambda_prior - length(PCRcounts) * r_nb[S + j]
    
    lambda[S + j] <- log(GIGrvg::rgig(1, lambda = lambda_gig, chi = 2 * chi_gig, psi = 2 * psi_gig ))
    
  }
  
  list_SP_cpp = convertCPtoSP_cpp(beta_bar,
                                  lambda, mu_bar,
                                  logz_bar, v_bar,
                                  delta,
                                  gamma,
                                  beta_theta_bar,
                                  M_site, S_star,
                                  emptyTubes)
  
  beta0 <- list_SP_cpp$beta0
  mu <- list_SP_cpp$mu
  logz <- list_SP_cpp$logz
  v[,1:S] <- list_SP_cpp$v[,1:S]
  
  list("lambda" = lambda,
       "beta0" = beta0,
       "mu" = mu,
       "v" = v,
       "logz" = logz)
  
}


update_lambda_CP_beta0 <- function(beta0, beta_z, X_z, logz, 
                             mu, lambda, v, u, lambda_ijk, r_nb,
                             c_imk, delta, gamma, X_w, beta_theta, 
                             M_site, tau, sigma_mu,
                             lambda_prior, sigma_lambda,
                             S_star, emptyTubes){
  
  df_t  = 3
  
  S <- length(mu)
  
  list_CP_cpp <- convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, 
                                   delta, gamma, beta_theta, M_site, S_star,
                                   emptyTubes)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  beta_theta_bar <- list_CP_cpp$beta_theta_bar
  
  # update paramters
  
  for (j in 1:S) {
    
    # nonPCRcounts <- as.vector(y[,,j])[as.vector(c_imk[,,j]) == 2]
    # nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
    # a <- length(nonPCRcounts) * lambdatilde[j]
    
    Xwbetatheta <- beta_theta[j,1] + X_w %*% beta_theta[j,-c(1,2)]
    
    Xbetalogz <- X_z %*% beta_z[,j]
    
    X_l <- rep(logz_bar[,j], M_site)
    
    optim_fun <- function(x){
      -logdpost_cpp_beta0(x, X_l, 
                    beta_theta[j,2], Xwbetatheta, Xbetalogz, delta[,j], 
                    logz_bar[,j], tau[j], mu_bar[j], sigma_mu, lambda_prior[j], 
                    sigma_lambda)
    }
    
    # lambda_star <- optimize(optim_fun, c(lambda[j] - 10, lambda[j] + 10), tol = .001)$minimum
    # try two mins
    # optim1 <- 
    # lambda_star_1 <- optimize(optim_fun, c(lambda[j] - 1, lambda[j] + 1), tol = .001)$minimum
    # lambda_obj_1 <- optimize(optim_fun, c(lambda[j] - 1, lambda[j] + 1), tol = .001)$objective
    # lambda_star_2 <- optimize(optim_fun, c(lambda[j] - 4, lambda[j] + 4), tol = .01)$minimum
    # lambda_obj_2 <- optimize(optim_fun, c(lambda[j] - 4, lambda[j] + 4), tol = .01)$objective
    # lambda_star_3 <- optimize(optim_fun, c(lambda[j] - 10, lambda[j] + 10), tol = .01)$minimum
    # lambda_obj_3 <- optimize(optim_fun, c(lambda[j] - 10, lambda[j] + 10), tol = .01)$objective
    # if(lambda_obj_1 < lambda_obj_2 & lambda_obj_1 < lambda_obj_3){
    #   lambda_star <- lambda_star_1
    # } else if(lambda_obj_2 < lambda_obj_1 & lambda_obj_2 < lambda_obj_3){
    #   lambda_star <- lambda_star_2
    # } else {
    #   lambda_star <- lambda_star_3
    # }
    lambda_star <- optimize(optim_fun, c(lambda[j] - 10, lambda[j] + 10))$minimum
    
    sd_star <-  1 / sqrt(-der2_logdpost_cpp_beta0(lambda_star, X_l, 
                                            beta_theta[j,2], Xwbetatheta, Xbetalogz, delta[,j], 
                                            logz_bar[,j], tau[j], mu_bar[j], sigma_mu, lambda_prior[j], 
                                            sigma_lambda))
    
    # lambda_grid <- seq(lambda[j] - 3, lambda[j] + 3, length.out = 100)
    # # lambda_grid <- seq(lambda_star - 01, lambda_star + 01, length.out = 100)
    # y_grid <- sapply(lambda_grid, function(x){
    #   # log(dt2(x, lambda_star, sd_star, df_t))
    #   # optim_fun(x)
    #   logdpost_cpp_beta0(x, X_l,
    #                      beta_theta[j,2], Xwbetatheta, Xbetalogz, delta[,j], 
    #                      logz_bar[,j], sigma_beta, mu_bar[j], sigma_mu, lambda_prior[j], 
    #                      sigma_lambda)
    # })
    # # 
    # qplot(lambda_grid, y_grid) +
    # qplot(lambda_grid, exp(y_grid - max(y_grid)) / sum(exp(y_grid - max(y_grid)))) +
    #   geom_vline(aes(xintercept = lambda_true[j])) +
    #   geom_vline(aes(xintercept = lambda[j]), color = "red") +
    #   geom_vline(aes(xintercept = lambda_star), color = "green")
    
    lambda_new <- rt2(lambda_star, sd_star, 3)
    
    logproposal_ratio <-  log(dt2(lambda[j], lambda_star, sd_star, df_t)) - 
      log(dt2(lambda_new, lambda_star, sd_star, df_t))
    
    logposterior <- logdpost_cpp_beta0(lambda[j], X_l, 
                                 beta_theta[j,2], Xwbetatheta, Xbetalogz, delta[,j], 
                                 logz_bar[,j], tau[j], mu_bar[j], sigma_mu, lambda_prior[j], 
                                 sigma_lambda)
    logposterior_star <- logdpost_cpp_beta0(lambda_new, X_l, 
                                      beta_theta[j,2], Xwbetatheta, Xbetalogz, delta[,j], 
                                      logz_bar[,j], tau[j], mu_bar[j], sigma_mu, lambda_prior[j], 
                                      sigma_lambda)
    
    (mh_ratio <- exp(logposterior_star - logposterior + logproposal_ratio))
    
    if(runif(1) < mh_ratio){
      lambda[j] <- lambda_new
    }
    
  }
  
  for (j in seq_len(S_star)) {
    
    a_lambda_prior <- lambda_prior[S + j]^2 / sigma_lambda^2
    b_lambda_prior <- lambda_prior[S + j] / sigma_lambda^2
    
    # nonPCRcounts <- as.vector(y[,,j])[as.vector(c_imk[,,j]) == 2]
    # nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
    
    # psi_gig <- length(nonPCRcounts) * lambdatilde[S + j] + b_lambda_prior
    psi_gig <- b_lambda_prior
    
    PCRcounts <- as.vector(lambda_ijk[,,S + j])[as.vector(c_imk[,,S + j]) == 1]
    PCRcounts <- PCRcounts[!is.na(PCRcounts)]
    
    vpu <- u 
    vpu_PCR <- as.vector(vpu)[as.vector(c_imk[,,S + j]) == 1]
    vpu_PCR <- vpu_PCR[!is.na(vpu_PCR)]
    
    chi_gig <- sum(PCRcounts * r_nb[S + j] / exp(vpu_PCR)) 
    
    lambda_gig <- a_lambda_prior - length(PCRcounts) * r_nb[S + j]
    
    lambda[S + j] <- log(GIGrvg::rgig(1, lambda = lambda_gig, chi = 2 * chi_gig, psi = 2 * psi_gig ))
    
  }
  
  list_SP_cpp = convertCPtoSP_cpp(beta_bar,
                                  lambda, mu_bar,
                                  logz_bar, v_bar,
                                  delta,
                                  gamma,
                                  beta_theta_bar,
                                  M_site, S_star,
                                  emptyTubes)
  
  beta0 <- rep(0, S)
  mu <- list_SP_cpp$mu
  logz <- list_SP_cpp$logz
  v[,1:S] <- list_SP_cpp$v[,1:S]
  
  list("lambda" = lambda,
       "beta0" = beta0,
       "mu" = mu,
       "v" = v,
       "logz" = logz)
  
}

update_lambda_NP <- function(lambda_ijk, c_imk, mu,
                             r_nb, v, u,
                             lambda_prior, sigma_lambda){
  
  S <- length(mu)
  
  # update paramters
  
  for (j in seq_len(S)) {
    
    mean_lambda <- lambda_prior[j]
    mean_explambda <- exp(mean_lambda + sigma_lambda^2 / 2)
    var_explambda <- exp(sigma_lambda^2 / 2 - 1) *  exp(2 * mean_lambda + sigma_lambda^2)
    a_lambda_prior <- mean_explambda^2 / var_explambda
    b_lambda_prior <- mean_explambda / var_explambda
    # a_lambda_prior <- exp(lambda_prior)[j]^2 / sigma_lambda^2
    # b_lambda_prior <- exp(lambda_prior)[j] / sigma_lambda^2
    
    # nonPCRcounts <- as.vector(y[,,j])[as.vector(c_imk[,,j]) == 2]
    # nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
    # lengthnonPCRcounts <- sum(as.vector(c_imk[,,j]) == 2)
    
    # psi_gig <- length(nonPCRcounts) + b_lambda_prior
    psi_gig <- b_lambda_prior
    
    PCRcounts <- as.vector(lambda_ijk[,,j])[as.vector(c_imk[,,j]) == 1]
    PCRcounts <- PCRcounts[!is.na(PCRcounts)]
    
    vpu <- u + matrix(v[,j], nrow = sum(M_site) + emptyTubes, ncol = max(K), byrow = F) 
    vpu_PCR <- as.vector(vpu)[as.vector(c_imk[,,j]) == 1]
    vpu_PCR <- vpu_PCR[!is.na(vpu_PCR)]
    
    chi_gig <- sum(PCRcounts * r_nb[j] / exp(vpu_PCR)) 
    
    lambda_gig <- a_lambda_prior - length(PCRcounts) * r_nb[j]
    
    lambda[j] <- log(rgig(lambda = lambda_gig, chi = 2 * chi_gig, psi = 2 * psi_gig))
    
  }
  
  lambda
}

# LAMBDA 0 ----------------------------------------------------------------

update_lambda0_mixt <- function(y, c_imk, mu0, n0, sd_mu0 = .05, sd_n0 = .025){
  
  nonPCRcounts <- as.vector(y)[as.vector(c_imk) == 0]
  nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
  
  numZeroCounts <- sum(nonPCRcounts > 0)
  
  pi0 <- rbeta(1, 1 + length(nonPCRcounts) - numZeroCounts, 1 + numZeroCounts)
  
  # propose new sets of parameters
  
  nonzeroPCRcounts <- nonPCRcounts[nonPCRcounts > 0]
  
  mu0_star <- rnorm(1, mu0, sd_mu0)
  n0_star <- rnorm(1, n0, sd_n0)
  
  if(n0_star > 1 & mu0_star > 0){
    lik_star <- sum(dnbinom(nonzeroPCRcounts, mu = mu0_star, size = n0_star, log = T))
    lik_current <- sum(dnbinom(nonzeroPCRcounts, mu = mu0, size = n0, log = T))
    
    if(runif(1) < exp(lik_star - lik_current)){
      mu0 <- mu0_star
      n0 <- n0_star
    }
    
  }
  
  list("mu0" = mu0,
       "n0" = n0,
       "pi0" = pi0)
}

# LAMBDA TILDE ------------------------------------------------------------

update_lambda_tilde_NB <- function(y, c_imk, mu_tilde, n_tilde, sd_mu0 = 1, sd_n0 = 1){
  
  nonPCRcounts <- as.vector(y)[as.vector(c_imk) == 2]
  nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
  
  # propose new sets of parameters
  mu_tilde_star <- rnorm(1, mu_tilde, sd_mu0)
  n_tilde_star <- exp(rnorm(1, log(n_tilde), sd_n0))
  
  if(n_tilde_star > 1 & mu_tilde_star > 0){
    lik_star <- sum(dnbinom(nonPCRcounts, mu = mu_tilde_star, size = n_tilde_star, log = T)) + 
      n_tilde_star
    lik_current <- sum(dnbinom(nonPCRcounts, mu = mu_tilde, size = n_tilde, log = T)) + 
      n_tilde
    
    logproposal <- log(n_tilde_star) - log(n_tilde)
    
    if(runif(1) < exp(lik_star - lik_current + logproposal)){
      mu_tilde <- mu_tilde_star
      n_tilde <- n_tilde_star
    }
    
  }
  
  list("mu_tilde" = mu_tilde,
       "n_tilde" = n_tilde)
}

# TAU ---------------------------------------------------------------------

update_Tau <- function(X_z, logz, beta0, beta_z,
                       Tau_params, Tau_priors,
                       jointSpecies){
  
  
  if(jointSpecies){
    
    logz_tilde <- logz - cbind(1, X_z) %*% rbind(t(beta0), beta_z)
    
    n <- nrow(logz_tilde)
    d <- ncol(logz_tilde)
    
    # sum_S <- 0
    # S_matrices <- sapply(1:n, function(i){
    #   sum_S <<- sum_S + logz_tilde[i,] %*% t(logz_tilde[i,]) 
    # })
    
    sum_S <- t(logz_tilde) %*% logz_tilde
    
    lambda_Y <- Tau_priors$lambda_Y
    # Tau_params <- updateGHParameters(n, sum_S, Tau_params, lambda_Y)
    Tau_params <- sample_GraphHorseshoe(n, sum_S, Tau_params, lambda_Y)
    # list_Tau <- update_Tau(invTau, Tau, tau_NG, lambda_NG, gamma_NG, lambda_Y,
    #                        logz, beta0, beta_z, X_z)
    
    # invTau <- list_Tau$Omega
    # Tau <- list_Tau$Sigma
    # tau_NG <- list_Tau$tau
    # lambda_NG <- list_Tau$lambda
    # gamma_NG <- list_Tau$gamma
  } else {
    a_tau <- Tau_priors$a_tau
    b_tau <- Tau_priors$b_tau
    
    tau <- Tau_params$tau
    S <- length(tau)
    tau <- update_tau_cpp(tau, logz, X_z, beta_z, beta0,
                          a_tau, rep(b_tau, S)
                          # a_tau = .5, b_tau = 1 / a_tau
    )
    
    Tau_params <- list("tau" = tau)
  }
  
  
  Tau_params
}

# BETA --------------------------------------------------------------------

update_betaz_CP <- function(beta0, beta_z, logz, tau, X_z, sigma_beta, updatebeta0){
  
  ncov_z <- ncol(X_z)
  
  if(ncov_z > 0 | updatebeta0){
    
    S <- length(beta0)
    
    if(updatebeta0){
      X_beta <- cbind(1, X_z)
    } else {
      X_beta <- X_z
    }
    
    {
      tXX <- t(X_beta) %*% X_beta
      
      for (j in 1:S) {
        
        Lambda_beta <- (tXX / tau[j]^2) + (diag(1, nrow = ncov_z + updatebeta0) / sigma_beta^2)
        
        mu_beta_x <- matrix(apply(X_beta, 2, function(x){
          x * logz[,j]
        }) / tau[j]^2, nrow(X_beta), updatebeta0 + ncov_z)
        
        mu_beta <- apply(mu_beta_x, 2, sum)
        b_prior <- rep(0, ncov_z)
        if(updatebeta0){
          b_prior <- c(0, b_prior)
        }
        
        mu_beta <- mu_beta + b_prior
        
        beta_bar_beta <- mvrnorm(1, solve(Lambda_beta) %*% mu_beta, solve(Lambda_beta))
        
        if(updatebeta0){
          beta0[j] <- beta_bar_beta[1]
          beta_z[,j] <- beta_bar_beta[-1]  
        } else {
          beta_z[,j] <- beta_bar_beta
        }
        
      }
    }
    
  }
  
  list("beta0" = beta0,
       "beta_z" = beta_z)
}

update_betaz_CP_corr <- function(beta0, beta_z, logz, Tau, X_z, sigma_beta, updatebeta0){
  
  ncov_z <- ncol(X_z)
  S <- ncol(Tau)
  n <- nrow(X_z)
  
  if(ncov_z > 0 | updatebeta0){
    
    if(updatebeta0){
      X_beta <- cbind(1, X_z)
    } else {
      X_beta <- X_z
    }
    
    Tau <- Matrix(Tau)
    X_beta <- Matrix(X_beta)
    Id_n <- Matrix(diag(1, nrow = n))
    Id_s <- Matrix(diag(1, nrow = S))
    
    invSigma_tilde <- kronecker(Id_n, solve(Tau))
    
    X_tilde <-  kronecker(X_beta, Id_s)
    
    Lambda_beta <- t(X_tilde) %*% invSigma_tilde %*% X_tilde + diag(1 / sigma_beta^2, S * ncol(X_beta))
    
    term1 <- invSigma_tilde %*% as.vector(t(logz))
    txsigmalogz <- t(X_tilde) %*% term1
    mu_beta <- solve(Lambda_beta) %*% txsigmalogz
    # mu_beta2 <- solve(Lambda_beta) %*% t(X_tilde) %*% invSigma_tilde %*% as.vector(t(logz))
    
    betavec <- t(mvrnorm(1, as.vector(mu_beta), as.matrix(solve(Lambda_beta))))
    betamat <- matrix(betavec, ncol(X_beta), S, byrow = F)
    
    #
    
    # tXX <- t(X_beta) %*% X_beta
    # tXl <- t(X_beta) %*% logz
    # 
    # prior_mean <- matrix(0, nrow = ncov_z + updatebeta0, ncol = S)
    # # prior_mean[1,] <- rep(0, S)
    # 
    # M_term <- tXl + prior_mean
    # U_term <- tXX + diag(sigma_beta^2, nrow = (ncov_z + updatebeta0))
    # post_U <- solve(U_term)
    # post_M <- post_U %*% M_term
    # 
    # beta_bar_beta <- rmtrnorm(post_M, post_U, Tau)
    
    if(updatebeta0){
      beta0 <- as.matrix(betamat[1,])
      beta_z <- as.matrix(betamat[-1,,drop = F])  
    } else {
      beta_z <- as.matrix(betamat)
    }
    
  }
  
  list("beta0" = beta0,
       "beta_z" = beta_z)
}
