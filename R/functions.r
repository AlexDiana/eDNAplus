# GENERAL -----------------------------------------------------------------

fitModel <- function(data,
                     priors,
                     jointSpecies,
                     trueParams,
                     paramsUpdate = list(updateAll = T,
                                         params = NULL,
                                         correct = NULL),
                     MCMCparams,
                     paramsToSave = list(
                       "lambda" = T,
                       "beta0" = T,
                       "mu" = T,
                       "beta_z" = T,
                       "logz" = T,
                       "u" = F,
                       "v" = F,
                       "beta_w" = T,
                       "Tau" = T,
                       "sigma" = T,
                       "beta_theta" = T,
                       "theta" = F,
                       "csi" = T,
                       "delta" = F,
                       "gamma" = F,
                       "r" = T,
                       "p11" = T,
                       "p10" = T,
                       "p10" = T,
                       "mutilde" = T,
                       "c_imk" = F,
                       "mu0" = T,
                       "n0" = T,
                       "pi0" = T,
                       "eta" = F
                     )){
  
  # mcmc params
  {
    nchain <- MCMCparams$nchain
    nburn <- MCMCparams$nburn
    niter <- MCMCparams$niter
    nthin <- MCMCparams$nthin
    iterToAdapt <- MCMCparams$iterToAdapt  
  }
  
  # clean data
  {
    y <- data$y
    M_site <- data$M_site
    K <- data$K
    emptyTubes <- data$emptyTubes
    S_star <- data$S_star
    X_w <- data$X_w
    X_z <- data$X_z
    v_spikes <- data$v_spikes
    spikedSample <- data$spikedSample
    
    ncov_z <- ncol(X_z)
    ncov_w <- ncol(X_w)
    n <- length(M_site)
    S <- dim(y)[3] - S_star
  }
  
  # priors
  {
    lambda_prior <- priors$lambda_prior
    sigma_beta <- priors$sigma_beta
    sigma_lambda <- priors$sigma_lambda
    sigma_gamma <- priors$sigma_gamma
    sigma_beta_theta <- priors$sigma_beta_theta
    sigma_u <- priors$sigma_u
    a_p11 <- priors$a_p11
    b_p11 <- priors$b_p11
    a_p10 <- priors$a_p10
    b_p10 <- priors$b_p10
    a_theta0 <- priors$a_theta0
    b_theta0 <- priors$b_theta0
    
    # set lambda prior
    {
      if(is.null(lambda_prior)){
        lambda_prior <- apply(y, 3, function(x){
          log(mean(x[x > quantile(x[x > 0], probs = c(0.1), na.rm = T)], na.rm = T))
        })
      }
    }
    
    if(jointSpecies){
      lambda_Y <- 1
      Tau_Priors <- list("lambda_Y" = lambda_Y)
    } else {
      
      a_tau <- 1
      b_tau <- 1
      
      Tau_Priors <- list("a_tau" = a_tau,
                         "b_tau" = rep(1, S))
    }
  }
  
  # params to update
  {
    if(paramsUpdate$updateAll){
      
      updateLambda_CP <- T; correctLambda <- F
      updateLambda_NP <- F; correctLambda <- F
      
      updateBeta_w <- T; correctBeta_w <- F
      updateBeta_z <- T; correctBeta_z <- F
      updateMu <- T; correctMu <- F
      
      updateL <- T; correctL <- F
      
      updateV <- T; correctV <- F
      updateU <- T; correctU <- F
      updateUV <- T;
      
      updateLambdaijk <- T; correctLambdaijk <- F
      
      updateSigma <- T; correctSigma <- F
      updateTau <- T; correctTau <- F
      
      updateR <- T; correctR <- F
      
      updateDeltaGammaC <- T; correctDeltaGammaC <- F
      
      updateLambdaTilde <- T; correctLambdaTilde <- F
      updateLambda0 <- T; correctLambda0 <- F
      
      updateP11 <- T; correctP11 <- F
      updateP10 <- T; correctP10 <- F
      
      updateBetaTheta <- T; correctBetaTheta <- F
      updateTheta10 <- T; correctTheta10 <- F
    } else {
      
      paramsToUpdate <- paramsUpdate$params
      
      updateLambda_CP <- paramsToUpdate$lambda; correctLambda <- !updateLambda_CP
      updateLambda_NP <- F
      
      updateL <- paramsToUpdate$logz; correctL <- !updateL
      updateMu <- paramsToUpdate$mu; correctMu <- !updateMu
      
      updateV <- paramsToUpdate$v; correctV <- !updateV
      updateU <- paramsToUpdate$u; correctU <- !updateU
      updateUV <- paramsToUpdate$uv; correctU <- !updateUV; correctV <- !updateUV
      
      updateBetaTheta <- paramsToUpdate$beta_theta; correctBetaTheta <- !updateBetaTheta
      updateTheta10 <- paramsToUpdate$csi; correctTheta10 <- !updateTheta10
      
      updateBeta_w <- paramsToUpdate$beta_w; correctBeta_w <- !updateBeta_w
      updateBeta_z <- paramsToUpdate$beta_z; correctBeta_z <- !updateBeta_z
      
      updateSigma <- paramsToUpdate$sigma; correctSigma <- !updateSigma
      updateTau <- paramsToUpdate$tau; correctTau <- !updateTau
      
      updateDeltaGammaC <- paramsToUpdate$deltagammac; correctDeltaGammaC <- !updateDeltaGammaC#F
      
      updateLambdaijk <- paramsToUpdate$eta; correctLambdaijk <- !updateLambdaijk
      updateR <- paramsToUpdate$r; correctR <- !updateR
      
      updateLambdaTilde <- paramsToUpdate$lambdatilde; correctLambdaTilde <- !updateLambdaTilde
      updateLambda0 <- paramsToUpdate$lambda0; correctLambda0 <- !updateLambda0
      
      updateP11 <- paramsToUpdate$p11; correctP11 <- !updateP11
      updateP10 <- paramsToUpdate$p10; correctP10 <- !updateP10
      
      if(!is.null(paramsUpdate$correct)){
        
        paramsCorrect <- paramsUpdate$correct
        
        correctLambda <- paramsCorrect$lambda
        
        correctL <- paramsCorrect$logz
        correctMu <- paramsCorrect$mu
        
        correctV <- paramsCorrect$v
        correctU <- paramsCorrect$u
        
        
        correctBetaTheta <- paramsCorrect$beta_theta
        correctTheta10 <- paramsCorrect$csi
        
        correctBeta_w <- paramsCorrect$beta_w
        correctBeta_z <- paramsCorrect$beta_z
        
        correctSigma <- paramsCorrect$sigma
        correctTau <- paramsCorrect$tau
        
        correctDeltaGammaC <- paramsCorrect$deltagammac
        
        correctLambdaijk <- paramsCorrect$eta
        correctR <- paramsCorrect$r
        
        correctLambdaTilde <- paramsCorrect$lambdatilde
        correctLambda0 <- paramsCorrect$lambda0
        
        correctP11 <- paramsCorrect$p11
        correctP10 <- paramsCorrect$p10
        
      }
      
    }
    
  }
  
  # output
  {
    
    if(paramsToSave$lambda){
      lambda_output <- array(NA, c(nchain, S + S_star, niter))
    } else {
      lambda_output <- NULL
    }
    
    if(paramsToSave$beta0){
      beta0_output <- array(NA, dim = c(nchain, S, niter))
    } else {
      beta0_output <- NULL
    }
    
    if(paramsToSave$mu){
      mu_output <- array(NA, dim = c(nchain, S, niter))
    } else {
      mu_output <- NULL
    }
    
    if(paramsToSave$beta_z){
      beta_z_output <- array(NA, dim = c(nchain, ncov_z, S, niter))
    } else {
      beta_z_output <- NULL
    }
    
    if(paramsToSave$logz){
      logz_output <- array(NA, dim = c(nchain, n, S, niter))
    } else {
      logz_output <- NULL
    }
    
    if(paramsToSave$u){
      u_output <- array(NA, c(nchain, sum(M_site), max(K), niter))
    } else {
      u_output <- NULL
    }
    
    if(paramsToSave$v){
      v_output <- array(NA, dim = c(nchain, sum(M_site), S + S_star, niter))
    } else {
      v_output <- NULL
    }
    
    if(paramsToSave$beta_w){
      beta_w_output <- array(NA, c(nchain, ncov_w, S, niter))
    } else {
      beta_w_output <- NULL
    }
    
    if(paramsToSave$Tau){
      if(jointSpecies){
        Tau_output <- array(NA, dim = c(nchain, niter, S, S))
      } else {
        Tau_output <- array(NA, c(nchain, S, niter))
      }
    } else {
      Tau_output <- NULL
    }
    
    if(paramsToSave$sigma){
      sigma_output <- array(NA, c(nchain, S, niter))
    } else {
      sigma_output <- NULL
    }
    
    if(paramsToSave$beta_theta){
      beta_theta_output <- array(NA, dim = c(nchain, S, 2 + ncov_w, niter))
    } else {
      beta_theta_output <- NULL
    }
    
    if(paramsToSave$theta){
      theta_output <- array(NA, dim = c(nchain, sum(M_site), S, niter))
    } else {
      theta_output <- NULL
    }
    
    if(paramsToSave$csi){
      csi_output <- array(NA, dim = c(nchain, S, niter))
    } else {
      csi_output <- NULL
    }
    
    if(paramsToSave$delta){
      delta_output <- array(NA, dim = c(nchain, sum(M_site) + emptyTubes, S + S_star, niter))
    } else {
      delta_output <- NULL
    }
    
    if(paramsToSave$gamma){
      gamma_output <- array(NA, dim = c(nchain, sum(M_site) + emptyTubes, S + S_star, niter))
    } else {
      gamma_output <- NULL
    }
    
    if(paramsToSave$r){
      r_output <- array(NA, dim = c(nchain, S + S_star, niter))
    } else {
      r_output <- NULL
    }
    
    if(paramsToSave$p11){
      p_11_output <- array(NA, dim = c(nchain, S + S_star, niter))
    } else {
      p_11_output <- NULL
    }
    
    if(paramsToSave$p10){
      p_10_output <- array(NA, dim = c(nchain, S + S_star, niter))
    } else {
      p_10_output <- NULL
    }
    
    if(paramsToSave$mutilde){
      mutilde_output <- array(NA, c(nchain, niter))
    } else {
      mutilde_output <- NULL
    }
    
    if(paramsToSave$c_imk){
      cimk_output <- array(NA, dim = c(nchain, sum(M_site) + emptyTubes, max(K), S + S_star, niter))
    } else {
      cimk_output <- NULL
    }
    
    if(paramsToSave$mu0){
      mu0_output <- matrix(NA, nchain, niter)  
    } else {
      mu0_output <- NULL
    }
    
    if(paramsToSave$n0){
      n0_output <- matrix(NA, nchain, niter) 
    } else {
      n0_output <- NULL
    }
    
    if(paramsToSave$pi0){
      pi0_output <- matrix(NA, nchain, niter) 
    } else {
      pi0_output <- NULL
    }
    
    if(paramsToSave$eta){
      eta_output <- array(NA, dim = c(nchain, sum(M_site) + emptyTubes, max(K), S + S_star, niter))
    } else {
      eta_output <- NULL
    }
    
  }
  
  for (chain in 1:nchain) {
    
    # chain output
    {
      if(paramsToSave$lambda){
        lambda_output_iter <- matrix(NA, nrow = S + S_star, ncol = niter)
      } else {
        lambda_output_iter <- NULL
      }
      
      if(paramsToSave$beta0){
        beta0_output_iter <- array(NA, dim = c(S, niter))
      } else {
        beta0_output_iter <- NULL
      }
      
      if(paramsToSave$mu){
        mu_output_iter <- array(NA, dim = c(S, niter))
      } else {
        mu_output_iter <- NULL
      }
      
      if(paramsToSave$beta_z){
        beta_z_output_iter <- array(NA, dim = c(ncov_z, S, niter))
      } else {
        beta_z_output_iter <- NULL
      }
      
      if(paramsToSave$logz){
        logz_output_iter <- array(NA, dim = c(n, S, niter))
      } else {
        logz_output_iter <- NULL
      }
      
      if(paramsToSave$u){
        u_output_iter <- array(NA, dim = c(sum(M_site), max(K), niter))
      } else {
        u_output_iter <- NULL
      }
      
      if(paramsToSave$v){
        v_output_iter <- array(NA, dim = c(sum(M_site), S + S_star, niter))
      } else {
        v_output_iter <- NULL
      }
      
      if(paramsToSave$beta_w){
        beta_w_output_iter <- array(NA, c(ncov_w, S, niter))
      } else {
        beta_w_output_iter <- NULL
      }
      
      if(paramsToSave$Tau){
        if(jointSpecies){
          Tau_output_iter <- array(NA, dim = c(niter, S, S))
        } else {
          Tau_output_iter <- matrix(NA, S, niter)
        }
      } else {
        Tau_output_iter <- NULL
      }
      
      if(paramsToSave$sigma){
        sigma_output_iter <- matrix(NA, S, niter)
      } else {
        sigma_output_iter <- NULL
      }
      
      if(paramsToSave$beta_theta){
        beta_theta_output_iter <- array(NA, dim = c(S, 2 + ncov_w, niter))
      } else {
        beta_theta_output_iter <- NULL
      }
      
      if(paramsToSave$theta){
        theta_output_iter <- array(NA, dim = c(sum(M_site), S, niter))
      } else {
        theta_output_iter <- NULL
      }
      
      if(paramsToSave$csi){
        csi_output_iter <- array(NA, dim = c(S, niter))
      } else {
        csi_output_iter <- NULL
      }
      
      if(paramsToSave$delta){
        delta_output_iter <- array(NA, dim = c(sum(M_site) + emptyTubes, S + S_star, niter))
      } else {
        delta_output_iter <- NULL
      }
      
      if(paramsToSave$gamma){
        gamma_output_iter <- array(NA, dim = c(sum(M_site) + emptyTubes, S + S_star, niter))
      } else {
        gamma_output_iter <- NULL
      }
      
      if(paramsToSave$r){
        r_output_iter <- matrix(NA, S + S_star, niter)
      } else {
        r_output_iter <- NULL
      }
      
      if(paramsToSave$p11){
        p_11_output_iter <- array(NA, dim = c(S + S_star, niter))
      } else {
        p_11_output_iter <- NULL
      }
      
      if(paramsToSave$p10){
        p_10_output_iter <- array(NA, dim = c(S + S_star, niter))
      } else {
        p_10_output_iter <- NULL
      }
      
      if(paramsToSave$mutilde){
        mutilde_output_iter <- rep(NA, niter)
      } else {
        mutilde_output_iter <- NULL
      }
      
      if(paramsToSave$c_imk){
        cimk_output_iter <- array(NA, dim = c(sum(M_site) + emptyTubes, max(K), S  + S_star, niter))
      } else {
        cimk_output_iter <- NULL
      }
      
      if(paramsToSave$mu0){
        mu0_output_iter <- rep(NA, niter)
      } else {
        mu0_output_iter <- NULL
      }
      
      if(paramsToSave$n0){
        n0_output_iter <- rep(NA, niter)
      } else {
        n0_output_iter <- NULL
      }
      
      if(paramsToSave$pi0){
        pi0_output_iter <- rep(NA, niter)
      } else {
        pi0_output_iter <- NULL
      }
      
      if(paramsToSave$eta){
        eta_output_iter <- array(NA, dim = c(sum(M_site), max(K), S + S_star, niter))
      } else {
        eta_output_iter <- NULL
      }
      
    }
    
    # starting values
    {
      
      if(correctLambda){
        
        lambda <- lambda_true
        
      } else
      {
        lambda <- lambda_prior
        
      }
      
      if(correctL){
        
        logz <- logz_true
        
      } else
      {
        logz <- matrix(0, n, S)
      }
      
      if(correctMu){
        mu <- mu_true
      } else
      {
        mu <- rep(0, S)
      }
      
      if(correctV){
        v <- v_true
      } else
      {
        v <- matrix(0, sum(M_site) + emptyTubes, S)
        if(S_star > 0){
          v <- cbind(v, v_spikes)
        }
      }
      
      if(correctU){
        u <- u_true
      } else
      {
        u <- matrix(0, sum(M_site) + emptyTubes, max(K))
      }
      
      if(correctR){
        r_nb <- r_true
      } else
      {
        r_nb <- rpois(S + S_star, 100)
      }
      
      if(correctBeta_w){
        beta_w <- beta_w_true
      } else
      {
        beta_w <- matrix(0, ncov_w, S)
      }
      
      if(correctBeta_z){
        beta0 <- beta0_true
        beta_z <- beta_z_true
      } else
      {
        beta0 <- rep(0, S)
        beta_z <- matrix(0, nrow = ncov_z, ncol = S)
      }
      
      if(correctDeltaGammaC){
        c_imk <- c_imk_true
        
        delta <- delta_true
        
        gamma <- gamma_true
      } else
      {
        c_imk <- array(NA, dim = c(sum(M_site) + emptyTubes, max(K), S + S_star))
        for (i in 1:n) {
          for (m in 1:M_site[i]) {
            numRep <- K[m + sum(M_site[seq_len(i-1)])]
            for (k in 1:numRep) {
              for (j in 1:(S + S_star)) {
                
                # if(y[m + sum(M_site[seq_len(i-1)]),k,j] > 100){
                if(y[m + sum(M_site[seq_len(i-1)]),k,j] > sample(c(100,500), 1)){
                  c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- 1
                  # } else if(y[m + sum(M_site[seq_len(i-1)]),k,j] > 10){
                } else if(y[m + sum(M_site[seq_len(i-1)]),k,j] > sample(c(10,50), 1)){
                  c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- 2
                } else {
                  c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- 0
                }
                
              }
            }
          }
        }
        for (m in seq_len(emptyTubes)) {
          if(y[sum(M_site) + m,k,j] > sample(c(10,50), 1)){
            c_imk[sum(M_site) + m,k,j] <- 2
          } else {
            c_imk[sum(M_site) + m,k,j] <- 0
          }
        }
        
        delta <- array(NA, dim = c(sum(M_site) + emptyTubes, S + S_star))
        for (i in 1:n) {
          for (m in 1:M_site[i]) {
            for (j in 1:S) {
              if(sum(c_imk[m + sum(M_site[seq_len(i-1)]),
                           1:K[m + sum(M_site[seq_len(i-1)])],j]) > 0){
                delta[m + sum(M_site[seq_len(i-1)]),j] <- 1
              } else {
                delta[m + sum(M_site[seq_len(i-1)]),j] <- 0
              }
              
            }
            for (j in seq_len(S_star)) {
              delta[m + sum(M_site[seq_len(i-1)]),S + j] <- 1
            }
          }
        }
        for (m in seq_len(emptyTubes)) {
          delta[sum(M_site) + m,] <- 0
        }
        
        gamma <- array(NA, dim = c(sum(M_site) + emptyTubes, S + S_star))
        for (i in 1:n) {
          for (m in 1:M_site[i]) {
            for (j in 1:(S + S_star)) {
              gamma[m + sum(M_site[seq_len(i-1)]),j] <- 0
            }
          }
        }
        for (m in seq_len(emptyTubes)) {
          gamma[sum(M_site) + m,] <- 0
        }
        
      }
      
      if(correctLambda0){
        mu0 <- mu0_true
        n0 <- n0_true
        pi0 <- pi0_true
      } else
      {
        mu0 <- 2
        n0 <- 1
        pi0 <- .99
      }
      
      if(correctLambdaTilde){
        mu_tilde <- mu_tilde_true
        n_tilde <- n_tilde_true
      } else
      {
        mu_tilde <- 100
        n_tilde <- 100000
      }
      
      if(correctSigma){
        sigma <- sigma_true
      } else
      {
        sigma <- rep(1, S)
      }
      
      if(correctTau){
        if(jointSpecies){
          # invTau <- solve(Tau_true) # Omega_Y
          # Tau <- Tau_true # Sigma_y #solve(Omega)
          # tau_NG <- matrix(1, nrow = S, ncol = S)
          # lambda_NG <- 1
          # gamma_NG <- 1
          
          Tau_params <- list("Omega" = Tau_true,
                             "lambdasq" = matrix(1, nrow = S, ncol = S),
                             "tausq" = 1,
                             "nu" = matrix(1, nrow = S, ncol = S),
                             "csi" = 1,
                             "Sigma" = solve(Tau_true))
          
        } else {
          Tau_params <- list("tau" = tau_true)
        }
      } else
      {
        if(jointSpecies){
          # invTau <-  # Omega_Y
          # Tau <- diag(1, nrow = S) # Sigma_y #solve(Omega)
          # tau_NG <- matrix(1, nrow = S, ncol = S)
          # lambda_NG <- 1
          # gamma_NG <- 1
          
          Tau_params <- list("Omega" = diag(1, nrow = S),
                             "lambdasq" = matrix(1, nrow = S, ncol = S),
                             "tausq" = 1,
                             "nu" = matrix(1, nrow = S, ncol = S),
                             "csi" = 1,
                             "Sigma" = diag(1, nrow = S))
        } else {
          tau <- rep(1, S)
          Tau_params <- list("tau" = tau)
        }
      }
      
      if(correctP11){
        p_11 <- p11_true
      } else
      {
        p_11 <- rep(0.9, S + S_star)
      }
      
      if(correctP10){
        p_10 <- p10_true
      } else
      {
        p_10 <- rep(0.01, S + S_star)
      }
      
      if(correctBetaTheta){
        beta_theta <- beta_theta_true
        theta11 <- theta11_true
      } else
      {
        beta_theta <- cbind(rep(0, S),pmax(rnorm(S, 1), 0), matrix(0, nrow = S, ncol = ncov_w))
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
      
      if(correctTheta10){
        theta10 <- theta10_true
      } else
      {
        theta10 <- rep(.01, S)
      }
      
      if(correctLambdaijk){
        lambda_ijk <- y
      } else 
      {
        lambda_ijk <- array(NA, c(sum(M_site) + emptyTubes, max(K), S + S_star))
        for (i in 1:n) {
          # print(i)
          for (m in 1:M_site[i]) {
            numRep <- K[m + sum(M_site[seq_len(i-1)])]
            for (k in 1:numRep) {
              u_imk <- u[m + sum(M_site[seq_len(i-1)]), k]
              for (j in 1:(S + S_star)) {
                lambda_ijk[m + sum(M_site[seq_len(i - 1)]), k, j] <-  exp(v[m + sum(M_site[seq_len(i - 1)]), j] +
                                                                            lambda[j] +
                                                                            u_imk)
              }
            }
          }
        }
        for (m in seq_len(emptyTubes)) {
          numRep <- K[m + sum(M_site)]
          for (k in 1:numRep) {
            u_imk <- u[m + sum(M_site), k]
            for (j in 1:(S + S_star)) {
              lambda_ijk[m + sum(M_site), k, j] <-  exp(v[m + sum(M_site), j] +
                                                          lambda[j] +
                                                          u_imk)
            }
          }
        }
      }
      
    }
    
    for (iter in 1:(nburn + nthin * niter)) {
      print(apply(gamma, 2, mean))
      if(iter <= nburn){
        print(paste0("Chain = ",chain," - Burn-in Iteration = ",iter))
      } else {
        print(paste0("Chain = ",chain," - Iteration = ",iter - nburn))
      }
      
      # LAMBDA ----------------------------------------------
      
      if(updateLambda_CP){
        
        list_lambda <- update_lambda_CP(beta0, beta_z, logz, 
                                        mu, lambda, v, u, lambda_ijk, r_nb,
                                        c_imk, delta, gamma, beta_theta, 
                                        M_site, 
                                        sigma_beta, sigma_mu,
                                        lambda_prior, sigma_lambda,
                                        S_star)
        lambda <- list_lambda$lambda
        beta0 <- list_lambda$beta0
        mu <- list_lambda$mu
        v <- list_lambda$v
        logz <- list_lambda$logz
        
      }
      
      # LAMBDA ----------------------------------------------
      
      if(updateLambda_NP){
        # print("Update lambda")
        
        # lambda <- update_lambda(beta0, mu, lambda,
        #                              sigma_beta, sigma_mu,
        #                              exp(lambda_prior), sigma_lambda,
        #                              S_star, betaThetaEqual1)
        
        lambda <- update_lambda_NP(lambda_ijk, c_imk, mu,
                                   r_nb, v, u,
                                   lambda_prior, sigma_lambda)
        
      }
      
      # MU -----------------------------------
      
      if(updateMu){
        
        mu <- update_mu_cpp(mu, lambda, delta, gamma, sigma, sigma_gamma, beta0,
                            beta_z, logz, v, beta_theta, M_site, sigma_mu, S_star, emptyTubes)
        
      }
      
      
      # LAMBDA IJK ----
      
      if(updateLambdaijk){
        lambda_ijk <- update_lambdaijk(lambda, lambda_ijk, v, u, r_nb, c_imk, M_site, y, K,
                                       S_star, emptyTubes)
      }
      
      # VU ---------
      
      if(updateUV){
        list_uv <- update_uv_poisgamma_cpp(u, v, logz, lambda, X_z, beta_theta, beta_z, beta0,
                                           r_nb, mu, lambda_ijk, c_imk, delta, gamma, sigma, sigma_gamma,
                                           sigma_u, M_site, X_w, beta_w, K, S_star, emptyTubes)
        u <- list_uv$u
        v <- list_uv$v
        lambda <- list_uv$lambda
      }
      
      # V ------------------------------------------
      
      if(updateV){
        # print("Update v")
        
        v <- update_v_poisgamma_cpp(v, logz,
                                    lambda,  X_z,
                                    beta_theta, u, beta_z,
                                    beta0, r_nb, mu, lambda_ijk,
                                    c_imk, delta, gamma, sigma,
                                    sigma_gamma, M_site,
                                    X_w, beta_w,
                                    K, S_star, emptyTubes)
      }
      
      # U --------------------------------------
      
      if(updateU){
        # print("Update u")
        
        list_u <- update_u_poisgamma_cpp(v, u, lambda, beta0, beta_z, logz,
                                         mu, lambda_ijk, r_nb, X_w, beta_w, c_imk,
                                         delta, gamma, sigma_u, beta_theta, sigma,
                                         sigma_gamma, M_site,
                                         K, S_star, emptyTubes)
        u <- list_u$u
        lambda <- list_u$lambda
        
      }
      
      # LOG_Z ------------------------------
      
      if(updateL){
        
        if(jointSpecies){
          Tau <- Tau_params$Omega
          logz <- update_logz_corr_cpp(logz, beta0, X_z, beta_z, mu,
                                       v, lambda, beta_theta, X_w, beta_w,
                                       Tau, delta, gamma, sigma,
                                       M_site, S_star, emptyTubes)
        } else {
          tau <- Tau_params$tau
          logz <- update_logz_cpp(logz, beta0, X_z, beta_z,
                                  mu, v, lambda, beta_theta,
                                  X_w, beta_w, tau, delta, gamma, sigma,
                                  M_site, S_star, emptyTubes)
        }
        
      }
      
      # BETA_Z --------------
      
      if(updateBeta_z){
        if(jointSpecies){
          Tau <- Tau_params$Sigma
          list_beta_z <- update_betaz_CP_corr(beta0, beta_z, logz, Tau, X_z, sigma_beta)
        } else {
          tau <- Tau_params$tau
          list_beta_z <- update_betaz_CP(beta0, beta_z, logz, tau, X_z, sigma_beta)
        }
        
        beta0 <- list_beta_z$beta0
        beta_z <- list_beta_z$beta_z
        
      }
      
      # BETA_W ---------
      
      if(updateBeta_w){
        # print("Update betaw")
        if(ncov_w > 0){
          beta_w <- update_betaw_cpp(beta_w, v, delta, logz, X_w, sigma, sigma_beta, M_site)
        }
      }
      
      # TAU ----------------------------------------------------------
      
      Tau_params <- update_Tau(X_z, logz, beta0, beta_z,
                               Tau_params, Tau_Priors,
                               jointSpecies)
      
      # if(updateTau){
      #   
      #   if(jointSpecies){
      #     # list_Tau <- update_Tau_Wish(nu0, Sigma0, logz, beta0, betaz, X_z)
      #     
      #     # invTau <- list_Tau$invSigma
      #     # Tau <- list_Tau$Sigma
      #     list_Tau <- update_Tau(invTau, Tau, tau_NG, lambda_NG, gamma_NG, lambda_Y,
      #                            logz, beta0, beta_z, X_z)
      #     
      #     invTau <- list_Tau$Omega
      #     Tau <- list_Tau$Sigma
      #     tau_NG <- list_Tau$tau
      #     lambda_NG <- list_Tau$lambda
      #     gamma_NG <- list_Tau$gamma
      #   } else {
      #     tau <- update_tau_cpp(tau, logz, X_z, beta_z, beta0,
      #                           a_tau, rep(b_tau, S)
      #                           # a_tau = .5, b_tau = 1 / a_tau
      #     )
      #     # list_Tau <- update_Tau_Wish(nu0, Sigma0, logz, beta0, betaz, X_z)
      #     #
      #     # invTau <- list_Tau$invSigma
      #     # Tau <- list_Tau$Sigma
      #     # tau <- sqrt(diag(Tau))
      #   }
      #   
      # }
      
      # UPDATE a TAU ----------------------------------------------
      
      # a_tau <- sapply(1:S, function(j) {
      #   rinvgamma_cpp(.5, 1 / A_tau + 1 / tau[j]^2)
      # })
      
      # SIGMA ----------------------------------------------------------
      
      if(updateSigma){
        # print("Update sigma")
        # sigma <- rep(1, S)
        sigma <- update_sigma_cpp(sigma, lambda, beta_z, beta0,
                                  mu, logz, v, X_w, beta_w,
                                  delta, gamma, beta_theta,
                                  a_sigma, rep(b_sigma, S),
                                  # a_sigma = .5, b_sigma = 1 / a_sigma,
                                  M_site, S_star, emptyTubes)
        
      }
      
      # UPDATE a SIGMA ----------------------------------------------
      
      # a_sigma <- sapply(1:S, function(j) {
      #   rinvgamma_cpp(.5, 1 / A_sigma + 1 / sigma[j]^2)
      # })
      
      # R - NB ---------------------
      
      if(updateR){
        # print("Update r")
        
        r_nb <- update_r_nb_cpp(r_nb, lambda, u,
                                v, y, delta, gamma, c_imk, M_site, K,
                                optimStep =  F,
                                sd_r_proposal = .05)
        
      }
      
      # DELTA/C/D --------------------------------------------------------
      
      if(updateDeltaGammaC){
        
        v_pres <- (delta == 1) | (gamma == 1)
        list_deltagammac <- update_delta_c_d_rjmcmc(v_pres, y, v, lambda, r_nb,
                                                    M_site, K, 
                                                    mu0, n0, pi0, mu_tilde,
                                                    n_tilde, u, logz, X_w, beta_w,
                                                    sigma, mu, sigma_gamma, v_sd = .5,
                                                    p_11, p_10, theta11, theta10, spikedSample,
                                                    emptyTubes, S_star)
        delta <- list_deltagammac$delta
        gamma <- list_deltagammac$gamma
        c_imk <- list_deltagammac$c_imk
        v <- list_deltagammac$v
        
      }
      
      # LAMBDA IJK ----
      
      if(updateLambdaijk){
        
        lambda_ijk <- update_lambdaijk(lambda, lambda_ijk, v, u, r_nb, c_imk, M_site, y, K,
                                       S_star, emptyTubes)
      }
      
      # LAMBDA 0 ---------------------------------------------------------
      
      if(updateLambda0){
        # print("Update lambda0")
        
        list_lambda0 <- update_lambda0_mixt(y, c_imk, mu0, n0, sd_mu0 = .05, sd_n0 = .05)
        mu0 <- list_lambda0$mu0
        n0 <- list_lambda0$n0
        pi0 <- list_lambda0$pi0
      }
      
      # LAMBDA TILDE ---------------------------------------------------------
      
      if(updateLambdaTilde){
        # print("Update lambda tilde")
        
        # lambdatilde <- update_lambdatilde(y, c_imk, lambda)
        list_lambda_tilde <- update_lambda_tilde_NB(y, c_imk, mu_tilde, n_tilde, sd_mu0 = 1, sd_n0 = 1)
        mu_tilde <- list_lambda_tilde$mu_tilde
        n_tilde <- list_lambda_tilde$n_tilde
      }
      
      # P11 --------------------------------------------------------------
      
      if(updateP11){
        # print("Update p11")
        
        p_11 <- update_p_11_cpp(p_11, delta, gamma, c_imk, M_site,
                                K, a_p11, b_p11)
      }
      
      # P10 --------------------------------------------------------------
      
      if(updateP10){
        
        p_10 <- update_p_10_cpp(p_10, delta, gamma, c_imk, M_site, K, a_p10, b_p10, emptyTubes)
      }
      
      # UPDATE BETA THETA11 --------------------------------------------------------------
      
      if(updateBetaTheta){
        
        list_beta_theta <- update_betatheta11_cpp(logz,
                                                  beta_theta,
                                                  theta11,
                                                  delta[1:sum(M_site),],
                                                  X_w, M_site,
                                                  b_theta11 = c(1,rep(0, ncov_w)),
                                                  B_theta11 = diag(sigma_beta_theta, nrow = 1 + ncov_w),
                                                  F)
        beta_theta <- list_beta_theta$beta_theta
        theta11 <- list_beta_theta$theta11
        
      }
      
      # UPDATE THETA10 --------------------------------------------------------------
      
      if(updateTheta10){
        
        theta10 <- update_theta10_cpp(theta10, delta, gamma, M_site,
                                      a_theta0, b_theta0)
        
      }
      
      # WRITE RESULTS -----------------------------------------------------------
      
      if(iter > nburn & (iter - nburn) %% nthin == 0){
        
        trueIter <- (iter - nburn) / nthin
        
        if(paramsToSave$lambda){
          lambda_output_iter[,trueIter] <- lambda
        } 
        
        if(paramsToSave$beta0){
          beta0_output_iter[,trueIter] <- beta0
        } 
        
        if(paramsToSave$mu){
          mu_output_iter[,trueIter] <- mu
        } 
        
        if(paramsToSave$beta_z){
          beta_z_output_iter[,,trueIter] <- beta_z
        } 
        
        if(paramsToSave$logz){
          logz_output_iter[,,trueIter] <- logz
        } 
        
        if(paramsToSave$u){
          u_output_iter[,,trueIter] <- u
        } 
        
        if(paramsToSave$v){
          v_output_iter[,,trueIter] <- v
        } 
        
        if(paramsToSave$beta_w){
          beta_w_output_iter[,,trueIter] <- beta_w
        }  
        
        if(paramsToSave$Tau){
          if(jointSpecies){
            Tau_output_iter[trueIter,,] <- Tau
          } else {
            Tau_output_iter[,trueIter] <- tau
          }
        } 
        
        if(paramsToSave$sigma){
          sigma_output_iter[,trueIter] <- sigma
        } 
        
        if(paramsToSave$beta_theta){
          beta_theta_output_iter[,,trueIter] <- beta_theta
        } 
        
        if(paramsToSave$theta){
          theta_output_iter[,,trueIter] <- theta11
        } 
        
        if(paramsToSave$csi){
          csi_output_iter[,trueIter] <- theta10
        } 
        
        if(paramsToSave$delta){
          delta_output_iter[,,trueIter] <- delta
        } 
        
        if(paramsToSave$gamma){
          gamma_output_iter[,,trueIter] <- gamma
        }
        
        if(paramsToSave$r){
          r_output_iter[,trueIter] <- r_nb
        }
        
        if(paramsToSave$p11){
          p_11_output_iter[,trueIter] <- p_11
        } 
        
        if(paramsToSave$p10){
          p_10_output_iter[,trueIter] <- p_10
        } 
        
        if(paramsToSave$mutilde){
          mutilde_output_iter[trueIter] <- mu_tilde
        } 
        
        if(paramsToSave$c_imk){
          cimk_output_iter[,,,trueIter] <- c_imk
        } 
        
        if(paramsToSave$mu0){
          mu0_output_iter[trueIter] <- mu0
        } 
        
        if(paramsToSave$n0){
          n0_output_iter[trueIter] <- n0
        } 
        
        if(paramsToSave$pi0){
          pi0_output[chain,] <- pi0_output_iter
        } 
        
        if(paramsToSave$eta){
          eta_output_iter[,,,trueIter] <- lambda_ijk
        } 
        
      }
      
    }
    
    # copy chains output
    {
      if(paramsToSave$lambda){
        lambda_output[chain,,] <- lambda_output_iter
      } 
      
      if(paramsToSave$beta0){
        beta0_output[chain,,] <- beta0_output_iter
      } 
      
      if(paramsToSave$mu){
        mu_output[chain,,] <- mu_output_iter
      } 
      
      if(paramsToSave$beta_z){
        beta_z_output[chain,,,] <- beta_z_output_iter
      } 
      
      if(paramsToSave$logz){
        logz_output[chain,,,] <- logz_output_iter
      } 
      
      if(paramsToSave$u){
        u_output[chain,,,] <- u_output_iter
      } 
      
      if(paramsToSave$v){
        v_output[chain,,,] <- v_output_iter
      } 
      
      if(paramsToSave$beta_w){
        beta_w_output[chain,,,] <- beta_w_output_iter
      }  
      
      if(paramsToSave$Tau){
        if(jointSpecies) {
          Tau_output[chain,,,] <- Tau_output_iter
        } else {
          Tau_output[chain,,] <- Tau_output_iter
        }
      } 
      
      if(paramsToSave$sigma){
        sigma_output[chain,,] <- sigma_output_iter
      } 
      
      if(paramsToSave$beta_theta){
        beta_theta_output[chain,,,] <- beta_theta_output_iter
      } 
      
      if(paramsToSave$theta){
        theta_output[chain,,,] <- theta_output_iter
      } 
      
      if(paramsToSave$csi){
        csi_output[chain,,] <- csi_output_iter
      } 
      
      if(paramsToSave$delta){
        delta_output[chain,,,] <- delta_output_iter
      } 
      
      if(paramsToSave$gamma){
        gamma_output[chain,,,] <- gamma_output_iter
      }
      
      if(paramsToSave$r){
        r_output[chain,,] <- r_output_iter
      }
      
      if(paramsToSave$p11){
        p_11_output[chain,,] <- p_11_output_iter
      } 
      
      if(paramsToSave$p10){
        p_10_output[chain,,] <- p_10_output_iter
      } 
      
      if(paramsToSave$mutilde){
        mutilde_output[chain,] <- mutilde_output_iter
      } 
      
      if(paramsToSave$c_imk){
        cimk_output[chain,,,,] <- cimk_output_iter
      } 
      
      if(paramsToSave$mu0){
        mu0_output[chain,] <- mu0_output_iter
      } 
      
      if(paramsToSave$n0){
        n0_output[chain,] <- n0_output_iter
      } 
      
      if(paramsToSave$pi0){
        pi0_output[chain,] <- pi0_output_iter
      } 
      
      if(paramsToSave$eta){
        eta_output[chain,,,,] <- eta_output_iter
      } 
      
    }
    
  }
  
  params_output <- list(
    "lambda_output" = lambda_output,
    "beta0_output" = beta0_output,
    "mu_output" = mu_output,
    "beta_z_output" = beta_z_output,
    "logz_output" = logz_output,
    "v_output" = v_output,
    "u_output" = u_output,
    "beta_w_output" = beta_w_output,
    "Tau_output" = Tau_output,
    "sigma_output" = sigma_output,
    "beta_theta_output" = beta_theta_output,
    "theta_output" = theta_output,
    "csi_output" = csi_output,
    "delta_output" = delta_output,
    "gamma_output" = gamma_output,
    "p_11_output" = p_11_output,
    "p_10_output" = p_10_output,
    "r_output" = r_output,
    "mutilde_output" = mutilde_output,
    "cimk_output" = cimk_output,
    "mu0_output" = mu0_output,
    "n0_output" = n0_output,
    "pi0_output" = pi0_output,
    "eta_output" = eta_output
  )
  
  output <- list("params_output" = params_output)
}


# UTILITY --------

logistic <- function(x){
  1 / (1 + exp(-x))
}

# PLOT FUNCTION --------------

diagnosticsPlot <- function(chain_output){
  
  chain_output_all <-
    matrix(chain_output, nrow = nchain, ncol = niter)
  
  chain_output_long <- reshape2::melt(chain_output_all)
  
  ggplot2::ggplot(data = chain_output_long, ggplot2::aes(
    x = Var2,
    y = value,
    group = Var1,
    color = factor(Var1)
  )) + ggplot2::geom_line() +
    ggplot2::xlab("Iterations") + ggplot2::ylab("Value") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 17),
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 11, face = "bold"),
      axis.text.x = ggplot2::element_text(
        size = 11,
        face = "bold",
        hjust = 1
      ),
      axis.line = ggplot2::element_line(colour = "black", size = 0.15),
      # panel.grid.minor = element_line(colour="grey", size=0.15),
      panel.grid.major = ggplot2::element_line(colour = "grey", size = 0.15),
      panel.background = ggplot2::element_rect(fill = "white", color = "black"),
      legend.position = "none"
    )
  
}

# INTERWEAVING FUNCTIONS --------

convertSPtoCP <- function(lambda, beta_z, mu, logz, v,
                          delta, gamma){
  
  beta_bar <- lambda + beta_z[1,]
  mu_bar <- lambda + mu
  logz_bar <- logz
  for (j in 1:S) {
    logz_bar[,j] <- logz[,j] + beta_bar[j] - beta_z[1,j]
  }
  v_bar <- v
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        if(delta[m + sum(M_site[seq_len(i-1)]), j] == 1){
          v_bar[m + sum(M_site[seq_len(i-1)]), j] <- 
            v[m + sum(M_site[seq_len(i-1)]), j] + logz_bar[i,j] -
            logz[i,j]
        } else if (gamma[m + sum(M_site[seq_len(i-1)]), j] == 1){
          v_bar[m + sum(M_site[seq_len(i-1)]), j] <- 
            v[m + sum(M_site[seq_len(i-1)]), j] + mu_bar[j] - mu[j]
          
        } 
        
      }
    }
  }
  
  list("beta_bar" = beta_bar,
       "mu_bar" = mu_bar,
       "logz_bar" = logz_bar,
       "v_bar" = v_bar)
}

convertCPtoSP <- function(beta_bar, beta_z, lambda, mu_bar, logz_bar,
                          v_bar, delta, gamma){
  
  beta_z[1,] <- beta_bar - lambda 
  mu <- mu_bar - lambda
  for (j in 1:S) {
    logz[,j] <- logz_bar[,j] + beta_z[1,j] - beta_bar[j] 
  }
  
  v <- v_bar
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        if(delta[m + sum(M_site[seq_len(i-1)]), j] == 1){
          v[m + sum(M_site[seq_len(i-1)]), j] <- 
            v_bar[m + sum(M_site[seq_len(i-1)]), j] + logz[i,j] - logz_bar[i,j] 
        } else if (gamma[m + sum(M_site[seq_len(i-1)]), j] == 1){
          v[m + sum(M_site[seq_len(i-1)]), j] <- 
            v_bar[m + sum(M_site[seq_len(i-1)]), j] + mu[j] - mu_bar[j]
        } 
        
      }
    }
  }
  
  list("beta_z" = beta_z,
       "mu" = mu,
       "logz" = logz,
       "v" = v)
}

convertSPtoNP <- function(logz, beta_z, v, mu, 
                          delta, gamma){
  
  logz_tilde <- logz
  for (j in 1:S) {
    logz_tilde[,j] <- logz[,j] - beta_z[1,j] 
  }
  v_tilde <- v
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        if(delta[m + sum(M_site[seq_len(i-1)]), j] == 1){
          v_tilde[m + sum(M_site[seq_len(i-1)]), j] <- 
            v[m + sum(M_site[seq_len(i-1)]), j] - 
            logz[i,j] 
        } else if (gamma[m + sum(M_site[seq_len(i-1)]), j] == 1){
          v_tilde[m + sum(M_site[seq_len(i-1)]), j] <- 
            v[m + sum(M_site[seq_len(i-1)]), j] - 
            mu[j] 
        } 
      }
    }
  }
  
  list("logz_tilde" = logz_tilde,
       "v_tilde" = v_tilde)
}

convertNPtoSP <- function(beta_z, logz_tilde, v_tilde, 
                          mu, delta, gamma){
  
  logz <- logz_tilde
  for (j in 1:S) {
    logz[,j] <- logz_tilde[,j] + beta_z[1,j] 
  }
  
  v <- v_tilde
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        if(delta[m + sum(M_site[seq_len(i-1)]), j] == 1){
          v[m + sum(M_site[seq_len(i-1)]), j] <- 
            v_tilde[m + sum(M_site[seq_len(i-1)]), j] + 
            logz[i,j] 
        } else if (gamma[m + sum(M_site[seq_len(i-1)]), j] == 1){
          v[m + sum(M_site[seq_len(i-1)]), j] <- 
            v_tilde[m + sum(M_site[seq_len(i-1)]), j] + 
            mu[j] 
        } 
      }
    }
  }
  
  list("logz" = logz,
       "v" = v)
}

update_lambda_iw <- function(beta_z, mu, lambda, log_z, v, u,
                             sigma_beta, sigma_mu){
  
  S <- length(mu)
  
  # UPDATE LAMBDA CP ------
  
  # list_CP <- convertSPtoCP(lambda, beta_z, mu, logz, v, delta, gamma)
  # beta_bar <- list_CP$beta_bar
  # mu_bar <- list_CP$mu_bar
  # logz_bar <- list_CP$logz_bar
  # v_bar <- list_CP$v_bar
  
  list_CP_cpp <- convertSPtoCP_cpp(lambda, beta_z, mu, logz, v, delta, gamma, M_site)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  
  # update paramters
  
  for (j in 1:S) {
    
    sigma_inv <- 1 / sigma_mu^2 + 1 / sigma_beta^2
    
    mu_inv <- mu_bar[j] / sigma_mu^2 + beta_bar[j] / sigma_beta^2
    
    posterior_var <- 1 /  sigma_inv
    posterior_mean <- mu_inv * posterior_var
    
    lambda[j] <- rnorm(1, posterior_mean, sqrt(posterior_var))
    
  }
  
  # list_SP <- convertCPtoSP(beta_bar, beta_z, lambda, mu_bar, logz_bar,
  #                          v_bar,
  #               delta, gamma)
  # beta_z <- list_SP$beta_z
  # mu <- list_SP$mu
  # logz <- list_SP$logz
  # v <- list_SP$v
  
  list_SP <- convertCPtoSP_cpp(beta_bar, beta_z, lambda, mu_bar, logz_bar,
                               v_bar, delta, gamma, M_site)
  beta_z <- list_SP$beta_z
  mu <- list_SP$mu
  logz <- list_SP$logz
  v <- list_SP$v
  
  # UPDATE LAMBDA NP ------
  
  # list_NP <- convertSPtoNP(logz, beta_z, v, mu, 
  #                          delta, gamma)
  # logz_tilde <- list_NP$logz_tilde
  # v_tilde <- list_NP$v_tilde
  
  list_NP <- convertSPtoNP_cpp(logz, beta_z, v, mu, 
                               delta, gamma, M_site)
  logz_tilde <- list_NP$logz_tilde
  v_tilde <- list_NP$v_tilde
  
  # update parameters
  
  {
    for (j in 1:S) {
      
      sum_v <- 0
      
      availableCounts <- sum(c_imk[,,j] == 1 & !is.na(c_imk[,,j]))
      
      PCR_counts <- rep(NA, availableCounts)
      PCR_v <- rep(NA, availableCounts)
      
      l <- 1
      for (i in 1:n) {
        for (m in 1:M_site[i]) {
          for (k in 1:K[m + sum(M_site[seq_len(i-1)])]) {
            
            if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1 &
               c_imk[m + sum(M_site[seq_len(i-1)]),k,j] == 1){
              PCR_counts[l] <- y[m + sum(M_site[seq_len(i-1)]),k,j]
              PCR_v[l] <-  beta_z[1,j] +
                logz_tilde[i,j] +
                v_tilde[m + sum(M_site[seq_len(i-1)]),j] +
                u[m + sum(M_site[seq_len(i-1)]),k]
              l <- l + 1
            } else if (gamma[m + sum(M_site[seq_len(i-1)]),j] == 1 &
                       c_imk[m + sum(M_site[seq_len(i-1)]),k,j] == 1){
              PCR_counts[l] <- y[m + sum(M_site[seq_len(i-1)]),k,j]
              PCR_v[l] <- mu[j] +
                v_tilde[m + sum(M_site[seq_len(i-1)]),j] +
                u[m + sum(M_site[seq_len(i-1)]),k]
              l <- l + 1
            }
            
          }
        }
      }
      
      mean_likelihood <- log(sum(PCR_counts)) - log(sum(exp(PCR_v)))
      sigma_likelihood <- 1 / sum(PCR_counts)
      
      mean_prior <- 0
      prior_var <- 1000
      
      posterior_var <- 1 / (1 / prior_var + 1 / sigma_likelihood)
      posterior_mean <- ((mean_prior / prior_var) + (mean_likelihood / sigma_likelihood)) * posterior_var
      lambda_star <- rnorm(1, posterior_mean, sqrt(posterior_var))
      
      lambda_current <- lambda[j]
      
      logposterior <- logposterior_lambda(PCR_counts, PCR_v,
                                          lambda_current, lambda_star)
      
      if(runif(1) < logposterior){
        lambda[j] <- lambda_star
      }
      
    }
  }
  
  # list_SP <- convertNPtoSP(beta_z, logz_tilde, v_tilde, 
  #                          mu, delta, gamma)
  # logz <- list_SP$logz
  # v <- list_SP$v
  
  list_SP <- convertNPtoSP_cpp(beta_z, logz_tilde, v_tilde, 
                               mu, delta, gamma, M_site)
  logz <- list_SP$logz
  v <- list_SP$v
  
  list("lambda" = lambda,
       "mu" = mu,
       "beta_z" = beta_z,
       "logz" = logz,
       "v" = v)
  
}

logposterior_lambda <- function(PCR_counts, PCR_v,
                                x_current, x_star){
  sum(dpois(PCR_counts, exp(x_star + PCR_v), log = T)) -
    sum(dpois(PCR_counts, exp(x_current + PCR_v), log = T))
}

update_beta_iw <- function(beta_z, X_z, lambda, logz, v, u,
                           tau, sigma_beta){
  
  # UPDATE BETA CP --------------------------------------------
  
  list_CP_cpp <- convertSPtoCP_cpp(lambda, beta_z, mu, logz, v, delta, gamma, M_site)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  
  # update parameters
  
  {
    X_beta_noempty <- X_z[emptySites == 0,]
    l_noempty <- logz_bar[emptySites == 0,]
    
    ncov_z <- ncol(X_z)
    tXX <- t(X_beta_noempty) %*% X_beta_noempty
    for (j in 1:S) {
      
      Lambda_beta <- (tXX / tau[j]^2) + (diag(1, nrow = ncov_z) / sigma_beta^2)
      mu_beta <- apply(t(sapply(1:nrow(X_beta_noempty), function(i){
        X_beta_noempty[i,] * l_noempty[i,j] / tau[j]^2
      })), 2, sum) +
        c(lambda[j] / sigma_beta^2, rep(0, ncov_z - 1))
      
      beta_bar_beta <- mvrnorm(1, solve(Lambda_beta) %*% mu_beta, solve(Lambda_beta))
      
      beta_bar[j] <- beta_bar_beta[1]
      beta_z[-1,j] <- beta_bar_beta[-1]
      
    }
  }
  
  list_SP <- convertCPtoSP_cpp(beta_bar, beta_z, lambda, mu_bar, logz_bar,
                               v_bar, delta, gamma, M_site)
  beta_z <- list_SP$beta_z
  mu <- list_SP$mu
  logz <- list_SP$logz
  v <- list_SP$v
  
  # UPDATE BETA STANDARD PARAM --------------------------------------------
  
  # list_NP <- convertSPtoNP(logz, beta_z, v, mu, 
  #                          delta, gamma)
  # logz_tilde <- list_NP$logz_tilde
  # v_tilde <- list_NP$v_tilde
  
  # update parameters
  {
    
    X_beta_noempty <- X_z[emptySites == 0,]
    l_noempty <- logz[emptySites == 0,]
    
    ncov_z <- ncol(X_beta_noempty)
    tXX <- t(X_beta_noempty) %*% X_beta_noempty
    for (j in 1:S) {
      
      Lambda_beta <- (tXX / tau[j]^2) + (diag(1, nrow = ncov_z) / sigma_beta^2)
      mu_beta <- apply(t(sapply(1:nrow(X_beta_noempty), function(i){
        X_beta_noempty[i,] * l_noempty[i,j] / tau[j]^2
      })), 2, sum)
      
      beta_z[,j] <- mvrnorm(1, solve(Lambda_beta) %*% mu_beta, solve(Lambda_beta))
      
    }
    
  }
  
  list("lambda" = lambda,
       "beta_z" = beta_z,
       "logz" = logz,
       "v" = v)
  
}

logposterior_logz <- function(l, y_all, v_all, x_all, sigma_j,
                              v_delta, prior_mean, tau_j){
  
  x_all_l <- x_all * as.vector(l)
  
  sum(dnorm(v_delta, l, sigma_j)) + 
    sum(sapply(1:length(y_all), function(i){
      dbinom(y_all[i], 1, prob = 1 / (1 + exp(- x_all_l - v_all[i])), log = T)
    })) +
    dnorm(l, prior_mean, tau_j, log = T)
  
}

update_logz <- function(logz, lambda, beta_z, v, c_imk, delta,
                        mu, X_z, sigma, tau){
  
  df_t = 3
  
  # UPDATE LOGZ CP  --------------------------------------------
  
  list_CP_cpp <- convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, 
                                   gamma, beta_theta, M_site, S_star)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  beta_theta_bar <- list_CP_cpp$beta_theta_bar
  
  Xz_beta <- X_z %*% beta_z
  Xw_beta <- X_w %*% beta_w
  Xw_beta_theta <- X_w %*% t(beta_theta[,-c(1,2)] )
  
  # update parameters
  
  for (i in 1:n) {
    if(!(emptySites[i] == 1)){
      for (j in 1:S) {
        
        logz_current <- logz_bar[i,j]
        
        # data from delta
        y_all <- delta[1:M_site[i] + sum(M_site[seq_len(i-1)]),j]
        v_all <- beta_theta[j, 1] + 
          Xw_beta_theta[1:M_site[i] + sum(M_site[seq_len(i-1)]),j]
        x_all <- beta_theta[j,2] / exp(lambda[j])
        
        # data from vbar
        idxDelta1 <- delta[1:M_site[i] + sum(M_site[seq_len(i-1)]),j] == 1
        v_delta <- v_bar[1:M_site[i] + sum(M_site[seq_len(i-1)]),j][idxDelta1] -
          Xw_beta[1:M_site[i] + sum(M_site[seq_len(i-1)]), j][idxDelta1] 
        nsamples <- sum(idxDelta1)
        
        prior_mean <- beta_bar[j] + X_z[i,] %*% beta_z[,j]
        
        logf <- function(l){
          
          x_all_l <- l * x_all
          
          return(- (- 1 / (2 * sigma[j]^2) * sum((l - v_delta)^2) + 
                      sum(y_all * (v_all + x_all_l)) - 
                      sum(log(1 + exp(x_all_l + v_all))) - 
                      1 / (2 * tau[j]^2) * (l - prior_mean)^2  ))
          
        }
        
        H_f <- function(l){
          
          x_all_l <- l * x_all
          
          - 1 / (sigma[j]^2) * length(v_delta) + 
            -sum(x_all * x_all * exp(v_all + x_all_l) / (1 + exp(v_all + x_all_l))^2) +
            
            # -sum(x_all * exp(v_all + x_all_l)) + sum(y_all * x_all) -
            # sum(y_all * (v_all + x_all_l)) - sum(log(1 + exp(x_all_l + v_all))) - 
            (- 1 / (tau[j]^2)   )
          
          
        }
        
        # optim_res <- optimize(logf, c(logz_bar[i,j] - 30,
        #                               logz_bar[i,j] + 30), 
        #                       tol = exp(-7))
        # 
        # logz_star <- optim_res$minimum
        # sd_star <- sqrt(solve(-H_f(logz_star)))
        
        logz_star <- findzero_cpp(logz_bar[i,j] - 40,
                                  logz_bar[i,j] + 40,
                                  exp(-7),
                                  x_all,
                                  y_all,
                                  v_all,
                                  v_delta, tau[j],
                                  sigma[j], prior_mean)
        
        sd_star <- sqrt(solve(-h_f_cpp(logz_star,
                                       x_all,
                                       y_all,
                                       v_all,
                                       v_delta, tau[j],
                                       sigma[j], prior_mean)))
        
        
        
        # H_f(logz_star)
        # h_f_cpp(logz_star, 
        #         x_all[1], 
        #         y_all,
        #         v_all,
        #         v_delta,
        #         tau[j],
        #         sigma[j],
        #         prior_mean)
        # 
        # logfp(logz_bar[i,j])
        # 
        # logf_cpp(logz_bar[i,j],
        #          x_all[1], sigma[j], v_delta, v_all, y_all, tau[j], prior_mean)
        # 
        # findzero_cpp(logz_bar[i,j] - 4,
        #              logz_bar[i,j] + 4,
        #              tol = .01,
        #              x_all[1], y_all, v_all, v_delta, tau[j],
        #              sigma[j], prior_mean)
        
        logz_new <- rnorm(1, logz_star, sd_star)  
        
        loglikelihood <- logposterior_logz_cpp(v_delta, sigma[j],
                                               logz_current, logz_new,
                                               prior_mean, tau[j]^2)
        
        logposterior_ratio_logistic <- logposterior_logz_logistic_cpp(y_all,
                                                                      x_all,
                                                                      v_all,
                                                                      logz_current,
                                                                      logz_new);
        
        logproposal_ratio = dnorm(logz_current, logz_star, sd_star, 1) - 
          dnorm(logz_new, logz_star, sd_star, 1)
        
        logposterior <- loglikelihood + logposterior_ratio_logistic + logproposal_ratio
        
        logz_grid <- seq(0, 10, length.out = 100)
        true_logz_values <- sapply(1:length(logz_grid), function(l){
          
          # logf_cpp(logz_grid[l], x_all, sigma[j],
          #          v_delta,
          #          v_all, y_all, tau[j], prior_mean)
          
          logposterior_logz_cpp(v_delta, sigma[j],
                                10, logz_grid[l],
                                prior_mean, tau[j] * tau[j]) + 
            logposterior_logz_logistic_cpp(y_all,
                                           x_all,
                                           v_all,
                                           10,
                                           logz_grid[l])
          
        })
        
        # density_grid <- seq(5, 7, length.out = 100)
        explogz_values <- sapply(1:length(logz_grid), function(l){
          dnorm(logz_grid[l], logz_star, sd_star, log = F)
        })
        
        qplot(logz_grid, exp(true_logz_values - max(true_logz_values))) #+
        # qplot(logz_grid, true_logz_values) +
        # geom_vline(aes(xintercept = logz_star)) +
        geom_vline(aes(xintercept = logz_new)) +
          geom_vline(aes(xintercept = logz_current)) +
          geom_point(data = NULL, aes(x = logz_grid, y = explogz_values / max(explogz_values)), color = "red")
        # qplot(logz_grid, logz_values)
        
        if(runif(1) < exp(logposterior)){
          logz_bar[i,j] <- logz_new
        }
        
      }
    }
  }
  
  list_SP <- convertCPtoSP_cpp(beta_bar, beta_z, lambda, mu_bar, logz_bar,
                               v_bar, delta, gamma, M_site)
  beta_z <- list_SP$beta_z
  mu <- list_SP$mu
  logz <- list_SP$logz
  v <- list_SP$v
  
  list("lambda" = lambda,
       "mu" = mu,
       "beta_z" = beta_z,
       "logz" = logz,
       "v" = v)
}

logposterior_logz <- function(PCR_counts, PCR_v,
                              logz_current, logz_star,
                              prior_mean, prior_var){
  sum(dpois(PCR_counts, exp(logz_star + PCR_v), log = T)) -
    sum(dpois(PCR_counts, exp(logz_current + PCR_v), log = T)) + 
    dnorm(logz_star, prior_mean, sqrt(prior_var), log = T) + 
    dnorm(logz_current, prior_mean, sqrt(prior_var), log = T)
}

update_mu_iw <- function(mu, logz, lambda, beta_z, v, c_imk, delta,
                         X_z, sigma_gamma, sigma_mu){
  
  # UPDATE MU CP --------------------------------------------
  
  list_CP_cpp <- convertSPtoCP_cpp(lambda, beta_z, mu, logz, v, delta, gamma, M_site)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  
  # update parameters
  
  for (j in 1:S) {
    
    samples_vbar <- sum(gamma[,j] == 1 & !is.na(gamma[,j]))
    
    v_bar_mu <- v_bar[gamma[,j] == 1 & !is.na(gamma[,j]),j]
    
    if(samples_vbar > 0){
      lik_mean <- sum(v_bar_mu)
      lik_var <- sigma_gamma[j]^2  
    } else {
      lik_mean <- 0
      lik_var <- Inf
    }
    
    prior_mean <- lambda[j]
    prior_var <- tau[j]^2
    
    posterior_var <- 1 / (1 / prior_var + samples_vbar / lik_var)
    posterior_mean <- ((prior_mean / prior_var) + (lik_mean / lik_var)) * posterior_var
    mu_bar[j] <- rnorm(1, posterior_mean, sqrt(posterior_var))
  }
  
  list_SP <- convertCPtoSP_cpp(beta_bar, beta_z, lambda, mu_bar, logz_bar,
                               v_bar, delta, gamma, M_site)
  beta_z <- list_SP$beta_z
  mu <- list_SP$mu
  logz <- list_SP$logz
  v <- list_SP$v
  
  # UPDATE MU NP --------------------------------------------
  
  list_NP <- convertSPtoNP_cpp(logz, beta_z, v, mu, 
                               delta, gamma, M_site)
  logz_tilde <- list_NP$logz_tilde
  v_tilde <- list_NP$v_tilde
  
  # update parameters
  
  for (j in 1:S) {
    
    availableCounts <- 0#sum(c_imk[,,j] == 1 & !is.na(c_imk[,,j]))
    
    PCR_counts <- rep(NA, max(M_site) * max(K))
    PCR_v <- rep(NA, max(M_site) * max(K))
    
    l <- 1
    
    for (i in 1:n) {
      if(!(emptySites[i] == 1)){
        for (m in 1:M_site[i]) {
          
          if(gamma[m + sum(M_site[seq_len(i-1)]),j] == 1){
            
            for (k in 1:K[m + sum(M_site[seq_len(i-1)])]) {
              
              if(c_imk[m + sum(M_site[seq_len(i-1)]),k,j] == 1){
                
                PCR_counts[l] <- y[m + sum(M_site[seq_len(i-1)]),k,j]
                PCR_v[l] <-  lambda[j] +
                  v_tilde[m + sum(M_site[seq_len(i-1)]),j] + 
                  u[m + sum(M_site[seq_len(i-1)]),k]
                l <- l + 1
                
              }
              
            }
            
          }
          
        }
      }
    }
    
    PCR_counts <- PCR_counts[seq_len(l-1)]
    PCR_v <- PCR_v[seq_len(l-1)]
    
    if(l > 1){
      mean_likelihood <- log(sum(PCR_counts)) - log(sum(exp(PCR_v)))
      sigma_likelihood <- 1 / sum(PCR_counts)
    } else {
      mean_likelihood <- 0
      sigma_likelihood <- Inf
    }
    
    mean_prior <- 0
    prior_var <- sigma_mu^2
    
    posterior_var <- 1 / (1 / prior_var + 1 / sigma_likelihood)
    posterior_mean <- ((mean_prior / prior_var) + (mean_likelihood / sigma_likelihood)) * posterior_var
    mu_star <- rnorm(1, posterior_mean, sqrt(posterior_var))
    
    mu_current <- mu[j]
    
    logposterior <- logposterior_mu_iw(PCR_counts, PCR_v,
                                       mu_current, mu_star,
                                       mean_prior, prior_var)
    
    if(runif(1) < logposterior){
      mu[j] <- mu_star  
    }
    
  }
  
  list_SP <- convertNPtoSP_cpp(beta_z, logz_tilde, v_tilde, 
                               mu, delta, gamma, M_site)
  logz <- list_SP$logz
  v <- list_SP$v
  
  list("lambda" = lambda,
       "mu" = mu,
       "beta_z" = beta_z,
       "logz" = logz,
       "v" = v)
}

logposterior_mu_iw <- function(PCR_counts, PCR_v,
                               mu_current, mu_star,
                               prior_mean, prior_var){
  sum(dpois(PCR_counts, exp(mu_star + PCR_v), log = T)) -
    sum(dpois(PCR_counts, exp(mu_current + PCR_v), log = T)) + 
    dnorm(mu_star, prior_mean, sqrt(prior_var), log = T) + 
    dnorm(mu_current, prior_mean, sqrt(prior_var), log = T)
}

update_v_iw <- function(v, logz, lambda, u, beta_z, c_imk, delta,
                        mu, sigma, sigma_gamma){
  
  # UPDATE V CP --------------------------------------------
  
  list_CP_cpp <- convertSPtoCP_cpp(lambda, beta_z, mu, logz, v, delta, gamma, M_site)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  
  # update parameters
  
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      if(!(emptySites[i] == 1)){
        for (j in 1:S) {
          if(delta[m + sum(M_site[seq_len(i-1)]), j] == 1 | 
             gamma[m + sum(M_site[seq_len(i-1)]), j] == 1){
            
            currentK <- K[m + sum(M_site[seq_len(i-1)])]
            
            idxReplicatePositive <- 
              which(as.vector(c_imk[m + sum(M_site[seq_len(i-1)]),1:currentK,j]) == 1)
            
            PCR_counts <- 
              as.vector(y[m + sum(M_site[seq_len(i-1)]),1:currentK,j])[idxReplicatePositive]
            
            PCR_v <- u[m + sum(M_site[seq_len(i-1)]),idxReplicatePositive]
            
            if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1){
              
              prior_mean <- logz_bar[i,j]
              prior_var <- sigma[j]^2
              
            } else {
              
              prior_mean <- mu_bar[j]
              prior_var <- sigma_gamma[j]^2
              
            }
            
            if(length(idxReplicatePositive) > 0){
              
              mean_likelihood <- log(sum(PCR_counts)) - log(sum(exp(PCR_v)))
              sigma_likelihood <- 1 / sum(PCR_counts)
              
            } else {
              
              mean_likelihood <- 0
              sigma_likelihood <- Inf
              
            }
            
            posterior_var <- 1 / (1 / prior_var + 1 / sigma_likelihood)
            posterior_mean <- ((prior_mean / prior_var) + (mean_likelihood / sigma_likelihood)) * posterior_var
            
            v_current <- v_bar[m + sum(M_site[seq_len(i-1)]),j]
            v_star <- rnorm(1, posterior_mean, sqrt(posterior_var))
            
            logposterior_ratio <- logposterior_v_iw(PCR_counts, PCR_v,
                                                    v_current, v_star,
                                                    prior_mean, prior_var)
            
            if(runif(1) < exp(logposterior_ratio)){
              v_bar[m + sum(M_site[seq_len(i-1)]),j] <- v_star
            }
            
          } 
        }
      }  
    }
    
  }
  
  list_SP <- convertCPtoSP_cpp(beta_bar, beta_z, lambda, mu_bar, logz_bar,
                               v_bar, delta, gamma, M_site)
  beta_z <- list_SP$beta_z
  mu <- list_SP$mu
  logz <- list_SP$logz
  v <- list_SP$v
  
  # UPDATE V NP --------------------------------------------
  
  list_NP <- convertSPtoNP_cpp(logz, beta_z, v, mu, 
                               delta, gamma, M_site)
  logz_tilde <- list_NP$logz_tilde
  v_tilde <- list_NP$v_tilde
  
  # update parameters
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      if(!(emptySites[i] == 1)){
        for (j in 1:S) {
          if(delta[m + sum(M_site[seq_len(i-1)]), j] == 1 | 
             gamma[m + sum(M_site[seq_len(i-1)]), j] == 1){
            
            currentK <- K[m + sum(M_site[seq_len(i-1)])]
            
            idxReplicatePositive <- 
              which(as.vector(c_imk[m + sum(M_site[seq_len(i-1)]),1:currentK,j]) == 1)
            
            PCR_counts <- 
              as.vector(y[m + sum(M_site[seq_len(i-1)]),1:currentK,j])[idxReplicatePositive]
            
            if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1){
              
              PCR_v <- lambda[j] + 
                beta_z[1,j] + 
                logz_tilde[i,j] + 
                u[m + sum(M_site[seq_len(i-1)]),idxReplicatePositive]
              
              mean_prior <- 0
              prior_var <- sigma[j]^2
              
            } else {
              
              PCR_v <- lambda[j] + 
                mu[j] + 
                u[m + sum(M_site[seq_len(i-1)]),idxReplicatePositive]
              
              mean_prior <- 0
              prior_var <- sigma_gamma[j]^2
              
            }
            
            if(length(idxReplicatePositive) > 0){
              
              mean_likelihood <- log(sum(PCR_counts)) - log(sum(exp(PCR_v)))
              sigma_likelihood <- 1 / sum(PCR_counts)
              
            } else {
              
              mean_likelihood <- 0
              sigma_likelihood <- Inf
              
            }
            
            posterior_var <- 1 / (1 / prior_var + 1 / sigma_likelihood)
            posterior_mean <- ((mean_prior / prior_var) + (mean_likelihood / sigma_likelihood)) * posterior_var
            
            v_current <- v_tilde[m + sum(M_site[seq_len(i-1)]),j]
            v_star <- rnorm(1, posterior_mean, sqrt(posterior_var))
            
            logposterior_ratio <- logposterior_v_iw(PCR_counts, PCR_v,
                                                    v_current, v_star,
                                                    mean_prior, prior_var)
            if(runif(1) < exp(logposterior_ratio)){
              v_tilde[m + sum(M_site[seq_len(i-1)]),j] <- v_star
            }
            
          }
        }
      }  
    }
    
  }
  
  list_SP <- convertNPtoSP_cpp(beta_z, logz_tilde, v_tilde, 
                               mu, delta, gamma, M_site)
  logz <- list_SP$logz
  v <- list_SP$v
  
  list("lambda" = lambda,
       "mu" = mu,
       "beta_z" = beta_z,
       "logz" = logz,
       "v" = v)
}

update_v_poisgamma_r <- function(v, logz, lambda, u, beta_z, c_imk, delta,
                                 mu, sigma, sigma_gamma){
  
  Xw_beta <- X_w %*% beta_w
  
  df_t <- 3
  
  # UPDATE V CP --------------------------------------------
  
  list_CP_cpp <- convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, gamma, beta_theta, M_site, S_star)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  
  # update parameters
  
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      if(!(emptySites[i] == 1)){
        for (j in 1:S) {
          if(delta[m + sum(M_site[seq_len(i-1)]), j] == 1 | 
             gamma[m + sum(M_site[seq_len(i-1)]), j] == 1){
            
            v_current <- v_bar[m + sum(M_site[seq_len(i-1)]), j]
            
            a <- 0
            b <- 0
            rnb_current <- r_nb[j]
            
            lambdas <- lambda_ijk[m + sum(M_site[seq_len(i-1)]),
                                  c_imk[m + sum(M_site[seq_len(i-1)]),,j] == 1,j]
            u_present <- u[m + sum(M_site[seq_len(i-1)]),c_imk[m + sum(M_site[seq_len(i-1)]),,j] == 1]
            
            
            a <- sum(lambdas * rnb_current / exp(u_present))
            b <- length(u_present) * rnb_current
            
            if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1){
              
              prior_mean = logz_bar[i,j] + Xw_beta[m + sum(M_site[seq_len(i-1)]),j]
              prior_var = sigma[j] * sigma[j]
              
            } else {
              
              prior_mean = mu_bar[j];
              prior_var = sigma_gamma * sigma_gamma;
              
            }
            
            mu_v = prior_mean
            s_v = 1 / prior_var
            c_v = mu_v * s_v - b
            
            log_lam_argument = log(a)  - (c_v / s_v) - log(s_v);
            if(log_lam_argument > 40){
              v_star = c_v / s_v + log_lam_argument - log(log_lam_argument);
            } else {
              v_star = c_v / s_v + lambertW0_CS(exp(log_lam_argument));
            }
            
            var_star = - 1 / (- a * exp(- v_star) - s_v)
            
            v_new = rt2(v_star, sqrt(var_star), 3)
            
            logpost_new = logpost_v(v_new, lambdas, u_present, r_nb[j], 
                                    prior_mean, prior_var);
            logpost_current = logpost_v(v_current,
                                        lambdas, u_present, r_nb[j], 
                                        prior_mean, prior_var);
            
            log_posterior = logpost_new - logpost_current;
            
            log_proposal = log(dt2(v_current, v_star, sqrt(var_star), df_t)) - 
              log(dt2(v_new, v_star, sqrt(var_star), df_t));
            
            log_posterior - log_proposal
            
          } 
        }
      }  
    }
    
  }
  
  list_SP <- convertCPtoSP_cpp(beta_bar, beta_z, lambda, mu_bar, logz_bar,
                               v_bar, delta, gamma, M_site)
  beta_z <- list_SP$beta_z
  mu <- list_SP$mu
  logz <- list_SP$logz
  v <- list_SP$v
  
  # UPDATE V NP --------------------------------------------
  
  list_NP <- convertSPtoNP_cpp(logz, beta_z, v, mu, 
                               delta, gamma, M_site)
  logz_tilde <- list_NP$logz_tilde
  v_tilde <- list_NP$v_tilde
  
  # update parameters
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      if(!(emptySites[i] == 1)){
        for (j in 1:S) {
          if(delta[m + sum(M_site[seq_len(i-1)]), j] == 1 | 
             gamma[m + sum(M_site[seq_len(i-1)]), j] == 1){
            
            currentK <- K[m + sum(M_site[seq_len(i-1)])]
            
            idxReplicatePositive <- 
              which(as.vector(c_imk[m + sum(M_site[seq_len(i-1)]),1:currentK,j]) == 1)
            
            PCR_counts <- 
              as.vector(y[m + sum(M_site[seq_len(i-1)]),1:currentK,j])[idxReplicatePositive]
            
            if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1){
              
              PCR_v <- lambda[j] + 
                beta_z[1,j] + 
                logz_tilde[i,j] + 
                u[m + sum(M_site[seq_len(i-1)]),idxReplicatePositive]
              
              mean_prior <- 0
              prior_var <- sigma[j]^2
              
            } else {
              
              PCR_v <- lambda[j] + 
                mu[j] + 
                u[m + sum(M_site[seq_len(i-1)]),idxReplicatePositive]
              
              mean_prior <- 0
              prior_var <- sigma_gamma[j]^2
              
            }
            
            if(length(idxReplicatePositive) > 0){
              
              mean_likelihood <- log(sum(PCR_counts)) - log(sum(exp(PCR_v)))
              sigma_likelihood <- 1 / sum(PCR_counts)
              
            } else {
              
              mean_likelihood <- 0
              sigma_likelihood <- Inf
              
            }
            
            posterior_var <- 1 / (1 / prior_var + 1 / sigma_likelihood)
            posterior_mean <- ((mean_prior / prior_var) + (mean_likelihood / sigma_likelihood)) * posterior_var
            
            v_current <- v_tilde[m + sum(M_site[seq_len(i-1)]),j]
            v_star <- rnorm(1, posterior_mean, sqrt(posterior_var))
            
            logposterior_ratio <- logposterior_v_iw(PCR_counts, PCR_v,
                                                    v_current, v_star,
                                                    mean_prior, prior_var)
            if(runif(1) < exp(logposterior_ratio)){
              v_tilde[m + sum(M_site[seq_len(i-1)]),j] <- v_star
            }
            
          }
        }
      }  
    }
    
  }
  
  list_SP <- convertNPtoSP_cpp(beta_z, logz_tilde, v_tilde, 
                               mu, delta, gamma, M_site)
  logz <- list_SP$logz
  v <- list_SP$v
  
  list("lambda" = lambda,
       "mu" = mu,
       "beta_z" = beta_z,
       "logz" = logz,
       "v" = v)
}

loglik_other <- function(u, v, r, y_i){
  sum(y_i * (u + v)) - sum((y_i + r) * log(exp(u + v) + r))
}

d2loglik <- function(u, y, v, r){
  - r * sum((y + r) * exp(u + v) / (r + exp(u + v))^2 )
}

update_v_nb <- function(v, logz, lambda, r_nb, u, beta_z, c_imk, delta,
                        mu, sigma, sigma_gamma){
  
  # UPDATE V CP --------------------------------------------
  
  list_CP_cpp <- convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, 
                                   delta, gamma, beta_theta, M_site)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  beta_theta0_bar <- list_CP_cpp$beta_theta0_bar
  
  # update parameters
  
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      if(!(emptySites[i] == 1)){
        for (j in 1:S) {
          if(delta[m + sum(M_site[seq_len(i-1)]), j] == 1 | 
             gamma[m + sum(M_site[seq_len(i-1)]), j] == 1){
            
            currentK <- K[m + sum(M_site[seq_len(i-1)])]
            
            idxReplicatePositive <- 
              which(as.vector(c_imk[m + sum(M_site[seq_len(i-1)]),1:currentK,j]) == 1)
            
            PCR_counts <- 
              as.vector(y[m + sum(M_site[seq_len(i-1)]),1:currentK,j])[idxReplicatePositive]
            
            # PCR_v <- u[m + sum(M_site[seq_len(i-1)]),idxReplicatePositive]
            
            if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1){
              
              prior_mean <- logz_bar[i,j]
              prior_var <- sigma[j]^2
              
            } else {
              
              prior_mean <- mu_bar[j]
              prior_var <- sigma_gamma^2
              
            }
            
            if(length(idxReplicatePositive) > 0 & sum(PCR_counts) > 0){
              
              # v_present <- v_bar[m + sum(M_site[seq_len(i-1)]), j]
              u_present <- u[m + sum(M_site[seq_len(i-1)]), idxReplicatePositive]
              
              # a <- sum(PCR_counts)
              # b <- sum((PCR_counts + r_nb[j]) * exp(u_present))
              # c_ab <- a / b
              # v_star <- log( c_ab * r_nb[j] / (1 - c_ab * mean(exp(u_present))) )  
              
              # approx_new
              
              a <- sum(PCR_counts)
              b_i <- (PCR_counts + r_nb[j]) * exp(u_present)
              c_i <- exp(u_present)
              b <- sum(b_i / c_i)
              c_ab <- a / b
              r_bar <- mean(r_nb[j] / c_i)
              v_star <- log( c_ab * r_bar / (1 - c_ab ) )
              
              #
              
              var_star <- - 1 / d2loglik(v_star, PCR_counts, u_present, r_nb[j])
              # d2loglik_cpp(v_star, PCR_counts, u_present, r_nb[j])
              mean_likelihood <- v_star
              sigma_likelihood <- var_star
              
              posterior_var <- 1 / (1 / prior_var + 1 / var_star)
              posterior_mean <- ((prior_mean / prior_var) + (mean_likelihood / sigma_likelihood)) * posterior_var
              
              v_current <- v_bar[m + sum(M_site[seq_len(i-1)]), j]
              v_new <- rnorm(1, posterior_mean, sqrt(posterior_var))
              
              # logposterior_ratio <- logposterior_u(PCR_counts, PCR_v,
              #                                      u_current, u_star,
              #                                      mean_prior, prior_var)
              
              # u_grid <- seq(4 - 2, 4 + 2, length.out = 40)
              # 
              # betafun_grid <- sapply(seq_along(u_grid), function(i){
              #   loglik_other(u_grid[i], u_present, r_nb[j], PCR_counts)
              # })
              # qplot(u_grid, exp(betafun_grid - max(betafun_grid)))
              
              loglik_new <- loglik_other(v_new, u_present, r_nb[j], PCR_counts)
              loglik_current <- loglik_other(v_bar[m + sum(M_site[seq_len(i-1)]),j], 
                                             u_present, r_nb[j], PCR_counts)
              
              logprior_new <- dnorm(v_new, 
                                    prior_mean, sqrt(prior_var), log = T)
              logprior_current <- dnorm(v_bar[m + sum(M_site[seq_len(i-1)]), j], 
                                        prior_mean, sqrt(prior_var), log = T)
              
              logposterior_new <- loglik_new + logprior_new
              logposterior_current <- loglik_current + logprior_current
              exp(logposterior_new - logposterior_current)
              # 
              if(runif(1) < exp(logposterior_new - logposterior_current)){
                v_bar[m + sum(M_site[seq_len(i-1)]), j] <- v_new
              }
              
            } else {
              
              v_bar[m + sum(M_site[seq_len(i-1)]), j] <- rnorm(1, prior_mean, sqrt(prior_var))
              
            }
            
            
          } 
        }
      }  
    }
    
  }
  
  list_SP <- convertCPtoSP_cpp(beta_bar, lambda, mu_bar, logz_bar,
                               v_bar, delta, gamma, beta_theta0_bar, 
                               beta_theta,M_site)
  beta_z <- list_SP$beta_z
  mu <- list_SP$mu
  logz <- list_SP$logz
  v <- list_SP$v
  
  v
}

logposterior_v_iw <- function(PCR_counts, PCR_v,
                              v_current, v_star,
                              prior_mean, prior_var){
  sum(dpois(PCR_counts, exp(v_star + PCR_v), log = T)) -
    sum(dpois(PCR_counts, exp(v_current + PCR_v), log = T)) + 
    dnorm(v_star, prior_mean, sqrt(prior_var), log = T) + 
    dnorm(v_current, prior_mean, sqrt(prior_var), log = T)
}


update_u_iw <- function(u, zeta, v, logz, lambda, beta_z, c_imk, delta,
                        mu, sigma_u, sigma, sigma_gamma){
  
  # UPDATE U PX --------------------------------------------
  
  # define parameters in CP
  
  list_CP_cpp <- convertSPtoCP_cpp(lambda, beta_z, mu, logz, v, delta, gamma, M_site)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  
  # define parameters in PX
  u_hat <- u + zeta
  v_hat <- v_bar - zeta
  
  # update zeta
  
  sum_sigma <- 0
  sum_mu <- 0
  
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (k in 1:K[m + sum(M_site[seq_len(i-1)])]) {
        sum_sigma <- sum_sigma + (1 / sigma_u^2)
        sum_mu <- sum_mu + u_hat[m + sum(M_site[seq_len(i-1)]),k] / sigma_u^2
      }
    }
  }
  
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1){
          # sum_mu <- sum_mu + (- v_hat[m + sum(M_site[seq_len(i-1)]),j]) / sigma[j]^2
          sum_mu <- sum_mu + (logz_bar[i,j] - v_hat[m + sum(M_site[seq_len(i-1)]),j]) / sigma[j]^2
          sum_sigma <- sum_sigma + (1 / sigma[j]^2)
        } else if (gamma[m + sum(M_site[seq_len(i-1)]),j] == 1){
          # sum_mu <- sum_mu + (- v_hat[m + sum(M_site[seq_len(i-1)]),j]) / sigma_gamma[j]^2
          sum_mu <- sum_mu + (mu_bar[j] - v_hat[m + sum(M_site[seq_len(i-1)]),j]) / sigma_gamma[j]^2
          sum_sigma <- sum_sigma + (1 / sigma_gamma[j]^2)
        }
        
      }
    }
  }
  
  posterior_var <- 1 / sum_sigma
  posterior_mean <- sum_mu * posterior_var
  
  (zeta <- rnorm(1, posterior_mean, sqrt(posterior_var)))
  
  # reupdate parameters
  u <- u_hat - zeta
  v_bar <- v_hat + zeta
  
  # UPDATE U CP --------------------------------------------
  
  # update parameters
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (k in 1:K[m + sum(M_site[seq_len(i-1)])]) {
        
        species_present <- which(c_imk[m + sum(M_site[seq_len(i-1)]),k,] == 1)
        PCR_counts <- y[m + sum(M_site[seq_len(i-1)]),k,species_present]
        
        PCR_v <- v_bar[m + sum(M_site[seq_len(i-1)]),species_present]
        
        mean_prior <- 0
        prior_var <- sigma_u^2
        
        if(length(species_present) > 0){
          
          mean_likelihood <- log(sum(PCR_counts)) - log(sum(exp(PCR_v)))
          sigma_likelihood <- 1 / sum(PCR_counts)
          
        } else {
          
          mean_likelihood <- 0
          sigma_likelihood <- Inf
          
        }
        
        posterior_var <- 1 / (1 / prior_var + 1 / sigma_likelihood)
        posterior_mean <- ((mean_prior / prior_var) + (mean_likelihood / sigma_likelihood)) * posterior_var
        
        u_current <- u[m + sum(M_site[seq_len(i-1)]),k]
        u_star <- rnorm(1, posterior_mean, sqrt(posterior_var))
        
        logposterior_ratio <- logposterior_u(PCR_counts, PCR_v,
                                             u_current, u_star,
                                             mean_prior, prior_var)
        
        if(runif(1) < exp(logposterior_ratio)){
          u[m + sum(M_site[seq_len(i-1)]),k] <- u_star
        }
        
      }
    }
  }
  
  # redefine parameters from the CP
  {
    for (i in 1:n) {
      for (m in 1:M_site[i]) {
        for (j in 1:S) {
          if(delta[m + sum(M_site[seq_len(i-1)]), j] == 1){
            v[m + sum(M_site[seq_len(i-1)]), j] <-
              v_bar[m + sum(M_site[seq_len(i-1)]), j] - logz_bar[i,j] +
              logz[i,j]
          } else if (gamma[m + sum(M_site[seq_len(i-1)]), j] == 1){
            v[m + sum(M_site[seq_len(i-1)]), j] <- v_bar[m + sum(M_site[seq_len(i-1)]), j] -
              mu_bar[j] + mu[j]
          }
          
        }
      }
    }
    
  }
  
  list("lambda" = lambda,
       "mu" = mu,
       "beta_z" = beta_z,
       "logz" = logz,
       "v" = v,
       "u" = u,
       "zeta" = zeta)
}


logposterior_u <- function(PCR_counts, PCR_v,
                           u_current, u_star,
                           prior_mean, prior_var){
  sum(dpois(PCR_counts, exp(u_star + PCR_v), log = T)) -
    sum(dpois(PCR_counts, exp(u_current + PCR_v), log = T)) + 
    dnorm(u_star, prior_mean, sqrt(prior_var), log = T) + 
    dnorm(u_current, prior_mean, sqrt(prior_var), log = T)
}

loglikelihood_zeta_z <- function(zeta_zj, logz_tilde, beta_theta1_tilde){
  
  sum <- 0
  
  Xbeta <- X_z %*% rbind(beta0_true, beta_z_true)
  
  loglik_logztilde <- sum(sapply(1:n, function(i){
    dnorm(logz_tilde[i,j], zeta_zj * Xbeta[i,j], tau[j] * zeta_zj, log = T)
  }))
  
  loglik_betatheta <- dnorm(beta_theta1_tilde[j], 0, sigma_beta / zeta_zj, log = T)
  
  loglik_logztilde + loglik_betatheta
  
}

update_zeta_z <- function(zeta_z, logz, beta_theta,
                          sd_zetaj = .0005){
  
  list_CP_cpp <- convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, 
                                   delta, gamma, beta_theta, M_site)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  beta_theta0_bar <- list_CP_cpp$beta_theta0_bar
  
  zeta_z <- rep(1, S)
  
  logz_tilde <- sapply(1:S, function(j){
    logz_bar[,j] * zeta_z[j]
  })
  
  beta_theta1_tilde <- sapply(1:S, function(j){
    beta_theta[j,2] / zeta_z[j]
  })
  
  for (j in 1:S) {
    
    zeta_zj_current <- zeta_z[j]
    zeta_zj_star <- rnorm(1, zeta_z[j], sd_zetaj)
    
    loglik_zetastar<- loglikelihood_zeta_z(zeta_zj_star, logz_tilde, beta_theta1_tilde)
    loglik_zetacurrent <- loglikelihood_zeta_z(zeta_zj_current, logz_tilde, beta_theta1_tilde)
    
    exp(loglik_zetastar - loglik_zetacurrent)
    
    if(runif(1) < exp(loglik_zetastar - loglik_zetacurrent)){
      zeta_z[j] <- zeta_zj_star
    }
    
  }
  
  logz_bar <- sapply(1:S, function(j){
    logz_tilde[,j] / zeta_z[j]
  })
  
  beta_theta[,2] <- sapply(1:S, function(j){
    beta_theta1_tilde[j] / zeta_z[j]
  })
  
  list_SP <- convertCPtoSP_cpp(beta_bar, lambda, mu_bar, logz_bar, v_bar, delta, 
                               gamma, beta_theta0_bar, beta_theta, M_site)
  beta0 <- list_SP$beta0
  mu <- list_SP$mu
  logz <- list_SP$logz
  v <- list_SP$v
  beta_theta <- list_SP$beta_theta
  
  list("beta_theta" = beta_theta,
       "zeta_z" = zeta_z,
       "logz" = logz)
  
}

loglik_other_u <- function(u, v, r, y_i){
  - (sum(y_i * (u + v)) - sum((y_i + r) * log(exp(u + v) + r)))
}

d2loglik_u <- function(u, y, v, r){
  - sum(r * (y + r) * exp(u + v) / (r + exp(u + v))^2 )
}

update_u_nb <- function(u){
  
  list_CP_cpp <- convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, 
                                   delta, gamma, beta_theta, M_site)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  beta_theta0_bar <- list_CP_cpp$beta_theta0_bar
  
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (k in 1:K[m + sum(M_site[seq_len(i-1)])]) {
        
        speciesPresent <- which(c_imk[m + sum(M_site[seq_len(i-1)]), k,] == 1)
        
        if(length(speciesPresent) > 0){
          
          y_present <- y[m + sum(M_site[seq_len(i-1)]), k, speciesPresent]
          
          v_present <- v_bar[m + sum(M_site[seq_len(i-1)]), speciesPresent]
          # u_present <- u[m + sum(M_site[seq_len(i-1)]), k]
          r_present <- r_nb[speciesPresent]
          
          # a <- sum(PCR_counts)
          # b <- sum((PCR_counts + r_nb[j]) * exp(u_present))
          # c_ab <- a / b
          # v_star <- log( c_ab * r_nb[j] / (1 - c_ab * mean(exp(u_present))) )  
          
          # approx_new
          
          # a <- sum(y_present)
          # b_i <- (y_present + r_present) 
          # c_i <- exp(v_present)
          # b <- sum(b_i)
          # c_ab <- a / b
          # r_bar <- mean(r_present / c_i)
          # u_star <- log( c_ab * r_bar / (1 - c_ab ) )
          # 
          # u_grid <- seq(- 2,  6, length.out = 40)
          # 
          # betafun_grid <- sapply(seq_along(u_grid), function(i){
          #   loglik_other(u_grid[i], v_present, r_present, y_present)
          # })
          # qplot(u_grid, betafun_grid)
          
          optim_u_rcpp(u[m + sum(M_site[seq_len(i-1)]), k],
                       y_present,
                       v_present,
                       r_present)
          
          prior_mean <- 0
          prior_var <- 1
          
          u_star <- optimize(loglik_other_u, interval = c(u[m + sum(M_site[seq_len(i-1)]), k] - 2,
                                                          u[m + sum(M_site[seq_len(i-1)]), k] + 2),
                             v = v_present,
                             r = r_present,
                             y_i = y_present)$minimum
          
          #
          
          var_star <- - 1 / d2loglik_u(u_star, y_present, v_present, r_present)
          # d2loglik_cpp(v_star, PCR_counts, u_present, r_nb[j])
          mean_likelihood <- u_star
          sigma_likelihood <- var_star
          
          posterior_var <- 1 / (1 / prior_var + 1 / var_star)
          posterior_mean <- ((prior_mean / prior_var) + (mean_likelihood / sigma_likelihood)) * posterior_var
          
          u_current <- u[m + sum(M_site[seq_len(i-1)]), k]
          u_new <- rnorm(1, posterior_mean, sqrt(posterior_var))
          
          logposterior_unew <- loglik_other_u(u_new, v_present, r_present, y_present) + 
            dnorm(u_new, prior_mean, sqrt(prior_var), 1)
          
          logposterior_current <- loglik_other_u(u_current, v_present, r_present, y_present) + 
            dnorm(u_current, prior_mean, sqrt(prior_var), 1)
          
          if(runif(1) < exp(logposterior_unew - logposterior_current)){
            
            u[m + sum(M_site[seq_len(i-1)]), k] <- u_new
            
          }
          
        } 
        
      }
    }
  }
  
  u
}

update_r_nb <- function(r_nb){
  
  list_CP_cpp <- convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, 
                                   delta, gamma, beta_theta, M_site)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  beta_theta0_bar <- list_CP_cpp$beta_theta0_bar
  
  for (j in 1:S) {
    
    samples_present <- which(c_imk[, ,j] == 1, arr.ind = T)
    
    y_present <- mapply(function(k1, k2){ y[k1, k2, j]}, 
                        samples_present[,1], 
                        samples_present[,2])
    
    # samples_Present <- unique(samples_present[,1])
    # PCR_present <- unique(samples_present[,2])
    
    uv <- v_bar[samples_present[,1], j] + 
      mapply(function(k1, k2) u[k1, k2], 
             samples_present[,1], 
             samples_present[,2])
    mean_uv <- exp(uv)
    
    r_grid <- seq(.02, r_nb[j] + 2, length.out = 100)
    f_grid <- sapply(seq_along(r_grid), function(i){
      loglik_r(r_grid[i], y_present, mean_uv)
    })
    
    qplot(r_grid, f_grid)
    
    # pi <- mean_uv / (mean_uv + r_nb[j])
    # sqrt(pi * r_nb[j] / (1 - pi))
    
    r_star <- optim_rcpp(r_nb[j], y_present, mean_uv)
    sd_star <- - 1 / d2loglik_r_cpp(r_star, y_present, mean_uv)
    
    r_new <- rnorm(1, r_star, sqrt(sd_star))
    
    loglik_new <- loglik_star(r_new, y_present, mean_uv)
    loglik_current <- loglik_star(r_nb[j], y_present, mean_uv)
    # exp(loglik_new - loglik_current)
    if(runif(1) < exp(loglik_new - loglik_current)){
      r_nb[j] <- r_new
    }
    
  }
  
  r_nb
  
}

update_logz_corr <- function(beta0, beta_z, lambda, mu, 
                             logz, v, Tau, delta, X_z){
  
  list_CP_cpp = convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, delta, 
                                  gamma, beta_theta, M_site)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  
  Xb <-  X_z %*% beta_z
  
  invTau <- solve(Tau)
  
  # logz <- matrix(NA, n, S)
  for (i in 1:n) {
    
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
      v_all <- beta_theta0[j, 1] +
        r[1:M_site[i] + sum(M_site[seq_len(i-1)])] * beta_theta[j, 3]
      x_all <- beta_theta[j, 2] / exp(lambda[j])
      
      logz_star_mean[j] <- findzero_likonly_cpp_corr(logzbar_current_j - 4,
                                                     logzbar_current_j + 4,
                                                     .01,
                                                     x_all,
                                                     y_all, v_all, 
                                                     v_samples,
                                                     sigma[j])
      
      logz_star_sd[j] <- sqrt(1 / (-h_f_loglik_cpp(logz_star_mean[j],
                                                   x_all,
                                                   y_all, v_all, 
                                                   v_samples,
                                                   sigma[j])))
      
      # invSigma[j,j] <- sqrt(1 /  logz_star_sd[j]^2)
      
    }
    
    invSigma <- diag(1 / logz_star_sd^2)
    post_cov <- solve(invTau + invSigma)
    post_mu <- (invSigma %*% logz_star_mean) + invTau %*% t(Xb_i)
    
    post_mean <- post_cov %*% post_mu
    
    logz_star <- mvrnorm(1, post_mean, post_cov)
    
    logposterior_ratios <- 0
    
    for (j in 1:S) {
      
      logzbar_current_j = logz_bar_current[j]
      
      idxPresent <- which(delta[1:M_site[i] + sum(M_site[seq_len(i-1)]),j] == 1)
      
      v_samples <- v_bar[idxPresent + sum(M_site[seq_len(i-1)]),j] - 
        r[idxPresent + sum(M_site[seq_len(i-1)])] * alpha[j]
      
      y_all <- delta[1:M_site[i] + sum(M_site[seq_len(i-1)]),j]
      v_all <- beta_theta0_bar[j] +
        r[1:M_site[i] + sum(M_site[seq_len(i-1)])] * beta_theta[j, 3]
      x_all <- beta_theta[j, 2]
      
      {
        # logz_grid <- seq(logzbar_current_j - 1, logzbar_current_j + 1, length.out = 100)
        # fun_grid <- sapply(seq_along(logz_grid), function(i){
        #   logposterior_logz_both_cpp(logz_grid[i], v_samples, sigma[j],
        #                              y_all, x_all,
        #                              v_all)
        #   # logf_likonly_cpp(logz_grid[i], x_all, sigma[j], v_samples, v_all, y_all, 0, prior_mean)
        # })
        # 
        # approx_grid <- sapply(seq_along(logz_grid), function(i){
        #   dnorm(logz_grid[i], logz_star_mean[j], logz_star_sd[j], log = T)
        # })
        # 
        # qplot(logz_grid, fun_grid)
        # ggplot() + 
        #   geom_line(data = NULL, aes(x = logz_grid, y = exp(fun_grid - max(fun_grid)))) + 
        #   geom_line(data = NULL, aes(x = logz_grid, y = exp(approx_grid - max(approx_grid))), color = "red") 
      }
      
      logpost_diff <- 
        logposterior_logz_both_cpp(logz_star[j], v_samples, sigma[j],
                                   y_all, x_all,
                                   v_all) -
        logposterior_logz_both_cpp(logzbar_current_j, v_samples, sigma[j],
                                   y_all, x_all,
                                   v_all)
      # logposterior_logz_cpp(v_samples, sigma[j], logzbar_current_j, logz_star[j], Xb_i[j], 
      # Tau[j,j])  +
      # logposterior_logz_logistic_cpp(y_all, x_all,
      # v_all, logzbar_current_j, logz_star[j])
      
      logposterior_ratios <- logposterior_ratios + logpost_diff
      
      # print(logpost_diff)
      
    }
    
    # logz_star <- mvrnorm(1, Xb_i, Tau)
    
    logprior_ratio <- dmvnorm_cpp(logz_star, Xb_i, Tau, 1) - 
      dmvnorm_cpp(logz_bar_current, Xb_i, Tau, 1)
    
    logproposal_ratio <- dmvnorm_cpp(logz_bar_current, post_mean, post_cov, 1) - 
      dmvnorm_cpp(logz_star, post_mean, post_cov, 1) 
    
    (mh_ratio <- exp(logposterior_ratios +  logprior_ratio + logproposal_ratio))
    
    if(runif(1) < mh_ratio){
      
      logz_bar[i,] <- logz_star
      
    }
    
  }
  
  list_SP <- convertCPtoSP_cpp(beta_bar, lambda, mu_bar, logz_bar, v_bar, delta, 
                               gamma, beta_theta0_bar, beta_theta, M_site)
  beta0 <- list_SP$beta0
  mu <- list_SP$mu
  logz <- list_SP$logz
  v <- list_SP$v
  beta_theta <- list_SP$beta_theta
  
  list("logz" = logz,
       "v" = v)
}

# DELTA GAMMA C ---------

update_delta_foreach <- function(y, 
                                 v, 
                                 lambda, 
                                 r_nb,
                                 M_site, K, 
                                 lambdatilde,
                                 mu0, n0, pi0, u, 
                                 logz, X_w, 
                                 beta_w,
                                 sigma, 
                                 mu, sigma_gamma, v_sd = .5,
                                 p_11, 
                                 p_10, 
                                 theta11, 
                                 theta10, emptySites){
  
  list_delta_all <- foreach(m = 1:S, .combine=c, 
                            .packages = "functionsForForeach") %dopar% {
                              
                              list_delta <- functionsForForeach::update_delta_c_d_rjmcmc(y[,,m,drop=F], 
                                                                                         v[,m,drop=F], 
                                                                                         lambda[m], 
                                                                                         r_nb[m],
                                                                                         M_site, K, 
                                                                                         lambdatilde[m],
                                                                                         mu0, n0, pi0, u, 
                                                                                         logz[,m,drop=F], X_w, 
                                                                                         beta_w[,m,drop=F],
                                                                                         sigma[m], 
                                                                                         mu[m], sigma_gamma, v_sd,
                                                                                         p_11[m], 
                                                                                         p_10[m], 
                                                                                         theta11[,m,drop=F], 
                                                                                         theta10[m], emptySites)
                            }
  
  delta <- Reduce("cbind",list_delta_all[1 + 4 * 0:(S-1)])
  c_imk <- Reduce("abind",list_delta_all[2 + 4 * 0:(S-1)])
  gamma <- Reduce("cbind",list_delta_all[3 + 4 * 0:(S-1)])
  v <- Reduce("cbind",list_delta_all[4 + 4 * 0:(S-1)])
  
  list("delta" = delta,
       "c_imk" = c_imk,
       "gamma" = gamma,
       "v" = v)
}

# LAMBDA ----------

update_lambda_CP_old <- function(beta0, mu, lambda,
                                 sigma_beta, sigma_mu,
                                 explambda_prior, sigma_lambda,
                                 S_star,
                                 betaThetaEqual1,
                                 lambda_prop = .5){
  
  S <- length(mu)
  
  list_CP_cpp <- convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, 
                                   delta, gamma, beta_theta, M_site, S_star)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  beta_theta_bar <- list_CP_cpp$beta_theta_bar
  
  # update paramters
  
  b_lambda <- .001
  
  for (j in 1:S) {
    
    a_lambda_prior <- explambda_prior[j] / 1000
    
    nonPCRcounts <- as.vector(y[,,j])[as.vector(c_imk[,,j]) == 2]
    nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
    
    dpost <- function(lambda){
      
      Xbeta <- beta_theta[j,1] + (1 / exp(lambda)) * exp(rep(logz_bar[,j], M_site)) + 
        X_w %*% beta_theta[j,-c(1,2)]
      
      Xbeta <- pmax(pmin(Xbeta, 10),-10)
      
      sum1 <- sum(dpois(nonPCRcounts, exp(lambda) * lambdatilde[j], log = T)) + 
        dnorm(beta_bar[j], lambda, sigma_beta, log = T) +
        dnorm(mu_bar[j], lambda, sigma_mu, log = T) +
        dgamma(exp(lambda), a_lambda_prior, b_lambda, log = T)
      
      sumdbern <- sum(sapply(1:length(Xbeta), function(i){
        dbinom(delta[i,j], 1, logistic(Xbeta[i]), log = T)
      }))
      
      sum1 + sumdbern 
      
    }
    
    lambda_star <- rnorm(1, lambda[j], sd = lambda_prop)
    
    lambda_post <- dpost(lambda[j])
    lambdastar_post <- dpost(lambda_star)
    
    if(runif(1) < exp(lambdastar_post - lambda_post)){
      lambda[j] <- lambda_star
    }
    
  }
  
  for (j in seq_len(S_star)) {
    
    a_lambda_prior <- explambda_prior[S + j] / 1000
    
    nonPCRcounts <- as.vector(y[,,j])[as.vector(c_imk[,,j]) == 2]
    nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
    
    psi_gig <- length(nonPCRcounts) * lambdatilde[S + j] + b_lambda
    
    PCRcounts <- as.vector(lambda_ijk[,,S + j])[as.vector(c_imk[,,S + j]) == 1]
    PCRcounts <- PCRcounts[!is.na(PCRcounts)]
    
    vpu <- u 
    vpu_PCR <- as.vector(vpu)[as.vector(c_imk[,,S + j]) == 1]
    vpu_PCR <- vpu_PCR[!is.na(vpu_PCR)]
    
    chi_gig <- sum(PCRcounts * r_nb[S + j] / exp(vpu_PCR)) 
    
    lambda_gig <- a_lambda_prior + sum(nonPCRcounts) - length(PCRcounts) * r_nb[S + j]
    
    lambda[S + j] <- log(rgig(1, lambda = lambda_gig, chi = 2 * chi_gig, psi = 2 * psi_gig ))
    
  }
  
  lambda
  
}

logdpost <- function(lambda){
  
  Xbeta <- beta_theta[j,1] + (1 / exp(lambda)) * X_logz + 
    X_w %*% beta_theta[j,-c(1,2)]
  
  Xbeta <- pmax(pmin(Xbeta, 10),-10)
  
  sum1 <- sum(dpois(nonPCRcounts, exp(lambda) * lambdatilde[j], log = T)) + 
    dnorm(beta_bar[j], lambda, sigma_beta, log = T) +
    dnorm(mu_bar[j], lambda, sigma_mu, log = T) +
    dnorm(lambda_prior[j], lambda, sigma_lambda, log = T) 
  
  sumdbern <- sum(sapply(1:length(Xbeta), function(i){
    dbinom(delta[i,j], 1, logistic(Xbeta[i]), log = T)
  }))
  
  sum1 + sumdbern 
  
}

der_logdpost <- function(lambda){
  
  Xbeta <- beta_theta[j,1] + X_w %*% beta_theta[j,-c(1,2)]
  
  Xbeta <- pmax(pmin(Xbeta, 10),-10)
  
  sum1 <-  - a * exp(lambda) + sum(nonPCRcounts) -
    (1 / (sigma_mu^2 )) * (lambda - mu_bar[j]) -
    (1 / (sigma_beta^2 )) * (lambda - beta_bar[j]) -
    (1 / (sigma_lambda^2 )) * (lambda - lambda_prior[j]) 
  
  sumdbern <- sum(sapply(1:length(Xbeta), function(i){
    if(delta[i,j] == 1){
      - logistic(Xbeta[i]) * (X_logz[i] / exp(lambda)) * exp(- Xbeta[i] - X_logz[i] / exp(lambda))
    } else {
      logistic(Xbeta[i]) * X_logz[i] / exp(lambda) 
    }
  }))
  
  sum1 + sumdbern 
  
}

update_lambda_CP <- function(beta0, beta_z, logz, 
                             mu, lambda, v, u, lambda_ijk, r_nb,
                             c_imk, delta, gamma, beta_theta, 
                             M_site, 
                             sigma_beta, sigma_mu,
                             lambda_prior, sigma_lambda,
                             S_star){
  
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
    
    lambda[S + j] <- log(rgig(1, lambda = lambda_gig, chi = 2 * chi_gig, psi = 2 * psi_gig ))
    
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

rgig_trunc <- function(lambda, chi, psi, trunc){
  
  x <- rgig(1, lambda, chi, psi)
  
  while(x < trunc){
    x <- rgig(1, lambda, chi, psi)
  }
  
  x
}

update_lambda_NP <- function(lambda_ijk, c_imk, mu,
                             r_nb, v, u,
                             lambda_prior, sigma_lambda){
  
  S <- length(mu)
  
  # update paramters
  
  for (j in seq_len(S)) {
    
    a_lambda_prior <- exp(lambda_prior)[j]^2 / sigma_lambda^2
    b_lambda_prior <- exp(lambda_prior)[j] / sigma_lambda^2
    
    # nonPCRcounts <- as.vector(y[,,j])[as.vector(c_imk[,,j]) == 2]
    # nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
    # lengthnonPCRcounts <- sum(as.vector(c_imk[,,j]) == 2)
    
    # psi_gig <- length(nonPCRcounts) + b_lambda_prior
    psi_gig <- b_lambda_prior
    
    PCRcounts <- as.vector(lambda_ijk[,,j])[as.vector(c_imk[,,j]) == 1]
    PCRcounts <- PCRcounts[!is.na(PCRcounts)]
    
    vpu <- u + matrix(v[,j], nrow = sum(M_site), ncol = max(K), byrow = F) 
    vpu_PCR <- as.vector(vpu)[as.vector(c_imk[,,j]) == 1]
    vpu_PCR <- vpu_PCR[!is.na(vpu_PCR)]
    
    chi_gig <- sum(PCRcounts * r_nb[j] / exp(vpu_PCR)) 
    
    lambda_gig <- a_lambda_prior - length(PCRcounts) * r_nb[j]
    
    lambda[j] <- log(rgig(lambda = lambda_gig, chi = 2 * chi_gig, psi = 2 * psi_gig))
    
  }
  
  lambda
}

update_lambda_spikeins <- function(mu, 
                                   lambda,
                                   u,
                                   r_nb,
                                   lambda_prior, 
                                   sigma_lambda,
                                   S_star){
  
  S <- length(mu)
  
  b_lambda <- .001
  
  for (j in seq_len(S_star)) {
    
    a_lambda_prior <- exp(lambda_prior[S + j]) / 1000
    
    nonPCRcounts <- as.vector(y[,,j])[as.vector(c_imk[,,j]) == 2]
    nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
    
    psi_gig <- length(nonPCRcounts) * lambdatilde[S + j] + b_lambda
    
    PCRcounts <- as.vector(lambda_ijk[,,S + j])[as.vector(c_imk[,,S + j]) == 1]
    PCRcounts <- PCRcounts[!is.na(PCRcounts)]
    
    vpu <- u 
    vpu_PCR <- as.vector(vpu)[as.vector(c_imk[,,S + j]) == 1]
    vpu_PCR <- vpu_PCR[!is.na(vpu_PCR)]
    
    chi_gig <- sum(PCRcounts * r_nb[S + j] / exp(vpu_PCR)) 
    
    lambda_gig <- a_lambda_prior + sum(nonPCRcounts) - length(PCRcounts) * r_nb[S + j]
    
    lambda[S + j] <- log(rgig(1, lambda = lambda_gig, chi = 2 * chi_gig, psi = 2 * psi_gig ))
    
  }
  
  lambda
  
}

update_lambda <- function(beta0, mu, lambda,
                          sigma_beta, sigma_mu,
                          explambda_prior, sigma_lambda,
                          S_star,
                          betaThetaEqual1){
  
  S <- length(mu)
  
  # update paramters
  
  
  b_lambda <- .01
  
  for (j in 1:(S + S_star)) {
    
    a_lambda <- explambda_prior[j] / 100
    
    nonPCRcounts <- as.vector(y[,,j])[as.vector(c_imk[,,j]) == 2]
    nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
    
    psi_gig <- length(nonPCRcounts) * lambdatilde[j] + b_lambda
    
    PCRcounts <- as.vector(lambda_ijk[,,j])[as.vector(c_imk[,,j]) == 1]
    PCRcounts <- PCRcounts[!is.na(PCRcounts)]
    
    vpu <- u + matrix(v[,j], nrow = nrow(u), ncol = ncol(u), byrow = F)
    vpu_PCR <- as.vector(vpu)[as.vector(c_imk[,,j]) == 1]
    vpu_PCR <- vpu_PCR[!is.na(vpu_PCR)]
    
    chi_gig <- sum(PCRcounts * r_nb[j] / exp(vpu_PCR)) 
    
    lambda_gig <- a_lambda + sum(nonPCRcounts) - length(PCRcounts) * r_nb[j]
    
    lambda[j] <- log(rgig(1, lambda = lambda_gig, chi = 2 * chi_gig, psi = 2 * psi_gig ))
    
  }
  
  lambda
  
}

update_lambda_old <- function(beta0, mu, lambda,
                              sigma_beta, sigma_mu,
                              lambda_prior, sigma_lambda,
                              S_star,
                              betaThetaEqual1){
  
  S <- length(mu)
  
  list_CP_cpp <- convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, 
                                   delta, gamma, beta_theta, M_site, S_star)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  beta_theta_bar <- list_CP_cpp$beta_theta_bar
  
  # update paramters
  
  for (j in 1:S) {
    # print(j)
    
    # # l_differences <- logz_bar[,j] - X_z %*% beta_z[,j]
    # # n_obs <- length(l_differences)
    # 
    nonPCRcounts <- as.vector(y[,,j])[as.vector(c_imk[,,j]) == 2]
    nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
    # 
    a <- length(nonPCRcounts) * lambdatilde[j]
    # 
    # # laplace 
    # 
    # b <- - 1 / (2 * sigma_mu^2) - 1 / (2 * sigma_beta^2) - 1 / (2 * sigma_lambda^2) 
    # 
    # # b <- - 1 / (sigma_mu^2) - 
    # # n_obs / (tau[j]^2) - beta_theta[j,2]^2 / (sigma_beta^2)
    # 
    # c <- sum(nonPCRcounts) + 
    #   (beta_bar[j] / sigma_beta^2) + 
    #   (lambda_prior[j] / sigma_lambda^2) + 
    #   # (sum(l_differences) / tau[j]^2) +
    #   (mu_bar[j] / sigma_mu^2) #+ 
    #   # (- beta_theta[j,2] * beta_theta0_bar[j] / sigma_beta^2)
    # 
    # if(a > 0){
    #   
    #   log_lam_argument <- log(-a / (2*b) ) + (-c/(2*b))
    #   
    #   if(log_lam_argument > 100){
    #     lamW0 <- log_lam_argument - log(log_lam_argument)
    #     lambda_star <- (- 2 * b * lamW0 - c) / (2 * b)
    #   } else {
    #     lambda_star <- (- 2 * b * lambertW0(exp(log_lam_argument)) - c) / (2 * b)
    #   }
    #   
    # } else {
    #   
    #   lambda_star <- -c / (2*b)
    #   
    # }
    
    
    #
    if(!betaThetaEqual1){
      f_j <- function(lambda){
        
        - a * exp(lambda) + sum(nonPCRcounts) * lambda - 
          (1 / (2 * sigma_lambda^2 )) * (lambda - lambda_prior[j])^2 -
          (1 / (2 * sigma_mu^2 )) * (lambda - mu_bar[j])^2 -
          (1 / (2 * sigma_beta^2 )) * (lambda - beta_bar[j])^2 +
          2 * lambda - beta_theta_bar[j,2] * exp(2 * lambda) / (2 * sigma_beta^2)
      }
      
      df_j <- function(lambda){
        - a * exp(lambda) + sum(nonPCRcounts) -
          (1 / (sigma_lambda^2 )) * (lambda - lambda_prior[j]) -
          (1 / (sigma_mu^2 )) * (lambda - mu_bar[j]) -
          (1 / (sigma_beta^2 )) * (lambda - beta_bar[j]) +
          2 - 2 * beta_theta_bar[j,2] * exp(2 * lambda) / (2 * sigma_beta^2)
      }
    } else {
      f_j <- function(lambda){
        
        - a * exp(lambda) + sum(nonPCRcounts) * lambda - 
          (1 / (2 * sigma_lambda^2 )) * (lambda - lambda_prior[j])^2 -
          (1 / (2 * sigma_mu^2 )) * (lambda - mu_bar[j])^2 -
          (1 / (2 * sigma_beta^2 )) * (lambda - beta_bar[j])^2 
      }
      
      df_j <- function(lambda){
        - a * exp(lambda) + sum(nonPCRcounts) -
          (1 / (sigma_lambda^2 )) * (lambda - lambda_prior[j]) -
          (1 / (sigma_mu^2 )) * (lambda - mu_bar[j]) -
          (1 / (sigma_beta^2 )) * (lambda - beta_bar[j]) 
      }
    }
    
    
    
    # lambda_grid <- seq(lambda[j] - 20, lambda[j] + 20, length.out = 100)
    # y_grid <- rep(NA, length(lambda_grid))
    # for (i in seq_along(lambda_grid)) {
    #   y_grid[i] <- f_j(lambda_grid[i])
    # }
    # 
    # qplot(lambda_grid, exp(y_grid - max(y_grid)))
    
    lambda[j] <- ars::ars(1, f_j, df_j, x = c(lambda[j] - 20, 
                                              lambda[j], 
                                              lambda[j] + 20),
                          ns = 100)
    
  }
  
  list_SP <- convertCPtoSP_cpp(beta_bar, lambda, mu_bar, logz_bar, v_bar, delta, 
                               gamma, beta_theta, M_site, S_star)
  beta0 <- list_SP$beta0
  mu <- list_SP$mu
  logz <- list_SP$logz
  v <- list_SP$v
  beta_theta <- list_SP$beta_theta
  
  for (j in S + seq_len(S_star)) {
    
    PCRcounts <- as.vector(y[,,j])[as.vector(c_imk[,,j]) == 1]
    PCRcounts <- PCRcounts[!is.na(PCRcounts)]
    PCRlambdas <- as.vector(lambda_ijk[1:sum(M_site),,j])[as.vector(c_imk[,,j]) == 1]
    PCRlambdas <- PCRlambdas[!is.na(PCRlambdas)]
    uplusv <- matrix(NA, sum(M_site), max(K))
    for (i in 1:n) {
      for (m in 1:M_site[i]) {
        for (k in 1:K[m + sum(M_site[seq_len(i-1)])]) {
          uplusv[m + sum(M_site[seq_len(i-1)]),k] <- 
            v[m + sum(M_site[seq_len(i-1)]),j] +
            u[m + sum(M_site[seq_len(i-1)]),k]
        }
      }
    }
    uplusv_PCR <- as.vector(uplusv)[as.vector(c_imk[,,j]) == 1]
    uplusv_PCR <- uplusv_PCR[!is.na(uplusv_PCR)]
    
    nonPCRcounts <- as.vector(y[,,j])[as.vector(c_imk[,,j]) == 2]
    nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
    # 
    a <- length(nonPCRcounts) * lambdatilde[j]
    
    b <- length(PCRcounts) * r_nb[j]
    
    c_term <- - r_nb[j] * sum(PCRlambdas / exp(uplusv_PCR))  
    
    f_j <- function(lambda){
      
      - a * exp(lambda) + sum(nonPCRcounts) * lambda +
        - b * lambda + exp(-lambda) * c_term + 
        (1 / (2 * sigma_lambda^2 )) * (lambda - lambda_prior[j])^2 
    }
    
    df_j <- function(lambda){
      - a * exp(lambda) + sum(nonPCRcounts) -
        b - exp(-lambda) * c_term +
        (1 / (sigma_lambda^2 )) * (lambda - lambda_prior[j]) 
    }
    
    
    lambda[j] <- ars::ars(1, f_j, df_j, x = c(lambda[j] - 20, 
                                              lambda[j], 
                                              lambda[j] + 20),
                          ns = 100)
    
    
  }
  
  list("beta0" = beta0,
       "mu" = mu,
       "logz" = logz,
       "v" = v,
       "lambda" = lambda,
       "beta_theta" = beta_theta)
  
}

# LAMBDA 0 ----------------------------------------------------------------

update_lambda0 <- function(y, c_imk, gamma_0, gamma_1){
  
  nonPCRcounts <- as.vector(y)[as.vector(c_imk) == 0]
  nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
  
  rgamma(1, gamma_0 + sum(nonPCRcounts), gamma_1 + length(nonPCRcounts))
}

update_lambda0_NB <- function(y, c_imk, mu0, n0, sd_mu0 = .05, sd_n0 = .025){
  
  nonPCRcounts <- as.vector(y)[as.vector(c_imk) == 0]
  nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
  
  # propose new sets of parameters
  mu0_star <- rnorm(1, mu0, sd_mu0)
  n0_star <- rnorm(1, n0, sd_n0)
  
  if(n0_star > 0 & mu0_star > 0){
    lik_star <- sum(dnbinom(nonPCRcounts, mu = mu0_star, size = n0_star, log = T))
    lik_current <- sum(dnbinom(nonPCRcounts, mu = mu0, size = n0, log = T))
    
    if(runif(1) < exp(lik_star - lik_current)){
      mu0 <- mu0_star
      n0 <- n0_star
    }
    
  }
  
  list("mu0" = mu0,
       "n0" = n0)
}

update_lambda0_mixt <- function(y, c_imk, mu0, n0, sd_mu0 = .05, sd_n0 = .025){
  
  nonPCRcounts <- as.vector(y)[as.vector(c_imk) == 0]
  nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
  
  numZeroCounts <- sum(nonPCRcounts > 0)
  
  pi0 <- rbeta(1, 1 + length(nonPCRcounts) - numZeroCounts, 1 + numZeroCounts)
  
  # propose new sets of parameters
  
  nonzeroPCRcounts <- nonPCRcounts[nonPCRcounts > 0]
  
  mu0_star <- rnorm(1, mu0, sd_mu0)
  n0_star <- rnorm(1, n0, sd_n0)
  
  if(n0_star > 0 & mu0_star > 0){
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

rgamma_trunc <- function(alpha, beta){
  
  x <- rgamma(1, shape = alpha, 
              rate = beta)
  
  while(x > 1){
    x <- rgamma(1, shape = alpha, 
                rate = beta)
  }
  
  return(x)
}

# update_lambdatilde <- function(y, c_imk, lambda){
#   
#   lambdatilde <- sapply(1:(S + S_star), function(j){
#     
#     nonPCRcounts <- as.vector(y[,,j])[as.vector(c_imk[,,j]) == 2]
#     nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
#     
#     if(length(nonPCRcounts) > 0){
#       if((sum(nonPCRcounts) + 1) / (length(nonPCRcounts) * exp(lambda[j])) > 1){
#         lambdatilde_j <- 1
#       } else {
#         lambdatilde_j <- rgamma_trunc(sum(nonPCRcounts) + 1, 
#                                       length(nonPCRcounts) * exp(lambda[j]))
#       }
#       # rgamma(1, shape = sum(nonPCRcounts) + 1, 
#       #              rate = length(nonPCRcounts) * exp(lambda[j]))
#     } else {
#       lambdatilde_j <- runif(1)
#     }
#     
#     return(lambdatilde_j)
#   })
#   
#   lambdatilde
# }
# 
# update_lambdatilde <- function(y, c_imk){
#   
#   nonPCRcounts <- as.vector(y)[as.vector(c_imk) == 2]
#   nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
#   
#   rexp_trun
#   
#   lambdatilde
# }

update_lambda_tilde_NB <- function(y, c_imk, mu_tilde, n_tilde, sd_mu0 = 1, sd_n0 = 1){
  
  nonPCRcounts <- as.vector(y)[as.vector(c_imk) == 2]
  nonPCRcounts <- nonPCRcounts[!is.na(nonPCRcounts)]
  
  # propose new sets of parameters
  mu_tilde_star <- rnorm(1, mu_tilde, sd_mu0)
  n_tilde_star <- n_tilde#rnorm(1, n_tilde, sd_n0)
  
  if(n_tilde_star > 0 & mu_tilde_star > 0){
    lik_star <- sum(dnbinom(nonPCRcounts, mu = mu_tilde_star, size = n_tilde_star, log = T)) + 
      n_tilde_star
    lik_current <- sum(dnbinom(nonPCRcounts, mu = mu_tilde, size = n_tilde, log = T)) + 
      n_tilde
    
    if(runif(1) < exp(lik_star - lik_current)){
      mu_tilde <- mu_tilde_star
      n_tilde <- n_tilde_star
    }
    
  }
  
  list("mu_tilde" = mu_tilde,
       "n_tilde" = n_tilde)
}


# TAU ---------------------------------------------------------------------

rinvgamma <- function(a, b){
  1 / stats::rgamma(1, shape = a, rate = b)
}

update_tau <- function(tau, logz, X_z, beta_z, beta0, a_tau, b_tau){
  
  for (j in 1:S) {
    
    n_samples <- 0
    sumsq <- 0
    
    for (i in 1:n) {
      if(emptySites[i] == 0){
        
        sumsq <- sumsq + ((X_z[i,-1,drop=T] %*% beta_z[,j] + beta0[j]) - logz[i,j])^2
        n_samples <- n_samples + 1
      }
    }
    
    tau[j] <- sqrt(rinvgamma(a_tau + n_samples / 2, b_tau[j] + sumsq / 2))
  }
  
  tau
}

updateGHParameters <- function(n, S, GH_params, lambda_Y){
  
  Omega <- GH_params$Omega
  lambdasq <- GH_params$lambdasq
  tausq <- GH_params$tausq
  xi <- GH_params$csi
  nu <- GH_params$nu
  Sigma <- GH_params$Sigma
  
  p <- ncol(Sigma)
  
  for (i in 1:p) {
    
    gamm <- rgamma(1, shape = n/2 + 1, rate = (S[i,i] + lambda_Y) / 2)
    
    invOmegaii <- Sigma[-i,-i] - Sigma[-i,i] %*% t(Sigma[-i,i]) / Sigma[i,i]
    
    C <- solve((S[i,i] + lambda_Y) * invOmegaii + solve(diag(lambdasq[-i,i] * tausq)))
    # C <- (1 / S[i,i]) * solve(invOmegaii + solve(S[i,i] * Lambda_star * tau^2))
    
    beta <- mvrnorm(1, - C %*% S[-i,i], C)
    # beta <- sampleNormFast(S[i,i] * invOmegaii + solve(S[i,i] * Lambda_star * tau^2), - S[-i,i])
    
    Omega[-i,i] <- beta
    Omega[i,-i] <- beta
    
    Omega[i,i] <- gamm + t(beta) %*% invOmegaii %*% beta
    
    Lambda_i <- sapply(1:(p-1), function(l){ 
      rinvgamma(1, (1 / nu[-i,i][l]) + beta[l]^2 / (2 * tausq))
      # 1 / rgamma(1, shape = 1, scale = (1 / N[-i,i][l]) + beta[l]^2 / (2 * tausq))
    })
    
    lambdasq[-i,i] <- Lambda_i
    
    lambdasq[i,-i] <- lambdasq[-i,i]
    
    nu[-i,i] <- sapply(1:(p-1), function(l){
      rinvgamma(1, 1 + 1 / Lambda_i[l])
      # 1 / rgamma(1, shape = 1, scale = 1 + 1 / Lambda_i[l])
    })
    
    nu[i,-i] <- nu[-i,i]
    
    # recompute Sigma
    {
      Sigma[-i,-i] <- invOmegaii + (invOmegaii %*% beta) %*% t(invOmegaii %*% beta) / gamm
      
      Sigma[-i,i] <- - (invOmegaii %*% beta) / gamm
      Sigma[i,-i] <- Sigma[-i,i]
      
      Sigma[i,i] <- 1 / gamm  
    }
    
  }
  
  # update tau
  tausq <- rinvgamma((choose(p, 2) + 1) / 2,
                     (1 / xi) + (0.5) * sum( (Omega[lower.tri(Omega)]^2) / lambdasq[lower.tri(lambdasq)] ))
  
  xi <- rinvgamma(1, 1 + (1 / tausq))
  
  list("Omega" = Omega,
       "lambdasq" = lambdasq,
       "tausq" = tausq,
       "nu" = nu,
       "csi" = xi,
       "Sigma" = Sigma)
}

update_Tau <- function(X_z, logz, beta0, beta_z,
                       Tau_params, Tau_priors,
                       jointSpecies){
  
  
  if(jointSpecies){
    
    logz_tilde <- logz - cbind(1, X_z) %*% rbind(t(beta0), beta_z)
    
    n <- nrow(logz_tilde)
    d <- ncol(logz_tilde)
    
    sum_S <- 0
    S_matrices <- sapply(1:n, function(i){
      sum_S <<- sum_S + logz_tilde[i,] %*% t(logz_tilde[i,]) 
    })
    
    lambda_Y <- Tau_priors$lambda_Y
    Tau_params <- updateGHParameters(n, sum_S, Tau_params, lambda_Y)
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
    
    tau <- update_tau_cpp(tau, logz, X_z, beta_z, beta0,
                          a_tau, rep(b_tau, S)
                          # a_tau = .5, b_tau = 1 / a_tau
    )
    
    Tau_params <- list("tau" = tau)
  }
  
  
  Tau_params
}

update_Tau_old <- function(Omega, Sigma, tau, lambda, gamma, lambda_Y,
                           logz, beta0, beta_z, X_z){
  
  logz_tilde <- logz - cbind(1, X_z) %*% rbind(t(beta0), beta_z)
  
  n <- nrow(logz_tilde)
  d <- ncol(logz_tilde)
  
  S <- 0
  S_matrices <- sapply(1:n, function(i){
    S <<- S + logz_tilde[i,] %*% t(logz_tilde[i,]) 
  })
  
  # update parameters
  for (i in 1:d) {
    
    gamm <- rgamma(1, shape = n/2 + 1, rate = (S[i,i] + lambda_Y) / 2)
    
    invOmegaii <- Sigma[-i,-i] - Sigma[-i,i] %*% t(Sigma[-i,i]) / Sigma[i,i]
    
    C <- solve((S[i,i] + lambda_Y) * invOmegaii + diag(1 / (tau[-i,i]), nrow = d - 1))
    
    beta <- mvrnorm(1, - C %*% S[-i,i], C)
    # beta <- sampleNormFast(S[i,i] * invOmegaii + solve(S[i,i] * Lambda_star * tau^2), - S[-i,i])
    
    Omega[-i,i] <- beta
    Omega[i,-i] <- beta
    
    Omega[i,i] <- gamm + t(beta) %*% invOmegaii %*% beta
    
    # recompute Sigma
    {
      Sigma[-i,-i] <- invOmegaii + (invOmegaii %*% beta) %*% t(invOmegaii %*% beta) / gamm
      
      Sigma[-i,i] <- - (invOmegaii %*% beta) / gamm
      Sigma[i,-i] <- Sigma[-i,i]
      
      Sigma[i,i] <- 1 / gamm  
    }
    
  }
  
  for (i in seq_len(d)) {
    for (j in seq_len(i-1)) {
      tau[i,j] <- rgig(1, lambda - 1/2, Omega[i,j]^2, 1 / gamma^2)
      tau[j,i] <- tau[i,j]
    }
  }
  
  # update lambda
  
  lambda <- update_lambda_GH(lambda, tau, gamma, sigma_z = .0002)
  
  # update tau
  
  gamma <- update_gamma_GH(gamma, M = 1, lambda)
  
  
  list("Omega" = Omega,
       "Sigma" = Sigma,
       "tau" = tau,
       "lambda" = lambda,
       "gamma" = gamma)
}

update_Tau_Wish <- function(nu0, Sigma0, logz, beta0, betaz, X_z){
  
  logz_tilde <- logz - cbind(1, X_z) %*% rbind(t(beta0), beta_z)
  
  n <- nrow(logz_tilde)
  d <- ncol(logz_tilde)
  
  S <- 0
  S_matrices <- sapply(1:n, function(i){
    S <<- S + logz_tilde[i,] %*% t(logz_tilde[i,]) 
  })
  
  Sigma <- rinvwishart(1, nu = nu0 + n, Omega = Sigma0 + S)[,,1]
  
  invSigma <- solve(Sigma)
  
  list("invSigma" = invSigma,
       "Sigma" = Sigma)
}

# OTHER PARAMETERS ---------

update_lambda_GH <- function(lambda, tau, gamma_NG, sigma_z){
  
  tau_all <- tau[lower.tri(tau)]
  p <- length(tau_all)
  
  lambda_star <- exp(rnorm(1, sd = sigma_z)) * lambda
  
  mh_ratio_term1 <- exp(dexp(lambda_star, rate = 1, T) - dexp(lambda, rate = 1, T))
  
  mh_ratio_term2 <- exp(p * (lgamma(lambda) - lgamma(lambda_star)))
  
  mh_ratio_term3 <- exp( (lambda_star - lambda) * ( -p * log(2 * gamma_NG^2) + sum(log(tau))  ) )
  
  mh_ratio <- mh_ratio_term1 * mh_ratio_term2 * mh_ratio_term3
  
  if(runif(1) < mh_ratio){
    lambda <- lambda_star
  }
  
  lambda
}

update_gamma_GH <- function(tau, M, lambda){
  
  tau_all <- tau[lower.tri(tau)]
  p <- length(tau_all)
  
  invgammasq <- rgamma(1, 2 + p * lambda, (M / (2 * lambda)) + (1/2) * sum(tau))
  
  sqrt(1 / invgammasq)
}

update_sigma <- function(sigma, v, delta, r, alpha, a_sigma, b_sigma){
  
  for (j in 1:S) {
    
    n_samples <- 0
    sumsq <- 0
    
    for (i in 1:n) {
      if(emptySites[i] == 0){
        for (m in 1:M_site[i]) {
          
          if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1){
            sumsq <- sumsq + (v[m + sum(M_site[seq_len(i-1)]),j] - 
                                (logz[i,j] + 
                                   r[m + sum(M_site[seq_len(i-1)])] * alpha[j]))^2
            
            n_samples <- n_samples + 1
          }
          
        }
        
      }
    }
    
    sigma[j] <- sqrt(rinvgamma(a_sigma + (n_samples / 2), b_sigma[j] + (sumsq / 2)))
  }
  
  sigma
}

update_beta_theta <- function(theta11, beta_theta, delta, r, logz){
  
  for (j in 1:S) {
    
    y_long <- c()
    X_long <- matrix(NA, nrow = 0, ncol = 3)
    for (i in 1:n) {
      if(!(emptySites[i] == 1)){
        for (m in 1:M_site[i]) {
          
          y_long <- c(y_long, delta[m + sum(M_site[seq_len(i-1)]),j])
          X_long <- rbind(X_long, 
                          c(1, logz[i,j], 
                            r[m + sum(M_site[seq_len(i-1)])]))
          
        }
      }
    }
    
    beta_theta[j,] <- sample_betaPG(beta_theta[j,], X_long, b = rep(0, 3),
                                    B = diag(1, nrow = 3), n = rep(1, nrow(X_long)),
                                    k = y_long - .5)
    
    theta11[,j] <- logistic(X_long %*% beta_theta[j,])
    
  }
  
  list("beta_theta" = beta_theta,
       "theta11" = theta11)
}

update_p_11 <- function(p_11, delta, gamma, c_imk){
  
  for (j in 1:S) {
    
    numPresents <- 0
    numSuccesses <- 0
    
    for (i in 1:n) {
      
      for (m in 1:M_site[i]) {
        
        if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1 | 
           gamma[m + sum(M_site[seq_len(i-1)]),j] == 1){
          
          
          currentK <- K[m + sum(M_site[seq_len(i-1)])]
          
          for (k in 1:currentK) {
            
            if(c_imk[m + sum(M_site[seq_len(i-1)]),k,j] == 1){
              
              numSuccesses <- numSuccesses + 1
              
            }
            
            numPresents <- numPresents + 1
            
          }
          
        }
        
      }
      
    }
    
    p_11[j] <- rbeta(1, 1 + numSuccesses, 1 + numPresents - numSuccesses)
    
  }
  
  p_11
}

update_p_10 <- function(p_10, delta, gamma, c_imk){
  
  for (j in 1:S) {
    
    numPresents <- 0
    numSuccesses <- 0
    
    for (i in 1:n) {
      
      for (m in 1:M_site[i]) {
        
        if(delta[m + sum(M_site[seq_len(i-1)]),j] == 0 & 
           gamma[m + sum(M_site[seq_len(i-1)]),j] == 0){
          
          currentK <- K[m + sum(M_site[seq_len(i-1)])]
          
          for (k in 1:currentK) {
            
            if(c_imk[m + sum(M_site[seq_len(i-1)]),k,j] == 2){
              
              numSuccesses <- numSuccesses + 1
              
            }
            
            numPresents <- numPresents + 1
            
          }
          
        }
        
      }
      
    }
    
    p_10[j] <- rbeta(1, 1 + numSuccesses, 1 + numPresents - numSuccesses)
    
  }
  
  p_10
}

update_theta10 <- function(theta10, delta, a0, b0){
  
  for (j in 1:S) {
    
    numPresents <- 0
    numSuccesses <- 0
    
    for (i in 1:n) {
      
      if(emptySites[i] == 0){
        
        for (m in 1:M_site[i]) {
          
          if(delta[m + sum(M_site[seq_len(i-1)]),j] == 0){
            
            numPresents <- numPresents + 1
            
            if(gamma[m + sum(M_site[seq_len(i-1)]),j] == 1){
              
              numSuccesses <- numSuccesses + 1
              
            }
            
          }
          
        }
        
      }
      
    }
    
    theta10[j] <- rbeta(1, a0 + numSuccesses, b0 + numPresents - numSuccesses)
    
  }
  
  theta10
}

update_lambdaijk_r <- function(lambda_ijk, lambda, v, u){
  
  for (i in 1:n) {
    # print(i)
    for (m in 1:M_site[i]) {
      numRep <- K[m + sum(M_site[seq_len(i-1)])]
      for (k in 1:numRep) {
        for (j in 1:S) {
          
          mean_lambdaijk <- exp(lambda[j] + v[m + sum(M_site[seq_len(i-1)]), j] +
                                  u[m + sum(M_site[seq_len(i-1)]), k])
          lambda_ijk[m + sum(M_site[seq_len(i-1)]),k,j] <- 
            rgamma(1, shape = r_nb[j] + y[m + sum(M_site[seq_len(i-1)]),k,j], 
                   rate = 1 + r_nb[j] / mean_lambdaijk )
          
        }
      }
    }
  }
  
  lambda_ijk
}

# BETA --------------------------------------------------------------------

update_betaz_CP_old <- function(logz, lambda, tau, X_z, sigma_beta, S_star){
  
  list_CP_cpp <- convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v,
                                   delta, gamma, beta_theta, M_site, S_star)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  beta_theta0_bar <- list_CP_cpp$beta_theta0_bar
  
  ncov_z <- ncol(X_z)
  
  # X_beta <- X_z
  X_beta <- cbind(1, X_z)
  
  {
    tXX <- t(X_beta) %*% X_beta
    
    # for (j in 1:S) {
    # 
    #   Lambda_beta <- (tXX / tau[j]^2) + (diag(1, nrow = ncov_z) / sigma_beta^2)
    #   
    #   mu_beta_x <- matrix(t(sapply(1:nrow(X_beta), function(i){
    #     X_beta[i,] * logz[i,j] / tau[j]^2
    #   })), nrow(X_beta), ncov_z)
    #   mu_beta <- apply(mu_beta_x, 2, sum)
    #   b_prior <- rep(0, ncov_z)
    #   
    #   mu_beta <- mu_beta + b_prior
    #   # mu_beta <- apply(t(sapply(1:nrow(X_beta), function(i){
    #   #   X_beta[i,] * logz_bar[i,j] / tau[j]^2
    #   # })), 2, sum) +
    #   #   rep(0, ncov_z)
    # 
    #   beta_bar_beta <- mvrnorm(1, solve(Lambda_beta) %*% mu_beta, solve(Lambda_beta))
    # 
    #   beta_z[,j] <- beta_bar_beta
    # 
    # }
    
    for (j in 1:S) {
      
      Lambda_beta <- (tXX / tau[j]^2) + (diag(1, nrow = ncov_z + 1) / sigma_beta^2)
      
      # mu_beta_x <- matrix(t(sapply(1:nrow(X_beta), function(i){
      #   X_beta[i,] * logz_bar[i,j] / tau[j]^2
      # })), nrow(X_beta), 1 + ncov_z)
      
      mu_beta_x <- matrix(apply(X_beta, 2, function(x){
        x * logz_bar[,j]
      }) / tau[j]^2, nrow(X_beta), 1 + ncov_z)
      
      mu_beta <- apply(mu_beta_x, 2, sum)
      b_prior <- c(lambda[j] / sigma_beta^2, rep(0, ncov_z))
      
      mu_beta <- mu_beta + b_prior
      
      beta_bar_beta <- mvrnorm(1, solve(Lambda_beta) %*% mu_beta, solve(Lambda_beta))
      
      beta_bar[j] <- beta_bar_beta[1]
      beta_z[,j] <- beta_bar_beta[-1]
      
    }
  }
  
  list_SP <- convertCPtoSP_cpp(beta_bar, lambda, mu_bar, logz_bar,
                               v_bar, delta, gamma,
                               beta_theta, M_site, S_star)
  beta0 <- list_SP$beta0
  
  list("beta0" = beta0,
       "beta_z" = beta_z)
}

update_betaz_CP <- function(beta0, beta_z, logz, tau, X_z, sigma_beta){
  
  ncov_z <- ncol(X_z)
  
  # X_beta <- X_z
  X_beta <- cbind(1, X_z)
  
  {
    tXX <- t(X_beta) %*% X_beta
    
    for (j in 1:S) {
      
      Lambda_beta <- (tXX / tau[j]^2) + (diag(1, nrow = ncov_z + 1) / sigma_beta^2)
      
      mu_beta_x <- matrix(apply(X_beta, 2, function(x){
        x * logz[,j]
      }) / tau[j]^2, nrow(X_beta), 1 + ncov_z)
      
      mu_beta <- apply(mu_beta_x, 2, sum)
      b_prior <- c(0, rep(0, ncov_z))
      
      mu_beta <- mu_beta + b_prior
      
      beta_bar_beta <- mvrnorm(1, solve(Lambda_beta) %*% mu_beta, solve(Lambda_beta))
      
      beta0[j] <- beta_bar_beta[1]
      beta_z[,j] <- beta_bar_beta[-1]
      
    }
  }
  
  list("beta0" = beta0,
       "beta_z" = beta_z)
}

update_betaz_CP_corr_old <- function(logz, Tau, X_z, sigma_beta){
  
  list_CP_cpp <- convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, 
                                   delta, gamma, beta_theta, M_site, S_star)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  beta_theta0_bar <- list_CP_cpp$beta_theta0_bar
  
  {
    ncov_z <- ncol(X_z)
    
    X_beta <- cbind(1, X_z)
    l_noempty <- logz_bar
    
    tXX <- t(X_beta) %*% X_beta
    tXl <- t(X_beta) %*% logz_bar
    
    prior_mean <- matrix(0, nrow = ncov_z + 1, ncol = S)
    prior_mean[1,] <- lambda[1:S]
    
    M_term <- tXl + prior_mean
    U_term <- tXX + diag(sigma_beta^2, nrow = (ncov_z + 1))
    post_U <- solve(U_term)
    post_M <- post_U %*% M_term
    
    beta_bar_beta <- rmtrnorm(post_M, post_U, Tau)
    
    beta_bar <- as.matrix(beta_bar_beta[1,])
    beta_z <- as.matrix(beta_bar_beta[-1,,drop = F])
    
  }
  
  list_SP <- convertCPtoSP_cpp(beta_bar, lambda, mu_bar, logz_bar,
                               v_bar, delta, gamma,
                               beta_theta, M_site, S_star)
  beta0 <- list_SP$beta0
  mu <- list_SP$mu
  logz <- list_SP$logz
  v <- list_SP$v
  
  list("beta0" = beta0,
       "beta_z" = beta_z)
}

update_betaz_CP_corr <- function(beta0, beta_z, logz, Tau, X_z, sigma_beta){
  
  {
    ncov_z <- ncol(X_z)
    S <- ncol(Tau)
    
    X_beta <- cbind(1, X_z)
    l_noempty <- logz
    
    tXX <- t(X_beta) %*% X_beta
    tXl <- t(X_beta) %*% logz
    
    prior_mean <- matrix(0, nrow = ncov_z + 1, ncol = S)
    prior_mean[1,] <- rep(0, S)
    
    M_term <- tXl + prior_mean
    U_term <- tXX + diag(sigma_beta^2, nrow = (ncov_z + 1))
    post_U <- solve(U_term)
    post_M <- post_U %*% M_term
    
    beta_bar_beta <- rmtrnorm(post_M, post_U, Tau)
    
    beta0 <- as.matrix(beta_bar_beta[1,])
    beta_z <- as.matrix(beta_bar_beta[-1,,drop = F])
    
  }
  
  list("beta0" = beta0,
       "beta_z" = beta_z)
}

logposterior_beta <- function(beta, 
                              X_all, v_all, y_all,
                              X_all2, v_all2, y_all2){
  
  Xbeta <- X_all %*% beta
  Xbeta2 <- X_all2 %*% beta
  
  sum(sapply(1:length(y_all), function(i){
    dbinom(y_all[i], 1, 
           prob = logistic(Xbeta[i] + v_all[i]), log = T)
  })) +
    sum(sapply(1:length(y_all2), function(i){
      dpois(y_all2[i], lambda = exp(Xbeta2[i] + v_all2[i]), log = T)
    }))
  
}

log_fp <- function(beta){
  
  Xbeta <- X_all %*% beta
  Xbeta2 <- X_all2 %*% beta
  
  sapply(1:ncov, function(k){
    
    -sum(X_all[,k] * exp(v_all + Xbeta)) + sumyallXall[k] + 
      (-sum(X_all2[,k] * exp(v_all2 + Xbeta2)) + sumyall2Xall2[k])
    # -sum(X_all[,k] * exp(v_all + Xbeta)) + sum(y_all * X_all[,k]) + 
    #   (-sum(X_all2[,k] * exp(v_all2 + Xbeta2)) + sum(y_all2 * X_all2[,k]))
    
  })
  
}

H_f <- function(beta){
  
  Xbeta <- X_all %*% beta
  Xbeta2 <- X_all2 %*% beta
  
  H <- matrix(NA, ncov, ncov)
  for (l1 in 1:ncov_z) {
    for (l2 in 1:ncov_z) {
      H[l1,l2] = -sum(X_all[,l1] * X_all[,l2] * exp(v_all + Xbeta) / (1 + exp(v_all + Xbeta))^2) + 
        (-sum(X_all2[,l1] * X_all2[,l2] * exp(v_all2 + Xbeta2)))
    }
  }
  
  H
  
}

# OLD STUFF ---------------

update_alpha <- function(alpha, v, delta, logz, r, sigma, sigma_beta){
  
  list_CP_cpp <- convertSPtoCP_cpp(lambda, beta_z, beta0, mu, logz, v, 
                                   delta, gamma, beta_theta, M_site)
  beta_bar <- list_CP_cpp$beta_bar
  mu_bar <- list_CP_cpp$mu_bar
  logz_bar <- list_CP_cpp$logz_bar
  v_bar <- list_CP_cpp$v_bar
  beta_theta0_bar <- list_CP_cpp$beta_theta0_bar
  
  for (j in 1:S) {
    
    r_long <- c()
    y_long <- c()
    r_idx <- 1
    for (i in 1:n) {
      for (m in 1:M_site[i]) {
        if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1){
          r_long[r_idx] <- r[m + sum(M_site[seq_len(i-1)])]
          y_long[r_idx] <- v_bar[m + sum(M_site[seq_len(i-1)]),j] - logz_bar[i,j]
          r_idx <- r_idx + 1
        }
      }
    }
    
    Lambda_beta <- (r_long %*% r_long / sigma[j]^2) + (diag(1, nrow = 1) / sigma_beta^2)
    mu_beta <- sum(y_long * r_long) / sigma[j]^2
    
    alpha[j] <- mvrnorm(1, solve(Lambda_beta) %*% mu_beta, solve(Lambda_beta))
    
  }
  
  list_SP <- convertCPtoSP_cpp(beta_bar, lambda, mu_bar, logz_bar,
                               v_bar, delta, gamma,
                               beta_theta0_bar, beta_theta, M_site)
  beta0 <- list_SP$beta0
  mu <- list_SP$mu
  logz <- list_SP$logz
  v <- list_SP$v
  
  return(alpha)
}

update_c <- function(c_imk, y, w, lambda, u, p11){
  
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        w_PCR <- w[m + sum(M_site[seq_len(i-1)]),j]
        
        if(w_PCR > 0){
          for (k in 1:K) {
            
            
            PCR_counts <- y[m + sum(M_site[seq_len(i-1)]),k,j]
            lambda_PCR <- lambda[j]
            
            ratio_prob <- p11[j] / (1 - p11[j])
            ratio_likelihood <- dpois(PCR_counts[k], lambda = lambda_PCR * w_PCR * u[i,m,k], log = T) -
              dpois(PCR_counts[k], lambda = lambda0, log = T)
            
            ratio <- ratio_prob * ratio_likelihood
            
            mh_ratio <- 1 / (1 + (1 / ratio))
            
            if(runif(1) < mh_ratio){
              c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- 1
            } else {
              c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- 0
            }
            
          }  
        } else {
          c_imk[m + sum(M_site[seq_len(i-1)]),k,j] <- 0
        }
        
      }
    }
  }
  
  c_imk
}

update_w <- function(y, c_imk, logz, sigmasq, u, lambda){
  
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        
        idxReplicatePositive <- which(as.vector(c_imk[m + sum(M_site[seq_len(i-1)]),,j]) == 1)
        PCR_counts <- as.vector(y[m + sum(M_site[seq_len(i-1)]),,j])[idxReplicatePositive]
        
        # if(sum(PCR_counts) == 0 & length(PCR_counts) != 0){
        # browser()
        # }
        # PCR_counts[PCR_counts == 0] <- .0001
        
        a_w <- sum(PCR_counts)
        if(a_w == 0) a_w <- .0001
        b_w <- lambda[j] * sum(u[i,m,idxReplicatePositive])
        
        nsamples <- length(PCR_counts)
        lik_mean <- a_w / b_w
        lik_var <- a_w / (b_w)^2
        
        prior_mean <- logz[i,j]
        prior_var <- sigmasq
        if(nsamples > 0){
          posterior_var <- 1 / (1 / prior_var + nsamples / lik_var)
          posterior_mean <- ((prior_mean / prior_var) + (lik_mean / lik_var)) * posterior_var  
        } else {
          posterior_var <- 1 / (1 / prior_var )
          posterior_mean <- (prior_mean / prior_var) * posterior_var
        }
        
        
        w_star <- rnorm(1, posterior_mean, posterior_var)
        w_current <- w[m + sum(M_site[seq_len(i-1)]),j]
        
        if(w_star > 0){
          
          logposterior_ratio <- logposterior_w(PCR_counts, w_star, w_current, 
                                               u[i,m,idxReplicatePositive], lambda[j], 
                                               prior_mean, prior_var)
          
          if(runif(1) < exp(logposterior_ratio)){
            w[m + sum(M_site[seq_len(i-1)]),j] <- w_star
          }  
          
        }
        
      }
    }
  }
  
  return(w)
}

logposterior_w <- function(PCR_counts, wstar, wcurrent, 
                           u_im, lambda, logz, sigmasq){
  sum(dpois(PCR_counts, lambda = lambda * u_im * wstar, log = T)) -
    sum(dpois(PCR_counts, lambda = lambda * u_im * wcurrent, log = T)) + 
    dnorm(w_star, logz, sigmasq, log = T) - 
    dnorm(wcurrent, logz, sigmasq, log = T)
}

update_z <- function(logz, X, beta, w, delta, gamma, sigmasq, tausq){
  
  for (i in 1:n) {
    for (j in 1:S) {
      idxDelta1 <- delta[i,seq_len(M_site[i]),j] == 1
      logtildew <- log(w[i,seq_len(M_site[i]),j] - 
                         gamma[i,seq_len(M_site[i]),j] * eta[i,seq_len(M_site[i]),j])
      logtildew <- logtildew[idxDelta1]
      nsamples <- sum(idxDelta1)
      prior_mean <- X[i,] %*% beta[,j]
      prior_var <- tausq
      lik_mean <- sum(logtildew)
      lik_var <- sigmasq
      posterior_var <- 1 / (1 / prior_var + nsamples / tausq)
      posterior_mean <- ((prior_mean / prior_var) + (lik_mean / lik_var)) * posterior_var
      logz[i,j] <- rnorm(1, posterior_mean, sqrt(posterior_var))
    }
  }
  
  logz
}


update_beta_nleq <- function(beta_z, X_z, y, delta, beta_theta,
                             v_bar, u, r, alpha, logz_tilde){
  
  ncov <- nrow(beta_z)
  
  for (j in 1:S) {
    
    {
      # data from delta
      X_all <- matrix(0, nrow = 0, ncol = ncov_z)
      y_all <- c()
      v_all <- c()
      for (i in 1:n) {
        for (m in 1:M_site[i]) {
          X_all <- rbind(X_all,beta_theta[j,2] * X_z[i,-1,drop=F])
          y_all[m + sum(M_site[seq_len(i-1)])] <- delta[m + sum(M_site[seq_len(i-1)]),j]
          v_all[m + sum(M_site[seq_len(i-1)])] <- beta_theta[j,1] +
            beta_theta[j,2] * logz_bar[i,j] +
            beta_theta[j,3] * r[m + sum(M_site[seq_len(i-1)])]
        }
      }
      
      # data from y
      X_all2 <- matrix(0, nrow = 0, ncol = ncov)
      y_all2 <- c()
      v_all2 <- c()
      l <- 1
      for (i in 1:n) {
        for (m in 1:M_site[i]) {
          for (k in 1:K[m + sum(M_site[seq_len(i-1)])]) {
            if(c_imk[m + sum(M_site[seq_len(i-1)]),k,j] == 1){
              X_all2 <- rbind(X_all2,X_z)
              y_all2[l] <- y[m + sum(M_site[seq_len(i-1)]),k,j]
              v_all2[l] <- v_bar[m + sum(M_site[seq_len(i-1)]),j] +
                alpha[j] * r[m + sum(M_site[seq_len(i-1)])] +
                u[m + sum(M_site[seq_len(i-1)]),k]
              l <- l + 1
            }
          }
        }
      }
      
    }
    
    list_matrices <- createMatricesNleq(j - 1, X_z, delta, 
                                        beta_theta, logz_tilde, r, c_imk, y, alpha, 
                                        u, v_bar, M_site, K)
    X_all <-  list_matrices$X_all
    y_all <-  list_matrices$y_all
    v_all <-  list_matrices$v_all
    X_all2 <-  list_matrices$X_all2
    y_all2 <-  list_matrices$y_all2
    v_all2 <-  list_matrices$v_all2
    
    #
    
    sumyallXall <- sapply(1:ncov, function(k){
      sum(y_all * X_all[,k])
    })
    
    sumyall2Xall2 <- sapply(1:ncov, function(k){
      sum(y_all2 * X_all2[,k])
    })
    
    log_fp <- function(beta){
      
      Xbeta <- X_all %*% beta
      Xbeta2 <- X_all2 %*% beta
      
      sapply(1:ncov, function(k){
        
        -sum(X_all[,k] * exp(v_all + Xbeta)) + sumyallXall[k] + 
          (-sum(X_all2[,k] * exp(v_all2 + Xbeta2)) + sumyall2Xall2[k])
        # -sum(X_all[,k] * exp(v_all + Xbeta)) + sum(y_all * X_all[,k]) + 
        #   (-sum(X_all2[,k] * exp(v_all2 + Xbeta2)) + sum(y_all2 * X_all2[,k]))
        
      })
      
    }
    
    modelOutput <- nleqslv(beta_z[,j], log_fp)
    
    # log_fp(beta_z[,j])
    # log_fp_cpp(beta_z[,j],
    #        X_all, v_all, y_all,
    #        X_all2, v_all2, y_all2)
    
    beta_star <- modelOutput$x
    
    # cov_star <- solve(-H_f(beta))
    
    cov_star <- solve(-H_f_cpp(beta_star, X_all, v_all, y_all,
                               X_all2, v_all2, y_all2))
    
    beta_new <- mvrnorm(1, beta_star, Sigma = cov_star)
    
    # logposterior_star <- logposterior_beta(beta_new, X_all, v_all, y_all,
    #                                        X_all2, v_all2, y_all2)
    # logposterior_current <- logposterior_beta(beta_z[,j], X_all, v_all, y_all,
    #                                        X_all2, v_all2, y_all2)
    logposterior_star <- logposterior_beta_cpp(beta_new, X_all, v_all, y_all,
                                               X_all2, v_all2, y_all2)
    logposterior_current <- logposterior_beta_cpp(beta_z[,j], X_all, v_all, y_all,
                                                  X_all2, v_all2, y_all2)
    
    # logprior_star <- 
    
    if(runif(1) < exp(logposterior_star - logposterior_current)){
      
      beta_z[,j] <- beta_new
      
    }
    
  }
  
  beta_z
  
}

BinToDec <- function(x) {
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1)) 
}

DecToBin <- function(x, M){
  x <- rev(intToBits(x))
  x <- x[length(x) - 0:(M-1)]
  x <- paste(as.integer(x), collapse = "")  
  as.numeric(strsplit(x,"")[[1]]) 
}

compute_logprob_y_delta0 <- function(y_counts, c_imk_current, lambda_j,
                                     mu0, n0, lambdatilde){
  
  # lambda_currents <- ifelse(c_imk_current == 2, exp(lambda_j) * lambdatilde, lambda0)
  # if(c_imk_current)
  # sum(dpois(y_counts, lambda_currents, log = T))
  
  sum(sapply(1:currentK, function(k){
    if(c_imk_current[k] == 2){
      sum(dpois(y_counts[k], exp(lambda) * lambdatilde, log = T))
    } else {
      sum(dnbinom(y_counts[k], size = n0, mu = mu0, log = T))
    }
  }))
  
}

compute_logprob_y_delta1 <- function(y_counts, c_imk_current,  
                                     lambda, mu0, n0, lambdatilde,
                                     u_im, v_im){
  
  currentK <- length(y_counts)
  
  sum(sapply(1:currentK, function(k){
    if(c_imk_current[k] == 1){
      sum(dpois(y_counts[k], exp(lambda + u_im[k] + v_im), log = T))
    } else if(c_imk_current[k] == 2){
      sum(dpois(y_counts[k], exp(lambda) * lambdatilde, log = T))
    } else {
      sum(dnbinom(y_counts[k], size = n0, mu = mu0, log = T))
    }
  }))
  
}

update_delta_c_d <- function(delta, gamma, c_imk, y, v, u, 
                             lambda, lambdatilde, lambda0, 
                             theta11, theta10, p11, p10){
  
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        
        currentK <- K[m + sum(M_site[seq_len(i-1)])]
        
        y_counts <- y[m + sum(M_site[seq_len(i-1)]), 1:currentK, j]
        
        u_im <- u[m + sum(M_site[seq_len(i-1)]), 1:currentK]
        
        if(v_pres[m + sum(M_site[seq_len(i-1)]), j] == 1){
          v_star <- v[m + sum(M_site[seq_len(i-1)]), j]
        } else {
          v_star <- 0
        }
        
        
        log_allProbs <- rep(NA, 3 * 2^currentK)
        mat_delta_c_d <- matrix(NA, 3 * 2^currentK, 2 + currentK)
        
        # delta = 0, gamma = 0
        for (l in 1:(2^currentK)) {
          
          c_imk_current <- 2 * DecToBin(l - 1, currentK)
          
          log_prob_y <- compute_logprob_y_delta0_cpp(y_counts,
                                                     c_imk_current,
                                                     currentK,
                                                     n0, mu0, pi0,
                                                     n_tilde, mu_tilde, 
                                                     lambda[j])
          
          prob_delta_gamma <- log(1 - theta11[m + sum(M_site[seq_len(i-1)]),j]) + 
            log(1 - theta10[j])
          
          prob_d <- sum(dbinom(c_imk_current / 2, 1, p_10[j], log = T))
          
          log_allProbs[l] <- log_prob_y + prob_delta_gamma + prob_d
          
          mat_delta_c_d[l,] <- c(0, 0, c_imk_current)
        }
        
        # delta = 1
        for (l in 1:(2^currentK)) {
          
          c_imk_current <- DecToBin(l - 1, currentK)
          
          log_prob_y <- compute_logprob_y_delta1_rnb_cpp(y_counts, c_imk_current, 
                                                         currentK, n0, mu0, pi0,
                                                         r_nb[j],
                                                         v_star,
                                                         n_tilde, mu_tilde, 
                                                         lambda[j],
                                                         u_im)
          
          prob_delta <- log(theta11[m + sum(M_site[seq_len(i-1)]),j]) 
          
          log_prior_v <- dnorm(v_star, logz[i,j], sigma[j], 1)
          
          prob_c <- sum(dbinom(c_imk_current, 1, p_11[j], log = T))
          
          mat_delta_c_d[2^currentK + l,] <- c(1,0,c_imk_current)
          log_allProbs[2^currentK + l] <- log_prob_y + prob_delta + prob_c + log_prior_v
          
          
        }
        
        # gamma = 1
        for (l in 1:(2^currentK)) {
          
          c_imk_current <- DecToBin(l - 1, currentK)
          
          log_prob_y <- compute_logprob_y_delta1_rnb_cpp(y_counts, c_imk_current, 
                                                         currentK, n0, mu0, pi0,
                                                         r_nb[j],
                                                         v_star,
                                                         n_tilde, mu_tilde, 
                                                         lambda[j],
                                                         u_im)
          
          prob_delta <- log(1 - theta11[m + sum(M_site[seq_len(i-1)]), j]) + 
            log(theta10[j]) 
          
          prob_c <- sum(dbinom(c_imk_current, 1, p_11[j], log = T))
          
          log_prior_v <- dnorm(v_star, mu[j], sigma_gamma, 1)
          
          mat_delta_c_d[2 * 2^currentK + l,] <- c(0,1,c_imk_current)
          log_allProbs[2 * 2^currentK + l] <- log_prob_y + prob_delta + prob_c + log_prior_v
          
          
        }
        
        
        log_allProbs <- log_allProbs - max(log_allProbs)
        allProbs <- exp(log_allProbs) / sum(exp(log_allProbs))
        
        sampledComb_idx <- sample(1:length(allProbs), 1, prob = allProbs)
        
        sampledComb <- mat_delta_c_d[sampledComb_idx,]
        
        delta[m + sum(M_site[seq_len(i-1)]),j] <- sampledComb[1]
        gamma[m + sum(M_site[seq_len(i-1)]),j] <- sampledComb[2]
        c_imk[m + sum(M_site[seq_len(i-1)]),1:currentK,j] <- sampledComb[2 + 1:currentK]
      }
    }
  }
  
  list("delta" = delta,
       "gamma" = gamma,
       "c_imk" = c_imk)
}

compute_logprob_y_delta1_rnb_cpp_r <- function(y_counts, c_imk_current){
  
  sapply(1:currentK, function(k){
    if(c_imk_current[k] == 1){
      
    }
  })
  
}

update_delta_c_d_rjmcmc_r <- function(delta, gamma, c_imk, y, v, u, 
                                      lambda, lambdatilde, lambda0, 
                                      theta11, theta10, p11, p10){
  
  p0 <- n0 / (n0 + mu0)
  
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        
        currentK <- K[m + sum(M_site[seq_len(i-1)])]
        
        y_counts <- y[m + sum(M_site[seq_len(i-1)]), 1:currentK, j]
        
        u_im <- u[m + sum(M_site[seq_len(i-1)]), 1:currentK]
        v_im <- v[m + sum(M_site[seq_len(i-1)]), j]
        
        if(emptySites[i] == 0){ # not an empty tube
          
          log_allProbs <- rep(NA, 3 * 2^currentK)
          mat_delta_c_d <- matrix(NA, 3 * 2^currentK, 2 + currentK)
          
          if(delta[m + sum(M_site[seq_len(i-1)]), j] == 1 | 
             gamma[m + sum(M_site[seq_len(i-1)]), j] == 1){
            
            v_star <- v_im
            
            # delta = 0, gamma = 0
            for (l in 1:(2^currentK)) {
              
              c_imk_current <- 2 * DecToBin(l - 1, currentK)
              
              log_prob_y <- compute_logprob_y_delta0(y_counts, c_imk_current, 
                                                     lambda[j], mu0, n0, lambdatilde[j])
              
              prob_delta_gamma <- log(1 - theta11[m + sum(M_site[seq_len(i-1)]),j]) + 
                log(1 - theta10[j])
              
              prob_d <- sum(dbinom(c_imk_current / 2, 1, p_10[j], log = T))
              
              log_allProbs[l] <- log_prob_y + prob_delta_gamma + prob_d
              
              mat_delta_c_d[l,] <- c(0, 0, c_imk_current)
            }
            
            # delta = 1
            for (l in 1:(2^currentK)) {
              
              c_imk_current <- DecToBin(l - 1, currentK)
              
              log_prob_y <- compute_logprob_y_delta1_rnb_cpp(y_counts, c_imk_current, 
                                                             currentK, n0, p0, pi0,
                                                             r_nb[j],
                                                             v_star,
                                                             lambdatilde[j],
                                                             lambda[j],
                                                             u_im)
              
              prob_delta <- log(theta11[m + sum(M_site[seq_len(i-1)]),j]) 
              
              prob_c <- sum(dbinom(c_imk_current, 1, p_11[j], log = T))
              
              mat_delta_c_d[2^currentK + l,] <- c(1,0,c_imk_current)
              log_allProbs[2^currentK + l] <- log_prob_y + prob_delta + prob_c 
              
              
            }
            
            # gamma = 1
            for (l in 1:(2^currentK)) {
              
              c_imk_current <- DecToBin(l - 1, currentK)
              
              log_prob_y <- compute_logprob_y_delta1(y_counts, c_imk_current, 
                                                     lambda[j], mu0, n0, lambdatilde[j],
                                                     u_im, v_im)
              
              prob_delta <- log(1 - theta11[m + sum(M_site[seq_len(i-1)]), j]) + 
                log(theta10[j]) 
              
              prob_c <- sum(dbinom(c_imk_current, 1, p_11[j], log = T))
              
              mat_delta_c_d[2 * 2^currentK + l,] <- c(0,1,c_imk_current)
              log_allProbs[2 * 2^currentK + l] <- log_prob_y + prob_delta + prob_c 
              
              
            }
            
          }
          
          
        } 
        
        log_allProbs <- log_allProbs - max(log_allProbs)
        allProbs <- exp(log_allProbs) / sum(exp(log_allProbs))
        
        sampledComb_idx <- sample(1:length(allProbs), 1, prob = allProbs)
        
        sampledComb <- mat_delta_c_d[sampledComb_idx,]
        
        delta[m + sum(M_site[seq_len(i-1)]),j] <- sampledComb[1]
        gamma[m + sum(M_site[seq_len(i-1)]),j] <- sampledComb[2]
        c_imk[m + sum(M_site[seq_len(i-1)]),1:currentK,j] <- sampledComb[2 + 1:currentK]
      }
    }
  }
  
  list("delta" = delta,
       "gamma" = gamma,
       "c_imk" = c_imk)
}


simulate_data <- function(beta_true, tau_true, sigma_true, 
                          X_beta, X_p11, theta11_true, 
                          beta_p11_true, lambda_true,
                          lambda0_true, lambdatilde_true){
  
  logz_true <- X_beta %*% beta_true + sapply(1:S, function(j) rnorm(n, 0, sd = tau_true[j]))
  
  delta_true <- array(NA, dim = c(sum(M_site), S))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      for (j in 1:S) {
        delta_true[m + sum(M_site[seq_len(i-1)]),j] <- rbinom(1, 1, theta11_true[j])
      }
    }
  }
  
  gamma_true <- array(0, dim = c(sum(M_site), S))
  
  v_true <- matrix(NA, sum(M_site), S)
  for (i in 1:n) {
    if(emptySites[i] == 0){
      for (m in 1:M_site[i]) {
        for (j in 1:S) {
          v_true[m + sum(M_site[seq_len(i-1)]),j] <- 
            rnorm(1, logz_true[i,j] + eta_true[j], sigma_true[j])
        }
      }
    }
  }
  
  w_true <- delta_true * exp(v_true)
  
  p11_true <- logistic(X_p11 %*% beta_p11_true)
  
  c_imk_true <- array(NA, dim = c(sum(M_site), max(K), S))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      numRep <- K[m + sum(M_site[seq_len(i-1)])]
      for (k in 1:numRep) {
        for (j in 1:S) {
          if(emptySites[i] == 1){
            c_imk_true[m + sum(M_site[seq_len(i-1)]),k,j] <- 0
          } else if (delta_true[m + sum(M_site[seq_len(i-1)]),j] == 1){
            c_imk_true[m + sum(M_site[seq_len(i-1)]),k,j] <- rbinom(1, 1, p11_true[i])
          } else {
            c_imk_true[m + sum(M_site[seq_len(i-1)]),k,j] <- 0
          }
        }
      }
    }
  }
  
  d_imk_true <- array(NA, dim = c(sum(M_site), max(K), S))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      numRep <- K[m + sum(M_site[seq_len(i-1)])]
      for (k in 1:numRep) {
        for (j in 1:S) {
          if(c_imk_true[m + sum(M_site[seq_len(i-1)]),k,j] == 0){
            d_imk_true[m + sum(M_site[seq_len(i-1)]),k,j] <- rbinom(1, 1, p10_true[j])
          } else {
            d_imk_true[m + sum(M_site[seq_len(i-1)]),k,j] <- 0
          }
        }
      }
    }
  }
  
  u_true <- matrix(NA, sum(M_site), max(K))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      numRep <- K[m + sum(M_site[seq_len(i-1)])]
      for (k in 1:numRep) {
        u_true[m + sum(M_site[seq_len(i-1)]),k] <- rgamma(1, a_u_true, a_u_true)
      }
    }
  }
  
  y <- array(NA, dim = c(sum(M_site), max(K), S))
  for (i in 1:n) {
    for (m in 1:M_site[i]) {
      numRep <- K[m + sum(M_site[seq_len(i-1)])]
      for (k in 1:numRep) {
        u_imk <- u_true[m + sum(M_site[seq_len(i-1)]), k]
        for (j in 1:S) {
          if(c_imk_true[m + sum(M_site[seq_len(i-1)]),k,j] == 0 &
             d_imk_true[m + sum(M_site[seq_len(i-1)]),k,j] == 0){
            y[m + sum(M_site[seq_len(i-1)]),k,j] <- rpois(1, lambda0_true)
          } else if(c_imk_true[m + sum(M_site[seq_len(i-1)]),k,j] == 0 &
                    d_imk_true[m + sum(M_site[seq_len(i-1)]),k,j] == 1){
            y[m + sum(M_site[seq_len(i-1)]),k,j] <- rpois(1, lambdatilde_true[j])
          } else {
            y[m + sum(M_site[seq_len(i-1)]),k,j] <- 
              rpois(1, lambda_true *
                      u_imk * w_true[m + sum(M_site[seq_len(i-1)]), j])
          }
        }
      }
    }
  }
  
  y
}
