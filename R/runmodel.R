# GENERAL -----------------------------------------------------------------

#' Model fitting
#' 
#' @description This function fits the occupancy model of Diana et al. (2021). 
#' Note that in the following the parameters are described with the notations used in the paper.
#' 
#' @param data A list containing the data. Detailed informations are given in the vignette. The elements are:
#'  \itemize{
#' 
#'  \item \code{PCR_table}: A data.frame with the reads of the wildlife species. 
#'  \item \code{infos}: A data.frame with the samples informations.
#'  \item \code{X_z}: Covariate matrix for the biomass abundance
#'  \item \code{X_w}: Covariate matrix for the biomass detection
#'  \item \code{PCR_spike}: A data.frame with the reads of the spike-in species. 
#'  }
#'  
#' @param priors List of prior parameters. If not specified, default values are used.
#' @param jointSpecies Does the model include species correlation?
#' @param MCMCparams A list with components:
#' \itemize{
#'    \item \code{nchain}: Number of chains
#'    \item \code{nburn}: Number of burn-in iterations
#'    \item \code{niter}: Number of iterations
#'    \item \code{nthin}: Number of iterations to thin
#' }
#' @param paramsToSave A list with components the elements of which the MCMC output need to be saved. 

#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return 
#' A list with components:
#' 
#' \itemize{
#' 
#'  \item \code{params_output} A list with the MCMC output of each parameter specified in the \code{paramesToSave} parameter
#'      
#'  
#'  }
#' 
#' @examples
#' 
#' modelResults <- runEDNA(data, 
#'                         priors,
#'                         jointSpecies = T,
#'                         MCMCparams)  
#' 
runEDNA <- function(data,
                    priors = NULL,
                    jointSpecies = F,
                    paramsUpdate = NULL,
                    MCMCparams,
                    paramsToSave = NULL){
  
  data <- cleanData(data)
  
  if(missing(paramsUpdate)){
    paramsUpdate <- list(updateAll = T,
                         params = NULL,
                         correct = NULL,
                         trueParams = NULL)
  }
  
  if(missing(paramsToSave)){
    paramsToSave <- list(
      "lambda" = T,
      "beta0" = T,
      "mu" = T,
      "beta_z" = T,
      "logz" = T,
      "u" = T,
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
    )
  }
  
  modelResults <- fitModel(data,
                           priors,
                           jointSpecies,
                           paramsUpdate,
                           MCMCparams)
  
  modelResults
}

cleanData <- function(data){
  
  PCR_table <- data$PCR_table
  infos <- data$infos
  X_w0 <- data$X_w
  X_z0 <- data$X_z
  PCR_spike <- data$PCR_spike
  spikeInBiomass <- data$spikeInBiomass
  spikeInPresent <- data$spikeInPresent
  
  K <- infos$Replicates
  
  S <- ncol(PCR_table) / max(K)
  
  # sites <- unique(infos$Site)
  sites <- infos$Site[!duplicated(infos$Site)]
  sites <- setdiff(sites,"empty")
  n <- length(sites)
  # samples <- unique(infos$Sample[infos$Site %in% sites])
  nonEmptySamples <- infos$Sample[infos$Site %in% sites]
  samples <- nonEmptySamples[!duplicated(nonEmptySamples)]  
  M_site <- sapply(1:n, function(i){
    length(unique(infos$Sample[infos$Site == sites[i]]))
  })
  
  emptyTubes <- length(unique(infos$Sample[infos$Site == "empty"]))
  
  if(is.null(PCR_spike)){
    S_star <- 0  
  } else {
    S_star <- ncol(PCR_spike) / max(K)
  }
  
  if(is.null(spikeInBiomass)){
    v_spikes <- matrix(0, nrow = sum(M_site) + emptyTubes,
                       S_star)
  } else {
    v_spikes <- spikeInBiomass
  }
  
  if(is.null(spikeInPresent)){
    spikedSample <- rbind(matrix(1, sum(M_site), S_star),
                          matrix(0, emptyTubes, S_star))
  } else {
    spikedSample <- spikeInPresent
  }
  
  OTUnames <- colnames(PCR_table[,1:S])
  if(S_star > 0){
    OTUspike_names <- colnames(PCR_spike[,1:S_star])
    OTUnames <- c(OTUnames, OTUspike_names)
  }
  
  y <- array(NA, dim = c(sum(M_site) + emptyTubes, max(K), S + S_star),
             dimnames = list(NULL, NULL, OTUnames))
  for (k in 1:max(K)) {
    y[,k,1:S] <- as.matrix(PCR_table[,1:S + (k-1)*S])
    y[,k,S + seq_len(S_star)] <- as.matrix(PCR_spike[,seq_len(S_star) + (k-1)*S_star])
  }
  
  if(is.null(X_w0)){
    namesCovW <- NULL
    X_w <- matrix(0, nrow = length(samples), ncol = 0)
  } else {
    namesCovW <- colnames(X_w0)[-1]
    X_w <- as.matrix(X_w0[match(samples, X_w0[,1]), -1])
    X_w <- scale(X_w)
  }
  
  if(is.null(X_z0)){
    namesCovZ <- NULL
    X_z <- matrix(0, nrow = length(sites), ncol = 0)
  } else {
    namesCovZ <- colnames(X_z0)[-1]
    X_z <- as.matrix(X_z0[match(sites, X_z0[,1]), -1])
    X_z <- scale(X_z)
  }
  
  data <- list("y" = y,
               "infos" = infos,
               "M_site" = M_site,
               "K" = K,
               "emptyTubes" = emptyTubes,
               "S_star" = S_star,
               "spikedSample" = spikedSample,
               "v_spikes" = v_spikes,
               "X_w" = X_w,
               "X_z" = X_z,
               "namesCovZ" = namesCovZ,
               "namesCovW" = namesCovW,
               "sites" = sites)
  
}

fitModel <- function(data,
                     priors = NULL,
                     jointSpecies,
                     paramsUpdate = list(updateAll = T,
                                         params = NULL,
                                         correct = NULL,
                                         trueParams = NULL),
                     MCMCparams,
                     idConstraints = list(beta0equal0 = F,
                                          betathetaequal0 = T),
                     paramsToSave = list(
                       "lambda" = T,
                       "beta0" = T,
                       "mu" = T,
                       "beta_z" = T,
                       "logz" = T,
                       "u" = T,
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
  
  # read data
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
    sites <- data$sites
    infos <- data$infos
    
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
    
    # standard priors
    {
      a_theta0 <- 1
      b_theta0 <- 50
      
      a_p11 <- 20
      b_p11 <- 1
      
      a_p10 <- 1
      b_p10 <- 100
      
      lambda_Y <- 1
      
      a_sigma <- 3
      b_sigma <- 2
      
      mean_r <- 100
      sd_r <- 100
      
      mean_ntilde <- 100
      sd_ntilde <- 100
      
      sigma_beta <- 1
      
      sigma_mu <- 1
      sigma_gamma <- 1
      
      sigma_u <- 1
      
      sigma_beta_theta <- 1
      
      sigma_lambda <- 3
      
      lambda_prior <- apply(y, 3, function(x){
        log(mean(x[x > quantile(x[x > 0], probs = c(0.1), na.rm = T)], na.rm = T))
      })
      
      lambda_Y <- 1
      a_tau <- 1
      b_tau <- 1
      
      
    }
    
    if(!is.null(priors)){
      
      sapply(seq_along(priors), function(i){
        val <- priors[[i]]
        namevar <- names(priors)[i]
        
        eval(parse(text = paste0(namevar," <<- ",val)))
      })
      
    }
    
    if(jointSpecies){
      
      Tau_Priors <- list("lambda_Y" = lambda_Y)
    } else {
      
      Tau_Priors <- list("a_tau" = a_tau,
                         "b_tau" = rep(1, S))
    }
    
  }
  
  # params to update
  {
    beta0equal0 = idConstraints$beta0equal0
    betathetaequal0 = idConstraints$betathetaequal0
    
    if(paramsUpdate$updateAll){
      
      updateLambda_CP <- T; correctLambda <- F
      updateLambda_NP <- T; correctLambda <- F
      
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
      updateLambda_NP <- paramsToUpdate$lambda; correctLambda <- !updateLambda_CP
      
      updateL <- paramsToUpdate$logz; correctL <- !updateL
      updateMu <- paramsToUpdate$mu; correctMu <- !updateMu
      
      updateV <- paramsToUpdate$v; correctV <- !updateV
      updateU <- paramsToUpdate$u; correctU <- !updateU
      updateUV <- paramsToUpdate$uv; #correctU <- !updateUV; correctV <- !updateUV
      
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
  
  # mcmc params
  {
    nchain <- MCMCparams$nchain
    nburn <- MCMCparams$nburn
    niter <- MCMCparams$niter
    nthin <- MCMCparams$nthin
    iterToAdapt <- Inf#MCMCparams$iterToAdapt  
  }
  
  # output
  {
    
    if(paramsToSave$lambda){
      lambda_output <- array(NA, c(nchain, niter, S + S_star))
    } else {
      lambda_output <- NULL
    }
    
    if(paramsToSave$beta0){
      beta0_output <- array(NA, dim = c(nchain, niter, S))
    } else {
      beta0_output <- NULL
    }
    
    if(paramsToSave$mu){
      mu_output <- array(NA, dim = c(nchain, niter, S))
    } else {
      mu_output <- NULL
    }
    
    if(paramsToSave$beta_z){
      beta_z_output <- array(NA, dim = c(nchain,  niter, ncov_z, S))
    } else {
      beta_z_output <- NULL
    }
    
    if(paramsToSave$logz){
      logz_output <- array(NA, dim = c(nchain, niter, n, S))
    } else {
      logz_output <- NULL
    }
    
    if(paramsToSave$u){
      u_output <- array(NA, c(nchain,  niter, sum(M_site) + emptyTubes, max(K)))
    } else {
      u_output <- NULL
    }
    
    if(paramsToSave$v){
      v_output <- array(NA, dim = c(nchain, niter, sum(M_site), S + S_star))
    } else {
      v_output <- NULL
    }
    
    if(paramsToSave$beta_w){
      beta_w_output <- array(NA, c(nchain, niter, ncov_w, S))
    } else {
      beta_w_output <- NULL
    }
    
    if(paramsToSave$Tau){
      if(jointSpecies){
        Tau_output <- array(NA, dim = c(nchain, niter, S, S))
      } else {
        Tau_output <- array(NA, c(nchain, niter, S))
      }
    } else {
      Tau_output <- NULL
    }
    
    if(paramsToSave$sigma){
      sigma_output <- array(NA, c(nchain, niter, S))
    } else {
      sigma_output <- NULL
    }
    
    if(paramsToSave$beta_theta){
      beta_theta_output <- array(NA, dim = c(nchain, niter, S, 2 + ncov_w))
    } else {
      beta_theta_output <- NULL
    }
    
    if(paramsToSave$theta){
      theta_output <- array(NA, dim = c(nchain, niter, sum(M_site), S))
    } else {
      theta_output <- NULL
    }
    
    if(paramsToSave$csi){
      csi_output <- array(NA, dim = c(nchain, niter, S))
    } else {
      csi_output <- NULL
    }
    
    if(paramsToSave$delta){
      delta_output <- array(NA, dim = c(nchain, niter, sum(M_site) + emptyTubes, S + S_star))
    } else {
      delta_output <- NULL
    }
    
    if(paramsToSave$gamma){
      gamma_output <- array(NA, dim = c(nchain, niter, sum(M_site) + emptyTubes, S + S_star))
    } else {
      gamma_output <- NULL
    }
    
    if(paramsToSave$r){
      r_output <- array(NA, dim = c(nchain, niter, S + S_star))
    } else {
      r_output <- NULL
    }
    
    if(paramsToSave$p11){
      p_11_output <- array(NA, dim = c(nchain, niter, S + S_star))
    } else {
      p_11_output <- NULL
    }
    
    if(paramsToSave$p10){
      p_10_output <- array(NA, dim = c(nchain, niter, S + S_star))
    } else {
      p_10_output <- NULL
    }
    
    if(paramsToSave$mutilde){
      mutilde_output <- array(NA, c(nchain, niter))
    } else {
      mutilde_output <- NULL
    }
    
    if(paramsToSave$c_imk){
      cimk_output <- array(NA, dim = c(nchain, niter, sum(M_site) + emptyTubes, max(K), S + S_star))
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
      eta_output <- array(NA, dim = c(nchain, niter, sum(M_site) + emptyTubes, max(K), S + S_star))
    } else {
      eta_output <- NULL
    }
    
    # latent variables means
    {
      theta11_mean <- array(0, dim = c(nchain, sum(M_site), S))
      delta_mean <- array(0, dim = c(nchain, sum(M_site) + emptyTubes, S))
      gamma_mean <- array(0, dim = c(nchain, sum(M_site) + emptyTubes, S))
      c_imk1_mean <- array(0, dim = c(nchain, sum(M_site) + emptyTubes, max(K), S))
      c_imk2_mean <- array(0, dim = c(nchain, sum(M_site) + emptyTubes, max(K), S))
    }
    
    
  }
  
  for (chain in 1:nchain) {
    
    # chain output
    {
      if(paramsToSave$lambda){
        lambda_output_iter <- matrix(NA, nrow = niter, ncol = S + S_star)
      } else {
        lambda_output_iter <- NULL
      }
      
      if(paramsToSave$beta0){
        beta0_output_iter <- array(NA, dim = c(niter, S))
      } else {
        beta0_output_iter <- NULL
      }
      
      if(paramsToSave$mu){
        mu_output_iter <- array(NA, dim = c(niter, S))
      } else {
        mu_output_iter <- NULL
      }
      
      if(paramsToSave$beta_z){
        beta_z_output_iter <- array(NA, dim = c(niter, ncov_z, S))
      } else {
        beta_z_output_iter <- NULL
      }
      
      if(paramsToSave$logz){
        logz_output_iter <- array(NA, dim = c(niter, n, S))
      } else {
        logz_output_iter <- NULL
      }
      
      if(paramsToSave$u){
        u_output_iter <- array(NA, dim = c(niter, sum(M_site) + emptyTubes, max(K)))
      } else {
        u_output_iter <- NULL
      }
      
      if(paramsToSave$v){
        v_output_iter <- array(NA, dim = c(niter, sum(M_site), S + S_star))
      } else {
        v_output_iter <- NULL
      }
      
      if(paramsToSave$beta_w){
        beta_w_output_iter <- array(NA, c(niter, ncov_w, S))
      } else {
        beta_w_output_iter <- NULL
      }
      
      if(paramsToSave$Tau){
        if(jointSpecies){
          Tau_output_iter <- array(NA, dim = c(niter, S, S))
        } else {
          Tau_output_iter <- matrix(NA, niter, S)
        }
      } else {
        Tau_output_iter <- NULL
      }
      
      if(paramsToSave$sigma){
        sigma_output_iter <- matrix(NA, niter, S)
      } else {
        sigma_output_iter <- NULL
      }
      
      if(paramsToSave$beta_theta){
        beta_theta_output_iter <- array(NA, dim = c(niter, S, 2 + ncov_w))
      } else {
        beta_theta_output_iter <- NULL
      }
      
      if(paramsToSave$theta){
        theta_output_iter <- array(NA, dim = c(niter, sum(M_site), S))
      } else {
        theta_output_iter <- NULL
      }
      
      if(paramsToSave$csi){
        csi_output_iter <- array(NA, dim = c(niter, S))
      } else {
        csi_output_iter <- NULL
      }
      
      if(paramsToSave$delta){
        delta_output_iter <- array(NA, dim = c(niter, sum(M_site) + emptyTubes, S + S_star))
      } else {
        delta_output_iter <- NULL
      }
      
      if(paramsToSave$gamma){
        gamma_output_iter <- array(NA, dim = c(niter, sum(M_site) + emptyTubes, S + S_star))
      } else {
        gamma_output_iter <- NULL
      }
      
      if(paramsToSave$r){
        r_output_iter <- matrix(NA, niter, S + S_star)
      } else {
        r_output_iter <- NULL
      }
      
      if(paramsToSave$p11){
        p_11_output_iter <- array(NA, dim = c(niter, S + S_star))
      } else {
        p_11_output_iter <- NULL
      }
      
      if(paramsToSave$p10){
        p_10_output_iter <- array(NA, dim = c(niter, S + S_star))
      } else {
        p_10_output_iter <- NULL
      }
      
      if(paramsToSave$mutilde){
        mutilde_output_iter <- rep(NA, niter)
      } else {
        mutilde_output_iter <- NULL
      }
      
      if(paramsToSave$c_imk){
        cimk_output_iter <- array(NA, dim = c(niter, sum(M_site) + emptyTubes, max(K), S  + S_star))
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
        eta_output_iter <- array(NA, dim = c(niter, sum(M_site), max(K), S + S_star))
      } else {
        eta_output_iter <- NULL
      }
      
     
    }
    
    # starting values
    {
      trueParams <- paramsUpdate$trueParams
      
      if(correctLambda){
        
        lambda <- trueParams$lambda_true
        
      } else
      {
        lambda <- lambda_prior
        
      }
      
      if(correctL){
        
        logz <- trueParams$logz_true
        
      } else
      {
        logz <- matrix(0, n, S)
      }
      
      if(correctMu){
        mu <- trueParams$mu_true
      } else
      {
        mu <- rep(0, S)
      }
      
      if(correctDeltaGammaC){
        c_imk <- trueParams$c_imk_true
        
        delta <- trueParams$delta_true
        
        gamma <- trueParams$gamma_true
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
      
      if(correctV){
        v <- trueParams$v_true
      } else
      {
        v <- matrix(NA, sum(M_site) + emptyTubes, S)
        for (i in 1:n) {
          for (m in 1:M_site[i]) {
            for (j in 1:S) {
              if(delta[m + sum(M_site[seq_len(i-1)]),j] == 1 | 
                 gamma[m + sum(M_site[seq_len(i-1)]),j] == 1){
                v[m + sum(M_site[seq_len(i-1)]),j] <- 0
              } 
            }
          }
        }
        if(S_star > 0){
          v <- cbind(v, v_spikes)
        }
      }
      
      if(correctU){
        u <- trueParams$u_true
      } else
      {
        u <- matrix(0, sum(M_site) + emptyTubes, max(K))
      }
      
      if(correctR){
        r_nb <- trueParams$r_true
      } else
      {
        r_nb <- rpois(S + S_star, 100)
      }
      
      if(correctBeta_w){
        beta_w <- trueParams$beta_w_true
      } else
      {
        beta_w <- matrix(0, ncov_w, S)
      }
      
      if(correctBeta_z){
        beta0 <- trueParams$beta0_true
        beta_z <- trueParams$beta_z_true
      } else
      {
        beta0 <- rep(0, S)
        beta_z <- matrix(0, nrow = ncov_z, ncol = S)
      }
      
      if(correctLambda0){
        mu0 <- trueParams$mu0_true
        n0 <- trueParams$n0_true
        pi0 <- trueParams$pi0_true
      } else
      {
        mu0 <- 2
        n0 <- 1
        pi0 <- .99
      }
      
      if(correctLambdaTilde){
        mu_tilde <- trueParams$mu_tilde_true
        n_tilde <- 1000
      } else
      {
        mu_tilde <- 100
        n_tilde <- 100000
      }
      
      if(correctSigma){
        sigma <- trueParams$sigma_true
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
          
          Tau_params <- list("Omega" = trueParams$Tau_true,
                             "lambdasq" = matrix(1, nrow = S, ncol = S),
                             "tausq" = 1,
                             "nu" = matrix(1, nrow = S, ncol = S),
                             "csi" = 1,
                             "Sigma" = solve(Tau_true))
          
        } else {
          Tau_params <- list("tau" = trueParams$tau_true)
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
        p_11 <- trueParams$p11_true
      } else
      {
        p_11 <- rep(0.9, S + S_star)
      }
      
      if(correctP10){
        p_10 <- trueParams$p10_true
      } else
      {
        p_10 <- rep(0.01, S + S_star)
      }
      
      if(correctBetaTheta){
        beta_theta <- trueParams$beta_theta_true
        theta11 <- trueParams$theta11_true
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
        theta10 <- trueParams$theta10_true
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
      
      if(iter <= nburn){
        print(paste0("Chain = ",chain," - Burn-in Iteration = ",iter))
      } else {
        print(paste0("Chain = ",chain," - Iteration = ",iter - nburn))
      }
      
      # LAMBDA ----------------------------------------------
      
      if(updateLambda_CP){
        
        if(beta0equal0){
          if(jointSpecies){
            tau <- sqrt(diag(Tau_params$Sigma))
          } else {
            tau <- Tau_params$tau
          }
          list_lambda <- update_lambda_CP_beta0(beta0, beta_z, X_z, logz, 
                                                mu, lambda, v, u, lambda_ijk, r_nb,
                                                c_imk, delta, gamma, X_w, beta_theta, 
                                                M_site, tau, sigma_mu,
                                                lambda_prior, sigma_lambda,
                                                S_star, emptyTubes)
        } else {
          list_lambda <- update_lambda_CP(beta0, beta_z, logz, 
                                          mu, lambda, v, u, lambda_ijk, r_nb,
                                          c_imk, delta, gamma, X_w, beta_theta, 
                                          M_site, sigma_beta, sigma_mu,
                                          lambda_prior, sigma_lambda,
                                          S_star, emptyTubes)
          
        }
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
        
        # lambda <- update_lambda_NP(lambda_ijk, c_imk, mu,
        #                            r_nb, v, u,
        #                            lambda_prior, sigma_lambda)
        
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
          Tau <- Tau_params$Sigma
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
          list_beta_z <- update_betaz_CP_corr(beta0, beta_z, logz, Tau, X_z, sigma_beta, !beta0equal0)
        } else {
          tau <- Tau_params$tau
          list_beta_z <- update_betaz_CP(beta0, beta_z, logz, tau, X_z, sigma_beta, !beta0equal0)
        }
        
        beta0 <- as.vector(list_beta_z$beta0)
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
      
      if(updateTau){
        Tau_params <- update_Tau(X_z, logz, beta0, beta_z,
                                 Tau_params, Tau_Priors,
                                 jointSpecies)  
      }
      
      
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
                                mean_r, sd_r,
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
        list_lambda_tilde <- update_lambda_tilde_NB(y, c_imk, mu_tilde, n_tilde, mean_ntilde,
                                                    sd_ntilde, sd_mu0 = 1, sd_n0 = 1)
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
        
        b_theta11 <- c(1,rep(0, ncov_w))
        B_theta11 <- diag(sigma_beta_theta, nrow = as.numeric(!betathetaequal0) + 1 + ncov_w)
        if(!betathetaequal0){
          b_theta11 <- c(0, b_theta11)
        }
        
        list_beta_theta <- update_betatheta11_cpp(logz,
                                                  beta_theta,
                                                  theta11,
                                                  delta[1:sum(M_site),],
                                                  X_w, M_site,
                                                  b_theta11,
                                                  B_theta11,
                                                  !betathetaequal0,
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
          lambda_output_iter[trueIter,] <- lambda
        } 
        
        if(paramsToSave$beta0){
          beta0_output_iter[trueIter,] <- beta0
        } 
        
        if(paramsToSave$mu){
          mu_output_iter[trueIter,] <- mu
        } 
        
        if(paramsToSave$beta_z){
          beta_z_output_iter[trueIter,,] <- beta_z
        } 
        
        if(paramsToSave$logz){
          logz_output_iter[trueIter,,] <- logz
        } 
        
        if(paramsToSave$u){
          u_output_iter[trueIter,,] <- u
        } 
        
        if(paramsToSave$v){
          v_output_iter[trueIter,,] <- v
        } 
        
        if(paramsToSave$beta_w){
          beta_w_output_iter[trueIter,,] <- beta_w
        }  
        
        if(paramsToSave$Tau){
          if(jointSpecies){
            Tau_output_iter[trueIter,,] <- Tau_params$Sigma
          } else {
            Tau_output_iter[trueIter,] <- tau
          }
        } 
        
        if(paramsToSave$sigma){
          sigma_output_iter[trueIter,] <- sigma
        } 
        
        if(paramsToSave$beta_theta){
          beta_theta_output_iter[trueIter,,] <- beta_theta
        } 
        
        if(paramsToSave$theta){
          theta_output_iter[trueIter,,] <- theta11
        } 
        
        if(paramsToSave$csi){
          csi_output_iter[trueIter,] <- theta10
        } 
        
        if(paramsToSave$delta){
          delta_output_iter[trueIter,,] <- delta
        } 
        
        if(paramsToSave$gamma){
          gamma_output_iter[trueIter,,] <- gamma
        }
        
        if(paramsToSave$r){
          r_output_iter[trueIter,] <- r_nb
        }
        
        if(paramsToSave$p11){
          p_11_output_iter[trueIter,] <- p_11
        } 
        
        if(paramsToSave$p10){
          p_10_output_iter[trueIter,] <- p_10
        } 
        
        if(paramsToSave$mutilde){
          mutilde_output_iter[trueIter] <- mu_tilde
        } 
        
        if(paramsToSave$c_imk){
          cimk_output_iter[trueIter,,,] <- c_imk
        } 
        
        if(paramsToSave$mu0){
          mu0_output_iter[trueIter] <- mu0
        } 
        
        if(paramsToSave$n0){
          n0_output_iter[trueIter] <- n0
        } 
        
        if(paramsToSave$pi0){
          pi0_output_iter[trueIter] <- pi0
        } 
        
        if(paramsToSave$eta){
          eta_output_iter[trueIter,,,] <- lambda_ijk
        } 
        
        theta11_mean[chain,,] <- theta11_mean[chain,,] +
          (theta11 / niter) 
        delta_mean[chain,,] <- delta_mean[chain,,] + 
          (delta[,1:S] / niter)
        gamma_mean[chain,,] <- gamma_mean[chain,,] + 
          (gamma[,1:S] / niter)
        c_imk1_mean[chain,,,] <- (c_imk1_mean[chain,,,]) + 
          (c_imk[,,1:S] == 1) / niter
        c_imk2_mean[chain,,,] <- (c_imk2_mean[chain,,,]) + 
          (c_imk[,,1:S] == 2) / niter
         
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
    "eta_output" = eta_output,
    "theta11_mean" = theta11_mean,
    "delta_mean" = delta_mean,
    "gamma_mean" = gamma_mean,
    "c_imk1_mean" = c_imk1_mean,
    "c_imk2_mean" = c_imk2_mean
  )
  
  data_infos = list(
    "y" = y,
    "infos" = data$infos,
    "sites" = sites,
    "OTUnames" = dimnames(y)[[3]],
    "jointSpecies" = jointSpecies,
    "namesCovZ" = data$namesCovZ,
    "namesCovW" = data$namesCovW
  )
  
  output <- list("params_output" = params_output,
                 "data_infos" = data_infos)
}
