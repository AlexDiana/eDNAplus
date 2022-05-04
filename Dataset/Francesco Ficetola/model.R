library(MASS); library(Rcpp); library(RcppArmadillo)
setwd("C:/Users/alexd/Dropbox/R Folder/PostDoc/Francesco Ficetola")
sourceCpp("code.cpp")

# sim data ----------

Y <- 45
S <- 30
M <- 2
K <- 4
P <- 7

M_y <- rep(M, Y)

library(gtools)

beta_true <- matrix(rnorm((P + 1) * S), nrow = P + 1, ncol = S)
beta_true[sample(1:((P+1)*S), 120)] <- 0
X <- cbind(1,matrix(runif(Y * P), nrow = Y, ncol = P))
N_imk <- array(sample(1:2000, size = Y * K * M, replace = TRUE), dim = c(Y, M, K))

gamma_true <- exp(X %*% beta_true)

p_jik <- array(NA, dim = c(Y, M, S))
for (i in 1:Y) {
  p_jik[i,,] <- rdirichlet(M, gamma_true[i,])
}

y_jimk <- array(NA, dim = c(Y, M, K, S))

for (i in 1:Y) {
  for (m in 1:M) {
    for (k in 1:K) {
      y_jimk[i,m,k,] <- rmultinom(1, N_imk[i,m,k], prob = p_jik[i,m,])
    }
  }
}

# real data ---------------------------------------------------------------

load("C:/Users/alexd/Dropbox/R Folder/PostDoc/Francesco Ficetola/francescodata.rda")

P <- ncol(X)
S <- dim(y_jimk)[4]

X <- as.matrix(cbind(1, scale(X)))

# functions ---------

logbeta_fun_r <- function(alpha){
  
  sum(lgamma(alpha))- lgamma(sum(alpha))
  
}

# loglikelihood_single <- function(gamma_tilde_ik){
#   
#   
#   
# }

loglikelihood_y_r_old <- function(beta, y_jim, N_imk){
  
  gamma <- exp(X %*% beta)
  
  loglikelihood <- 0
  
  for (i in 1:Y) {
    
    for (m in 1:M) {
      
      loglikelihood <- loglikelihood - logbeta_fun(gamma[i,])
      
      gamma_tilde_ik <- y_jim[i,m,] + gamma[i,]
      
      loglikelihood <- loglikelihood + logbeta_fun(gamma_tilde_ik)
      
    }
    
  }
  
  loglikelihood
}

logsingleterm1_r <- function(gamma_i){
  
  lgamma(sum(gamma_i)) - sum(lgamma(gamma_i))
}

logsingleterm2_r <- function(gamma_i, eta_j_im){
  
  lgamma(sum(gamma_i + eta_j_im)) - sum(lgamma(gamma_i + eta_j_im))
  
}

loglikelihood_y_r <- function(beta, eta_jim){
  
  gamma <- exp(X %*% beta)
  
  loglikelihood1 <- sum(sapply(1:Y, function(i){
    M_y[i] * (lgamma(sum(gamma[i,])) - sum(lgamma(gamma[i,])))
  }))
  
  loglikelihood2 <- sum(sapply(1:Y, function(i){
    sum(sapply(1:M_y[i], function(m){
      - lgamma(sum(gamma[i,] + eta_jim[i,m,])) + sum(lgamma(gamma[i,] + eta_jim[i,m,]))
    }))
  }))
  
  loglikelihood1 + loglikelihood2
  
}

# mcmc ------------

eta_jim <- apply(y_jimk,c(1,2,4),sum)

# initialize params
{
  optim_result <- nlm(f = function(betaVec){
    beta <- matrix(betaVec, nrow = P + 1, ncol = S)
    -(loglikelihood_y_optim(beta, X, eta_jim, M_y))
  }, p = rep(0, (P+1)*S))
  
  beta_mle <- matrix(optim_result$estimate, nrow = P + 1, ncol = S)
  
  hessian_mle <- matrix(NA, nrow = (P + 1) * S, ncol = (P + 1) * S)
  for(p1 in 1:(P+1)){
    for (p2 in 1:(P+1)) {
      print(paste0("P1 = ",p1," - P2 = ",p2))
      for (j1 in 1:S) {
        for (j2 in 1:S) {
          hessian_mle[p1 + (j1 - 1) * (P+1), p2 + (j2 - 1) * (P+1)] <- hessian_loglikelihood_y(p1 - 1, j1 - 1,
                                                                                           p2 - 1, j2 - 1,
                                                                                           beta_mle, X, eta_jim, M_y)
        }  
      }
    }
  }
  
  Sigma_cov <- solve(-hessian_mle)
  Sigma_cov <- Sigma_cov
  cholSigma_P <- array(NA, dim = c(S, P+1, P+1))
  rooti_P <- array(NA, dim = c(S, P+1, P+1))
  for (i in 1:S) {
    Sigma_cov_i <- Sigma_cov[1:(P+1) + (i - 1)*(P+1), 1:(P+1) + (i - 1)*(P+1)]
    cholSigma_P[i,,] <- t(chol(Sigma_cov_i))
    rooti_P[i,,] <- t(solve(chol(Sigma_cov_i)))
  }
  # rooti <- t(solve(chol(Sigma_cov)))
  # cholSigma <- t(chol(Sigma_cov))
}

beta <- beta_mle
# oldLogLik <- loglikelihood_y_optim(beta, X, eta_jim, M_y) 

nburn <- 0
niter <- 100000

iterAfterAdapting <- 100
beta_output <- array(NA, dim = c(niter, P + 1, S))

params_values <- array(NA, dim = c(niter, P + 1, S))

sumMatrices <- array(0, dim = c(S, P + 1, P + 1))
sumMeans <- matrix(0, nrow = S, ncol =  P + 1)

epsilon_params <- 0.2
beta_proposal <- .05

current_loglikelihood <- loglikelihood_y_optim(beta, X, eta_jim, M_y) 

for (iter in 1:(niter+nburn)) {
  
  if(iter %% 100 == 0){
    print(iter)  
  }

  # SPECIES-WISE ADAPTIVE MH -----
  
  for (j in 1:S) {
    
    if(iter > iterAfterAdapting){
      # Sigma_n <- cov(params_values[1:(iter-1),,j])
      Sigma_n <- sumMatrices[j,,] / (iter) - ((sumMeans[j,]) %*% t(sumMeans[j,])) / (iter^2)
      Sigma_proposal <- (1 - beta_proposal) * (2.38) * Sigma_n / 3 +
        beta_proposal * diag(epsilon_params, nrow = P + 1)
    } else {
      Sigma_proposal <- diag(epsilon_params, nrow = P + 1) 
    }
    
    beta_star <- beta
    beta_star_i_vec <- mvrnormArmaQuick(as.vector(beta[,j]), Sigma_proposal)
    beta_star[,j] <- beta_star_i_vec
    
    loglikelihood_current <- current_loglikelihood#loglikelihood_y_optim(beta, X, eta_jim, M_y)
    loglikelihood_star <- loglikelihood_y_optim(beta_star, X, eta_jim, M_y) 
    exp(loglikelihood_star - loglikelihood_current)
    
    # logproposal_current <- sum(sapply(1:nrow(Sigma_cov), function(i){
    #   dnorm(beta[i], beta_mle[i], diag(Sigma_cov)[j])
    # }))
    
    mhratio <- exp(loglikelihood_star - loglikelihood_current)
    mhratio# print(mhratio)
    
    if(runif(1) < mhratio){
      beta <- beta_star
      current_loglikelihood <- loglikelihood_star
      # print("accepted")
    }
      
  }
  
  params_values[iter,,] <- beta
  
  for (j in 1:S) {
      
    sumMatrices[j,,] <- sumMatrices[j,,] + beta[,j] %*% t(beta[,j])
    sumMeans[j,] <- sumMeans[j,] + beta[,j]
    
  }

  
  # SPECIES-WISE LAPLACE MH ----------
  
  # for (i in 1:S) {
  #   
  #   beta_star <- beta
  #   beta_star_i_vec <- mvrnormArmaQuick(as.vector(beta_mle[,i]), cholSigma_P[i,,])
  #   beta_star[,i] <- beta_star_i_vec
  #   
  #   loglikelihood_current <- loglikelihood_y_optim(beta, X, eta_jim, M_y) 
  #   loglikelihood_star <- loglikelihood_y_optim(beta_star, X, eta_jim, M_y) 
  #   exp(loglikelihood_star - loglikelihood_current)
  #   
  #   logproposal_current <- dmvnorm_cpp_fast(as.vector(beta[,i]), as.vector(beta_mle[,i]), rooti_P[i,,], returnLog = T)
  #   logproposal_star <- dmvnorm_cpp_fast(as.vector(beta_star[,i]), as.vector(beta_mle[,i]), rooti_P[i,,], returnLog = T)
  #   exp(logproposal_current - logproposal_star)
  #   
  #   # logproposal_current <- sum(sapply(1:nrow(Sigma_cov), function(i){
  #   #   dnorm(beta[i], beta_mle[i], diag(Sigma_cov)[j])
  #   # }))
  #   
  #   mhratio <- exp(logproposal_current + loglikelihood_star - loglikelihood_current - logproposal_star)
  #   mhratio# print(mhratio)
  #   
  #   if(runif(1) < mhratio){
  #     beta <- beta_star
  #     # print("accepted")
  #   }
  #     
  # }
 
  # LAPLACE APPROXIMATION COMPLETE  ------------
  
  # beta_star_vec <- mvrnormArmaQuick(as.vector(beta_mle), cholSigma)
  # beta_star <- matrix(beta_star_vec, nrow = P + 1, ncol = S, byrow = F)
  # 
  # loglikelihood_current <- loglikelihood_y_optim(beta, X, eta_jim, M_y) 
  # loglikelihood_star <- loglikelihood_y_optim(beta_star, X, eta_jim, M_y) 
  # exp(loglikelihood_star - loglikelihood_current)
  # 
  # logproposal_current <- dmvnorm_cpp_fast(as.vector(beta), as.vector(beta_mle), rooti, returnLog = T)
  # logproposal_star <- dmvnorm_cpp_fast(as.vector(beta_star), as.vector(beta_mle), rooti, returnLog = T)
  # exp(logproposal_current - logproposal_star)
  # 
  # mhratio <- exp(logproposal_current + loglikelihood_star - loglikelihood_current - logproposal_star)
  # mhratio# print(mhratio)
  # 
  # if(runif(1) < mhratio){
  #   beta <- beta_star
  #   # print("accepted")
  # }
  
  # COMPONENT-WISE MH UPDATE ----------
  
  #  for (p in 1:(P+1)) {
  #   
  #   for (s in 1:S) {
  #     
  #     beta_star <- beta
  #     beta_star[p,s] <- rnorm(1, beta[p,s], sd = 0.05)
  #     
  #     loglikelihood_current <-  loglikelihood_y(beta, X, y_jim, M_y) 
  #     loglikelihood_star <-  loglikelihood_y(beta_star, X, y_jim, M_y) 
  #     exp(loglikelihood_star - loglikelihood_current)
  #     if(runif(1) < exp(loglikelihood_star - loglikelihood_current)){
  #       beta <- beta_star
  #     }
  #     
  #   }
  #   
  # }
  
  if(iter > nburn){
    beta_output[iter - nburn,,] <- beta
  }
  # beta_output[accepted,,] <- beta
  
}

library(ggplot2)
library(coda)

j <- 1
p <- 5
qplot(1:niter, beta_output[,p,j], geom = "line") + geom_hline(aes(yintercept = beta_true[p,j])) #+ ylim(c(-2,2))
  
  
qplot(1:niter, beta_output[,p,j], geom = "line") #+ geom_hline(aes(yintercept = beta_mle[p,j]))

CI_intervals <- apply(beta_output, c(2,3), function(x){
  quantile(x, probs = c(0.025,0.975))
})

ESS <- apply(beta_output, c(2,3), function(x){
  effectiveSize(as.mcmc(x))
})

CI_intervals[,2,]

ggplot() + geom_point(data = NULL, aes(x = 1:niter, y = beta_output[,p,j])) + geom_hline(aes(yintercept = beta_true[p,j])) + 
  ylim(c(min(beta_true), max(beta_true)))

# speed ----------

loglikelihood_y(beta_true, X, eta_jim, M_y)
loglikelihood_y_old(beta_true, X, eta_jim, M_y)

library(microbenchmark)

microbenchmark(logbeta_fun(alpha),
               logbeta_fun_r(alpha))

microbenchmark({
  loglikelihood_y(beta, X, eta_jim, M_y)
},{
  loglikelihood_y_r(beta, eta_jim)
}, times = 10)

microbenchmark({
  beta_star_vec <- mvrnorm(1, as.vector(beta_mle), Sigma_cov)
  beta_star <- matrix(beta_star_vec, nrow = P + 1, ncol = S, byrow = F)
},{
  
  loglikelihood_current <- oldLogLik
  loglikelihood_star <- loglikelihood_y_optim(beta_star, X, eta_jim, M_y) 
  # exp(loglikelihood_star - loglikelihood_current)
  
},{
  logproposal_current <- dmvnorm_cpp_fast(as.vector(beta), as.vector(beta_mle), rooti, returnLog = T)
  logproposal_star <- dmvnorm_cpp_fast(as.vector(beta_star), as.vector(beta_mle), rooti, returnLog = T)
  # exp(logproposal_current - logproposal_star)
}, times = 5)

# LAPLACE APPROXIMATION -----------


CIs <- array(NA, dim = c(P + 1, S, 3))
for (p in 1:(P+1)) {
  for (j in 1:S) {
    CIs[p,j,1:2] <- beta_mle[p,j] + c(-1,1) * Sigma_cov[p + (j - 1) * (P + 1), p + (j - 1) * (P + 1)]
    CIs[p,j,3] <- beta_true[p,j]
  }
}

CI_intervals <- apply(beta_output, c(2,3), function(x){
  quantile(x, probs = c(0.025,0.975))
})

CI_intervals <- CI_intervals[,-1,]

dimnames(CI_intervals)[[2]] <- colnames(X)
dimnames(CI_intervals)[[3]] <- speciesNames

# ggplot(data = NULL, aes(x = 1:(P+1), y = CIs[,j,3], ymin = CI_intervals[1,,j], ymax = CI_intervals[2,,j])) + 
  # geom_point(color = "red", size = 2) + geom_errorbar()


j <- 1

plotSpecies <- ggplot(data = NULL, aes(x = dimnames(CI_intervals)[[2]], ymin = CI_intervals[1,,j], ymax = CI_intervals[2,,j])) + 
   geom_errorbar() + ylim(c(-6,6)) + geom_hline(aes(yintercept = 0), size = 1, color = "red") +  
  theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,5,0), size = 14, face = "bold"),
        panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(fill = NA, colour = "grey20"),
        panel.grid.major = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(colour = "grey92", size = 0.25), 
        strip.background = element_rect(fill = "grey85", colour = "grey20"), 
        legend.key = element_rect(fill = "white", colour = NA)) + 
  ggtitle(speciesNames[j]) + xlab("Covariates")

ggsave(filename = paste(speciesNames[j],".jpg"), plotSpecies, device = "jpg")
j <- j + 1

# significant variables
apply(CI_intervals, c(2,3), function(x){
  (x[2] < 0) | (x[1] > 0)
})
