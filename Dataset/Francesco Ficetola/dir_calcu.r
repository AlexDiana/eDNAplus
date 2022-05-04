ddirichlet(gamma_true[1,] / (sum(gamma_true[1,])), gamma_true[1,])

ddirichlet(gamma_true[1,] / (sum(gamma_true[1,])), gamma_true[1,] / (sum(gamma_true[1,])))

beta_fun <- function(alpha){
  
  exp(sum(lgamma(alpha)) - lgamma(sum(alpha)))
  
}

alpha <- gamma_true[1,]

1 / beta_fun(gamma_true[1,])

prod((gamma_true[1,] / (sum(gamma_true[1,])))^(gamma_true[1,] -1)) / beta_fun(gamma_true[1,])

library(HMP)

eta_current <- eta_jim[,1,]

optimiz <- dirmult(eta_current, init = rep(1, S))

alpha <- optimiz$gamma
alpha

#

data <- data[, colSums(data) != 0, drop = FALSE]
alphap <- exp(X %*% beta_opt)[1,]
data <- eta_current[,,drop = F]
(ll <- sum(lgamma(rowSums(data) + 1) + lgamma(sum(alphap)) - 
            lgamma(sum(alphap) + rowSums(data))) + 
  sum(rowSums(lgamma(sweep(data, 
                           2, alphap, "+")) - lgamma(data + 1) - lgamma(t(replicate(nrow(data), 
                                                                                    alphap))))))

ddrmultinom(exp(beta_true)[1,], eta_current[1,]) + ddrmultinom(exp(beta_true)[1,], eta_current[2,])
ddrmultinom(exp(beta_opt)[1,], eta_current[1,])

ddrmultinom(alphap, data[1,]) + ddrmultinom(alphap, data[2,])

ddrmultinom <- function(alpha, eta_current_i){
  lfactorial(sum(eta_current_i)) + lgamma(sum(alpha)) -
    lgamma(sum(eta_current_i + alpha)) +
    sum(lgamma(eta_current_i + alpha) - lgamma(alpha)) - sum(lfactorial(eta_current_i)) 
}

sum(sapply(1:nrow(eta_current), function(j){
  ddrmultinom(exp(beta_true), eta_current[j,])
}))

sum(sapply(1:nrow(eta_current), function(j){
  ddrmultinom(exp(beta_opt), eta_current[j,])
}))

#

HMP::loglikDM(eta_current, exp(beta_true))
HMP::loglikDM(eta_current, exp(beta_opt))

ddrmultinom1 <- function(alpha){
  lgamma(sum(alpha)) - sum(lgamma(alpha))
}

ddrmultinom2 <- function(eta_ijm, alpha){
  sum(lgamma(alpha + eta_ijm)) - lgamma(sum(alpha + eta_ijm)) 
}

sum(sapply(1:nrow(eta_current), function(j){
  ddrmultinom1(exp(beta_true))
})) + sum(sapply(1:nrow(eta_current), function(j){
  ddrmultinom2(eta_current[j,], exp(beta_true))
}))

sum(sapply(1:nrow(eta_current), function(j){
  ddrmultinom1(exp(beta_opt))
})) + sum(sapply(1:nrow(eta_current), function(j){
  ddrmultinom2(eta_current[j,], exp(beta_opt))
}))

loglikelihood_y_r(beta_true, eta_jim) 
loglikelihood_y_r(beta_opt, eta_jim) 

gamma <- exp(X %*% beta_opt)[1,]

logsingleterm1(gamma)
ddirmultinom1(eta_current[1,], gamma)
