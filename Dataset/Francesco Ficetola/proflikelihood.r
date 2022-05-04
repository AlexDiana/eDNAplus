optim_result <- nlm(f = function(betaVec){
  beta <- matrix(betaVec, nrow = P + 1, ncol = S)
  -(loglikelihood_y_optim(beta, X, eta_jim, M_y))
}, p = rep(0, (P+1)*S))

beta_mle <- matrix(optim_result$estimate, nrow = P + 1, ncol = S)

# univariate profiles

p <- 5
j <- 2

x <- seq(beta_mle[p,j] - .1, beta_mle[p,j] + .1, by = .001)
y <- rep(NA, length(x))

beta <- beta_mle

for (i in 1:length(x)) {
  beta[p,j] <- x[i]
  y[i] <- loglikelihood_y_optim(beta, X, eta_jim, M_y)
}

y <- y - max(y)
qplot(x, exp(y))

beta_mle[p,j] + c(-1,1) * 1.96 * Sigma_cov[p + (j - 1) * (P + 1), p + (j - 1) * (P + 1)]


persp(x1, x2, exp(z))

# 


p1 <- 4; j1 <- 5
p2 <- 3; j2 <- 5

x1 <- seq(beta_mle[p1,j1] - .5, beta_mle[p1,j1] + .5, by = .04)
x2 <- seq(beta_mle[p2,j2] - .5, beta_mle[p2,j2] + .5, by = .04)

z <- matrix(NA, nrow = length(x1), ncol = length(x2))

beta <- beta_mle

for (i in 1:length(x1)) {
  print(i)
  for (j in 1:length(x2)) {
    beta[p1,j1] <- x1[i]
    beta[p2,j2] <- x2[j]
    z[i,j] <- loglikelihood_y_optim(beta, X, eta_jim, M_y)
  } 
}

z <- z - max(z)

persp(x1, x2, exp(z))
surface3d(x1, x2, exp(z), color = "grey")

# maximization --------

funToOptim <- function(betaVec){
  beta <- matrix(betaVec, nrow = P + 1, ncol = S)
  -(loglikelihood_y_optim(beta, X, eta_jim, M_y))
}

beta_start <- rep(0, (P+1)*S)
optim_result <- nlm(f = funToOptim, p = beta_start)

beta_opt <- matrix(optim_result$estimate, nrow = P + 1, ncol = S)

beta_opt
beta_true

#

grad_loglikelihood_y(1, 1, beta_opt, X, eta_jim, M_y)
grad_loglikelihood_y(1, 1, beta_true, X, eta_jim, M_y)

grad_loglikelihood_y(0, 1, beta, X, eta_jim, M_y)

gradient_loglikelihood <- matrix(NA, nrow = P + 1, ncol = S)
for(i in 1:(P+1)){
  for(j in 1:S){
    gradient_loglikelihood[i,j] <- grad_loglikelihood_y(i - 1, j -1, beta_opt, X, eta_jim, M_y)
  }
}

ddirichlet(p_jik[1,,], gamma[i,])
ddirichlet(p_jik[1,,], gamma_opt[i,])

# 

hessian <- matrix(NA, nrow = (P + 1) * S, ncol = (P + 1) * S)
for(p1 in 1:(P+1)){
  for (p2 in 1:(P+1)) {
    for (j1 in 1:S) {
      for (j2 in 1:S) {
        hessian[p1 + (j1 - 1) * (P+1), p2 + (j2 - 1) * (P+1)] <- hessian_loglikelihood_y(p1 - 1, j1 - 1,
                                                                                         p2 - 1, j2 - 1,
                                                                                         beta_opt, X, eta_jim, M_y)
      }  
    }
  }
}

p1 <- 1; j1 <- 1
p1 <- 1; j1 <- 2
hessian_loglikelihood_y(p1 - 1, j1 - 1,
                        p2 - 1, j2 - 1,
                        beta_opt, X, eta_jim, M_y)

hessian_loglikelihood_y(p2 - 1, j2 - 1,
                        p1 - 1, j1 - 1,
                        beta_opt, X, eta_jim, M_y)
