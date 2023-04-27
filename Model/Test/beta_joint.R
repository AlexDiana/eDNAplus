library(MASS); library(Matrix)

n <- 10
p <- 2
S <- 10

sigma_beta <- 2

# Sigma_n <- diag(1, nrow = n)
# Sigma_S <- diag(1, nrow = S)
Sigma_n <- rWishart(1, n + 1000, diag(1, nrow = n) / (n + 1000))[,,1]
Sigma_S <- rWishart(1, S + 100, diag(1, nrow = S) / (S + 10))[,,1]

invSigma_S <- solve(Sigma_S)
invSigma_n <- solve(Sigma_n)

X <- matrix(runif(n * p), n, p)
B <- matrix(runif(p * S), p, S)

Y <- rmtrnorm(X %*% B, Sigma_n, Sigma_S)

Yt <- as.vector(Y) 

# Sigma <- Matrix(Sigma)
# X <- Matrix(X)

# old 
{
  Id_n <- Matrix(diag(1, nrow = n))
  Id_s <- Matrix(diag(1, nrow = S))
  
  # invSigma_tilde <- kronecker(Id_n, Sigma)
  # 
  # X_tilde <- kronecker(X, Id_s)
  # 
  # Lambda_beta <- t(X_tilde) %*% invSigma_tilde %*% X_tilde
  # 
  # mu_beta <- solve(Lambda_beta) %*% t(X_tilde) %*% invSigma_tilde %*% Yt
  # 
  # cbind(mu_beta, as.vector(t(B)))
  
  invSigma_S <- Matrix(invSigma_S)
  invSigma_n <- Matrix(invSigma_n)
  X <- Matrix(X)
  
  # invSigma_tilde <- kronecker(invSigma_S, invSigma_n)
  
  X_tilde <-  kronecker(Id_s, X)
  
  # Lambda_beta <- Matrix::t(X_tilde) %*% invSigma_tilde %*% X_tilde + diag(1 / sigma_beta^2, S * ncol(X))
  Lambda_beta <- kronecker(invSigma_S, t(X) %*% invSigma_n %*% X) + diag(1 / sigma_beta^2, S * ncol(X)) 
  
  # all.equal(as.matrix(Lambda_beta), as.matrix(Lambda_beta2))
  
  # term1 <- invSigma_tilde %*% Yt
  term1 <- as.vector(invSigma_n %*% Y %*% invSigma_S)
  mu_beta <- Matrix::t(X_tilde) %*% term1
  
  # mean_beta <- solve(Lambda_beta) %*% mu_beta
  
  chol_Lambda_beta <- t(chol(as.matrix(Lambda_beta)))
  u <- backsolve(chol_Lambda_beta, mu_beta, upper.tri = F)
  z <- rnorm(ncol(X) * S)
  z <- rep(0, ncol(X) * S)
  betavec <- backsolve(t(chol_Lambda_beta), z + u, upper.tri = T)
  
  # betavec <- Matrix::t(mvrnorm(1, as.vector(mean_beta), as.matrix(solve(Lambda_beta))))
  matrix(betavec, p, S, byrow = F)
}

# new 
{
  inv_chol_Sigma_n <- solve(t(chol(Sigma_n)))
  inv_chol_Sigma_S <- solve(t(chol(Sigma_S)))
  
  U <- inv_chol_Sigma_n
  V <- inv_chol_Sigma_S
  
  Ltilde <- U %*% Y %*% t(V)
  
  # Xtilde <- kronecker(t(V), U %*% X)
  # XtildetXtilde <- t(Xtilde) %*% Xtilde
  
  XtildetXtilde <- kronecker(t(V) %*% V, t(X) %*% t(U) %*% U %*% X)
  
  # XtildeLtilde <- t(Xtilde) %*% as.vector(Ltilde)
  
  XtildeLtilde <- as.vector(t(X) %*% t(U) %*% Ltilde %*% V)
  
  # XtildeLtilde <- as.vector(t(X) %*% solve(Sigma_n) %*% Y %*% solve(Sigma_S))
  
  Lambda_beta <- XtildetXtilde + diag(1 / sigma_beta^2, S * ncol(X)) 
  mu_beta <- XtildeLtilde
  
  chol_Lambda_beta <- t(chol(as.matrix(Lambda_beta)))
  u <- backsolve(chol_Lambda_beta, mu_beta, upper.tri = F)
  z <- rnorm(ncol(X) * S)
  z <- rep(0, ncol(X) * S)
  betavec <- backsolve(t(chol_Lambda_beta), z + u, upper.tri = T)
  
  # betavec <- Matrix::t(mvrnorm(1, as.vector(mean_beta), as.matrix(solve(Lambda_beta))))
  matrix(betavec, p, S, byrow = F)
    
}

B
as.vector(B)

#
Xt <- as.matrix(X_tilde)
(lm(as.vector(Y) ~ Xt - 1))
as.vector(B)


library(microbenchmark)
microbenchmark({
  invSigma_tilde <- kronecker(invSigma_S, invSigma_n)
},{
  X_tilde <-  kronecker(Id_s, X)
},{  
  Lambda_beta <- kronecker(invSigma_S, t(X) %*% invSigma_n %*% X) + diag(1 / sigma_beta^2, S * ncol(X)) 
},{
  term1 <- as.vector(invSigma_n %*% Y %*% invSigma_S)
  mu_beta <- Matrix::t(X_tilde) %*% term1
},{
  chol_Lambda_beta <- t(chol(as.matrix(Lambda_beta)))
  u <- backsolve(chol_Lambda_beta, mu_beta, upper.tri = F)
  z <- rnorm(ncol(X) * S)
  z <- rep(0, ncol(X) * S)
  betavec <- backsolve(t(chol_Lambda_beta), z + u, upper.tri = T)
  
}, times = 5)
