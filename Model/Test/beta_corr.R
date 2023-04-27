library(MASS); library(Matrix)

n <- 5000
p <- 2
S <- 3

df_t <- S + 100
Sigma <- rWishart(1, df_t, diag(1, nrow = S) / df_t)[,,1]

X <- matrix(runif(n * p), n, p)
B <- matrix(runif(p * S), p, S)

Y <- matrix(NA, n, S)
for (i in 1:n) {
  Y[i,] <- mvrnorm(1, mu = X[i,] %*% B, Sigma)
}

Yt <- as.vector(t(Y)) 

Sigma <- Matrix(Sigma)
X <- Matrix(X)

Id_n <- Matrix(diag(1, nrow = n))
Id_s <- Matrix(diag(1, nrow = S))

invSigma_tilde <- kronecker(Id_n, Sigma)

X_tilde <- kronecker(X, Id_s)

Lambda_beta <- t(X_tilde) %*% invSigma_tilde %*% X_tilde

mu_beta <- solve(Lambda_beta) %*% t(X_tilde) %*% invSigma_tilde %*% Yt

cbind(mu_beta, as.vector(t(B)))

#
Xt <- as.matrix(X_tilde)
(lm(Yt ~ Xt - 1))



