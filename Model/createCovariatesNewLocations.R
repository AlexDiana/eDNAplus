# create prediction grid

x_min <- min(X_s[,1]) - .01
x_max <- max(X_s[,1]) + .01
y_min <- min(X_s[,2]) - .01
y_max <- max(X_s[,2]) + .01

X_s_star <- as.matrix(expand.grid(seq(x_min, x_max, length.out = 20),
                                  seq(y_min, y_max, length.out = 20)))

n_star <- nrow(X_s_star)

l_gp <- .1

Sigma_XsXs <- K2(X_s, X_s, 1, l_gp)
Sigma_XsXstar <- K2(X_s, X_s_star, 1, l_gp)

if(ncov_z > 0){
  X_z_star <- t(Sigma_XsXstar) %*% solve(Sigma_XsXs) %*% X_z[,1]
  
  ggplot(data = NULL, aes(x = X_s[,1],
                          y = X_s[,2], 
                          color = X_z[,1])) + geom_point()
  
  ggplot() + 
    # geom_point(data = NULL, aes(x = X_s[,1],
    #                                      y = X_s[,2], 
    #                                      color = X_z[,1]), size = 3) + 
    geom_point(data = NULL, aes(x = X_s_star[,1],
                                y = X_s_star[,2], 
                                color = X_z_star), 
               shape = 15, size = 5)
  
  # X_z_star
  
} else {
  X_z_star <- matrix(0, nrow = n_star, ncol = 0)
}
