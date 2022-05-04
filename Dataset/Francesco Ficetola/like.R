vstar <- log(aa / bb - cc)
sigmav <- aa / (aa - bb*cc)^2

u_im <- u[1,idxReplicatePositive]

v_grid <- seq(2.14, 2.19, length.out = 100)
posterior_w <- seq_along(v_grid)
for (j in seq_along(v_grid)) {
  posterior_w[j] <- sum(dpois(PCR_counts, 
                              lambda = lambda * u_im * exp(v_grid[j]), 
                              log = T)) 
}

library(ggplot2)
qplot(v_grid, exp(posterior_w))
