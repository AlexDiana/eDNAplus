var_biom <- function(S_star, S, n, M, K, tau, sigma, sigma_u, sigma_y){
  
  sigma_ratio <- sigma_u^2 / sigma_y^2
  term2 <- 1 + (sigma_ratio / (sigma_ratio * S_star + 1))
  (1 / M) * (sigma^2 + (sigma_y^2 / K) * term2)
  
}

var_cov <- function(S_star, S, n, M, K, tau, sigma, sigma_u, sigma_y){
  
  sigma_ratio <- sigma_u^2 / sigma_y^2
  term1 <- (1 / (n - 1)) * (tau^2 + (1 / M) * (sigma^2 + (sigma_y^2 / K) ))
  term2 <- 1 + (sigma_u^2 / (sigma_y^2 + (M * tau^2 + sigma^2) * K * (sigma_ratio * S_star + 1) +
                               sigma_u^2 * (S + S_star + 1)))
  
  term1 * term2
  
}

S <- 10
n <- 1000
sigma_y <- .1

betavar_true <- sapply(1:nrow(data_betabias), function(i){
  x <- data_betabias[i,]
  var_cov(x$S_star, S, n, x$M, x$K, x$tau, x$sigma, x$sigma_u, sigma_y) /
    var_cov(0, S, n, x$M, x$K, x$tau, x$sigma, x$sigma_u, sigma_y)
})

biomvar_true <- sapply(1:nrow(data_betabias), function(i){
  x <- data_betabias[i,]
  var_biom(x$S_star, S, n, x$M, x$K, x$tau, x$sigma, x$sigma_u, sigma_y) / 
    var_biom(0, S, n, x$M, x$K, x$tau, x$sigma, x$sigma_u, sigma_y)
})

data_betavars$vals <- betavar_true
data_biomdiff_vars$vals <- biomvar_true

# PLOTS ------

setwd(here("Model/Spike-in Sims/Plots True"))

generalPlot <- function(data, title){
  
  data <- data[data$S_star <= 5,]
  data$tau <- factor(data$tau)
  data$K <- factor(data$K)
  data$sigma <- factor(data$sigma)
  
  currentPlot <- ggplot(data = data, 
                        aes(x = S_star, 
                            y = vals,
                            shape = K,
                            color = tau,
                            linetype = sigma)) + 
    geom_line() + geom_point(size = 3) + facet_grid(cols = vars(M)) +
    coord_cartesian(ylim = c(min(data$vals), 
                             max(data$vals))) +
    scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
    scale_y_continuous(name = "") + 
    # scale_color_manual(values = c("#1C3FFD","#020873","#FF2D00","#900B0A")) +
    scale_color_manual(values = c("#1C3FFD","#FF2D00")) +
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(size = 14, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black")) + 
    ggtitle(title)
  
  currentPlot
}

(biomdiff_vars_plot <- generalPlot(data_biomdiff_vars,
                                   "Relative standard error of estimates of differences of biomasses"))
# ggsave(filename = paste0("biomdiff_vars",".jpg"), biomdiff_vars_plot)

(beta_vars_plot <- generalPlot(data_betavars,
                               "Relative standard error of estimates of covariate coefficients"))
# ggsave(filename = paste0("beta_vars",".jpg"), beta_vars_plot)