
# CREATE DATA ------

library(here); library(ggplot2)
load("~/eDNAPlus/Model/Spike-in Sims/results_new.rda")
# load("~/eDNAPlus/Model/Spike-in Sims/results_new_coeff.rda")
# load("~/eDNAPlus/Model/Spike-in Sims/results_new_coeff2.rda")

var_ratio_all <- apply(var_ratio_all, c(1,2), function(x){
  x / x[length(x)]
})
var_ratio_all <- aperm(var_ratio_all, c(2,3,1))
mae_ratio_all <- apply(mae_ratio_all, c(1,2), function(x){
  x / x[length(x)]
})
mae_ratio_all <- aperm(mae_ratio_all, c(2,3,1))


S_star_grid <- c(3, 2, 1, 0)

nsims <- 5
settings <- expand.grid(M = c(1,2,3),
                        K = c(1,3),
                        tau = c(.5, 1),
                        sigma = c(.5, 1))

# nsims <- 10
# settings <- expand.grid(M = c(1,2,3),
#                         K = c(1,3),
#                         tau = c(.2, .5),
#                         sigma = c(.2, .5))

var_ratio_all <- apply(var_ratio_all, c(1,3), mean)
mae_ratio_all <- apply(mae_ratio_all, c(1,3), mean)

data <- cbind(settings, mae_ratio_all)
# data <- cbind(settings, var_ratio_all)

library(reshape2)
data <- melt(data, id=c("M","K","tau","sigma"))

S_star_grid <- 3:0

data$variable <- S_star_grid[data$variable]
colnames(data)[5] <- "S_star"
data$S_star <- as.numeric(data$S_star)

# PLOTS ------

setwd(here("Model/Spike-in Sims/Plots"))

generalPlot <- function(data, title){
  
  # data <- data[data$S_star <= 5,]
  data$tau <- factor(data$tau)
  data$K <- factor(data$K)
  data$sigma <- factor(data$sigma)
  
  currentPlot <- ggplot(data = data, 
                        aes(x = S_star, 
                            y = value,
                            shape = K,
                            color = tau,
                            linetype = sigma)) + 
    geom_line() + geom_point(size = 3) + facet_grid(cols = vars(M)) +
    coord_cartesian(ylim = c(
      # min(data$value),
      0,
      max(data$value))) +
    scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
    scale_y_continuous(name = "Relative Variance") + 
    scale_color_manual(values= c("#1C3FFD","#FF2D00")) +
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 16, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          panel.background = element_rect(fill = "white", color = "black")) + 
    ggtitle(title)
  
  currentPlot
}

(beta_MAE_plot <- generalPlot(data, ""))
ggsave(filename = paste0("biom_MAE",".jpg"), beta_MAE_plot)

(beta_VAR_plot <- generalPlot(data, ""))
ggsave(filename = paste0("biom_VAR",".jpg"), beta_VAR_plot)

(beta_MAE_plot <- generalPlot(data, ""))
ggsave(filename = paste0("coeff_MAE",".jpg"), beta_MAE_plot)

(beta_VAR_plot <- generalPlot(data, ""))
ggsave(filename = paste0("coeff_VAR",".jpg"), beta_VAR_plot)

# -----------

(beta_MAE_plot <- generalPlot(data_betabias,
                              "Mean absolute error of estimates of covariate coefficients"))
ggsave(filename = paste0("beta_MAE",".jpg"), beta_MAE_plot)

(beta_vars_plot <- generalPlot(data_betavars,
                               "Relative standard error of estimates of covariate coefficients"))
ggsave(filename = paste0("beta_vars",".jpg"), beta_vars_plot)

{
  u_MAE_plot <- ggplot(data = data_ubias, aes(x = S_star, 
                                              y = vals,
                                              color = VarSettings)) + 
    geom_line() + geom_point() + facet_grid(cols = vars(M)) +
    coord_cartesian(ylim = c(0, max(data_ubias$vals))) +
    scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
    scale_y_continuous(name = "Variance") + 
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black")) + 
    ggtitle("Mean absolute error of estimate of pipeline noise")
  
  ggsave(filename = paste0("u_MAE",".jpg"), u_MAE_plot)
}

{
  
  # u var
  u_vars_plot <- ggplot(data = data_uvars, aes(x = S_star, 
                                               y = vals,
                                               color = Settings)) + 
    geom_line() + geom_point() + facet_grid(cols = vars(M)) +
    coord_cartesian(ylim = c(0, max(data_uvars$vals))) +
    scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
    scale_y_continuous(name = "Variance") + 
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black")) + 
    ggtitle("Standard error of estimate of pipeline noise")
  
  ggsave(filename = paste0("u_vars",".jpg"), u_vars_plot)
}

{
  # biomdiff bias
  biomdiff_MAE_plot <- ggplot(data = data_biomdiff_bias, aes(x = S_star, 
                                                             y = vals,
                                                             color = VarSettings)) + 
    geom_line() + geom_point() + facet_grid(cols = vars(M)) +
    coord_cartesian(ylim = c(0, max(data_biomdiff_bias$vals))) +
    scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
    scale_y_continuous(name = "Variance") + 
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black")) + 
    ggtitle("Mean absolute error of differences of biomasses")
  
  ggsave(filename = paste0("biomdiff_MAE",".jpg"), biomdiff_MAE_plot)
}

{
  # biomdiff vars
  biomdiff_vars_plot <- ggplot(data = data_biomdiff_vars[data_biomdiff_vars$S_star <= 5,], 
                               aes(x = S_star, 
                                   y = vals,
                                   shape = factor(K),
                                   color = factor(tau),
                                   linetype = factor(sigma))) + 
    geom_line() + geom_point(size = 4) + facet_grid(cols = vars(M)) +
    coord_cartesian(ylim = c(min(data_biomdiff_vars$vals), 
                             max(data_biomdiff_vars$vals))) +
    scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
    # scale_color_manual(values = c("#1C3FFD","#020873","#FF2D00","#900B0A")) +
    scale_color_manual(values = c("#1C3FFD","#FF2D00")) +
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black")) + 
    ggtitle("Relative standard error of estimates of differences of biomasses")
  
  ggsave(filename = paste0("biomdiff_vars",".jpg"), biomdiff_vars_plot)
}

{
  # beta bias
  beta_MAE_plot <- ggplot(data = data_betabias, aes(x = S_star, 
                                                    y = vals,
                                                    color = VarSettings)) + 
    geom_line() + geom_point() + facet_grid(cols = vars(M)) +
    coord_cartesian(ylim = c(0, max(data_betabias$vals))) +
    scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
    scale_y_continuous(name = "Variance") + 
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black")) + 
    ggtitle("Mean absolute error of estimates of covariate coefficients")
  
  ggsave(filename = paste0("beta_MAE",".jpg"), beta_MAE_plot)
}

{
  # beta vars
  ( beta_vars_plot <- ggplot(data = data_betavars[data_betavars$S_star <= 5,], 
                             aes(x = S_star, 
                                 y = vals,
                                 shape = factor(K),
                                 color = VarSettings,
                                 linetype = factor(sigma))) + 
      geom_line() + geom_point(size = 4) + facet_grid(cols = vars(M)) +
      coord_cartesian(ylim = c(min(data_betavars$vals), 
                               max(data_betavars$vals))) +
      scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
      scale_color_manual(values = c("#1C3FFD","#020873","#FF2D00","#900B0A")) +
      scale_y_continuous(name = "SD") + 
      theme(plot.title = element_text(hjust = 0.5, size = 17),
            axis.title = element_text(size = 12, face = "bold"),
            axis.text.y = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
            axis.line = element_line(colour="black", size=0.15),
            # panel.grid.minor = element_line(colour="grey", size=0.15),
            panel.grid.major = element_line(colour="grey", size=0.15),
            panel.background = element_rect(fill = "white", color = "black")) + 
      ggtitle("Relative standard error of estimates of covariate coefficients"))
  
  ggsave(filename = paste0("beta_vars",".jpg"), beta_vars_plot)
}
