library(here); library(ggplot2)
library(reshape2)
setwd(here("Model/Spike-in Sims/Plots"))

# BIOMASS DIFFERNECES ---

load("~/eDNAPlus/Model/Spike-in Sims/results_new.rda")

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

var_ratio_all <- apply(var_ratio_all, c(1,3), mean)
mae_ratio_all <- apply(mae_ratio_all, c(1,3), mean)

# mae
data <- cbind(settings, mae_ratio_all)

data <- melt(data, id=c("M","K","tau","sigma"))

data$variable <- S_star_grid[data$variable]
colnames(data)[5] <- "S_star"
data$S_star <- as.numeric(data$S_star)

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
    scale_y_continuous(name = "Relative Error") + 
    scale_color_manual(values= c("#1C3FFD","#FF2D00")) +
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          strip.text.x = element_text(size = 15, face = "bold", angle = 0),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(size = 18, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18),
          panel.background = element_rect(fill = "white", color = "black")) + 
    ggtitle(title)
  
  currentPlot
}

(beta_MAE_plot <- generalPlot(data, ""))
ggsave(filename = paste0("biom_MAE",".jpg"), beta_MAE_plot)

# var
data <- cbind(settings, var_ratio_all)

data <- melt(data, id=c("M","K","tau","sigma"))

data$variable <- S_star_grid[data$variable]
colnames(data)[5] <- "S_star"
data$S_star <- as.numeric(data$S_star)

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
          strip.text.x = element_text(size = 15, face = "bold", angle = 0),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(size = 18, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18),
          panel.background = element_rect(fill = "white", color = "black")) + 
    ggtitle(title)
  
  currentPlot
}

(beta_VAR_plot <- generalPlot(data, ""))
ggsave(filename = paste0("biom_VAR",".jpg"), beta_VAR_plot)

# COVARIATE COEFFICIENTS --------

load("~/eDNAPlus/Model/Spike-in Sims/results_new_coeff2.rda")

nsims <- 10
settings <- expand.grid(M = c(1,2,3),
                        K = c(1,3),
                        tau = c(.2, .5),
                        sigma = c(.2, .5))

var_ratio_all <- apply(var_ratio_all, c(1,2), function(x){
  x / x[length(x)]
})
var_ratio_all <- aperm(var_ratio_all, c(2,3,1))
mae_ratio_all <- apply(mae_ratio_all, c(1,2), function(x){
  x / x[length(x)]
})
mae_ratio_all <- aperm(mae_ratio_all, c(2,3,1))

var_ratio_all <- apply(var_ratio_all, c(1,3), mean)
mae_ratio_all <- apply(mae_ratio_all, c(1,3), mean)

# mae
data <- cbind(settings, mae_ratio_all)

data <- melt(data, id=c("M","K","tau","sigma"))

data$variable <- S_star_grid[data$variable]
colnames(data)[5] <- "S_star"
data$S_star <- as.numeric(data$S_star)

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
    scale_y_continuous(name = "Relative Error") + 
    scale_color_manual(values= c("#1C3FFD","#FF2D00")) +
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          strip.text.x = element_text(size = 15, face = "bold", angle = 0),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(size = 18, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18),
          panel.background = element_rect(fill = "white", color = "black")) + 
    ggtitle(title)
  
  currentPlot
}

(beta_MAE_plot <- generalPlot(data, ""))
ggsave(filename = paste0("coeff_MAE",".jpg"), beta_MAE_plot)

# var
data <- cbind(settings, var_ratio_all)

data <- melt(data, id=c("M","K","tau","sigma"))

data$variable <- S_star_grid[data$variable]
colnames(data)[5] <- "S_star"
data$S_star <- as.numeric(data$S_star)

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
          strip.text.x = element_text(size = 15, face = "bold", angle = 0),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(size = 18, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18),
          panel.background = element_rect(fill = "white", color = "black")) + 
    ggtitle(title)
  
  currentPlot
}

(beta_VAR_plot <- generalPlot(data, ""))
ggsave(filename = paste0("coeff_VAR",".jpg"), beta_VAR_plot)
