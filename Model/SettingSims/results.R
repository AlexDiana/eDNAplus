library(here); library(ggplot2)
setwd(here("Model/SettingSims"))


n_grid <- c(50, 100, 200)
M_grid <- c(1, 2, 3)
K_grid <- c(1, 2, 3)
nMK_grid <- expand.grid(n_grid, M_grid, K_grid)
settings_df <- as.data.frame(nMK_grid)
colnames(settings_df) <- c("n","M","K")

beta_vars_mean <- apply(beta_vars, 2, mean)
beta_biases_mean <- apply(beta_biases, 2, mean)

# PLOTS -------------------

generalPlot_beta <- function(vals, title, y_name){
  
  plot_data <- cbind(settings_df, data.frame(vals = vals))
  
  plot_data$K <- as.factor(plot_data$K)
  plot_data$M <- as.factor(plot_data$M)
  
  var_settings <- paste0(
    plot_data$M," - ",plot_data$K
  )
  
  plot_data <- cbind(plot_data, var_settings)
  
  ggplot(data = 
           plot_data,
         # plot_data,
         aes(x = n, 
             y = vals,
             alpha = M,
             group = factor(var_settings),
             linetype = K)) + 
    geom_line() + geom_point() + #facet_grid(cols = vars(K)) +
    coord_cartesian(ylim = c(0, max(plot_data$vals))) +
    scale_x_continuous(name = "n", breaks = n_grid) + 
    scale_y_continuous(name = y_name) + 
    scale_color_manual(values = c("#2962FF","#1510F0","#0003C7","#020873")) +
    scale_alpha_manual(values = c(.3,.65,1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 20, face = "bold"),
          axis.text.y = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(size = 18, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          legend.text= element_text(size=16),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black")) 
  # + ggtitle(title)
  
}

(beta_vars_plot <- generalPlot_beta(beta_vars_mean, 
                                    "Relative standard error of estimates of covariate coefficients",
                                    "SD"))
ggsave(filename = paste0("Beta_vars",".jpg"), beta_vars_plot)

(beta_MAE_plot <- generalPlot_beta(beta_biases_mean, 
                                   "Average absolute error of estimates of covariate coefficients",
                                   "MAE"))
ggsave(filename = paste0("Beta_MAE",".jpg"), beta_MAE_plot)

generalPlot_biomass <- function(vals, title){
  
  plot_data <- cbind(settings_df, data.frame(vals = vals))
  
  plot_data$K <- as.factor(plot_data$K)
  # plot_data$M <- as.factor(plot_data$M)
  
  var_settings <- paste0(
    plot_data$M," - ",plot_data$K
  )
  
  plot_data <- cbind(plot_data, var_settings)
  
  plot_data <- plot_data[plot_data$n == max(plot_data$n),]
  
  ggplot(data = 
           plot_data,
         aes(x = M, 
             y = vals,
             group = K,
             linetype = K)) + 
    geom_line() + geom_point() + #facet_grid(cols = vars(K)) +
    coord_cartesian(ylim = c(0, max(plot_data$vals))) +
    scale_x_continuous(name = "M", breaks = unique(plot_data$M)) + 
    scale_y_continuous(name = "SD") + 
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black")) + 
    ggtitle(title)
  
}

(Biom_Acc_plot <- generalPlot_biomass(ldiff_accuracy, "Brier score of biomasses differences"))
ggsave(filename = paste0("Biom_acc",".jpg"), Biom_Acc_plot)

(Biom_MAE_plot <- generalPlot_biomass(ldiff_biases, "Average absolute error of estimates of differences of biomasses"))
ggsave(filename = paste0("Biom_MAE",".jpg"), Biom_MAE_plot)

(Biom_Vars_plot <- generalPlot_biomass(ldiff_vars, "Relative standard error of estimates of differences of biomasses"))
ggsave(filename = paste0("BIom_Vars",".jpg"), Biom_Vars_plot)
