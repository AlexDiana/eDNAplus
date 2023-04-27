
# CREATE DATA ------

library(here); library(ggplot2)
files <- list.files(here("Model/Spike-in Sims"))

files <- files[grep(".rda",files)]

{
  
  data_betabias <- data.frame(vals = c(),
                              S_star = c(),
                              M = c(),
                              K = c(),
                              tau = c(),
                              sigma = c(),
                              sigma_u = c(),
                              Settings = c(),
                              VarSettings = c())
  
  data_betavars <- data.frame(vals = c(),
                              S_star = c(),
                              M = c(),
                              K = c(),
                              tau = c(),
                              sigma = c(),
                              sigma_u = c(),
                              Settings = c(),
                              VarSettings = c())
  
  data_biomdiff_bias <- data.frame(vals = c(),
                                   S_star = c(),
                                   M = c(),
                                   K = c(),
                                   tau = c(),
                                   sigma = c(),
                                   sigma_u = c(),
                                   Settings = c(),
                                   VarSettings = c())
  
  data_biomdiff_vars <- data.frame(vals = c(),
                                   S_star = c(),
                                   M = c(),
                                   K = c(),
                                   tau = c(),
                                   sigma = c(),
                                   sigma_u = c(),
                                   Settings = c(),
                                   VarSettings = c())
  
  data_ubias <- data.frame(vals = c(),
                            S_star = c(),
                            M = c(),
                            K = c(),
                            tau = c(),
                            sigma = c(),
                            sigma_u = c(),
                           Settings = c(),
                           VarSettings = c())
  
  
  data_uvars <- data.frame(vals = c(),
                            S_star = c(),
                            M = c(),
                            K = c(),
                            tau = c(),
                            sigma = c(),
                            sigma_u = c(),
                           Settings = c(),
                           VarSettings = c())
}

for (file in files) {
  
  load(here("Model/Spike-in Sims",file))
  
  if(length(beta_biases) == 5){
    S_star_grid <- c(0, 1, 2, 5, 10)  
  } else {
    S_star_grid <- c(0, 1, 2, 5)
  }
  
  
  settings <- paste0("tau=",simsSettings$tau,"_",
                    "sigma=",simsSettings$sigma,"_",
                    "sigma_u=",simsSettings$sigma_u,"_",
                    "M=",simsSettings$M,"_",
                    "n=",simsSettings$n,"_",
                    "K=",simsSettings$K)
  
  var_settings <- paste0(
    # expression(tau),
    "Tau",
    "=",simsSettings$tau," - ",
                         # "Sigma=",simsSettings$sigma," - ",
                         "K=",simsSettings$K)
  
  # Beta bias
  {
    data_betabias_current <- data.frame(vals = beta_biases / max(beta_biases),
                                        S_star = S_star_grid,
                                        M = simsSettings$M,
                                        K = simsSettings$K,
                                        tau = simsSettings$tau,
                                        sigma = simsSettings$sigma,
                                        sigma_u = simsSettings$sigma_u,
                                        Settings = settings,
                                        VarSettings = var_settings)
    
    data_betabias <- rbind(data_betabias, data_betabias_current)
    
  }

  # Beta vars
  {
    data_betavars_current <- data.frame(vals = beta_vars / max(beta_vars),
                                        S_star = S_star_grid,
                                        M = simsSettings$M,
                                        K = simsSettings$K,
                                        tau = simsSettings$tau,
                                        sigma = simsSettings$sigma,
                                        sigma_u = simsSettings$sigma_u,
                                        Settings = settings,
                                        VarSettings = var_settings)
    
    data_betavars <- rbind(data_betavars, data_betavars_current)
    
  }
  
  # Biomdiff bias
  {
    data_biombias_current <- data.frame(vals = ldiff_biases / max(ldiff_biases),
                                        S_star = S_star_grid,
                                        M = simsSettings$M,
                                        K = simsSettings$K,
                                        tau = simsSettings$tau,
                                        sigma = simsSettings$sigma,
                                        sigma_u = simsSettings$sigma_u,
                                        Settings = settings,
                                        VarSettings = var_settings)
    
    data_biomdiff_bias <- rbind(data_biomdiff_bias, data_biombias_current)
    
  }
  
  # Biomdiff vars
  {
    data_biomvars_current <- data.frame(vals = ldiff_vars / max(ldiff_vars),
                                        S_star = S_star_grid,
                                        M = simsSettings$M,
                                        K = simsSettings$K,
                                        tau = simsSettings$tau,
                                        sigma = simsSettings$sigma,
                                        sigma_u = simsSettings$sigma_u,
                                        Settings = settings,
                                        VarSettings = var_settings)
    
    data_biomdiff_vars <- rbind(data_biomdiff_vars, data_biomvars_current)
    
  }
  
  # u bias
  {
    data_ubias_current <- data.frame(vals = u_biases,
                                     S_star = S_star_grid,
                                     M = simsSettings$M,
                                     K = simsSettings$K,
                                     tau = simsSettings$tau,
                                     sigma = simsSettings$sigma,
                                     sigma_u = simsSettings$sigma_u,
                                     Settings = settings,
                                     VarSettings = var_settings)
    
    data_ubias <- rbind(data_ubias, data_ubias_current)
    
  }
  
  # u vars
  {
    data_uvars_current <- data.frame(vals = u_vars,
                                     S_star = S_star_grid,
                                     M = simsSettings$M,
                                     K = simsSettings$K,
                                     tau = simsSettings$tau,
                                     sigma = simsSettings$sigma,
                                     sigma_u = simsSettings$sigma_u,
                                     Settings = settings,
                                     VarSettings = var_settings)
    
    data_uvars <- rbind(data_uvars, data_uvars_current)
    
  }
  
  
}

# PLOTS ------

setwd(here("Model/Spike-in Sims/Plots"))

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

(u_MAE_plot <- generalPlot(data_ubias,
                          "Mean absolute error of estimate of pipeline noise"))
ggsave(filename = paste0("u_MAE",".jpg"), u_MAE_plot)

(u_vars_plot <- generalPlot(data_uvars,
                            "Standard error of estimate of pipeline noise"))
ggsave(filename = paste0("u_vars",".jpg"), u_vars_plot)

(biomdiff_MAE_plot <- generalPlot(data_biomdiff_bias,
            "Mean absolute error of differences of biomasses"))
ggsave(filename = paste0("biomdiff_MAE",".jpg"), biomdiff_MAE_plot)

(biomdiff_vars_plot <- generalPlot(data_biomdiff_vars,
            "Relative standard error of estimates of differences of biomasses"))
ggsave(filename = paste0("biomdiff_vars",".jpg"), biomdiff_vars_plot)

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
