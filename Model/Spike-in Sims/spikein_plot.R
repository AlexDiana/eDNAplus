library(ggplot2)

load("~/MetabarcodingModels/Model/Spike-in Sims/results_tau1_sigma02_n1000_m1_k1.rda")
beta_biases_tau1_sigma_02_sigma_u1 <- beta_biases
beta_vars_tau1_sigma_02_sigma_u1 <- beta_vars
ldiff_biases_tau1_sigma_02_sigma_u1 <- ldiff_biases
ldiff_vars_tau1_sigma_02_sigma_u1 <- ldiff_vars
u_biases_tau1_sigma_02_sigma_u1 <- u_biases
u_vars_tau1_sigma_02_sigma_u1 <- u_vars
load("~/MetabarcodingModels/Model/Spike-in Sims/results_tau1_sigma02_n1000_m1_k1_sigmau_2.rda")
beta_biases_tau1_sigma_02_sigma_u2 <- beta_biases
beta_vars_tau1_sigma_02_sigma_u2 <- beta_vars
ldiff_biases_tau1_sigma_02_sigma_u2 <- ldiff_biases
ldiff_vars_tau1_sigma_02_sigma_u2 <- ldiff_vars
u_biases_tau1_sigma_02_sigma_u2 <- u_biases
u_vars_tau1_sigma_02_sigma_u2 <- u_vars
load("~/MetabarcodingModels/Model/Spike-in Sims/results_tau1_sigma1_n1000_m1_k1.rda")
beta_biases_tau1_sigma_1_sigma_u2 <- beta_biases
beta_vars_tau1_sigma_1_sigma_u2 <- beta_vars
ldiff_biases_tau1_sigma_1_sigma_u2 <- ldiff_biases
ldiff_vars_tau1_sigma_1_sigma_u2 <- ldiff_vars
u_biases_tau1_sigma_1_sigma_u2 <- u_biases
u_vars_tau1_sigma_1_sigma_u2 <- u_vars

S_star_grid <- c(0,1,2,5,10)

# beta bias
{
  beta_biases_data <- data.frame(MAS = c(beta_biases_tau1_sigma_02_sigma_u1,
                                         beta_biases_tau1_sigma_02_sigma_u2,
                                         beta_biases_tau1_sigma_1_sigma_u2),
                                 Spikein = rep(S_star_grid, times = 3),
                                 Setting = rep(c("Setting 1", "Setting 2", "Setting 3"), each = 5)
  )
  
  
  # beta bias
  ggplot(data = beta_biases_data, aes(x = Spikein, 
                                      y = MAS,
                                      color = Setting)) + 
    geom_line() + geom_point() + 
    coord_cartesian(ylim = c(0, max(beta_biases))) + 
    scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black"))
  
}

# beta vars
{
  beta_vars_data <- data.frame(MAS = c(beta_vars_tau1_sigma_02_sigma_u1,
                                       beta_vars_tau1_sigma_02_sigma_u2,
                                       beta_vars_tau1_sigma_1_sigma_u2),
                                 Spikein = rep(S_star_grid, times = 3),
                                 Setting = rep(c("Setting 1", "Setting 2", "Setting 3"), each = 5)
  )
  
  
  # beta bias
  ggplot(data = beta_vars_data, aes(x = Spikein, 
                                      y = MAS,
                                      color = Setting)) + 
    geom_line() + geom_point() + 
    coord_cartesian(ylim = c(0, max(beta_vars_data$MAS))) + 
    scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black"))
  
}

# ldiff bias
{
  ldiff_biases_data <- data.frame(MAS = c(ldiff_biases_tau1_sigma_02_sigma_u1,
                                          ldiff_biases_tau1_sigma_02_sigma_u2,
                                         ldiff_biases_tau1_sigma_1_sigma_u2),
                                 Spikein = rep(S_star_grid, times = 3),
                                 Setting = rep(c("Setting 1", "Setting 2", "Setting 3"), each = 5)
  )
  
  
  # beta bias
  ggplot(data = ldiff_biases_data, aes(x = Spikein, 
                                      y = MAS,
                                      color = Setting)) + 
    geom_line() + geom_point() + 
    coord_cartesian(ylim = c(0, max(ldiff_biases_data$MAS))) + 
    scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black"))
  
}

# ldiff vars
{
  ldiff_vars_data <- data.frame(MAS = c(ldiff_vars_tau1_sigma_02_sigma_u1,
                                        ldiff_vars_tau1_sigma_02_sigma_u2,
                                        ldiff_vars_tau1_sigma_1_sigma_u2),
                                 Spikein = rep(S_star_grid, times = 3),
                                 Setting = rep(c("Setting 1", "Setting 2", "Setting 3"), each = 5)
  )
  
  
  # beta bias
  ggplot(data = ldiff_vars_data, aes(x = Spikein, 
                                      y = MAS,
                                      color = Setting)) + 
    geom_line() + geom_point() + 
    coord_cartesian(ylim = c(0, max(ldiff_vars_data$MAS))) + 
    scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black"))
  
}

# u vars
{
  u_vars_data <- data.frame(MAS = c(u_vars_tau1_sigma_02_sigma_u1,
                                    u_vars_tau1_sigma_02_sigma_u2,
                                    u_vars_tau1_sigma_1_sigma_u2),
                                 Spikein = rep(S_star_grid, times = 3),
                                 Setting = rep(c("Setting 1", "Setting 2", "Setting 3"), each = 5)
  )
  
  
  # u vars
  ggplot(data = u_vars_data, aes(x = Spikein, 
                                      y = MAS,
                                      color = Setting)) + 
    geom_line() + geom_point() + 
    coord_cartesian(ylim = c(0, max(u_vars_data$MAS))) + 
    scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black"))
  
}

# u biases
{
  u_biases_data <- data.frame(MAS = c(u_biases_tau1_sigma_02_sigma_u1,
                                    u_biases_tau1_sigma_02_sigma_u2,
                                    u_biases_tau1_sigma_1_sigma_u2),
                                 Spikein = rep(S_star_grid, times = 3),
                                 Setting = rep(c("Setting 1", "Setting 2", "Setting 3"), each = 5)
  )
  
  
  # u vars
  ggplot(data = u_biases_data, aes(x = Spikein, 
                                      y = MAS,
                                      color = Setting)) + 
    geom_line() + geom_point() + 
    coord_cartesian(ylim = c(0, max(u_biases_data$MAS))) + 
    scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
    theme(plot.title = element_text(hjust = 0.5, size = 17),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
          axis.line = element_line(colour="black", size=0.15),
          # panel.grid.minor = element_line(colour="grey", size=0.15),
          panel.grid.major = element_line(colour="grey", size=0.15),
          panel.background = element_rect(fill = "white", color = "black"))
  
}

#

# beta bias
ggplot(data = NULL, aes(x = S_star_grid, 
                        y = beta_biases)) + 
  geom_line() + geom_point() + 
  coord_cartesian(ylim = c(0, max(beta_biases))) + 
  scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
        axis.line = element_line(colour="black", size=0.15),
        # panel.grid.minor = element_line(colour="grey", size=0.15),
        panel.grid.major = element_line(colour="grey", size=0.15),
        panel.background = element_rect(fill = "white", color = "black"))

# beta bias
ggplot(data = NULL, aes(x = S_star_grid, 
                        y = beta_biases)) + 
  geom_line() + geom_point() + 
  coord_cartesian(ylim = c(0, max(beta_biases))) + 
  scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
        axis.line = element_line(colour="black", size=0.15),
        # panel.grid.minor = element_line(colour="grey", size=0.15),
        panel.grid.major = element_line(colour="grey", size=0.15),
        panel.background = element_rect(fill = "white", color = "black"))

# beta vars
ggplot(data = NULL, aes(x = S_star_grid, 
                        y = beta_vars)) + 
  geom_line() + geom_point() + 
  coord_cartesian(ylim = c(0, max(beta_vars))) + 
  scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
        axis.line = element_line(colour="black", size=0.15),
        # panel.grid.minor = element_line(colour="grey", size=0.15),
        panel.grid.major = element_line(colour="grey", size=0.15),
        panel.background = element_rect(fill = "white", color = "black"))

# l diff biases
ggplot(data = NULL, aes(x = S_star_grid, 
                        y = ldiff_biases)) + 
  geom_line() + geom_point() + 
  coord_cartesian(ylim = c(0, max(ldiff_biases))) + 
  scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
        axis.line = element_line(colour="black", size=0.15),
        # panel.grid.minor = element_line(colour="grey", size=0.15),
        panel.grid.major = element_line(colour="grey", size=0.15),
        panel.background = element_rect(fill = "white", color = "black"))

# l diff vars
ggplot(data = NULL, aes(x = S_star_grid, 
                        y = ldiff_vars)) + 
  geom_line() + geom_point() + 
  coord_cartesian(ylim = c(0, max(ldiff_vars))) + 
  scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
        axis.line = element_line(colour="black", size=0.15),
        # panel.grid.minor = element_line(colour="grey", size=0.15),
        panel.grid.major = element_line(colour="grey", size=0.15),
        panel.background = element_rect(fill = "white", color = "black"))

# u biases
ggplot(data = NULL, aes(x = S_star_grid, 
                        y = u_biases)) + 
  geom_line() + geom_point() + 
  coord_cartesian(ylim = c(0, max(u_biases))) + 
  scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
        axis.line = element_line(colour="black", size=0.15),
        # panel.grid.minor = element_line(colour="grey", size=0.15),
        panel.grid.major = element_line(colour="grey", size=0.15),
        panel.background = element_rect(fill = "white", color = "black"))

# u vars
ggplot(data = NULL, aes(x = S_star_grid, 
                        y = u_vars)) + 
  geom_line() + geom_point() + 
  coord_cartesian(ylim = c(0, max(u_vars))) + 
  scale_x_continuous(breaks = S_star_grid, name = "Spike-ins") + 
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold", angle = 0, hjust = 1),
        axis.line = element_line(colour="black", size=0.15),
        # panel.grid.minor = element_line(colour="grey", size=0.15),
        panel.grid.major = element_line(colour="grey", size=0.15),
        panel.background = element_rect(fill = "white", color = "black"))
