library(dplyr); library(tidyverse); library(reshape2)

dimnames(beta_z_output)[[3]] <- colnames(OTU)
dimnames(beta_z_output)[[2]] <- c("Eleveation","Distance")

beta_CI <- apply(beta_z_output, c(2,3), function(x){
  quantile(x, probs = c(.025,.25,.5,.75,.975))
})

beta_CI <- beta_CI[,1,] %>% t 
beta_CI <- as.data.frame(beta_CI)
beta_CI$Species <- rownames(beta_CI)

# colnames(beta_CI) <- c("CI","Species","Value")

beta_CI <- beta_CI %>% 
  mutate(Species = as.factor(Species)) %>% 
mutate(Species = fct_reorder(Species, `50%`)) 

beta_CI$Family <- sapply(1:nrow(beta_CI), function(i){
  strsplit(as.character(beta_CI$Species[i]),",")[[1]][1]
})


ggplot(data = beta_CI, aes(x = Species,
                           color = Family)) + 
                        geom_errorbar(aes(ymin = `2.5%`,
                                          ymax = `97.5%`), show.legend = F, 
                                      size = .15, fatten = 1) +
                        geom_point(aes(y = `50%`), show.legend = F) + 
  geom_linerange(aes(ymin = `25%`, ymax = `75%`), size = 1, show.legend = FALSE) +
  facet_grid(~ Family, scales = "free_x", space = "free") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                 axis.title = ggplot2::element_text(size = 20, face = "bold"),
                 axis.text.y = ggplot2::element_text(size = 13, face = "bold", angle = 90),
                 axis.text.x = ggplot2::element_text(size = 9, face = "bold", angle = 90),
                 panel.grid.major = ggplot2::element_line(colour="grey", size=0.15),
                 panel.background = ggplot2::element_rect(fill = "white", color = "black")) +
  geom_hline(aes(yintercept = 0), color = "red") + 
  # scale_x_continuous(breaks = subsetSpecies, name = "Species",
                     # labels = colnames(OTU)[subsetSpecies]) +
  scale_y_continuous(name = "Elevation")
