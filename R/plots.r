
# PLOT FUNCTION --------------

diagnosticsPlot <- function(chain_output){
  
  chain_output_all <-
    matrix(chain_output, nrow = nchain, ncol = niter)
  
  chain_output_long <- reshape2::melt(chain_output_all)
  
  ggplot2::ggplot(data = chain_output_long, ggplot2::aes(
    x = Var2,
    y = value,
    group = Var1,
    color = factor(Var1)
  )) + ggplot2::geom_line() +
    ggplot2::xlab("Iterations") + ggplot2::ylab("Value") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 17),
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 11, face = "bold"),
      axis.text.x = ggplot2::element_text(
        size = 11,
        face = "bold",
        hjust = 1
      ),
      axis.line = ggplot2::element_line(colour = "black", size = 0.15),
      # panel.grid.minor = element_line(colour="grey", size=0.15),
      panel.grid.major = ggplot2::element_line(colour = "grey", size = 0.15),
      panel.background = ggplot2::element_rect(fill = "white", color = "black"),
      legend.position = "none"
    )
  
}

