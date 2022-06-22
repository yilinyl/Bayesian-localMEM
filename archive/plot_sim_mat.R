# Visualize similarity matrix

require(ggplot2)
require(ggpubr)
require(reshape2)
# source("./utils.R")


#---------------------------------------
# partitions <- basket_part(x=x, n=n, a0=1, b0=1, 
#                           max_cl = length(x))
# 
# clusters <- basket_cluster(partitions, a0=1, b0=1, 
#                            x=x, n=n)

plot_sim_mat <- function(clusters, ...) {
  member <- clusters$member
  n <- clusters$n
  x <- clusters$x
  
  up_tri <- round(get_upper_tri(clusters$Smat), 2)
  mel_sym <- melt(up_tri, na.rm = T)
  info <- cbind.data.frame(Var2 = LETTERS[1:6], n = n, x = x)
  mel_sym <- merge(mel_sym, info, by = "Var2")
  
  ggheatmap <- ggplot()+
    geom_tile(data = mel_sym, 
              aes(Var2, Var1, fill = value), color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(0,1), space = "Lab", 
                         name="") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 0, vjust = 25, 
                                     size = 12, hjust = 0.5))+
    coord_fixed() + 
    geom_text(data = mel_sym, aes(Var2, Var1, label = value), 
              color = "black", size = 4) +
    geom_text(data = mel_sym, aes(Var2, -0.2, label = n),
              color = "black", size = 4, vjust = -0.5) +
    geom_text(data = mel_sym, aes(Var2, -0.6, label = x),
              color = "black", size = 4, vjust = -0.5) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          legend.justification = c(1, 0),
          legend.position = c(0.6, 0.7),
          legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", 
                                 title.hjust = 0.5)) +
    annotation_custom(grob = text_grob(expression(N[b]), face = "italic"),
                      xmin = 0.5, xmax = 0.65, ymin = -0.1, ymax = -0.1) +
    annotation_custom(grob = text_grob(expression(x[b])), xmin = 0.5, xmax = 0.65, ymin = -0.5, ymax = -0.5) 
  return(ggheatmap)
}


# plot_sim_mat(clusters)

