#' Title
#'
#' @return
#' @export
#'
#' @examples
plot_affinity<-function(bulk_knn){
  bulk_knn<-reshape2::melt(bulk_knn)
  ggheatmap <- ggplot(bulk_knn, aes(Var2, Var1, fill = value), vjust = 17.5, hjust = 5)+
    geom_tile()+
    scale_fill_gradient2(low = "white", high = "blue",midpoint = 0.008, limit = c(0,0.016), space = "Lab",name="Correlation") +
    labs(x = "",y = "", title = "kNN affinity")+
    theme_minimal()+
    scale_x_continuous(position = "top",expand = c(0.019,0,0,0))+
    scale_y_continuous(trans = 'reverse',expand = c(0.019,0,0,0))+
    coord_cartesian(xlim = c(2,2982),ylim=c(0,2980))+
    theme(panel.background=element_rect(fill='transparent',color='gray'),
          panel.grid=element_blank(),
          axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          plot.title = element_text(size = 32,face = "bold", vjust = 10.5, hjust = 0.5)
          # panel.margin = unit(2, "cm")
    )
  ggheatmap +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor=element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line=element_line(size=9,colour="#10B3AE"),
      plot.margin = unit(c(6, 2,4,  2), "lines"),
      legend.position = 'none'
    )
}
