#' Title
#'
#' @return
#' @export
#'
#' @examples
plot_atlas<-function(atlasinf,spatialposition,genename){
  if(!genename%in%colnames(atlasinf))
  {
    string<-paste('This gene "',genename,'" don\'t exist, input appropriate name',sep='')
    stop(string)
  }
  spatialposition$inf<-atlasinf[,genename]

  dist<-range(spatialposition$xcoord)[2]-range(spatialposition$xcoord)[1]
  start1<-range(spatialposition$xcoord)[1]+0.3*dist
  end1<-range(spatialposition$xcoord)[1]+0.7*dist
  start2<-range(spatialposition$zcoord)[1]
  end2<-range(spatialposition$zcoord)[2]
  p=ggplot(spatialposition,aes(x=xcoord,y=zcoord))+geom_point(size=5,aes(colour=inf))
  p=p+scale_colour_gradient(low="blue",high="red")
  p=p+theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none")+scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks = NULL)
  p=p+geom_rect(aes(xmin = start1, xmax =end1, ymin = start2, ymax = end2), fill = "white", alpha = 0,linetype = "dashed",color="grey",size=1.5)+
    annotate("text", x=-98, y=0, label= "DV", fontface =1,angle = 90,color="grey",size=8)+
    annotate("text", x=-104, y=start2, label= "100%", fontface =1,color="grey",size=8)+
    annotate("text", x=-98, y=end2, label= "0%", fontface =1,color="grey",size=8)+
    annotate("text", x=(start1+end1)/2, y=start2-5, label= "30%                  AP                  70%", fontface =1,color="grey",size=8)+
    annotate("text", x=end1+6, y=start2-5, label= "EL", fontface =1,color="grey",size=8)

  p=p+geom_segment(aes(x = -94, y = -11.8, xend =-94, yend = end2-5),arrow = arrow(length = unit(0.3, "cm")),color="grey")+
    geom_segment(aes(x = -94, y = -11.8, xend =-94, yend = start2+5),arrow = arrow(length = unit(0.3, "cm")),color="grey")+
    geom_segment(aes(x = -22.9, y = start2-2, xend = end1-25, yend = start2-2),arrow = arrow(length = unit(0.3, "cm")),color="grey")+
    geom_segment(aes(x = -22.9, y = start2-2, xend = start1+25, yend =start2-2),arrow = arrow(length = unit(0.3, "cm")),color="grey")
  p=p+coord_fixed(1.33)
  p
}
