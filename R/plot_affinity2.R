#此函数暂时不能用
plot_affinity2<-function(){
  #根据欧氏距离+各向异性高斯核算得的k近邻
  aff<-getAffinity(10,5,single_cell)
  aff[aff>0]<-1
  rotate <- function(x) t(apply(x, 2, rev))
  image(rotate(aff),col=c("black","white"))

  #根据皮尔逊相关系数算得的k近邻
  pearson<-dge_raw
  pearson_knn<-cor(pearson,method = "pearson")
  k_nearest_neighbors<-matrix(nrow=1297,ncol=10)
  for (i in 1: 1297)
  {
    pearson_knn[i,i]<-0
    k_nearest_neighbors[i,] <- order(pearson_knn[i,],decreasing = T)[1:10]
    pearson_knn[i,-k_nearest_neighbors[i,]] <- 0
  }
  pearson_knn[pearson_knn>0]<-1
  rotate <- function(x) t(apply(x, 2, rev))
  image(rotate(pearson_knn),col=c("black","white"))
}
