#构造87*2002 metacell矩阵
#' Title
#'
#' @param single_cell
#' @param sc_knn
#' @param lab
#'
#' @return
#' @export
#'
#' @examples
create_metacell<-function(single_cell,sc_knn,lab){
  cluster<- as.data.frame(table(lab))
  metacell<-matrix(0,nrow=dim(cluster)[1],ncol=ncol(single_cell))
  for(i in 1:dim(cluster)[1])
  {
    row_index<-which(lab==cluster[i,1])
    neighbors_num<-rowSums(sc_knn[row_index,row_index]>0)#计算在类的范围内，每个细胞拥有的邻居个数（每个细胞有k个近邻，但在类中，近邻个数小于k）
    if(sum(neighbors_num>0))
    {
      neighbors_num<-neighbors_num/sum(neighbors_num)#归一化得到权重向量
      metacell[i,]<-colSums(neighbors_num*single_cell[row_index,])#加权平均得到metacell
    }
    else
    {
      metacell[i,]<-colMeans(single_cell[row_index,])#加权平均得到metacell
    }
  }
  colnames(metacell)<-colnames(single_cell)
  metacell
}
