#' Title
#'
#' @param sc_knn
#' @param mc_knn
#' @param bulk_knn
#' @param mut_inf_mat1
#' @param mut_inf_mat2
#'
#' @return
#' @export
#' @importFrom Matrix Matrix t rowSums
#' @examples
network_integrate<-function(sc_knn,mc_knn,bulk_knn,mut_inf_mat1,mut_inf_mat2){

  dimen<-nrow(mc_knn)+nrow(bulk_knn)#扩散核矩阵的阶数

  SMatrix <- Matrix(0, nrow = dimen,ncol = dimen, sparse = TRUE)#稀疏方阵
  #sc_knn<-(sc_knn-min(sc_knn))/max(sc_knn)#归一化
  #SMatrix[1:nrow(sc_knn),1:nrow(sc_knn)]<-sc_knn#把小矩阵的值填入大矩阵

  mc_knn<-(mc_knn-min(mc_knn))/max(mc_knn)#归一化
  SMatrix[1:nrow(mc_knn),1:nrow(mc_knn)]<-mc_knn#把小矩阵的值填入大矩阵

  bulk_knn<-(bulk_knn-min(bulk_knn))/max(bulk_knn)#归一化
  SMatrix[(1+nrow(mc_knn)):dimen,(1+nrow(mc_knn)):dimen]<-bulk_knn#把小矩阵的值填入大矩阵



  #mut_inf_mat1<-(mut_inf_mat1-min(mut_inf_mat1))/(max(mut_inf_mat1)-min(mut_inf_mat1))#归一化
  #SMatrix[1:nrow(sc_knn),(1+nrow(sc_knn)):(nrow(sc_knn)+nrow(mc_knn))]<-mut_inf_mat1#把小矩阵的值填入大矩阵
  mut_inf_mat2<-(mut_inf_mat2-min(mut_inf_mat2))/(max(mut_inf_mat2)-min(mut_inf_mat2))#归一化
  SMatrix[1:nrow(mc_knn),(1+nrow(mc_knn)):dimen]<-mut_inf_mat2#把小矩阵的值填入大矩阵

  SMatrix<-(SMatrix+t(SMatrix))/2 #对称化

  D<-rowSums(SMatrix)+.Machine$double.eps #行和向量
  D<-1/sqrt(D)
  diagonal<-Matrix(0,dimen,dimen, sparse = TRUE)
  diag(diagonal) <-D #构造稀疏对角阵
  graph<-diagonal%*%SMatrix%*%diagonal #markov型规范化,G对称
  as.matrix(graph)
}
