#计算互信息熵
#' Title
#'
#' @param mat1
#' @param mat2
#'
#' @return
#' @export
#'
#' @examples
compute_mut_inf <- function(mat1,mat2){
  if(is.vector(mat2))
  {
    mut_inf_vector<-vector('integer',nrow(mat1))
    for(i in 1:nrow(mat1))
      mut_inf_vector[i]<-mpmi::cminjk.pw(mat1[i,],mat2)
    mut_inf_vector
  }
  else
  {
    mut_inf_mat <- matrix(0,nrow(mat1),nrow(mat2))
    for(i in 1:nrow(mat1)){
      for(j in 1:nrow(mat2))
      {
        mut_inf_mat[i,j]<-mpmi::cminjk.pw(mat1[i,],mat2[j,])
        warnings('off')
      }
    }
    mut_inf_mat
  }
}
