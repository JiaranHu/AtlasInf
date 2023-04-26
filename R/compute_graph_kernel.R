#' Title
#'
#' @param eps
#' @param count
#' @param graph
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
compute_graph_kernel<-function(eps=1e-5,graph,count,lambda){
  v<-rowSums(graph)
  v<-v+lambda
  initG<-Matrix::Matrix(0,nrow(graph),ncol(graph), sparse = TRUE)
  diag(initG) <- 1/v
  I<-Matrix::Matrix(0,nrow(graph),ncol(graph), sparse = TRUE)
  diag(I) <- 1
  H<-initG %*% graph
  ref<-I-H
  err<-0
  mat<-list(I,I)
  for(i in 1:count){
    if(i%%2)
    {
      mat[[2]]<-I+H%*%mat[[1]]
      refp<-ref%*%mat[[2]]%*%ref
      err<-max(abs(refp-ref))/max(max(abs(refp)),max(abs(ref)))
      if(err<eps)
        break
    }
    else
    {
      mat[[1]]<-I+H %*% (mat[[2]])
      refp<-ref%*%mat[[1]]%*%ref
      err<-max(abs(refp-ref))/max(max(abs(refp)),max(abs(ref)))
      if(err<eps)
        break
    }
  }
  if(err>=eps)
    print("can't achieve convergence wichin maximum iteration times")
  else
    print(paste("achieved convergence through iteration times:",as.character(count)))
  mat[[1+i%%2]]%*%initG
}
