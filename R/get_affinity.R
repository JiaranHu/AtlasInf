# k is a parameter representing the number of neighbors in knn graph
# l is a papameter using in determining the scaling factor of each cell
# bdtnp is cell*genes data and "data.frame" type
#构造knn图，图中元素是相似度
#' Title
#'
#' @param k
#' @param l
#' @param bdtnp
#'
#' @return
#' @export
#'
#' @examples
get_affinity<- function(k,l,bdtnp) {

  # distances between cells are measured by the euclidean distance
  adjacency_matrix = as.matrix(dist(bdtnp))

  scaling_factor<-1:nrow(bdtnp)

  knn_graph<-matrix(0,nrow=nrow(bdtnp),ncol=nrow(bdtnp))

  k_nearest_neighbors<-matrix(nrow=nrow(bdtnp),ncol=k)

  for (i in 1:nrow(bdtnp)) {
    k_nearest_neighbors[i,] <- order(adjacency_matrix[i,])[2:(k+1)]
    #scaling_factor[i]<-k_nearest_neighbors[i,l]
    scaling_factor[i]<-adjacency_matrix[i,k_nearest_neighbors[i,l]]
    # scaling_factor[i]<-1
  }

  for (i in 1:nrow(bdtnp)) {
    knn_graph[i,k_nearest_neighbors[i,]] <- adjacency_matrix[i,k_nearest_neighbors[i,]]
  }

  affinity <- matrix(nrow=nrow(bdtnp),ncol=nrow(bdtnp))

  # the similarity measure between two cells i and j is given by
  for (i in 1:nrow(bdtnp)){
    for (j in 1:nrow(bdtnp)){
      affinity[i,j] <- exp(-0.5*adjacency_matrix[i,j]*adjacency_matrix[i,j]/(scaling_factor[i]+scaling_factor[j]))/sqrt(2*pi*(scaling_factor[i]+scaling_factor[j]))
    }
  }

  for (i in 1:nrow(bdtnp)) {
    affinity[i,i]<-0
    k_nearest_neighbors[i,] <- order(affinity[i,],decreasing = T)[1:k]
    affinity[i,-k_nearest_neighbors[i,]] <- 0
  }
  # affinity is "matrix" type
  return(affinity)
}
