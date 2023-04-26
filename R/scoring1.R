#' Title
#'
#' @param inference_mat
#' @param spatialposition
#' @param rownumber
#' @param k
#'
#' @return
#' @export
#'
#' @examples
scoring1 <- function(inference_mat, spatialposition, rownumber, k)#计算指标1的函数，参数k表示邻居个数
{
  score1 <- 0
  sum <- 0
  tem <- 0
  distance = as.matrix(dist(spatialposition))#欧氏距离矩阵
  k_nearest_neighbors <- matrix(nrow = rownumber, ncol = k)
  for (i in 1:rownumber) {
    k_nearest_neighbors[i, ] <- order(distance[i, ])[2:(k + 1)]
    sum <- 0
    for (j in 1:k) {
      tem <- var(inference_mat[i, ] - inference_mat[k_nearest_neighbors[i,j], ])
      tem <- tem / distance[i, k_nearest_neighbors[i,j]]
      sum <- sum + tem
    }
    score1 <- score1 + sum
  }
  score1 <- score1 / k
}
