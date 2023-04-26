#平衡聚类
#' Title
#'
#' @param mat
#' @param clsize
#' @param method
#'
#' @return
#' @export
#'
#' @examples
nnit <- function(mat, clsize, method=c('random','maxd', 'mind')){
  clsize.rle = rle(as.numeric(cut(1:nrow(mat), ceiling(nrow(mat)/clsize))))
  clsize = clsize.rle$lengths
  lab = rep(NA, nrow(mat))
  dmat = as.matrix(dist(mat))
  cpt = 1
  while(sum(is.na(lab)) > 0){
    lab.ii = which(is.na(lab))
    dmat.m = dmat[lab.ii,lab.ii]
    if(method[1]=='random'){
      ii = sample.int(nrow(dmat.m),1)
    } else if(method[1]=='maxd'){
      ii = which.max(rowSums(dmat.m))
    } else if(method[1]=='mind'){
      ii = which.min(rowSums(dmat.m))
    } else {
      stop('unknown method')
    }
    lab.m = rep(NA, length(lab.ii))
    lab.m[head(order(dmat.m[ii,]), clsize[cpt])] = cpt
    lab[lab.ii] = lab.m
    cpt = cpt + 1
  }
  if(any(is.na(lab))){
    lab[which(is.na(lab))] = cpt
  }
  lab
}
