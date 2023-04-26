#' Title
#'
#' @return
#' @export
#'
#' @examples
#' 用science论文的方法算vISH，这个函数中的语句：
distmap_atlas<-function(rnaseq,insituexpr,spatialposition){

  if(!require(DistMap))
    stop("Lack of indispensable package \"DistMap\", please install it by devtools::install_github(\"rajewsky-lab/DistMap\")")

  dge_normalized <- system.file("extdata", "dge_normalized.txt", package = "AtlasInf")
  dge_normalized <- read.table(dge_normalized, header = T)
  #stopifnot(all(rownames(dge_raw) == rownames(dge_normalized)))
  binarized_bdtnp <- system.file("extdata", "binarized_bdtnp.csv", package = "AtlasInf")
  binarized_bdtnp <- read.csv(binarized_bdtnp, header = T)
  binarized_bdtnp <- fix_names(binarized_bdtnp)

  colnames(spatialposition) = c("x","y","z")
  dm = new("DistMap",
           raw.data=as.matrix(rnaseq),
           data=as.matrix(dge_normalized),
           insitu.matrix=as.matrix(binarized_bdtnp),
           geometry=as.matrix(spatialposition))

  dm <- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))

  dm <- mapCells(dm)

  atlas<-matrix(0,nrow(insituexpr),nrow(rnaseq))

  colnames(atlas)<-rownames(rnaseq)

  for(i in colnames(atlas)){
    atlas[,i] = computeVISH(dm, i, threshold=0)
  }
  atlas
}
