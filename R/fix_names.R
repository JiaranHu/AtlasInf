#' @title fix gene names to make them consistent between files
#' @description some genes have different names in different files of our demo data, make the names of these genes unique
#' @param data A express matrix with gene names as column names
#' @return a data frame with gene names fixed
#'
#' @export
#' @examples data<-fix_names(data)
#'
# fix gene names to make them consistent between files
fix_names <- function(data) {
  data <- as.matrix(data)
  genes <- colnames(data)
  genes <- gsub(".","-",genes,fixed = T)
  genes <- gsub("-spl-","(spl)",genes,fixed = T)
  colnames(data) <- genes
  return(as.data.frame(data))
}
