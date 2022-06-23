#' @name taxa.matrix
#' @title tax_table to matrix
#' @description Extract tax_table from a phyloseq object and turn it into a matrix
#' @param ps a phyloseq object
#' @seealso \code{\link{as.matrix}}, \code{\link{phyloseq}}, \code{\link{otu_table}}
#' @export
#' @examples
#' data(example_phyloseq)
#' taxa.mat <- taxa.matrix(example_phyloseq)
#' taxa.mat

taxa.matrix <- function(ps) {
  return(
    as(phyloseq::tax_table(ps), "matrix")
  )
}