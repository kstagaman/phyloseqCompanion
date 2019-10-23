#' otu_table to matrix
#'
#' Extract otu_table from a phyloseq object and turn it into a matrix
#' @param ps a phyloseq object
#' @seealso \code{\link{as}}, \code{\link{phyloseq}}, \code{\link{otu_table}}
#' @export
#' @examples
#' data(example_phyloseq)
#' otu.mat <- otu.matrix(example_phyloseq)
#' View(otu.mat)

otu.matrix <- function(ps) {
  return(as(phyloseq::otu_table(ps), "matrix"))
}