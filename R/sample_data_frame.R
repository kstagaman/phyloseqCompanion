#' sample_data to data.frame
#'
#' Extract sample_data from a phyloseq object and turn it into a data.frame
#' @param ps a phyloseq object
#' @param sample.column.name the name of the column to keep sample names (because data.table doesn't use row names). Defaults to "Sample". Can set to FALSE to drop sample names.
#' @seealso \code{\link{as}}, \code{\link{phyloseq}}, \code{\link{sample_data}}
#' @export
#' @examples
#' data(example_phyloseq)
#' smpl.DF <- sample.data.frame(example_phyloseq)
#' smpl.DF

sample.data.frame <- function(ps) {
  return(as(phyloseq::sample_data(ps), "data.frame"))
}
