#' sample_data to data.table
#'
#' Extract sample_data from a phyloseq object and turn it into a data.table
#' @param ps a phyloseq object
#' @param sample.column.name the name of the column to keep sample names (because data.table doesn't use row names). Defaults to "Sample". Can set to FALSE to drop sample names.
#' @seealso \code{\link{data.table}}, \code{\link{phyloseq}}, \code{\link{sample_data}}
#' @export
#' @examples
#' data(example_phyloseq)
#' smpl.DT <- sample.data.table(example_phyloseq)
#' smpl.DT

sample.data.table <- function(ps, sample.column.name = "Sample") {
  return(
    data.table::as.data.table(
      as(phyloseq::sample_data(ps), "data.frame"),
      keep.rownames = sample.column.name
      )
    )
}
