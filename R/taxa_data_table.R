#' tax_table to data.table
#'
#' Extract tax_table from a phyloseq object and turn it into a data.table
#' @param ps a phyloseq object
#' @param sample.column.name the name of the column to keep sample names (because data.table doesn't use row names). Defaults to "Sample". Can set to FALSE to drop sample names.
#' @seealso \code{\link{data.table}}, \code{\link{phyloseq}}, \code{\link{otu_table}}
#' @export
#' @examples
#' data(example_phyloseq)
#' taxa.DT <- taxa.data.table(example_phyloseq)
#' taxa.DT

taxa.data.table <- function(ps, sample.column.name = "Sample") {
  return(
    data.table::as.data.table(
      as(phyloseq::tax_table(ps), "data.frame"),
      keep.rownames = sample.column.name
    )
  )
}