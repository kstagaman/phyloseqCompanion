#' @name otu.data.table
#' @title otu_table to data.table
#' @description Extract otu_table from a phyloseq object and turn it into a dat.table
#' @param ps a phyloseq object
#' @param sample.column.name the name of the column to keep sample names (because data.table doesn't use row names). Defaults to "Sample". Can set to FALSE to drop sample names.
#' @seealso \code{\link{data.table}}, \code{\link{phyloseq}}, \code{\link{otu_table}}
#' @export
#' @examples
#' data(example_phyloseq)
#' otu.DT <- otu.data.table(example_phyloseq)
#' otu.DT

otu.data.table <- function(ps, sample.column.name = "Sample") {
  return(
    data.table::as.data.table(
      as(phyloseq::otu_table(ps), "matrix"),
      keep.rownames = sample.column.name
    )
  )
}
