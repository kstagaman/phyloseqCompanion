#' tax_table to data.table
#'
#' Extract tax_table from a phyloseq object and turn it into a data.table
#' @param ps a phyloseq object
#' @param taxon.column.name the name of the column to keep taxon names (because data.table doesn't use row names). Defaults to "Taxon". Can set to FALSE to drop sample names.
#' @seealso \code{\link{data.table}}, \code{\link{phyloseq}}, \code{\link{otu_table}}
#' @export
#' @examples
#' data(example_phyloseq)
#' taxa.DT <- taxa.data.table(example_phyloseq)
#' taxa.DT

taxa.data.table <- function(ps, taxon.column.name = "Taxon") {
  if (taxon.column.name %in% names(as(tax_table(ps), "data.frame"))) {
    dt <- data.table::as.data.table(
      as(phyloseq::tax_table(ps), "matrix")
    )
  } else {
    dt <- data.table::as.data.table(
      as(phyloseq::tax_table(ps), "matrix"),
      keep.rownames = taxon.column.name
    )
  }

  return(dt)
}