#' @name taxa.data.table
#' @title tax_table to data.table
#' @description Extract tax_table from a phyloseq object and turn it into a data.table
#' @param ps a phyloseq object
#' @param taxon.column.name the name of the column to keep taxon names (because data.table doesn't use row names). Defaults to "Taxon". Can set to FALSE to drop sample names.
#' @seealso \code{\link{data.table}}, \code{\link{phyloseq}}, \code{\link{otu_table}}
#' @export
#' @examples
#' data(example_phyloseq)
#' taxa.DT <- taxa.data.table(example_phyloseq)
#' taxa.DT

taxa.data.table <- function(ps, taxon.column.name = "Taxon") {
  tax.mat <- as(phyloseq::tax_table(ps), "matrix")
  if (taxon.column.name %in% colnames(tax.mat)) {
    dt <- data.table::as.data.table(tax.mat)
  } else {
    dt <- data.table::as.data.table(
      tax.mat,
      keep.rownames = taxon.column.name
    )
  }

  return(dt)
}