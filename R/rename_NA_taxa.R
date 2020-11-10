#' rename NA taxonomic levels
#'
#' Rename taxonomic levels so that NAs don't all get lumped together (i.e. uses higher-level assignments in lower-level names to preserve lineages).
#' @param ps a phyloseq object
#' @seealso \code{\link{tax_table}}, \code{\link{phyloseq}}
#' @export
#' @examples
#' data(example_phyloseq)
#' new.ps <- rename.NA.taxa(example_phyloseq)
#' head(tax_table(new.ps))

rename.NA.taxa <- function(ps) {
  tax.tbl <- taxa.data.table(ps)
  taxa.names <- taxa_names(ps)
  reordered <- c(names(tax.tbl)[-1], names(tax.tbl)[1])
  tax.tbl <- tax.tbl[, ..reordered]
  for (col in 2:ncol(tax.tbl)) {
    taxlevel <- names(tax.tbl)[col]
    prev.col <- tax.tbl[[col - 1]]
    curr.col <- tax.tbl[[col]]
    to.replace <- is.na(curr.col) |
      grepl("Unknown_Family", curr.col) |
      grepl("Incertae_Sedis", curr.col)
    curr.col[to.replace] <- paste0(prev.col[to.replace], "_", taxlevel)
    tax.tbl[[col]] <- curr.col
  }
  tax.mat <- as.matrix(tax.tbl)
  rownames(tax.mat) <- taxa.names
  tax_table(ps) <- tax_table(tax.mat)
  return(ps)
}
