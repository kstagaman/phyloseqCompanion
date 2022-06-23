#' @name taxa.prevalence
#' @title get prevalence (presence in number of samples) for each taxon
#' @description This function uses the abundance table in a `phyloseq` object to determine how many samples a given taxon (OTU/ASV) occurs in and returns either a named vector or a `data.table` with the results
#' @param physeq a phyloseq object with an `otu_table`
#' @param sorted character; how do you want the results sorted? "from.table" returns the results in the same order as the column names in the `otu_table`. "by.taxon returns the results sorted alphanumerically by taxon name. "decreasing" returns the results sorted by abundance, highest to lowest. "increasing" returns the results sorted by abundance, lowest to highest. Default is "from.table".
#' @param return.table logical; if TRUE, returns a 2-column `data.table` with results. If FALSE, returns a named vector
#' @seealso \code{\link{phyloseq}}, \code{\link{otu_table}}, \code{\link{otu.matrix}}, \code{\link{data.table}}
#' @export
#' @examples
#' data(example_phyloseq)
#' head(sort(taxa.prevalence(example_phyloseq), decreasing = TRUE))

taxa.prevalence <- function(
    physeq,
    sorted = c("from.table", "by.taxon", "decreasing", "increasing"),
    return.table = FALSE
) {
  sorted <- rlang::arg_match(sorted)
  abund.mat <- otu.matrix(physeq)
  abund.mat[abund.mat > 0] <- 1
  prev <- colSums(abund.mat)
  if (return.table) {
    prev <- as.data.table(prev, keep.rownames = "Taxon")
    names(prev)[2] <- "Prevalence"
    setkey(prev, Taxon)
    if (sorted == "from.table") {
      prev <- prev[colnames(abund.mat)]
    } else if (sorted == "decreasing") {
      prev <- prev[order(Prevalence, decreasing = TRUE)]
    } else if (sorted == "increasing") {
      prev <- prev[order(Prevalence, decreasing = FALSE)]
    }
  } else {
    if (sorted == "by.taxon") {
      prev <- prev[sort(names(prev))]
    } else if (sorted == "decreasing") {
      prev <- sort(prev, decreasing = TRUE)
    } else if (sorted == "increasing") {
      prev <- sort(prev, decreasing = FALSE)
    }
  }
  return(prev)
}