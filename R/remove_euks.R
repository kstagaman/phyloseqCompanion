#' @name remove.euks
#' @title Remove eukaryotic taxa from a phyloseq object
#' @description Removes taxa that match Kingdom "Eukaryota", Order "Chloroplast", or Family "Mitochondria".
#' @param physeq a phyloseq object.
#' @seealso \code{\link{phyloseq}}, \code{\link{subset_taxa}}
#' @export

remove.euks <- function(physeq) {
  subset_taxa(
    physeq,
    Kingdom != "Eukaryota" & Order != "Chloroplast" & Family != "Mitochondria"
  ) %>%
    return()
}