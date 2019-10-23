#' Numbered ASVs
#'
#' Replace ASV names (sequences) as output by dada2 with ASV numbers (e.g. ASV001, ASV002, ...). Returns a phyloseq object with the new names.
#' @param ps a phyloseq object
#' @param prefix The string you want to come before the id numbers. Defaults to "ASV". Cannot be "" or logical (TRUE/FALSE) as numbers can't be used as column names in data.frames.
#' @param save.dir directory in which to save the original ASV sequences, for reference. Defaults to "." (current directory).
#' @param save.file file name to save the original ASV sequences. Does NOT require an extension (saves as a .rds). Defaults to "ASV_sequences".
#' @seealso \code{\link{taxa_names}}, \code{\link{phyloseq}}
#' @export
#' @examples
#' data(example_phyloseq)
#' head(taxa_names(example_phyloseq))
#' new.phyloseq <- numbered.ASVs(example_phyloseq)
#' head(taxa_names(new.phyloseq))

numbered.ASVs <- function(ps, prefix = "ASV", save.dir = ".", save.file = "ASV_sequences") {
  if (prefix == "" | length(prefix) == 0 | is.logical(prefix)) {
    stop("prefix must be a non-numeric-only string of length greater than 0 and not logical (TRUE/FALSE)")
  }
  asv.seqs <- taxa_names(ps)
  saveRDS(asv.seqs, file = file.path(save.dir, paste0(save.file, ".rds")))
  n.digits <- nchar(length(asv.seqs))
  id.nums <- sapply(1:length(asv.seqs), function(d) {
    if (nchar(d) < 4) {
      zeroes <- paste(rep("0", n.digits - nchar(d)), collapse = "")
      paste0(zeroes, d)
    } else {
      paste0(d)
    }
  }) %>% unlist()
  taxa_names(ps) <- paste0(prefix, id.nums)
  return(ps)
}
