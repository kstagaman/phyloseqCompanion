#' Center Log Ratio tranform phyloseq
#'
#' CLR transform abundance counts in the otu_matrix of a phyloseq object and return whole phyloseq object.
#' @param ps a phyloseq object.
#' @param methods The names of the dissimilarity (or distance) indices. Defaults ("all") are all of Bray-Curtis, Canberra, Sørensen, W Unifrac, 0.5 Unifrac, and U Unifrac. Choosing "taxonomic" just runs Bray-Curtis, Canberra, and Sørensen. Choosing "phylogenetic" just runs W Unifrac, 0.5 Unifrac, and U Unifrac, which are implemented with \code{\link{GUniFrac}}. Can also supply a vector with any subset of these choices.
#' @param cores integer indicating how many cores to run in parallel. Default 1 does not run in parallel.
#' @param verbose TRUE/FALSE passed to any functions that can be verbose
#' @seealso \code{\link{phyloseq}}, \code{\link{capscale}}, \code{\link{vegdist}}, \code{\link{distance}}
#' @export
#' @examples
#' data(example_phyloseq)
#'
#' distance.mats <- gen.dist.matrices(example_phyloseq, methods = "taxonomic", cores = 2)
#'

clr.transform <- function(
  ps,
  min_reads = NULL,
  min_prop = 0.001,
  min_occur = 0,
  smpls_by_row = TRUE,
  method = "CZM",
  lab = 0
) {
  asv.mat <- otu.matrix(ps)
  if (is.null(min_reads)) {
    min_reads <- min(sample_sums(ps))
  }
  asv.mat.f <- codaSeq.filter(
    asv.mat,
    min.reads = min_reads,
    min.prop = min_prop,
    min.occurrence = min_occur,
    samples.by.row = smpls_by_row
  )
  # replace 0 values with an estimate
  asv.mat.fn0 <- cmultRepl(t(asv.mat.f), method = method, label = lab)
  otu_table(ps) <- otu_table(codaSeq.clr(asv.mat.fn0), taxa_are_rows = !smpls_by_row)
  return(ps)
}