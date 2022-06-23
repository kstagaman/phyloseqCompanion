#' @name clr.transform
#' @title Center Log Ratio tranform phyloseq
#' @description CLR transform abundance counts in the otu_matrix of a phyloseq object and return whole phyloseq object.
#' @param ps a phyloseq object.
#' @param min_reads The minimum reads per sample. Default=min(sample_sums(ps))
#' @param min_prop The minimum proportional abundance of a read in any sample. Default=0.001.
#' @param min_occur The minimum fraction of non-0 reads for each variable in all samples. Default=0
#' @param smpls_by_row True if rows contain samples, false if rows contain variables. Default=TRUE
#' @param method Geometric Bayesian multiplicative (GBM, default); square root BM (SQ); Bayes-Laplace BM (BL); count zero multiplicative (CZM); user-specified hyper-parameters (user). Default="CZM"
#' @param lab Unique label (numeric or character) used to denote count zeros in X (default label=0).
#' @seealso \code{\link{phyloseq}}, \code{\link{codaSeq.filter}}, \code{\link{cmultRepl}}, \code{\link{codaSeq.clr}}
#' @export
#' @examples
#' data(example_phyloseq)
#'
#' clr.ps <- clr.transform(example_phyloseq, min_reads = 10000)
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