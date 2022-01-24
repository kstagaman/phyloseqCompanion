#' A simple phyloseq object
#'
#' @format A phyloseq object with 20 samples and 13 taxa:
#' \describe{
#'   \item{OTU Table}{13 ASVs (sequences are completely random and all 290 bp) across 20 samples}
#'   \item{Sample Data}{20 samples with 3 samples variables: Age, Weight, Genotype}
#'   \item{Taxa Table}{the 13 ASVs have been randomly assigned to 6 taxonomic ranks: Kingdom, Phylum, Class, Order, Family, Genus}
#' }
"example_phyloseq"

#' A phyloseq object with numbered ASVs
#'
#' @format A phyloseq object with 20 samples and 13 taxa, where the ASV IDs have been changed from sequences to ASVXX (XX are digits):
#' \describe{
#'   \item{OTU Table}{13 ASVs (with new IDs) across 20 samples}
#'   \item{Sample Data}{20 samples with 3 samples variables: Age, Weight, Genotype}
#'   \item{Taxa Table}{the 13 ASVs have been randomly assigned to 6 taxonomic ranks: Kingdom, Phylum, Class, Order, Family, Genus}
#' }
"numberedASVs_phyloseq"
