#' Phyloseq to LefSe
#'
#' Convert phyloseq otu and tax tables into a file that can be fed into LefSe
#' @param ps A phyloseq object that contains a sample data table, an OTU (or ASV) table, and a taxonomy table.
#' @param covars A character vector with, the names of the variables in the sample data that will be used in LefSe. Max length: 2.
#' @param file.name Character string with name of file you want to write the output to. The convention is to use the .txt extension.
#' @param taxa.levels A character vector with the taxonomic levels you want considered by LefSe. Default is "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus".
#' @param transpose.otus Phyloseq generally makes the OTUs/ASVs the column names, but this function needs them as row names. If OTUs/ASVs are already row names, set this to FALSE. Default is TRUE.
#' @seealso \code{\link{phyloseq}}
#' @export
#' @examples
#' data(example_phyloseq)
#'
#' phyloseq2lefse(ps = example_phyloseq, covars = "Genotype", file.name = "example_lefse_data.txt")

phyloseq2lefse <- function(
  ps,
  covars,
  file.name = "lefse_data.txt",
  taxa.levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
  transpose.otus = TRUE
) {
  if (length(covars) > 2) {
    warning("The length of the `covars` vector is greater than 2. File must be edited manually for use in LEfSe, which throws and error when there are more rows than specific classes, subclasses, and subjects.")
  }
  smpl.data <- sample.data.frame(ps)
  smpl.data$Sample <- row.names(smpl.data)
  t.smpl.data <- t(smpl.data)
  t.smpl.data <- as.data.frame(t.smpl.data[c("Sample", covars), ])
  if (transpose.otus) {
    otu.tbl <- t(otu.matrix(ps)) # grab the otu table from the phyloseq object
  } else {
    otu.tbl <- otu.matrix(ps) # grab the otu table from the phyloseq object
  }
  tax.tbl <- taxa.matrix(ps) # grab the taxa table from the phyloseq object and coerce into a matrix
  tax.tbl <- tax.tbl[, taxa.levels]

  # The following loop goes through the taxa table starting from the highest level and gradually moving to the lowest taxonomic level. For each level it appends the unique taxonomic levels to the uniq.lvls vector.
  uniq.lvls <- c()
  for (i in c(1:length(tax.tbl[1, ]))) {
    lvls <- as.data.frame(do.call(paste, c(as.data.frame(tax.tbl[, 1:i]), sep = "|")))
    names(lvls) <- "tax.lvl"
    uniq.i <- as.character(unique(lvls$tax.lvl))
    uniq.lvls <- c(uniq.lvls, uniq.i)
  }
  tax.tbl.join <- as.data.frame(do.call(paste, c(as.data.frame(tax.tbl), sep = "|")))
  row.names(tax.tbl.join) <- row.names(tax.tbl)
  names(tax.tbl.join) <- "tax.lvl"

  # This loop goes through each sample (which are now column names for t.smpl.data), and calculates the relative abundance for each unique taxonomic level (from above). These abundances only sum to 1 for *each taxonomic level*
  uniq.tax.lvl.abunds <- data.frame(row.names = uniq.lvls)
  for (smpl in names(t.smpl.data)) {
    abunds <- as.data.frame(otu.tbl[row.names(otu.tbl), smpl])
    total.abund <- sum(abunds[, 1])
    smpl.tax.lvl.abunds <- cbind(abunds, tax.tbl.join)

    smpl.uniq.lvl.abunds <- data.frame()
    for (uniq.lvl in uniq.lvls) {
      uniq.sub <- subset(smpl.tax.lvl.abunds, grepl(uniq.lvl, smpl.tax.lvl.abunds$tax.lvl, fixed = TRUE))
      lvl.total <- as.factor(sum(uniq.sub[, 1]) / total.abund)
      uniq.lvl.smpl <- data.frame(row.names = uniq.lvl, "sample" = lvl.total)
      names(uniq.lvl.smpl) <- smpl
      smpl.uniq.lvl.abunds <- rbind(smpl.uniq.lvl.abunds, uniq.lvl.smpl)
    }

    uniq.tax.lvl.abunds <- cbind(uniq.tax.lvl.abunds, smpl.uniq.lvl.abunds)
  }

  final.data <- rbind(t.smpl.data, uniq.tax.lvl.abunds)
  write.table(final.data, file=file.name, col.names=FALSE, sep="\t", quote=FALSE)
}
