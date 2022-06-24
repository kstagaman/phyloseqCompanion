#' @name prep.physeq
#' @title Prepare a phyloseq with raw counts for analysis
#' @description Perform common, important pre-analysis steps such as renaming NA taxa, removing eukaryotic taxa, rarefaction, and estimation of alpha- and beta-diversities
#' @param min.prevalence numeric; the minimum acceptable number of samples a taxon must be present in to be included in subsequent analyses. Default is 1 (no prevalence filtering).
#' @param physeq a phyloseq object.
#' @param remove.euks logical; whether to remove taxa assigned to eukaryotic taxonomies with \code{\link{phyloseqCompanion::remove.euks}}. Default is TRUE.
#' @param min.rarefaction numeric; the minimum acceptable depth to rarefy to (actual number we be the smallest sample sum equal to or greater than this number). Set to 0 for no rarefaction. Default is 10000.
#' @param alpha.metrics character; a vector of alpha-diveristy metrics to estimate. Must match arguments for \code{\link{phyloseq::estimate_richness}}, or be "Phylogenetic" or "Richness" (both assessed with \code{\link{picante::pd}}). Set to NULL to skip alpha-diversity estimation. Default is c("Chao1", "Shannon", "Simpson", "Phylogenetic", "Richness").
#' @param beta.metrics; character; a vector of beta-diversity metrics (or groups of metrics) that match arugments for \code{\link{phyloseqCompanion::gen.dist.matrices}}. Set to NULL to skip beta-diversity estimation. Default is c("taxonomic", "phylogenetic").
#' @param sample.col.name character; name to use for sample column in extracted tables. Default is "Sample".
#' @param taxon.col.name character; name to use for taxon column in extracted tables. Default is "Taxon".
#' @param n.cores numeric; number of cores to use for beta-diversity estimation. Default is 1.
#' @param user.seed numeric; random seed for reproducibility. Default is 42.
#' @param ... additional arguments to pass to certain functions (like `clr_ps` to \code{\link{gen.dist.matrices}})
#' @seealso \code{\link{phyloseq}}, \code{\link{subset_taxa}}
#' @export

prep.physeq <- function(
    physeq,
    min.rarefaction = 1e4,
    min.prevalence = 1,
    alpha.metrics = c("Chao1", "Shannon", "Simpson", "Phylogenetic", "Richness"),
    beta.metrics = c("taxonomic", "phylogenetic"),
    sample.col.name = "Sample",
    taxon.col.name = "Taxon",
    n.cores = 1,
    user.seed = 42,
    ...
) {
  # Checks
  er.allowed <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
  okay.alphas <- c(er.allowed, "Phylogenetic", "Richness")
  if (!all(alpha.metrics %in% okay.alphas)) {
    rlang::abort(
      "One or more arguments provided to `alpha.metrics' is not supported by either phyloseq::estimate_richness() or picante::pd()"
    )
  } else {
    alpha.metrics <- setNames(alpha.metrics, alpha.metrics)
    er.alphas <- alpha.metrics[alpha.metrics %in% er.allowed]
    pd.alphas <- alpha.metrics[alpha.metrics %in% c("Phylogenetic", "Richness")]
  }
  gdm.allowed <- c(
    "all",
    "taxonomic",
    "taxonomic with Aitchison",
    "phylogenetic",
    "Bray-Curtis",
    "Canberra",
    "Sørensen",
    "Aitchison",
    "W Unifrac",
    "0.5 Unifrac",
    "U Unifrac"
  )
  if (!all(beta.metrics %in% gdm.allowed)) {
    rlang::abort(
      "One or more arguments provided to `beta.metrics' is not supported by phyloseqCompanion::gen.dist.matrices()"
    )
  }
  # Process
  if (min.prevalence == 1) {
    ps1 <- physeq
  } else {
    tax.prev <- taxa.prevalence(physeq, sorted = "decreasing")
    ps1 <- prune_taxa(taxa = names(tax.prev[tax.prev >= min.prevalence]), x = physeq)
  }
  if (remove.euks) {
    ps2 <- subset_taxa(
      ps1,
      Kingdom != "Eukaryota" & Order != "Chloroplast" & Family != "Mitochondria"
    ) %>%
      rename.NA.taxa()
  } else {
    ps2 <- rename.NA.taxa(ps1)
  }
  if (min.rarefaction == 0) {
    ps3 <- ps2
  } else {
    rar.depth <- min(sample_sums(ps2)[sample_sums(ps2) >= min.rarefaction])
    ps3 <- rarefy_even_depth(physeq = physeq, sample.size = rar.depth, rngseed = user.seed)
  }
  sample.dt <- sample.data.table(ps3) %>% setkeyv(sample.col.name)
  if (!is.null(alpha.metrics)) {
    if (length(er.alphas) > 0) {
      tax.alpha.dt <- estimate_richness(ps3, measures = er.alphas) %>%
        as.data.table(keep.rownames = sample.col.name) %>%
        .[, se.chao1 := NULL] %>%
        setkeyv(sample.col.name)
      sample.dt <- copy(sample.dt)[tax.alpha.dt]
    }
    if (length(pd.alphas) > 0) {
      require(picante)
      require(phytools)
      phy.alpha.dt <- pd(samp = otu.matrix(ps3), tree = midpoint.root(phy_tree(ps3))) %>%
        as.data.table(keep.rownames = sample.col.name)
      if (length(pd.alphas) < 2) {
        if (pd.alphas == "Phylogenetic") {
          phy.alpha.dt <- phy.alpha.dt[, 1:2]
        } else {
          phy.alpha.dt <- phy.alpha.dt[, c(1, 3)]
        }
        names(phy.alpha.dt)[2] <- pd.alphas
      } else {
        names(phy.alpha.dt)[2:3] <- pd.alphas
      }
      setkeyv(phy.alpha.dt, sample.col.name)
      sample.dt <- copy(sample.dt)[phy.alpha.dt]
    }
  }
  dist.list <- NULL
  if (!is.null(beta.metrics)) {
    dist.list <- gen.dist.matrices(
      methods = beta.metrics,
      ps = ps3,
      cores = n.cores
    )
    beta.metrics <- setNames(beta.metrics, beta.metrics)
  }
  list(
    Physeq = ps3,
    Samples = sample_names(ps3),
    Taxa = taxa_names(ps3),
    Metadata = sample.dt,
    Abund.mat = otu.matrix(ps3),
    Taxonomy = taxa.data.table(ps3),
    Dist.mats = dist.list,
    Alpha.metrics = alpha.metrics,
    Beta.metrics = beta.metrics
  ) %>% return()
}