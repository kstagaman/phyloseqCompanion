#' Capscale Phyloseq
#'
#' Run the function capscale on a phyloseq object. Returns a cca object.
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

gen.dist.matrices <- function(
  ps,
  methods = c("all", "taxonomic", "phylogenetic", "Bray-Curtis", "Canberra", "Sørensen", "W Unifrac", "0.5 Unifrac", "U Unifrac"),
  cores = 1,
  verbose = TRUE
){
  dists.dt <- data.table(
    Name = c("Bray-Curtis", "Canberra", "Sørensen", "W Unifrac", "0.5 Unifrac", "U Unifrac"),
    Method = c("bray", "canberra", "bray", NA, NA, NA),
    Param = c(0, NA, 1, 1, 0.5, 0)
  )
  if (methods == "taxonomic") {
    dists.dt <- dists.dt[!is.na(Method)]
  } else if (methods == "phylogenetic") {
    dists.dt <- dists.dt[is.na(Method)]
  } else if (methods != "all") {
    dists.dt <- dists.dt[Name %in% methods]
  }
  dist.names <- copy(dists.dt$Name)
  names(dist.names) <- dist.names
  setkey(dists.dt, Name)

  if (cores > 1) {
    cl <- makeCluster(cores, type = "FORK")
    registerDoParallel(cl, cores)
    dist.list <- foreach(
      n = dist.names,
      .final = function(x) setNames(x, names(dist.names)),
      .verbose = verbose
    ) %dopar% {
      if (is.na(dists.dt[n, "Method"])) {
        unit.tbl <- otu.matrix(ps)
        tree <- phy_tree(ps)
        dist <- GUniFrac(unit.tbl, tree, alpha = dists.dt[n]$Param)$unifracs[, , 1] %>% as.dist()
      } else {
        if (is.na(dists.dt[n, "Param"])) {
          dist <- phyloseq::distance(ps, method = dists.dt[n]$Method)
        } else {
          dist <- phyloseq::distance(
            ps,
            method = dists.dt[n]$Method,
            binary = dists.dt[n]$Param
          )
        }
      }
      return(dist)
    }
    stopCluster(cl)
  } else {
    dist.list <- lapply(dist.names, function(n) {
      if (is.na(dists.dt[n, "Method"])) {
        asv.tbl <- otu.matrix(ps)
        tree <- phy_tree(ps)
        dist <- GUniFrac(asv.tbl, tree, alpha = dists.dt[n]$Param)$unifracs[, , 1] %>% as.dist()
      } else {
        if (is.na(dists.dt[n, "Param"])) {
          dist <- phyloseq::distance(ps, method = dists.dt[n]$Method)
        } else {
          dist <- phyloseq::distance(
            ps,
            method = dists.dt[n]$Method,
            binary = dists.dt[n]$Param
          )
        }
      }
      return(dist)
    })
  }
  return(dist.list)
}