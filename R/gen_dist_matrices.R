#' Generate Distance Matrices
#'
#' Generate distance matrices for multiple beta-diversity metrics. Returns a list of matrices.
#' @param ps a phyloseq object.
#' @param methods The names of the dissimilarity (or distance) indices. Defaults ("all") are all of Bray-Curtis, Canberra, Sørensen, W Unifrac, 0.5 Unifrac, and U Unifrac. Choosing "taxonomic" just runs Bray-Curtis, Canberra, and Sørensen. Choosing "phylogenetic" just runs W Unifrac, 0.5 Unifrac, and U Unifrac, which are implemented with \code{\link{GUniFrac}}. Can also supply a vector with any subset of these choices.
#' @param cores integer indicating how many cores to run in parallel. Default 1 does not run in parallel.
#' @param verbose TRUE/FALSE passed to any functions that can be verbose
#' @param clr_ps required if methods that include "Aitchison" are called ("all", "taxonomic", "Aitchison")
#' @seealso \code{\link{phyloseq}}, \code{\link{capscale}}, \code{\link{vegdist}}, \code{\link{distance}}
#' @export
#' @examples
#' data(example_phyloseq)
#'
#' distance.mats <- gen.dist.matrices(example_phyloseq, methods = "taxonomic", cores = 2)

gen.dist.matrices <- function(
  ps,
  methods = c(
    "all",
    "taxonomic",
    "phylogenetic",
    "Bray-Curtis",
    "Canberra",
    "Sørensen",
    "Aitchison",
    "W Unifrac",
    "0.5 Unifrac",
    "U Unifrac"
    ),
  cores = 1,
  verbose = TRUE,
  clr_ps = NULL
){
  if (methods %in% c("all", "taxonomic", "Aitchison") & is.null(clr_ps)) {
    stop(
      '`clr_ps` must be supplied if methods that include "Aitchison" are called ("all", "taxonomic", "Aitchison")'
    )
  }
  dists.dt <- data.table(
    Name = c("Bray-Curtis", "Canberra", "Sørensen", "Aitchison", "W Unifrac", "0.5 Unifrac", "U Unifrac"),
    Method = c("bray", "canberra", "bray", "euclidean", NA, NA, NA),
    Param = c(0, NA, 1, NA, 1, 0.5, 0)
  )
  if (methods[1] == "taxonomic") {
    dists.dt <- dists.dt[!is.na(Method)]
  } else if (methods[1] == "phylogenetic") {
    dists.dt <- dists.dt[is.na(Method)]
  } else if (methods[1] != "all") {
    dists.dt <- dists.dt[Name %in% methods]
  }
  dist.names <- copy(dists.dt$Name) %>% magrittr::set_names(., .)
  names(dist.names) <- dist.names
  setkey(dists.dt, Name)
  unit.tbl <- otu.matrix(ps)
  if ("Aitchison" %in% dist.names) {
    clr.tbl <- otu.matrix(clr_ps)
  }

  if (cores > 1) {
    cl <- makeCluster(cores, type = "FORK")
    registerDoParallel(cl, cores)
    dist.list <- foreach(
      n = dist.names,
      .final = function(x) setNames(x, names(dist.names)),
      .verbose = verbose
    ) %dopar% {
      if (is.na(dists.dt[n, "Method"])) {
        tree <- phy_tree(ps)
        dist <- GUniFrac(unit.tbl, tree, alpha = dists.dt[n]$Param)$unifracs[, , 1] %>% as.dist()
      } else {
        if (is.na(dists.dt[n, "Param"])) {
          if (n == "Aitchison") {
            dist <- vegan::vegdist(clr.tbl, method = dists.dt[n]$Method)
          } else {
            dist <- vegan::vegdist(unit.tbl, method = dists.dt[n]$Method)
          }
        } else {
          dist <- vegan::vegdist(
            unit.tbl,
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
        tree <- phy_tree(ps)
        dist <- GUniFrac(unit.tbl, tree, alpha = dists.dt[n]$Param)$unifracs[, , 1] %>% as.dist()
      } else {
        if (is.na(dists.dt[n, "Param"])) {
          if (n == "Aitchison") {
            dist <- vegan::vegdist(clr.tbl, method = dists.dt[n]$Method)
          } else {
            dist <- vegan::vegdist(unit.tbl, method = dists.dt[n]$Method)
          }
        } else {
          dist <- vegan::vegdist(
            unit.tbl,
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
