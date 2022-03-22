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
    "taxonomic with Aitchison",
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
  method.ids <-  c(
    "bray.curtis", "canberra", "sorensen", "aitchison",
    "w.unifrac", "h.unifrac", "u.unifrac"
  )
  method.choices <- tolower(methods) %>%
    str_replace_all("ø", "o") %>%
    str_replace_all("weighted", "w") %>%
    str_replace_all("0.5|half-weighted|half", "h") %>%
    str_replace_all("unweighted", "u") %>%
    str_replace_all("[ -]", ".")
  aitch.match <- method.choices == "taxonomic.with.aitchison"
  if (any(aitch.match)) {
    method.choices[which(aitch.match)] <- "tax.with.aitch"
  }
  arg.res <- list(
    all            = list(Req.clrPS = TRUE,  Methods = method.ids),
    taxonomic      = list(Req.clrPS = FALSE, Methods = method.ids[1:3]),
    tax.with.aitch = list(Req.clrPS = TRUE,  Methods = method.ids[1:4]),
    phylogenetic   = list(Req.clrPS = FALSE, Methods = method.ids[5:7]),
    bray.curtis    = list(
      Req.clrPS = FALSE,
      Methods = "bray.curtis",
      Name = "Bray-Curtis"
    ),
    canberra       = list(
      Req.clrPS = FALSE,
      Methods = "canberra",
      Name = "Canberra"
    ),
    sorensen       = list(
      Req.clrPS = FALSE,
      Methods = "sorensen",
      Name = "Sørensen"
    ),
    aitchison      = list(
      Req.clrPS = TRUE,
      Methods = "aitchison",
      Name = "Aitchison"
    ),
    w.unifrac      = list(
      Req.clrPS = FALSE,
      Methods = "w.unifrac",
      Name = "Weighted UniFrac"
    ),
    h.unifrac      = list(
      Req.clrPS = FALSE,
      Methods = "h.unifrac",
      Name = "Half-weighted UniFrac"
    ),
    u.unifrac      = list(
      Req.clrPS = FALSE,
      Methods = "u.unifrac",
      Name = "Unweighted UniFrac"
    )
  )
  get.arg.meta <- function(choices, name) {
    sapply(choices, function(choice) { arg.res[[choice]][[name]] })
  }
  req.clrPS <- any(get.arg.meta(method.choices, "Req.clrPS")) &
    is.null(clr_ps)
  if (req.clrPS) {
    stop(
      '`clr_ps` must be supplied if methods that include "Aitchison" are called ("all", "taxonomic with Aitchison", "Aitchison")'
    )
  }
  dists.dt0 <- data.table(
    ID = c(
      "bray.curtis",
      "canberra",
      "sorensen",
      "aitchison",
      "w.unifrac",
      "h.unifrac",
      "u.unifrac"
    ),
    Arg = c("bray", "canberra", "bray", "euclidean", NA, NA, NA),
    Param = c(0, NA, 1, NA, 1, 0.5, 0)
  )
  setkey(dists.dt0, ID)
  dists.dt <- dists.dt0[
    unique(unname(unlist(get.arg.meta(method.choices, "Methods"))))
    ]
  setkey(dists.dt, ID)
  dist.ids <- copy(dists.dt$ID) %>% set_names(., get.arg.meta(., "Name"))
  unit.tbl <- otu.matrix(ps)
  if ("aitchison" %in% dist.ids) {
    clr.tbl <- otu.matrix(clr_ps)
  }

  if (cores > 1) {
    usable.cores <- min(c(cores, length(dist.ids)))
    cl <- makeCluster(usable.cores, type = "FORK")
    registerDoParallel(cl, usable.cores)
    dist.list <- foreach(
      n = dist.ids,
      .final = function(x) setNames(x, names(dist.ids)),
      .verbose = verbose
    ) %dopar% {
      if (is.na(dists.dt[n, "Arg"])) {
        tree <- phy_tree(ps)
        dist <- GUniFrac(
          unit.tbl,
          tree,
          alpha = dists.dt[n]$Param
          )$unifracs[, , 1] %>%
          as.dist()
      } else {
        if (is.na(dists.dt[n, "Param"])) {
          if (n == "aitchison") {
            dist <- vegan::vegdist(clr.tbl, method = dists.dt[n]$Arg)
          } else {
            dist <- vegan::vegdist(unit.tbl, method = dists.dt[n]$Arg)
          }
        } else {
          dist <- vegan::vegdist(
            unit.tbl,
            method = dists.dt[n]$Arg,
            binary = dists.dt[n]$Param
          )
        }
      }
      return(dist)
    }
    stopCluster(cl)
  } else {
    dist.list <- lapply(dist.ids, function(n) {
      if (is.na(dists.dt[n, "Arg"])) {
        tree <- phy_tree(ps)
        dist <- GUniFrac(
          unit.tbl,
          tree,
          alpha = dists.dt[n]$Param
          )$unifracs[, , 1] %>%
          as.dist()
      } else {
        if (is.na(dists.dt[n, "Param"])) {
          if (n == "Aitchison") {
            dist <- vegan::vegdist(clr.tbl, method = dists.dt[n]$Arg)
          } else {
            dist <- vegan::vegdist(unit.tbl, method = dists.dt[n]$Arg)
          }
        } else {
          dist <- vegan::vegdist(
            unit.tbl,
            method = dists.dt[n]$Arg,
            binary = dists.dt[n]$Param
          )
        }
      }
      return(dist)
    })
  }
  return(dist.list)
}
