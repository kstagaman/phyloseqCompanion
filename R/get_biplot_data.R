#' Extract CCA/RDA biplot data
#'
#' Get all the data required to plot a CCA or RDA biplot in ggplot. Returns a list of data.tables (and one character vector) that can easily be used to create a biplot in ggplot. (Or even in the base plot function.)
#' @param ps the phyloseq object that was used to create the CCA or RDA ordination
#' @param ord a cca or rda object (created by cca or capscale from vegan)
#' @param plot.axes the CCA/RDA axes you wish to plot (must be integers). Defaults to 1 and 2
#' @seealso \code{\link{cca}}, \code{\link{capscale}}
#' @export
#' @examples
#' data(example_phyloseq)


get.biplot.data <- function(ps, ord, plot.axes = c(1, 2)) {
  smpl.dt <- sample.data.table(ps)
  setkey(smpl.dt, Sample)
  sites0 <- as.data.table(scores(ord)$sites[, plot.axes], keep.rownames = "Sample")
  setkey(sites0, Sample)
  sites.dt <- sites0[smpl.dt, nomatch = 0]
  scale <- ordiArrowMul(scores(ord)$sites[, plot.axes])
  if (ncol(ord$CCA$biplot) < 2) {
    arrows.all <- as.data.table(
      ord$CCA$biplot / scale,
      keep.rownames = "Variable"
    )
    arrows.all[, MDS1 := 0 ]
  } else {
    arrows.all <- as.data.table(
      ord$CCA$biplot[, plot.axes] / scale,
      keep.rownames = "Variable"
    )
  }
  setkey(arrows.all, Variable)

  cntr.dt <- data.table(
    scores(ord)$centroids[, plot.axes],
    keep.rownames = "Variable"
  )
  setkey(cntr.dt, Variable)
  pc.exp <- eigenvals(ord) / sum(eigenvals(ord))
  var.exp <- unname(paste0(round(pc.exp[plot.axes] * 100, 1), "%"))
  var.exp.dt <- data.table(
    axis = names(eigenvals(ord)[plot.axes]),
    var.exp = var.exp
    )
  labs <- var.exp.dt[, .(paste0(axis, " (", var.exp, ")"))][[1]]
  return(
    list(
      sample.coords = sites.dt,
      vector.coords = arws.dt,
      centroid.coords = cntr.dt,
      axes.labs = labs
    )
  )
}
