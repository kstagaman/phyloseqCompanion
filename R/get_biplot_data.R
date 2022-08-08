#' @name get.biplot.data
#' @title Extract CCA/RDA biplot data
#' @description Get all the data required to plot a CCA or RDA biplot in ggplot. Returns a list of data.tables (and one character vector) that can easily be used to create a biplot in ggplot (or the base plot function).
#' @param smpls a data.frame/data.table or the phyloseq object that was used to create the CCA or RDA ordination
#' @param ord a cca or rda object (created by cca or capscale from vegan)
#' @param plot.axes the CCA/RDA axes you wish to plot (must be integers). Defaults to 1 and 2
#' @seealso \code{\link{cca}}, \code{\link{capscale}}
#' @export
#' @examples
#' data(example_phyloseq)

get.biplot.data <- function(smpls, ord, plot.axes = c(1, 2)) {
  if ("phyloseq" %in% class(smpls)) {
    smpl.dt <- sample.data.table(smpls)
  } else if ("data.table" %in% class(smpls)) {
    smpl.dt <- smpls
  } else if ("data.frame" %in% class(smpls)) {
    smpl.dt <- as.data.table(smpls, keep.rownames = "Sample")
  } else {
    rlang::abort("Argument `smpls' must be a data.frame, data.table, or phyloseq object.")
  }
  setkey(smpl.dt, Sample)
  ord.scores <- vegan::scores(ord, choices = plot.axes)
  sites0 <- as.data.table(ord.scores$sites, keep.rownames = "Sample")
  setkey(sites0, Sample)
  sites.dt <- sites0[smpl.dt, nomatch = 0]
  scale <- ordiArrowMul(ord.scores$sites)
  if (ncol(ord$CCA$biplot) < 2) {
    arws.dt <- as.data.table(
      ord$CCA$biplot / scale,
      keep.rownames = "Variable"
    )
    arws.dt[, MDS1 := 0]
  } else {
    arws.dt <- as.data.table(
      ord$CCA$biplot[, plot.axes] / scale,
      keep.rownames = "Variable"
    )
  }
  setkey(arws.dt, Variable)
  if (anyNA(ord.scores$centroids)) {
    cntr.dt <- NA
  } else if (is.null(ord.scores$centroids)) {
    cntr.dt <- NULL
  } else {
    cntr.dt <- data.table(
      ord.scores$centroids,
      keep.rownames = "Variable"
    )
    setkey(cntr.dt, Variable)
  }

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
      axes.labs = labs,
      coord.scale = scale
    )
  )
}

