#' Extract CCA/RDA biplot data
#'
#' Get all the data required to plot a CCA or RDA biplot in ggplot. Returns a list of data.tables (and one character vector) that can easily be used to create a biplot in ggplot. (Or even in the base plot function.)
#' @param ps the phyloseq object that was used to create the CCA or RDA ordination
#' @param ord a cca or rda object (created by cca or capscale from vegan)
#' @param anova an anova.cca object that has been run on the ord object
#' @param alpha the p-value cutoff for plotting significant vectors and centroids. Defaults to 0.05
#' @param plot.axes the CCA/RDA axes you wish to plot (must be integers). Defaults to 1 and 2
#' @seealso \code{\link{cca}}, \code{\link{capscale}}
#' @export
#' @examples
#' data(example_phyloseq)

get.biplot.data <- function(ps, ord, anova, alpha = 0.05, plot.axes = c(1, 2)) {
  smpl.dt <- sample.data.table(ps)
  setkey(smpl.dt, Sample)
  aov.dt <- data.table::as.data.table(anova, keep.rownames = "Term")
  names(aov.dt)[5] <- "Pr.F"
  sig.labs <- sort(aov.dt$Term[aov.dt$Pr.F < alpha])
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

  cntr.dt <- data.table()
  xlevs <- attributes(ord$terminfo$xlev)$names
  factor.terms <-  xlevs[xlevs %in% sig.labs]
  if (any(grepl(":", sig.labs))) {
    factors.intr.terms <- gsub(
      paste0(factor.terms, ":"), "",
      sig.labs[grepl(":", sig.labs)]
    )
    sig.arw.terms <- grep(
      paste(c(factor.terms, ":"), collapse = "|"),
      sig.labs,
      invert = T,
      value = T
    )
    arws.sig <- lapply(
      c(sig.arw.terms, factors.intr.terms),
      function(t) grep(t, arrows.all$Variable, value = T)
    ) %>% unlist()
    cntr.dt <- as.data.table(
      scores(ord)$centroids[, plot.axes],
      keep.rownames = "Variable"
    )
  } else {
    arws.sig <- sapply(
      sig.labs[!(sig.labs %in% factor.terms)],
      function(t) grep(t, arrows.all$Variable, value = T)
    ) %>% c()
    if (any(factor.terms %in% sig.labs)) {
      cntr.dt <- as.data.table(
        scores(ord)$centroids[, plot.axes],
        keep.rownames = "Variable"
      )
    }
  }
  for (term in factor.terms) {
    cntr.dt[, Group := gsub(term, "", Variable)]
    names(cntr.dt)[which(names(cntr.dt) == "Group")] <- term
  }


  arws.dt <- arrows.all[arws.sig]
  arws.dt[, Color := ifelse(grepl(":", Variable), "gray50", "black")]
  # eigenvals() == c(ord$CCA$eig, ord$CA$eig)
  pc.exp <- eigenvals(ord) / sum(eigenvals(ord))
  var.exp <- unname(paste0(round(pc.exp[plot.axes] * 100, 1), "%"))
  var.exp.dt <- data.table(
    axis = names(eigenvals(ord)[plot.axes]),
    var.exp = var.exp
    )
  labs <- var.exp.dt[, .(paste(axis, var.exp))][[1]]
  return(
    list(
      sample.coords = sites.dt,
      vector.coords = arws.dt,
      centroid.coords = cntr.dt,
      axes.labs = labs
    )
  )
}
