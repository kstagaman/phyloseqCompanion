#' Capscale Phyloseq
#'
#' Run the function capscale on a phyloseq object. Returns a cca object.
#' @param formula Model formula. The function can be called only with the formula interface. The LHS should be a phyloseq object, and the RHS gives the constraining variables, and conditioning variables can be given within a special function Condition.
#' @param dist.method The name of the dissimilarity (or distance) index. See \code{\link{distanceMethodList}} for a detailed list of the supported options, and links to accompanying documentation.
#' @param binary Perform presence/absence standardization before analysis.
#' @param ... Other parameters to pass to \code{\link{capscale}}, like \code{\link{sqrt.dist}}
#' @seealso \code{\link{phyloseq}}, \code{\link{capscale}}, \code{\link{vegdist}}, \code{\link{distance}}
#' @export
#' @examples
#' data(example_phyloseq)
#'
#' rda.data <- capscale.phyloseq(example_phyloseq ~ Genotype * Weight + Condition(Age))
#' print(vegan::anova.cca(rda.data$rda, by = "terms"))


capscale.phyloseq <- function(
  formula,
  dist.method = "euclidean",
  binary = FALSE,
  ...
){
  if (class(formula) != "formula") {
    stop("argument formula must be of class formula")
  }
  if (binary != TRUE & binary != FALSE) {
    stop("argument binary must be TRUE or FALSE")
  }
  frm.vec <- as.character(formula)
  ps <- eval(parse(text = frm.vec[2]))
  frm.str <- frm.vec[3]
  constrain.vars <- grep("Condition", strsplit(frm.str, " [+*|] ")[[1]], invert = T, value = T) %>%
    gsub("\\(", "", .) %>%
    gsub("\\)", "", .)
  if (grepl("Condition", frm.str)) {
    condition.vars <- {strsplit(frm.str, "Condition\\(")[[1]] %>%
        grep("\\)", ., value = T) %>%
        sub(")", "", .) %>%
        strsplit(., " [+*|] ")}[[1]]
    all.vars <- c(constrain.vars, condition.vars)
  } else {
    all.vars <- constrain.vars
  }

  smpl.df <- sample.data.frame(ps)
  smpl.df <- smpl.df[, all.vars]
  good.smpls <- row.names(smpl.df)[complete.cases(smpl.df)]
  if (length(good.smpls) < length(nrow(smpl.df))) {
    cat(
      paste(
        "Samples dropped due to NAs:",
        paste(row.names(smpl.df)[!complete.cases(smpl.df)], collapse = ", ")
      )
    )
  }
  good.ps <- prune_samples(good.smpls, ps)

  if (binary) {
    dist.mat <- phyloseq::distance(good.ps, method = dist.method, binary = TRUE)
  } else {
    dist.mat <- phyloseq::distance(good.ps, method = dist.method)
  }

  rda.frm <- paste("dist.mat ~", frm.vec[3])
  rda <- vegan::capscale(as.formula(rda.frm), data = smpl.df[good.smpls, ], ...)
  return(list(rda = rda, dist.matrix = dist.mat, metadata = smpl.df[good.smpls, ]))
}

