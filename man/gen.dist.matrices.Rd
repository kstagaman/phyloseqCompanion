% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_dist_matrices.R
\name{gen.dist.matrices}
\alias{gen.dist.matrices}
\title{Generate Distance Matrices}
\usage{
gen.dist.matrices(
  ps,
  methods = c("all", "taxonomic", "taxonomic with Aitchison", "phylogenetic",
    "Bray-Curtis", "Canberra", "Sørensen", "Aitchison", "W Unifrac", "0.5 Unifrac",
    "U Unifrac"),
  cores = 1,
  verbose = TRUE,
  clr_ps = NULL
)
}
\arguments{
\item{ps}{a phyloseq object.}

\item{methods}{The names of the dissimilarity (or distance) indices. Defaults ("all") are all of Bray-Curtis, Canberra, Sørensen, W Unifrac, 0.5 Unifrac, and U Unifrac. Choosing "taxonomic" just runs Bray-Curtis, Canberra, and Sørensen. Choosing "phylogenetic" just runs W Unifrac, 0.5 Unifrac, and U Unifrac, which are implemented with \code{\link{GUniFrac}}. Can also supply a vector with any subset of these choices.}

\item{cores}{integer indicating how many cores to run in parallel. Default 1 does not run in parallel.}

\item{verbose}{TRUE/FALSE passed to any functions that can be verbose}

\item{clr_ps}{required if methods that include "Aitchison" are called ("all", "taxonomic", "Aitchison")}
}
\description{
Generate distance matrices for multiple beta-diversity metrics. Returns a list of matrices.
}
\examples{
data(example_phyloseq)

distance.mats <- gen.dist.matrices(example_phyloseq, methods = "taxonomic", cores = 2)
}
\seealso{
\code{\link{phyloseq}}, \code{\link{capscale}}, \code{\link{vegdist}}, \code{\link{distance}}
}
