#' Sum Unique Sequences
#'
#' Get the number of unique sequences for each sample, an auxilliary function primarily called by other functions in the pipeline.
#' @param x A sample with associated sequences
#' @export

getN <- function(x) sum(getUniques(x))
