#' My Concatenate
#'
#' A formated version of `cat` for printing progress during the pipeline. Will rarely, if ever, be called directly by the user
#' @param x a string of text
#' @seealso \code{\link{cat}}
#' @export


my.cat <- function(x) { cat(paste0("\n### ", x, "\n"), sep = "\n") }
