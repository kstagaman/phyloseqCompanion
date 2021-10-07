#' Read File
#'
#' An auxilliary function to read in metadata (or other) files in common format types
#' @param file A file in a common format: csv, tsv/txt, or rds
#' @seealso \code{\link{read.csv}}, \code{\link{read.table}}, \code{\link{readRDS}}
#' @export

read.file <- function(file) {
  if (str_detect(file, "\\.csv$")) {
    return(read.csv(file, row.names = 1))
  } else if (str_detect(file, "\\.tsv$|\\.txt$")) {
    return(
      read.table(
        file,
        sep = "\t",
        header = TRUE,
        row.names = 1
      )
    )
  } else if (str_detect(file, "\\.rds$")) {
    return(readRDS(file))

  } else {
    stop("Filetype not recognizable from extension (csv/tsv/txt/rds)")
  }
}
