#' Initiate Pipeline
#'
#' Start the Shaprton Lab dada2 pipeline: create output directories, create the run environment, and set ceratain variables
#' @export
#' @examples
#' initiate.pipeline()
#' ls()
#' list.files()

initiate.pipeline <- function() {
  dada.pkg.ver <- paste0("dada2_", packageVersion("dada2"))
  process.date <- Sys.Date()
  process.id <- paste(dada.pkg.ver, process.date, sep = "_")
  run.env$output.path <- paste0(process.id, "_output")
  if (!dir.exists(run.env$output.path)) {
    dir.create(run.env$output.path)
  }
}
