#' Generate Pipeline Template
#'
#' Generate an R script template for running the Sharpton Lab dada2 pipeline
#' @param filename A string to name the resulting template. If the .R extension is not included in the filename, it will be appended. Defaut is "dada2_processing.R"
#' @export
#' @examples
#' generate.template()
#' file.show("dada2_processing.Rmd")

generate.template <- function(filename = "dada2_processing.R") {
  if (!grepl("\\.R$", filename)) { filename <- paste0(filename, ".R")}
  writeLines(template.code, con = filename)
}
