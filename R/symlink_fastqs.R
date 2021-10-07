#' Symlink Fastqs
#'
#' A function to create symlinks with easier to parse names that the typical names that files have coming off the sequencing machine
#' @param seq.dir The path to the directory in which the raw sequence files are contained.
#' @param ids.tbl A data.frame or data.table containing sample names in one column and their corresponding file IDs in another (e.g. barcodes or sample IDs as used in the raw sequence files). Sample names and file IDs *may* be identical in some instances, please still provide a two column table.
#' @param smpl.id.col The name/number of the column in the `ids.tbl` that contains the name of the samples, which will be used in naming the symlinks. Default is 'Sample'.
#' @param file.id.col The name/number of the column in the `ids.tbl` that contains the corresponding file IDs for matching raw sequence files. Default is 2.
#' @param split.pattern A character string containg the pattern you want to use to separate the sample name and R1/R2 in the symlink names. Default is '--'.
#' @param quiet Logical, if TRUE, there will be no printing of progress. Default is FALSE
#' @export

symlink.fastqs <- function(
  seq.dir,
  ids.tbl,
  smpl.id.col = "Sample",
  file.id.col = 2,
  split.pattern = "--",
  quiet = FALSE
) {
  if ("data.table" %in% class(ids.tbl)) {
    ids.dt <- ids.tbl
  } else {
    ids.dt <- as.data.table(ids.dt)
  }
  setkeyv(ids.dt, smpl.id.col)
  for (smpl in ids.dt[[smpl.id.col]]) {
    file.id <- ids.dt[smpl][[file.id.col]]
    if (!quiet) {
      cat(paste0(file.id, " -> ", smpl, " ... "), sep = "")
    }
    seq.files <- list.files(
      path = seq.dir,
      pattern = file.id,
      full.names = T
    )
    if (length(seq.files) == 0) {
      stop("no files matched file ID")
    } else {
      if (!quiet) { cat("good", sep = "\n") }
    }
    read1.file <- seq.files[str_detect(seq.files, "_R1_")]
    read2.file <- seq.files[str_detect(seq.files, "_R2_")]
    lnName.r1 <- paste(smpl, "R1.fastq.gz", sep = split.pattern)
    lnName.r2 <- paste(smpl, "R2.fastq.gz", sep = split.pattern)
    cmd1 <- paste("ln -s", read1.file, lnName.r1)
    cmd2 <- paste("ln -s", read2.file, lnName.r2)
    system(cmd1)
    system(cmd2)
  }
}
