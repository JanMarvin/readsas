#' write.sas
#'
#'@author Jan Marvin Garbuszus \email{jan.garbuszus@@ruhr-uni-bochum.de}
#'
#'@param dat data frame to save
#'@param filepath path to save file to
#'@param compress option to compress file
#'@param debug print debug information
#'
#'@useDynLib readsas, .registration=TRUE
#'
#'@export
write.sas <- function(dat, filepath, compress = 0, debug = FALSE) {

  filepath <- path.expand(filepath)

  vartypes <- sapply(dat, is.character) + 1
  colwidth <- sapply(dat, function(x)max(nchar(x)))
  colwidth[vartypes == 1] <- 8

  labels <- "testlab";

  # for numerics
  formats <- rep("BEST", ncol(dat))

  attr(dat, "vartypes") <- as.integer(vartypes)
  attr(dat, "colwidth") <- as.integer(colwidth)
  attr(dat, "formats") <- formats
  attr(dat, "labels") <- labels


  writesas(filepath, dat, compress = 0, debug = debug)

}
