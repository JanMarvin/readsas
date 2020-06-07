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

  writesas(filepath, dat, compress = 0, debug = debug)

}
