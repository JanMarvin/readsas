#' read.sas
#'
#'@author Jan Marvin Garbuszus \email{jan.garbuszus@@ruhr-uni-bochum.de}
#'
#'@param file file to read
#'@param debug print debug information
#'
#'@useDynLib readsas, .registration=TRUE
#'@importFrom utils download.file
#'
#'@export
read.sas <- function(file, debug = FALSE) {

  # Check if path is a url
  if (length(grep("^(http|ftp|https)://", file))) {
    tmp <- tempfile()
    download.file(file, tmp, quiet = TRUE, mode = "wb")
    filepath <- tmp
    on.exit(unlink(filepath))
  } else {
    # construct filepath and read file
    filepath <- get.filepath(file)
  }
  if (!file.exists(filepath))
    return(message("File not found."))



  data <- readsas(filepath, debug)


  return(data)

}
