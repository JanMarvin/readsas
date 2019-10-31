#' read.sas
#'
#'@author Jan Marvin Garbuszus \email{jan.garbuszus@@ruhr-uni-bochum.de}
#'
#'@param file file to read
#'@param debug print debug information
#'@param convert.dates default is TRUE
#'@param recode default is TRUE
#'@param rowcount number of rows to import (from start). negative values are
#' ignored
#'
#'@useDynLib readsas, .registration=TRUE
#'@importFrom utils download.file
#'
#'@export
read.sas <- function(file, debug = FALSE, convert.dates = TRUE, recode = TRUE,
                     rowcount = -1) {

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



  data <- readsas(filepath, debug, rowcount)

  labels  <- attr(data, "labels")
  formats <- attr(data, "formats")
  encoding <- attr(data, "encoding")

  # override encoding argument if file contains no valid information
  if (encoding == "")
    recode = FALSE


  if (convert.dates) {

    vars <- which(formats == "MMDDYY" | formats == "DATE")

    for (var in vars) {
      data[[var]] <- as.Date(
        as.POSIXct( data[[var]] * 24 * 60 * 60, origin = "1960-01-01")
        )
    }

  }

  if (recode) {

    vars <- which(sapply(data, is.character))

    data[vars] <- mapply(iconv, data[vars], MoreArgs=list(from = encoding),
                         SIMPLIFY = FALSE)

    labels <- iconv(labels, from = encoding)
    attr(data, "labels") <- labels

    names(data) <- iconv(names(data), from = encoding)

  }

  created   <- attr(data, "created")
  modified  <- attr(data, "modified")
  created2   <- attr(data, "created2")
  modified2  <- attr(data, "modified2")
  thrdts    <- attr(data, "thrdts")

  created  <- as.POSIXct( created-created2,  origin = "1960-01-01")
  modified <- as.POSIXct( modified-modified2, origin = "1960-01-01")
  thrdts   <- as.POSIXct( thrdts,   origin = "1960-01-01", "GMT")

  attr(data, "created")  <- created
  attr(data, "modified") <- modified
  attr(data, "thrdts")   <- thrdts


  return(data)
}
