#' read.sas
#'
#'@author Jan Marvin Garbuszus \email{jan.garbuszus@@ruhr-uni-bochum.de}
#'
#'@param file file to read
#'@param debug print debug information
#'@param convert.dates default is TRUE
#'
#'@useDynLib readsas, .registration=TRUE
#'@importFrom utils download.file
#'
#'@export
read.sas <- function(file, debug = FALSE, convert.dates = TRUE) {

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

  labels  <- attr(data, "labels")[1:ncol(data)]
  formats <- attr(data, "formats")[1:ncol(data)]

  if (convert.dates) {

    vars <- which(formats == "MMDDYY" | formats == "DATE")
    # z <- 1472562988

    # no leapseconds applied (is it required?)
    for (var in vars) {
      data[[var]] <- as.Date(
        as.POSIXct( data[[var]] * 24 * 60 * 60, origin = "1960-01-01")
        )
    }

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
