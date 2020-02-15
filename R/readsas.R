#' read.sas
#'
#'@author Jan Marvin Garbuszus \email{jan.garbuszus@@ruhr-uni-bochum.de}
#'
#'@param file file to read
#'@param debug print debug information
#'@param convert.dates default is TRUE
#'@param recode default is TRUE
#'@param select.rows \emph{integer.} Vector of one or two numbers. If single
#' value rows from 1:val are selected. If two values of a range are selected
#' the rows in range will be selected.
#'@param select.cols \emph{character:} Vector of variables to select.
#'
#'@useDynLib readsas, .registration=TRUE
#'@importFrom utils download.file
#'
#'@export
read.sas <- function(file, debug = FALSE, convert.dates = TRUE, recode = TRUE,
                     select.rows = NULL, select.cols = NULL) {

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

  # some select.row checks
  if (!is.null(select.rows)) {
    # check that it is a numeric
    if (!is.numeric(select.rows)){
      return(message("select.rows must be of type numeric"))
    } else {
      # guard against negative values
      if (any(select.rows < 0) )
        select.rows <- abs(select.rows)

      # check that length is not > 2
      if (length(select.rows) > 2)
        return(message("select.rows must be of length 1 or 2."))

      # if length 1 start at row 1
      if (length(select.rows) == 1)
        select.rows <- c(1, select.rows)
    }
    # reorder if 2 is bigger than 1
    if (select.rows[2] < select.rows[1])
      select.rows <- c(select.rows[2], select.rows[1])

    # make sure to start at index position 1 if select.rows[2] > 0
    if (select.rows[2] > 0 & select.rows[1] == 0)
      select.rows[1] <- 1
  } else {
    # set a value
    select.rows <- c(0,0)
  }

  if (is.null(select.cols)){
    select.cols <- ""
  }

  data <- readsas(filepath, debug, select.rows, select.cols)

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

  created    <- attr(data, "created")
  modified   <- attr(data, "modified")
  created2   <- attr(data, "created2")
  modified2  <- attr(data, "modified2")
  thrdts     <- attr(data, "thrdts")

  created  <- as.POSIXct( created-created2,  origin = "1960-01-01")
  modified <- as.POSIXct( modified-modified2, origin = "1960-01-01")
  thrdts   <- as.POSIXct( thrdts,   origin = "1960-01-01", "GMT")

  attr(data, "created")  <- created
  attr(data, "modified") <- modified
  attr(data, "thrdts")   <- thrdts

  attr(data, "created2")  <- NULL
  attr(data, "modified2") <- NULL


  return(data)
}
