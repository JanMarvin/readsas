#' read.sas
#'
#'@author Jan Marvin Garbuszus \email{jan.garbuszus@@ruhr-uni-bochum.de}
#'
#'@param file file to read
#'@param debug print debug information
#'@param convert_dates default is TRUE
#'@param recode default is TRUE
#'@param rowcount number of rows to import (from start). negative values are
#' ignored
#'@param remove_deleted logical if deleted rows should be removed from data
#'
#'@useDynLib readsas, .registration=TRUE
#'@importFrom utils download.file
#'
#'@export
read.sas <- function(file, debug = FALSE, convert_dates = TRUE, recode = TRUE,
                     rowcount = -1, remove_deleted = TRUE) {

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
    recode <- FALSE


  if (convert_dates) {

    # TODO formula to create possible date formats?
    dates <- c(
      "b8601da", "e8601da", "date", "day", "ddmmyy", "ddmmyyb", "ddmmyyc",
      "ddmmyyd", "ddmmyyn", "ddmmyyp", "ddmmyys", "eurdfdd", "eurdfde",
      "eurdfdn", "eurdfdwn", "eurdfmy", "eurdfwdx", "eurdfmn", "eurdfwkx",
      "eurdfmn", "eurdfwkx", "weekdate", "weekdatx", "weekday", "downame",
      "worddate", "worddatx", "julday", "julian", "nengo", "pdjulg", "pdjuli",
      "yymm", "yymmc", "yymmd", "yymmn", "yymmp", "yymms", "yymmdd",
      "yymmddb", "yymmddc", "yymmddd", "yymmddn", "yymmddp", "yymmdds",
      "yymon", "yyq", "yyqc", "yyqd", "yyqp", "yyqs", "yyqn", "yyqr",
      "yyqrc", "yyqrd", "yyqrp", "yyqrs", "yyqrn", "year",
      "mmddyy", "mmddyyc", "mmddyyd", "mmddyyn", "mmddyyp",
      "mmddyys", "mmyy", "mmyyc", "mmyyd", "mmyyn", "mmyyp", "mmyys",
      "monname", "month", "monyy", "qtr", "qtrr"
    )

    vars <- which(toupper(formats) %in% toupper(dates))

    for (var in vars) {
      data[[var]] <- as.Date(
        as.POSIXct(data[[var]] * 24 * 60 * 60,
                   origin = "1960-01-01"),
      )
    }

    datetime <- c(
      "e8601dn", "e8601dt", "e8601dx", "e8601dz", "e8601lx", "b8601dn",
      "b8601dt", "b8601dx", "b8601dz", "b8601lx", "dateampm", "datetime",
      "dtdate", "dtmonyy", "dtwkdatx", "dtyear", "tod", "mdyampm"
    )

    vars <- which(toupper(formats) %in% toupper(datetime))

    for (var in vars) {
      data[[var]] <- as.POSIXct(data[[var]],
                                origin = "1960-01-01",
                                tz = "UTC")
    }


  }

  if (recode) {

    vars <- which(sapply(data, is.character))

    data[vars] <- mapply(iconv, data[vars], MoreArgs = list(from = encoding),
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

  created  <- as.POSIXct(created - created2,
                         origin = "1960-01-01", tz = "UTC")
  modified <- as.POSIXct(modified - modified2,
                         origin = "1960-01-01", tz = "UTC")
  thrdts   <- as.POSIXct(thrdts,
                         origin = "1960-01-01", tz = "UTC")

  attr(data, "created")  <- created
  attr(data, "modified") <- modified
  attr(data, "thrdts")   <- thrdts

  attr(data, "created2")  <- NULL
  attr(data, "modified2") <- NULL

  if (remove_deleted) {

    del_rows <- attr(data, "deleted_rows")

    val <- attr(data, "valid")
    del <- attr(data, "deleted")
    attr(data, "deleted") <- NULL
    attr(data, "valid") <- NULL

    if (del_rows > 0) {

      # deleted row in compressed data
      if (!all(val)) {

        if (del_rows != length(val[val == FALSE]))
          warning("number of deleted rows does not match the indicated number of deleted rows")

        data <- data[val, , drop = FALSE]

      } else {

        if (del_rows != length(del[del == TRUE]))
          warning("number of deleted rows does not match the indicated number of deleted rows")

        data <- data[!del, , drop = FALSE]
      }

    }

    # better save than sorry
    if (del_rows > 0 && (all(isFALSE(del)) || all(isTRUE(val))))
      warning("file indicated deleted rows, but none was found")

    if (del_rows == 0 && (any(del) || !all(val)))
      warning("file indicated no deleted rows, but some were found")
  }

  return(data)
}
