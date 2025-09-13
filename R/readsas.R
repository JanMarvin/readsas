#' read.sas
#'
#' @description `read.sas` is a general function for reading sas7bdat files.
#' It supports a variety of complex sas7bdat files, x86 and x64, big and small
#' endian, and both compression types. It has been tested with numeric and
#' character data, and provides helper functions for converting sas7bdat to R
#' types `Date` and `POSIXct`. Time variables are converted to characters.
#' Conversion to date variables is applied if a known date, datetime or time
#' format is found in the sas7bdat file. For user-defined formats, the package
#' provides functions to convert from sas7bdat to R.
#'
#' Input files may contain deleted rows that are marked as deleted instead of
#' being removed from the input data. These are removed on import, if you still
#' need them look at `remove_deleted`. Formats, labels and additional file
#' information are available with `attributes()`.
#'
#' @param file file to read
#' @param debug print debug information
#' @param convert_dates default is `TRUE`
#' @param recode default is `TRUE`
#' @param select.rows \emph{integer.} Vector of rows to import. Minimum 0. Rows
#' imported are sorted. If 0 is in `select.rows`, zero rows are returned.
#' @param select.cols \emph{character:} Vector of variables to select.
#' @param remove_deleted logical if deleted rows should be removed from data
#' @param rownames first column will be used as rowname and removed from data
#' @param empty_to_na logical. In SAS empty characters are missing. this option
#' allows to convert `""` to `NA_character_` when importing.
#' @param convert logical convert missings `.I` and `.M` to Inf and -Inf
#'
#' @useDynLib readsas, .registration=TRUE
#' @importFrom utils download.file
#' @importFrom stringi stri_encode
#'
#' @seealso \link[foreign]{read.xport}
#'
#' @examples
#' fl <- system.file("extdata", "cars.sas7bdat", package = "readsas")
#' read.sas(fl)
#'
#' @export
read.sas <- function(file, debug = FALSE, convert_dates = TRUE, recode = TRUE,
                     select.rows = NULL, select.cols = NULL, remove_deleted = TRUE,
                     rownames = FALSE, empty_to_na = FALSE, convert = FALSE) {

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
    if (!is.numeric(select.rows)) {
      return(message("select.rows must be of type numeric"))
    } else {
      if (any(select.rows < 0)) stop("select.rows must be >= 0")

      select.rows <- sort(select.rows - 1)
    }

  }

  if (!is.null(select.cols) && !is.character(select.cols)) {
    return(message("select.cols must be of type character"))
  }

  tempstr <- tempfile()
  on.exit(unlink(tempstr), add = TRUE)
  data <- readsas(filepath, debug, select.rows, select.cols, empty_to_na,
                  tempstr, convert)

  cvec <- ifelse(rownames, -1, substitute())

  # rownames start at 0
  row_names <- attr(data, "rvec") + 1

  # only if a row was returned
  if (nrow(data)) row.names(data) <- row_names

  if (rownames) {
    rownames(data) <- data[[1]]
    data[[1]] <- NULL
  }


  encoding <- attr(data, "encoding") <- attr(data, "encoding")

  ## shorten attributes and reassign
  labels  <- attr(data, "labels")  <- attr(data, "labels")[cvec]
  formats <- attr(data, "formats") <- attr(data, "formats")[cvec]

  attr(data, "varnames") <- attr(data, "varnames")[cvec]
  attr(data, "fmtkeys")  <- attr(data, "fmtkeys")[cvec]
  attr(data, "fmt32")    <- attr(data, "fmt32")[cvec]
  attr(data, "ifmt32")   <- attr(data, "ifmt32")[cvec]
  attr(data, "colwidth") <- attr(data, "colwidth")[cvec]
  attr(data, "vartyps")  <- attr(data, "vartyps")[cvec]
  attr(data, "c8vec")    <- attr(data, "c8vec")[cvec]


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
      data[[var]] <- convert_to_date(data[[var]])
    }

    datetime <- c(
      "e8601dn", "e8601dt", "e8601dx", "e8601dz", "e8601lx", "b8601dn",
      "b8601dt", "b8601dx", "b8601dz", "b8601lx", "dateampm", "datetime",
      "dtdate", "dtmonyy", "dtwkdatx", "dtyear", "tod", "mdyampm"
    )

    vars <- which(toupper(formats) %in% toupper(datetime))

    for (var in vars) {
      data[[var]] <- convert_to_datetime(data[[var]])
    }

    time <- c(
      "time", "timeampm", "tod", "hhmm", "hour", "mmss", "systime"
    )

    vars <- which(toupper(formats) %in% toupper(time))

    for (var in vars) {
      data[[var]] <- convert_to_time(data[[var]])
    }

  }

  if (recode) {

    vars <- which(sapply(data, is.character))

    data[vars] <- mapply(stringi::stri_encode, data[vars], MoreArgs = list(from = encoding),
                         SIMPLIFY = FALSE)

    labels <- stringi::stri_encode(labels, from = encoding)
    attr(data, "labels") <- labels

    names(data) <- stringi::stri_encode(names(data), from = encoding)

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

    sel <- row_names
    del_rows <- attr(data, "deleted_rows")

    val <- attr(data, "valid")[sel]
    del <- attr(data, "deleted")[sel]
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

    # better safe than sorry
    if (del_rows > 0 && (all(isFALSE(del)) || all(isTRUE(val))))
      warning("file indicated deleted rows, but none was found")

    if (del_rows == 0 && (any(del) || !all(val)))
      warning("file indicated no deleted rows, but some were found")
  }

  if (!debug) {
    attr(data, "cvec") <- attr(data, "rvec") <- NULL
  }

  return(data)
}

#' helper function to convert SAS date numeric to date
#' @param x date or datetime variable
#' @examples
#'  # 2000-03-17
#'  convert_to_date(14686)
#' @name converttimedate
#' @export
convert_to_date <- function(x) {

  # leap year bug in year 4000/8000
  feb29_4k <- 745154
  feb29_8k <- 2206123
  x[x >= feb29_4k] <- x[x >= feb29_4k] + 1
  x[x >= feb29_8k] <- x[x >= feb29_8k] + 1

  as.Date(
    as.POSIXct(x * 24 * 60 * 60,
               origin = "1960-01-01")
  )
}

#' @rdname converttimedate
#' @examples
#'  # 2012-11-10 03:49:19
#'  convert_to_datetime(1668138559)
#' @export
convert_to_datetime <- function(x) {

  # leap year bug in year 4000/8000
  feb29_4k <- 64381305600
  feb29_8k <- 190609027200
  x[x >= feb29_4k] <- x[x >= feb29_4k] + 24 * 60 * 60
  x[x >= feb29_8k] <- x[x >= feb29_8k] + 24 * 60 * 60

  as.POSIXct(x,
             origin = "1960-01-01",
             tz = "UTC")
}

#' @rdname converttimedate
#' @examples
#'  # 04:04:46
#'  convert_to_time(14686)
#' @export
convert_to_time <- function(x) {
  format(convert_to_datetime(x), format = "%H:%M:%S")
}
