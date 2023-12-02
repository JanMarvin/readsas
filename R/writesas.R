#' write.sas
#'
#'@author Jan Marvin Garbuszus \email{jan.garbuszus@@ruhr-uni-bochum.de}
#'
#'@param dat data frame to save
#'@param filepath path to save file to
#'@param compress option to compress file
#'@param debug print debug information
#'@param bit32 write 32bit file
#'
#'@useDynLib readsas, .registration=TRUE
#'
#'@export
write.sas <- function(dat, filepath, compress = 0, debug = FALSE, bit32 = FALSE,
                      varlabels, size) {

  filepath <- path.expand(filepath)

  if (missing(varlabels)) {
    varlabels <- rep("", ncol(dat))
  }

  if (missing(size)) {

    hsize <- 65536
    psize <- 65536

    if (bit32) {
      hsize <- 1024
      psize <- 8192
    }

    size = c(hsize, psize)
  } else if (length(size) != 2) {
    size = c(size[1], size[1])
  }

  if (!inherits(dat, "data.frame")) {
    dat <- as.data.frame(dat)
  }

  # convert from factor
  ff <- sapply(dat, is.factor)
  dat[ff] <- lapply(dat[ff], as.character)

  vartypes <- sapply(dat, is.character) + 1
  colwidth <- sapply(dat, function(x) max(nchar(x)))
  colwidth[vartypes == 1] <- 8

  labels <- "testlab"

  vartypen <- sapply(dat, is.numeric)

  formats <- NA
  formats[vartypen] <- "BEST"
  formats[!vartypen] <- ""

  width <- 0
  width[vartypen] <- 32 # fix for now
  width[!vartypen] <- 0 # colwidth[!vartypen]

  decim <- sapply(dat, is.integer)
  decim[!vartypen] <- TRUE

  is.Date <- function(x) inherits(x, "Date")
  is.POSIX <- function(x) inherits(x, "POSIXt")

  vartypen <- sapply(dat, is.Date)
  dat[vartypen] <- lapply(dat[vartypen], as_date)
  formats[vartypen] <- "DATE"
  width[vartypen] <- 9

  vartypen <- sapply(dat, is.POSIX)
  dat[vartypen] <- lapply(dat[vartypen], as_datetime)
  formats[vartypen] <- "DATETIME"
  width[vartypen] <- 22
  decim[vartypen] <- 3

  if (debug) {
    message("vartypes")
    print(vartypes)
    message("colwidth")
    print(colwidth)
    message("formats")
    print(formats)
    message("width")
    print(width)
    message("decim")
    print(decim)
    message("labels")
    print(labels)
  }

  # for numerics
  # formats <- rep("BEST", ncol(dat))

  attr(dat, "vartypes") <- as.integer(vartypes)
  attr(dat, "colwidth") <- as.integer(colwidth)
  attr(dat, "formats")  <- formats
  attr(dat, "width")    <- width
  attr(dat, "decim")    <- decim
  attr(dat, "labels")   <- labels
  attr(dat, "varlabels") <- varlabels


  writesas(filepath, dat, compress = 0, debug = debug, bit32 = bit32,
           headersize = size[1], pagesize = size[2], dateval = as_datetime(Sys.time()))

}
