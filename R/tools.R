#' Construct File Path
#'
#' @param path path to dta file
#' @author Jan Marvin Garbuszus \email{jan.garbuszus@@ruhr-uni-bochum.de}
#' @author Sebastian Jeworutzki \email{sebastian.jeworutzki@@ruhr-uni-bochum.de}
#' @keywords internal
#' @noRd
get.filepath <- function(path = "") {
  if (substring(path, 1, 1) == "~") {
    filepath <- path.expand(path)
  } else {
    filepath <- path
  }
  if (!file.exists(filepath)) {
    return("File does not exist.")
  }

  return(filepath)
}
