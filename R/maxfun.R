#' MaxMin Function
#'
#' Simply a wrapper for approxfun, returning the maximum value forward
#' and backwards from a lookup table.
#'
#' @param refx source x input.
#' @param refy source y input.
#' @param targetx target x interval i.e. hourly/daily.
#' @param option "max" or "min" to return interval max or interval min.
#'
#' @return  vector of quality codes or such
#'
#' @examples
#'
#' longts <- seq(0, 500, by = 1)
#' # create a random timeseries
#' randts <- runif(length(longts), min = 0, max = 1000)
#'
#' plot(longts, randts)
#' lines(longts, randts, col = "grey")
#'
#' # new interval timeseries
#' intervalts <- seq(0, 500, by = 20)
#'
#' intervalMax <- maxfun(longts, randts, intervalts)
#' f.max <- approxfun(intervalMax, method = "constant")
#' lines(f.max(longts), col = "red")
#'
#' intervalMin <- maxfun(longts, randts, intervalts, option = "min")
#' f.min <- approxfun(intervalMin, method = "constant")
#' lines(f.min(longts), col = "red")
#'
#' @export
maxfun <- function(t, y, x, option = "max") # look forward and back to pick the max (or min)
{
  # function looks to either side of point x and outputs the minimum
  f.qual0 <- approxfun(t, y, method = "constant", ties = option, f = 0, rule = 2)
  f.qual1 <- approxfun(t, y, method = "constant", ties = option, f = 1, rule = 2)

  if (option == "max") {
    return(pmax(f.qual0(x), f.qual1(x)))
  } else if (option == "min") {
    return(pmin(f.qual0(x), f.qual1(x)))
  }

}





