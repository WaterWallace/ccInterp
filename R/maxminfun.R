#' MaxMin Function
#'
#' function for looking up the maximum y-value in a reference dataset, for a time at time of newts
#'
#' newts will be a forward looking mean, i.e. 10:00, 11:00, 12:00.. 10:00 will be the mean between 10:00 and 11:00
#' And the function result will be the maximum in the reference period between 10:00 and 11:00
#' This is useful for quality coded data that has been averaged.
#' the maximum quality code can be extracted from the original data and applied to the averaged data.
#' also useful for a daily maximum value
#' or a daily minimum value
#'
#' @param refx source x input.
#' @param refy source y input.
#' @param newts target x interval i.e. hourly/daily.
#' @param option "max" or "min" to return interval max or interval min.
#' @param dt data type, 1 = inst, 2 = forward mean, 6 = insttotal
#'
#' @return  dataframe with newts xaxis, and a forward looking max or min of refy.
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
#' intervalMax <- maxminfun(longts, randts, intervalts)
#' f.max <- approxfun(intervalMax, method = "constant")
#' lines(f.max(longts), col = "red")
#'
#' intervalMin <- maxminfun(longts, randts, intervalts, option = "min")
#' f.min <- approxfun(intervalMin, method = "constant")
#' lines(f.min(longts), col = "red")
#'
#' @export
maxminfun <- function(refx, refy, newts, option = "max", dt = 1) {
  originalnewts <- newts
  newts <- as.data.frame(newts)
  oldts <- data.frame(refx, refy)

  if(dt == 3) # trailing mean
  {
    offs <- median(diff(as.numeric(newts[,1])))
    newts[,1] <- newts[,1] + offs
  }

  if(dt == 1 | dt == 6) {
    # max neighbor
    lagged <- oldts %>% mutate(lag = lag(refy,1)) %>% mutate(lead = lead(refy,1))
    if (option == "max"){
      oldts$refy <- pmax(lagged$lag, lagged$lead, na.rm = TRUE)
    }else if (option == "min"){
      oldts$refy <- pmin(lagged$lag, lagged$lead, na.rm = TRUE)
    }else{
      stop("option must be min or max")
    }
  }

  # Function of hourly hours
  f.hourlyHours <- approxfun(newts[, 1], newts[, 1], method = "constant", f = 0, rule = 2)
  # which hour does this belong
  oldts <- cbind(oldts, hour = f.hourlyHours(oldts[, 1]))
  # add max QC for period to df
  newts <- cbind(newts, QC = maxfun(oldts$hour, oldts[, 2], newts[, 1], option = option))

  if(dt == 3) # trailing mean
  {
    newts[,1] <- newts[,1] - offs
  }

  return(newts)
}



