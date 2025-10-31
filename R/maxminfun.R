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
#' @return  vector of interpolated ref data at time of newts
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
maxminfun <- function(refx, refy, newts, option = "max", dt = 2, inputdt = 1) {
  # dt is output datatrans, default 2 (forward mean)
  # assuming you're inputting instantaenous data
  originalnewts <- newts
  newts <- as.data.frame(newts)
  oldts <- data.frame(refx, refy)

  # outputting instantaneous data isn't really possible, it's just a forward
  #mean that has been offset in time.
  if (dt == 1) {
    offs <- median(diff(as.numeric(newts[,1])))
    newts[,1] <- newts[,1] - offs / 2
  }
  if(inputdt == 3) # trailing mean
  {
    offs <- median(diff(as.numeric(newts[,1])))
    newts[,1] <- newts[,1] + offs
  }

  if(inputdt == 1 | inputdt == 6) {
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
  targetfloor <- approx(newts[,1], lag(newts[,1],1), oldts[,1], f=1, method = "constant", rule = 2)$y %>% as.POSIXct() #floor of target timeseries
  targetceiling <- approx(newts[,1], lead(newts[,1],1), oldts[,1], f=0, method = "constant", rule = 2)$y %>% as.POSIXct() #ceiling of target timeseries

  if(dt == 1 | dt == 6 | dt == 2)
  {
    targethour <- targetfloor
  }else if(dt == 3)
  {
    targethour <- targetceiling
  }

  # which hour does this belong
  oldts$hour <- targethour

  # add max QC for period to df
  QC <- maxfun(oldts$hour, oldts[, 2], newts[, 1], option = option)

  if(inputdt == 3) # trailing mean
  {
    newts[,1] <- newts[,1] - offs
  }

  return(QC)
}



