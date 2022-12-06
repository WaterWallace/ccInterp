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
#'
#' @return  dataframe with newts xaxis, and a forward looking max or min of refy.
#'
#' @examples
#'
#' longts <- seq( 0, 500, by=1 )
#' # create a random timeseries
#' randts <- runif(length(longts), min = 0, max = 1000)
#'
#' plot(longts, randts)
#' lines(longts, randts, col="grey")
#'
#' # new interval timeseries
#' intervalts <- seq( 0, 500, by=20 )
#'
#' intervalMax <- maxminfun(longts, randts, intervalts)
#' f.max <- approxfun(intervalMax, method="constant")
#' lines(  f.max(longts), col="red" )
#'
#' intervalMin <- maxminfun(longts, randts, intervalts, option="min")
#' f.min <- approxfun(intervalMin, method="constant")
#' lines(  f.min(longts), col="red" )
#'
#'
#' @export
maxminfun <- function(refx, refy, newts, option="max", dt=2)
{

  originalnewts <- newts
  newts <- as.data.frame(newts)
  if (dt==1)
  {
    subtract <- (newts[2,1] -newts[1,1])/2
    #newts<- as.numeric(newts)
    #shift the target timeseries back by half an interval
    # function assumes default dt of 2, which is forward means.
    newts[,1] <- newts[,1] - (subtract )
  }

  #refx - reference x lookup
  #refy - reference y lookup i.e. a quality code
  #newts - a timeseries dataframe, with at least a column of timestamps
  #option - "max" or "min"

  maxfun <- function(t, y, x, option="max") #look forward and back to pick the max (or min)
  {

    # function looks to either side of point x and outputs the minimum
    f.qual0 <-  approxfun(t, y, method = "constant", ties=option,  f=0)
    f.qual1 <- approxfun(t, y, method = "constant", ties=option, f=1)

    if (option=="max") {
      return( pmax(f.qual0(x),  f.qual1(x)) )
    }else if(option=="min")
      return( pmin(f.qual0(x),  f.qual1(x)) )
  }

  oldts <- data.frame(refx, refy)
  #newts <- data.frame(Timestamp = newts ) # why did I add this?
  #colnames(newts)[1]

  # Function of hourly hours
  f.hourlyHours <- approxfun(newts[,1], newts[,1], method="constant")
  # which hour does this belong
  oldts <- cbind (oldts, hour = f.hourlyHours(oldts[,1]) )
  # add max QC for period to df
  newts <- cbind(newts, QC = maxfun(oldts[,3], oldts[,2], newts[,1], option=option))

  #f.newts <- approxfun(newts, method="constant")
  #return(f.newts)

  if (dt==1)
  {
    newts[,1] <-  newts[,1]  + subtract
  }

  return(newts)

}




