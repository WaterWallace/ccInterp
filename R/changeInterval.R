#' Change Interval
#'
#' Trapezoidal integration of a time series
#'
#' not appropriate for downsampling intervals from hourly back to daily
#' Use changeInterval() for that.
#' (even though theoretically spinterpConvert will work fine for downsampling, and in some cases technically preferred)
#' Using "spinterp" Interpolation over long gaps could have unpredicatable results.
#'
#' @param ts dataframe of posixct time (time in seconds) and an instantaneous value (per second i.e. cumecs), and optionally a quality code
#' @param dt datatype of the input data, 1 = inst; 2 = fmean
#' default is 1 for instantaenous, the inputs can be randomly spaced in time.
#' 2 for fmean when you are converting evenly spaced forward means such as daily averages.
#' @param start a timestamp written, that can be converted to posixct.
#' i.e. "2016-12-17 05:00" or numeric, seconds since "1970-01-01" i.e. 1481914800
#' @param end same as start
#'
#' @param offset minutes to offset the averaging window, not the result
#' @param option "fmean", "inst", "sum"
#' @param linearInterp adds linearly interpolated points at new interval timestep
#' leave as true, unpredictable result with false, particularly over long gaps
#' @param rounded round the time to the hour
#'
#' @return  dataframe with tiemstamp as posixct a trapezoidal integrated rate and quality code if included.
#'
#' @examples
#'
#'
#' data(xy_data)
#' D <- within(xy_data, {
#'   x <- as.POSIXct(x * 60 * 60 * 24, origin = "1970-01-01")
#' })
#'
#' plot(D[50:100, ])
#'
#' # daily average of forward means
#' daily <- changeInterval(D, Interval = "Daily")
#' head(daily)
#' lines(daily)
#'
#'
#' # the same interval, but without automatically adjusting the offset to output by using "Daily" as interval
#' twentyfourhour <- changeInterval(D, Interval = 24 * 60)
#' head(twentyfourhour)
#'
#'
#' # output an hourly forward mean (timestamp at beginning of period ),
#' hourly <- changeInterval(D, Interval = "Hourly", option = "fmean")
#' f.hourly <- approxfun(hourly$Date, hourly$FMean, method = "constant")
#' ts <- seq(hourly$Date[1], hourly$Date[nrow(hourly)], by = 10 * 60)
#' lines(ts, f.hourly(ts), col = "red")
#'
#' # output an hourly instanaenous (timestamp in middle of averaging period ),
#' # offset average window by 30 minutes to output datapoints on the hour.
#' # essentially it is the average of the period 30 mins before and 30 mins after the timestamp.
#' hourlyInst <- changeInterval(D, Interval = "Hourly", offset = 0, option = "inst")
#' head(hourlyInst)
#' lines(hourlyInst, col = "orange")
#' @export


changeInterval <- function(ts, dt = 1, Interval = "Daily", start = 0,
                           end = 0, offset = 0, option = "fmean",
                           linearInterp = TRUE, rounded = TRUE) {
  inputts <- ts

  # function omits na's, so will interpolate between
  ts <- na.omit(ts)
  ts <- ts[, 1:2]

  cullBefore <- ts[1, 1]
  if (start == 0) {
    if (Interval == "Daily") {
      # round to full day
      # if (ts[1,1] >  round(ts[1,1], units="days") ) cullBefore = ts[1,1]
      start <- round(ts[1, 1], units = "days")
    } else if (Interval == "Hourly") {
      # if (ts[1,1] >  round(ts[1,1], units="days") ) cullBefore = ts[1,1]
      # round to full hour
      start <- round(ts[1, 1], units = "hours")
    } else {
      # start at raw value

      if (rounded == TRUE) {
        # if (ts[1,1] >  round(ts[1,1], units="days") ) cullBefore = ts[1,1]
        start <- round(ts[1, 1], units = "hours")
      } else {
        start <- ts[1, 1]
      }
    }
  } else {
    start <- as.POSIXct(start, origin = "1970-01-01")
  }
  if (end == 0) {
    end <- ts[nrow(ts), 1]
  } else {
    end <- as.POSIXct(end, origin = "1970-01-01")
  }


  # set time step based on interval
  if (Interval == "Daily") {
    timestep <- 24 * 60 * 60
  } else if (Interval == "Hourly") {
    timestep <- 1 * 60 * 60
  } else {
    # assume an integer
    timestep <- Interval * 60
  }

  if (linearInterp) {
    # add new points at timestamp boundaries
    newintTS <- seq(start + offset * 60, end, by = timestep)
    f.timeties <- approxfun(ts[, 1], ts[, 2])
    lineardf <- data.frame(newintTS, f.timeties(newintTS))
    # plot(lineardf)
    # merge lineardf and ts
    names(lineardf) <- names(ts)
    merged <- rbind(lineardf, ts)
    merged <- merged[order(merged[, 1]), ]
    ts <- na.omit(distinct(merged))
  }

  # duration between this point and the previous point
  ts$dur <- c(0, diff(as.numeric(ts[, 1]))) # duration in seconds

  # calculate a volume

  if (dt == 1) {
    # same as above, except getting the average of
    # point data, or interval data (means, with value at centre point)
    averageRate <- (0.5 * (c(ts[, 2], 0) + c(0, ts[, 2])))
    kilolitres <- c(ts$dur, 0) * averageRate # duration in seconds times average of timeseries values in rate per second
    kilolitres <- head(kilolitres, -1)

    ts$megalitres <- kilolitres / 1000
  }
  if (dt == 2) {
    # for forward mean data
    kilolitres <- ts$dur * ts[, 2] # multiply duration in seconds by cubic metres per second, to equal total cubic metres (kilolitres)
    ts$megalitres <- kilolitres / 1000 # divide by 1000 to equal total megalitres
  }

  # cumulative sum
  ts$accum <- cumsum(ts$megalitres)


  f.linear <- approxfun(ts[, 1], ts$accum)
  f.accum <- splinefun(ts[, 1], ts$accum, method = "fmm")


  # new time sequence for new interval
  newintTS <- seq(start + offset * 60, end, by = timestep)

  # delta y
  if (linearInterp) {
    newTS <- f.linear(newintTS)
  } else {
    newTS <- f.accum(newintTS)
  }

  dy <- c(diff(newTS), 0)


  # delta y / delta time
  if (option == "fmean") {
    df <- data.frame(Date = newintTS, FMean = round(dy / (timestep / 1000), 3)) # convert back from ML/day to cumecs
  } else if (option == "sum") {
    # output total difference in cumulative
    # convert back to KL i.e. cubic metres
    df <- data.frame(Date = newintTS, Sum = round(newTS * 1000, 3))
  } else if (option == "inst") {
    # add half a timestep to report an instantaneous mean
    df <- data.frame(Date = newintTS + (timestep / 2), Inst = round(dy / (timestep / 1000), 3))
  }

  if (length(inputts) == 3) {
    df <- maxminfun(inputts[, 1], inputts[, 3], df, option = "max")
  }

  # tidy up
  df <- df[df[, 1] >= cullBefore, ]
  df <- head(df, -1) # remove redundant zero value added as dy
  df <- na.trim(df)
}
