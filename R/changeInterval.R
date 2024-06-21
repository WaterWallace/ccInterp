#' Change Interval
#'
#' Trapezoidal integration of a time series
#'
#' For downsampling intervals, for example hourly to daily.
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
#' @param option "fmean", "inst", "sum", "resample"
#'     fmean Calcaulte a forward mean, i.e. the daily mean for this day
#'     inst Calculate instantaenous values
#'     sum Accumulate daily i.e. to calculate a total or volume, units unchanged
#'     resample interpolate at intervals, no transformation
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
                           rounded = TRUE)
{

  stopifnot("duplicate timestamps" = length(ts[,1]) == length(unique(ts[,1])) )

  names(ts)[1] <- "Date"
  names(ts)[2] <- "value"

  ts$Date <- as.POSIXct(ts$Date)
  inputts <- ts

  if (start == 0)
  {
    if (Interval == "Daily")
    {
      start <- round.POSIXt(ts[1, 1], units = "days")
      Interval = 60*60*24
    }else if ( Interval == "Hourly")
    {
      start <- round.POSIXt(ts[1, 1], units = "hours")
      Interval = 60*60
    }else if (rounded == TRUE){
      start <- round.POSIXt(ts[1, 1], units = "hours")
      Interval = Interval*60
    }else{
      start <- ts[1, 1]
      Interval = Interval*60
    }
  }else{
    # start should already be posixct time
  }

  if (end == 0) {
    end <- ts[nrow(ts), 1]
  } else {
    end <- as.POSIXct(end, origin = "1970-01-01")
  }

  # set to numeric
  ts <- ts %>% mutate( numDate = as.numeric(Date))

  #range(ts$Date) %>% round("hours")

  # new time sequence
  newintTS <- seq(start+offset*60, end , by = Interval)

  if(dt == 1){
    # merge new timestamps into original dataset with linear interpolation
    f.timeties <- approxfun(ts$Date, ts$value)
    linearts <- data.frame(Date = newintTS, value = f.timeties(newintTS))
    merged <- rbind(linearts, dplyr::select(ts, c(Date, value)))
    merged <- merged[order(merged$Date), ]
    ts <- na.omit(distinct(merged))

    # convert rate to cumulative volume
    #ts <- ts %>% mutate( accum = cumsum(value * c(0, diff(numDate))))
    ts <- ts %>% mutate( accum = cumtrapz(Date %>% as.numeric, value) )

    # linear lookup accum
    f.accum <- approxfun(ts$Date, ts$accum)

    # new dataframe for output timestep
    newts <- data.frame(Date = newintTS,
                        accum = f.accum(newintTS))


  }else{
    stopifnot("other than dt of 1 must be interval data, i.e. daily forward mean" = max(diff(ts$numDate)) == mean(diff(ts$numDate)))
    stopifnot("select dt of 1(inst), 2(forward mean) or 3(trailing mean)" = max(diff(ts$numDate)) == mean(diff(ts$numDate)))

    # accumulate first
    ts <- ts %>% mutate( accum = cumsum(value * c(0,diff(numDate))))
    #ts <- ts %>% mutate( accum = cumtrapz(numDate, value) )

    if(dt == 2)
    {
      f.accum <- approxfun(ts$numDate + 1 * mean(diff(ts$numDate)), ts$accum)
    }else{
      f.accum <- approxfun(ts$numDate - 0 * mean(diff(ts$numDate)), ts$accum)
    }

    # new dataframe for output timestep
    newts <- data.frame(Date = newintTS,
                        accum = f.accum(newintTS))

  }

  # add half a timestep for instantaneous
  if(option == "inst"){
    #f.spline <- splinefun(newts$Date + 0.5 * mean(diff(newts$Date)), newts$accum  )
    f.spline <- splinefun(newts$Date, newts$accum )
    #newts$Inst <- f.spline(newintTS, deriv = 1)
    newts <- newts %>% mutate(Inst = f.spline(Date, deriv = 1))
    newts <- newts %>% dplyr::select(c(Date, Inst))

  }else if(option == "fmean")
  {
    newts$FMean <- c(diff(newts$accum)/Interval,NA )
    newts <- newts %>% dplyr::select(c(Date, FMean))
  }else if(option == "tmean")
  {
    newts$TMean <- c(NA,diff(newts$accum)/Interval )
    newts <- newts %>% dplyr::select(c(Date, TMean))
  }
  newts$Date <- as.POSIXct(newts$Date)

  if (length(inputts) == 3) {
    newts <- maxminfun(inputts[, 1], inputts[, 3], newts, option = "max")
  }
  newts[,2] <- round(newts[,2], 3)
  newts <- na.omit(newts)

  return(newts)

}

