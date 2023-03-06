#' Combined Cumulative Interpolation Filter
#'
#' Derivative of overlapping trapezoidally integrated datasets
#'
#' Useful for detiding or denoising of data
#' The overlapping averaging windows allow multiple recalculations of the same data
#' Which provides an estimate of the uncertainty in applying the filter
#'
#' @param ts dataframe of posixct time (time in seconds) and an instantaneous value (per second i.e. cumecs), and optionally a quality code
#' @param hours a positive integer of the number of hours to average over
#' @param discardbelowzero Where a number below zero is illogical, omit. Recommend leaving as default (FALSE)
#' @param centred bias the result towards the centre, not found to be useful for tides at least
#' @param type "spinterp", "spline"
#'
#' @return  dataframe with timestamp as posixct a trapezoidal integrated rate and quality code if included.
#'
#' @examples
#' data(xy_data)
#' D <- within(xy_data, {
#'   x <- as.POSIXct(x * 60 * 60 * 24, origin = "1970-01-01")
#' })
#'
#' plot(D)
#'
#' stopifnot(require(matrixStats))
#'
#' # as a tidal filter
#' spinterpData <- ccInterpFilter(D, 24, centred = FALSE)
#' lines(spinterpData$Date, spinterpData$avg, col = "red")
#'
#' # 95% confidence
#' stdev <- rowSds(as.matrix(spinterpData[2:(length(spinterpData) - 1)]))
#' spinterpConfi <- data.frame(spinterpData$Date, spinterpData$avg)
#' spinterpConfi <- cbind(spinterpConfi, upper = spinterpData$avg + 2*stdev)
#' spinterpConfi <- cbind(spinterpConfi, lower = spinterpData$avg - 2*stdev)
#'
#' lines(spinterpConfi$spinterpData.Date, spinterpConfi$upper)
#' lines(spinterpConfi$spinterpData.Date, spinterpConfi$lower)
#'
#' # as a noise removing filter
#' # keep tide but smooth
#' spinterpData <- ccInterpFilter(D, 2, centred = FALSE)
#' plot(D[50:100, ])
#' lines(spinterpData$Date, spinterpData$avg, col = "red")
#'
#'
#' # 95% confidence
#' stdev <- rowSds(as.matrix(spinterpData[2:(length(spinterpData) - 1)]))
#' spinterpConfi <- data.frame(spinterpData$Date, spinterpData$avg)
#' spinterpConfi <- cbind(spinterpConfi, upper = spinterpData$avg + 2*stdev)
#' spinterpConfi <- cbind(spinterpConfi, lower = spinterpData$avg - 2*stdev)
#'
#' lines(spinterpConfi$spinterpData.Date, spinterpConfi$upper)
#' lines(spinterpConfi$spinterpData.Date, spinterpConfi$lower)
#'
#' @export
#'
#'
#'
#'
ccInterpFilter <- function(ts, hours = 24, discardbelowzero = FALSE,
                           centred = FALSE, type = "spinterp") {

  # Trim times up, won't use part hours
  if (lubridate::ceiling_date(ts[1, 1], unit = "hours") > ts[1, 1]) {
    roundedtime <- lubridate::ceiling_date(ts[1, 1], unit = "hours") # round up
    f.round <- approxfun(ts[, 1], ts[, 2])
    startvalue <- f.round(as.POSIXct(roundedtime)) # just interpolates the first timestamp on the hour
    if (length(ts) > 2) {
      df <- data.frame(t = roundedtime, startvalue, ts[1, 3])
    } else {
      df <- data.frame(t = roundedtime, startvalue)
    }
    colnames(df) <- colnames(ts)
    ts <- ts[ts[,1] > roundedtime,]
    ts <- rbind(df, ts)
  }
  # converts time series to daily, and cumulates, interpolates, and then gets the derivative.
  # does this for each hour, i.e. for 24 hour averaging it will do it 24 times offset an hour each time
  # the result is the average of all hours interpolated.
  spinterpData <- data.frame()
  # Loop through all daily interpolations
  for (i in 0:(hours - 1))
  {
    # for each 24 hour interpolation, increment the offset by an hour each time
    daily <- changeInterval(ts, Interval = hours * 60, offset = i * 60)
    offset <- 0
    if (discardbelowzero == FALSE) {
      # spinterp can't do negative cumulative dishcharge, so add an offset to prevent any negative flows
      offset <- min(daily$FMean)
      daily$FMean <- daily$FMean - offset
    }

    daily$FMean[daily$FMean < 0] <- 0
    daily[is.na(daily)] <- 0

    # convert to numeric
    numdate <- as.numeric(as.POSIXct(daily$Date, format = format)) / 60 / 60 / 24
    # create hourly spinterp data ( interpolation of cumulative daily discharge )
    spinterpSegment <- spinterpConvert(numdate, daily$FMean, type = type)
    # remove offset again
    spinterpSegment$Data <- spinterpSegment$Data + offset
    names(spinterpSegment)[2] <- paste("Hour",i,sep="")

    if(i > 0)
    {
      spinterpData <- spinterpData %>% dplyr::inner_join(spinterpSegment, by="Date")
    }else{
      spinterpData <- spinterpSegment
    }

  }

  # centred == TRUE = bias towards inner averages, tested but not useful for tidal filter
  if (centred == TRUE) {
    avgS <- round(hours / 4, 0)
    avgF <- hours - avgS
    spinterpData <- cbind(spinterpData, avg = rowMeans(spinterpData[avgS + 2:avgF + 1]))
  } else {
    spinterpData <- cbind(spinterpData, avg = rowMeans(spinterpData[-1]))
  }
  spinterpData$Date <- as.POSIXct(spinterpData$Date * 60 * 60 * 24, origin = "1970-01-01")
  return(spinterpData)
}
