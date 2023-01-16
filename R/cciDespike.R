#' cci Despkie
#'
#' Cuts timeseries data outside 2 standard deviations from the 3 hour ccinterp filtered data
#'
#' Useful for detiding or denoising of data
#' The overlapping averaging windows allow multiple recalculations of the same data
#' Which provides an estimate of the uncertainty in applying the filter
#'
#' @param spiky dataframe of posixct time (time in seconds) and an instantaneous value (per second i.e. cumecs)
#'
#' @return  dataframe with timestamp as posixct, and despiked timeseries data
#'
#' @examples
#' data(xy_data)
#' D <- within(xy_data, {
#'   x <- as.POSIXct(x * 60 * 60 * 24, origin = "1970-01-01")
#' })
#'
#' # add some spikes
#' D[round(nrow(D)/2),2] <- 60
#' D[round(nrow(D)/4),2] <- 80
#'
#' plot(D, col="red", pch=16)
#'
#' # despike
#' despikedD <- cciDespike(D)
#' points(despikedD, col="black", pch=16)
#' @export
cciDespike <- function(spiky)
{

  #spiky <- mrdPlot
  spinterpData <- ccInterpFilter(spiky, 3, centred=FALSE)

  # 95% confidence
  stdev <- rowSds(as.matrix(spinterpData[2:(length(spinterpData) -1 )]  ))
  spinterpConfi <- data.frame(spinterpData$Date, spinterpData$avg)
  spinterpConfi <- cbind(spinterpConfi, upper = spinterpData$avg + 2*stdev ) # two standard deviations
  spinterpConfi <- cbind(spinterpConfi, lower = spinterpData$avg - 2*stdev ) # two standard deviations

  f.upper <- approxfun(spinterpConfi$spinterpData.Date, spinterpConfi$upper)
  f.lower <- approxfun(spinterpConfi$spinterpData.Date, spinterpConfi$lower)

  spinterpDespiked <- spiky[spiky[,2] < f.upper(spiky[,1]) & spiky[,2] > f.lower(spiky[,1]),]
  removedPoints <- spiky[spiky[,2] > f.upper(spiky[,1]) | spiky[,2] < f.lower(spiky[,1]),]

  f.spinterpDespiked <- approxfun(spinterpDespiked)
  removedPoints$resid <- pmax(
    removedPoints[,2] - f.upper(removedPoints[,1]),
    f.lower(removedPoints[,1]) - removedPoints[,2]
  )

  removedPoints <- na.omit(removedPoints)
  spikeSD <- sd(na.omit(removedPoints$resid))
  removedPoints <- removedPoints[ abs(removedPoints$resid) > spikeSD, ]

  return (spiky[!spiky[,1] %in% removedPoints[,1],])
}



