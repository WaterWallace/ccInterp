#' cci Despkie
#'
#' Cuts timeseries data outside 2 standard deviations from the 3 hour ccinterp filtered data
#'
#' Useful for detiding or denoising of data
#' The overlapping averaging windows allow multiple recalculations of the same data
#' Which provides an estimate of the uncertainty in applying the filter
#'
#' @param spiky dataframe of posixct time (time in seconds) and an instantaneous value (per second i.e. cumecs)
#' @param hoursAvg number of hours to smooth
#' @param stdevs number of standard deviations to include
#' @param doPlot do plots
#'
#' @return  dataframe with timestamp as posixct, and despiked timeseries data
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

cciDespike <- function(spiky, hoursAvg = 3, stdevs = 2, doPlot = FALSE)
{

  spinterpData <- ccInterpFilter(spiky, hoursAvg, centred=FALSE)

  # 95% confidence
  stdev <- rowSds(as.matrix(spinterpData[2:(length(spinterpData) -1 )]  ))
  spinterpConfi <- data.frame(spinterpData$Date, spinterpData$avg)
  spinterpConfi <- cbind(spinterpConfi, upper = spinterpData$avg + stdevs*stdev ) # two standard deviations
  spinterpConfi <- cbind(spinterpConfi, lower = spinterpData$avg - stdevs*stdev ) # two standard deviations

  f.upper <- approxfun(spinterpConfi$spinterpData.Date, spinterpConfi$upper)
  f.lower <- approxfun(spinterpConfi$spinterpData.Date, spinterpConfi$lower)

  spinterpDespiked <- spiky[spiky[,2] <= f.upper(spiky[,1]) & spiky[,2] >= f.lower(spiky[,1]),]
  removedPoints <- spiky[spiky[,2] > f.upper(spiky[,1]) | spiky[,2] < f.lower(spiky[,1]),]

  f.spinterpDespiked <- approxfun(spinterpDespiked)
  removedPoints$resid <- pmax(
    removedPoints[,2] - f.upper(removedPoints[,1]),
    f.lower(removedPoints[,1]) - removedPoints[,2]
  )

  removedPoints <- na.omit(removedPoints)

  SixHourlyResiduals <- changeInterval(data.frame(removedPoints[,1], removedPoints$resid), Interval = 6*60)

  SixHourlyResiduals$SD <- rollapply(SixHourlyResiduals$FMean,width=10,FUN=sd,fill=NA,align="c")

  f.SD <- approxfun(SixHourlyResiduals$Date, SixHourlyResiduals$SD, na.rm=TRUE, rule=2)
  removedPoints <- removedPoints[ abs(removedPoints$resid) > f.SD(removedPoints[,1]), ]

  if(doPlot)
  {
    plot(spiky)

    upperLine <- data.frame( x = rev(spiky[,1]),
                             y = rev(f.upper(spiky[,1])+f.SD(spiky[,1]) )
    )
    bottomLine <- data.frame(x = spiky[,1],
                             y = f.lower(spiky[,1])-f.SD(spiky[,1])
    )
    bigshape <- rbind(bottomLine, upperLine)

    #
    upperLine <- data.frame( x = rev(spiky[,1]),
                             y = rev(f.upper(spiky[,1]) )
    )
    bottomLine <- data.frame(x = spiky[,1],
                             y = f.lower(spiky[,1])
    )
    shape <- rbind(bottomLine, upperLine)

    shape <- na.omit(shape)
    bigshape <- na.omit(bigshape)

    polygon(bigshape, col="lightgrey" )
    polygon(shape, col="darkgrey")
    points(spinterpDespiked, col="black")
    points(removedPoints, col="red", pch=19)

  }

  return (spiky[!spiky[,1] %in% removedPoints[,1],])
}



