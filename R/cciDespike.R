#' cci Despkie
#'
#' Cuts timeseries data outside 2 standard deviations from the 3 hour ccinterp filtered data
#'
#' Worth experimenting with the number of standard deviations and averaging window for the particular dataset
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
#'
#' spiky <- StevesCoolRandomTS()
#' spiky <- spiky %>% mutate(Tided = Signal+Noise+TideInteraction) %>%
#' dplyr::select(Time, Tided)
#' dfoutliers <- sample_n(spiky, 20) %>%
#'   mutate( Tided = rnorm(20, 0, 200) )
#'
#' spiky <- spiky %>%
#' left_join(dfoutliers, by = "Time", suffix = c("", "_new")) %>%
#' mutate(Tided = dplyr::coalesce(Tided_new, Tided)) %>%
#' dplyr::select(Time, Tided)
#'
#' cciDespike(spiky, doPlot = TRUE, stdevs = 2)
#' points(dfoutliers, col = "blue", pch = 4)
#'
#' @export
cciDespike <- function(spiky, hoursAvg = 3, stdevs = 2, rocstdev = 3, doPlot = FALSE, lonelyPoints = 4, spikeaction = "omit")
{
  #savedspiky <- spiky
  # spiky <- savedspiky
  # spiky <- LRTDatasaved
  # spiky <- as.data.frame(spiky)
  spiky$tdiff <- c(diff(as.numeric(spiky[,1])),0) / 60 / 60
  #names(spiky)[1] <- "ts"
  spiky$Filter <- 0

  #nrow(spiky)

  # remove lonely points
  lonelypoints <- spiky %>% dplyr::mutate(tdifflag = lag(tdiff,1,0)) %>%
  mutate(mindiff = pmin(tdiff,tdifflag))

  numlonelypoints <- lonelypoints %>% dplyr::filter(mindiff >= lonelyPoints)
  if(nrow(numlonelypoints) > 0)
  {
    lonelypoints[lonelypoints$mindiff >= lonelyPoints,]$Filter <- 2
  }
  spiky <- lonelypoints %>% dplyr::select(-c(tdiff,tdifflag,mindiff))

  f.spline <- splinefun(spiky[,1], spiky[,2])
  spiky$dydx <- f.spline(spiky[,1], deriv = 1)
  #spiky$Filter <- 0

  if(!is.null(rocstdev))
  {
    spiky$Filter[ abs( spiky$dydx ) > sd(spiky$dydx) * rocstdev  ] <- 2
  }

  #cbind(xts(spiky$Value, spiky$ts),
  #      xts(spiky$Filter, spiky$ts)) %>% dygraph()


  #removedPoints <- spiky[spiky$Filter == 1,]
  #spiky[spiky$Filter == 0,]

  spinterpData <- ccInterpFilter(data.frame(spiky[spiky$Filter == 0,1], spiky[spiky$Filter == 0,2]), hoursAvg, centred=FALSE)

  # 95% confidence
  stdev <- rowSds(as.matrix(spinterpData[2:(length(spinterpData) -1 )]  ))
  spinterpConfi <- data.frame(spinterpData$Date, spinterpData$avg)
  spinterpConfi <- cbind(spinterpConfi, upper = spinterpData$avg + stdevs*stdev ) # two standard deviations
  spinterpConfi <- cbind(spinterpConfi, lower = spinterpData$avg - stdevs*stdev ) # two standard deviations

  f.upper <- approxfun(spinterpConfi$spinterpData.Date, spinterpConfi$upper, rule = 2)
  f.lower <- approxfun(spinterpConfi$spinterpData.Date, spinterpConfi$lower, rule = 2)

  spiky <- spiky %>% na.omit(spiky)

  #spikysaved <- spiky
  #spiky <- spikysaved
  spiky$outofrange <- (spiky[,2] > f.upper(spiky[,1]) | spiky[,2] < f.lower(spiky[,1]))
  spiky <- spiky %>% mutate(Filter = ifelse(Filter == 0 & outofrange, 1, Filter)) %>%
    dplyr::select(-outofrange)

  #cbind(xts(spiky$Value, spiky$ts),
  #      xts(spiky$Filter, spiky$ts)) %>% dygraph()
  #xts(spiky$Value, spiky$ts) %>% dygraph

  #removedPoints <- rbind( removedPoints,
  #                        spiky[(spiky[,2] > f.upper(spiky[,1]) | spiky[,2] < f.lower(spiky[,1])),])



  #nrow(removedPoints) + nrow(spinterpDespiked)

  #spinterpDespiked <- spinterpDespiked %>% select(-dydx)
  #removedPoints <- removedPoints %>% select(-dydx)
  #spiky <-

  # how far the "spike" is away from the limits
  removedPoints <- spiky %>% dplyr::filter(Filter == 1)
  removedPoints$resid <- pmax(
    removedPoints[,2] - f.upper(removedPoints[,1]),
    f.lower(removedPoints[,1]) - removedPoints[,2]
  )

  removedPoints <- na.omit(removedPoints)




  SixHourlyResiduals <- changeInterval(data.frame(removedPoints[,1], removedPoints$resid), Interval = hoursAvg*2*60, option = "tmean")
  SixHourlyResiduals$SD <- rollapply(SixHourlyResiduals$TMean,width=10,FUN=sd,fill=NA,align="c") * stdevs
  f.SD <- approxfun(SixHourlyResiduals$Date, SixHourlyResiduals$SD, na.rm=TRUE, rule=2)

  removedPoints <- removedPoints[removedPoints$resid > f.SD(removedPoints[,1]),]

  if(FALSE)
  {
    plotspiky <- spiky %>% dplyr::filter(ts > "2021-06-11 00:00" & ts < "2021-09-09 00:00")
    plotspiky <- spiky
    #lenspiky <- nrow(spiky)
    #plot(spiky[(lenspiky-200):lenspiky,1], spiky[(lenspiky-200):lenspiky,2])
    plot(plotspiky[,1], plotspiky[,2])
    #spiky

    upperLine <- data.frame( x = rev(plotspiky[,1]),
                             y = rev(f.upper(plotspiky[,1]) + f.SD(plotspiky[,1]) )
    )
    bottomLine <- data.frame(x = plotspiky[,1],
                             y = f.lower(plotspiky[,1]) - f.SD(plotspiky[,1])
    )

    cbind(
      xts(spiky$Value, spiky$ts),
      xts(removedPoints$Value, removedPoints$ts),
      xts(upperLine$y, upperLine$x),
      xts(bottomLine$y, bottomLine$x)
    ) %>% dygraph

    bigshape <- rbind(bottomLine, upperLine)

    #
    upperLine <- data.frame( x = rev(plotspiky[,1]),
                             y = rev(f.upper(plotspiky[,1]) )
    )
    bottomLine <- data.frame(x = plotspiky[,1],
                             y = f.lower(plotspiky[,1])
    )
    shape <- rbind(bottomLine, upperLine)

    shape <- na.omit(shape)
    bigshape <- na.omit(bigshape)

    polygon(bigshape, col="lightgrey" )
    polygon(shape, col="darkgrey")
    points(spinterpDespiked, col="black")
    points(removedPoints, col="red", pch=19)

  }

  spiky <- spiky %>% mutate(spike = ifelse(Filter > 1, "spike","OK")) %>%
    mutate(spike = ifelse(.[[1]] %in% removedPoints[,1], "spike",spike))


  #spiky$spike <- ifelse(spiky[,1] %in% removedPoints[,1], "spike",spike)
  spiky <- spiky %>% dplyr::select( -c(dydx, Filter) )


  if(spikeaction == "code")
  {
      return(spiky)
  }else if(spikeaction == "infill")
  {
    gooddata <- spiky %>% dplyr::filter(spike == "OK")
    f.gooddata <- approxfun(gooddata[,1],gooddata[,2])
    spiky[spiky$spike == "spike",][,2] <- f.gooddata(spiky[spiky$spike == "spike",][,1])
    return(spiky %>% dplyr::select(-spike))
  }else if(spikeaction == "omit")
  {
    return ( spiky %>% dplyr::filter(spike == "OK") %>% dplyr::select(-spike))
  }
}







