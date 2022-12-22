#as a spike removal tool

cciDespike <- function(D)
{

  #D <- mrdPlot
  spinterpData <- ccInterpFilter(D, 3, centred=FALSE)

  # 95% confidence
  stdev <- rowSds(as.matrix(spinterpData[2:(length(spinterpData) -1 )]  ))
  spinterpConfi <- data.frame(spinterpData$Date, spinterpData$avg)
  spinterpConfi <- cbind(spinterpConfi, upper = spinterpData$avg + 2*stdev ) # two standard deviations
  spinterpConfi <- cbind(spinterpConfi, lower = spinterpData$avg - 2*stdev ) # two standard deviations

  f.upper <- approxfun(spinterpConfi$spinterpData.Date, spinterpConfi$upper)
  f.lower <- approxfun(spinterpConfi$spinterpData.Date, spinterpConfi$lower)

  spinterpDespiked <- D[D[,2] < f.upper(D[,1]) & D[,2] > f.lower(D[,1]),]
  removedPoints <- D[D[,2] > f.upper(D[,1]) | D[,2] < f.lower(D[,1]),]

  f.spinterpDespiked <- approxfun(spinterpDespiked)
  removedPoints$resid <- pmax(
    removedPoints[,2] - f.upper(removedPoints[,1]),
    f.lower(removedPoints[,1]) - removedPoints[,2]
  )

  removedPoints <- na.omit(removedPoints)
  spikeSD <- sd(na.omit(removedPoints$resid))
  removedPoints <- removedPoints[ abs(removedPoints$resid) > spikeSD, ]

  return (D[!D[,1] %in% removedPoints[,1],])
}

library(ccInterp)
library(matrixStats)

data(xy_data)
D <- within(xy_data, {
  x <- as.POSIXct(x * 60 * 60 * 24, origin = "1970-01-01")
})

# add some spikes
D[round(nrow(D)/2),2] <- 60
D[round(nrow(D)/4),2] <- 80


plot(D, col="red", pch=16)

# despike
despikedD <- cciDespike(D)
points(despikedD, col="black", pch=16)
