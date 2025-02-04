#' Steves Cool Random Timeseries
#'
#' Generates a timeseries of random freshwater and tidal discharge
#'
#' Useful for examples and evaluating tide filters
#'
#' @param maxFlow The upper limit of generated freshwater data
#' @param maxNoise The upper limit of generated tidal data
#' @param obs Roughly how many hours to simulate
#' @param smoothed applies godin filter to simulated flow when TRUE
#'
#' @return  dataframe with timestamp as posixct a trapezoidal integrated rate and quality code if included.
#'
#' @examples
#'
#' par(mfrow=c(3,1))
#' set.seed(999)
#' smoothPlot <- StevesCoolRandomTS(obs=1000)
#' set.seed(999)
#' rough <- StevesCoolRandomTS(obs=1000, smoothed=FALSE)
#'  plot(smoothPlot$Time, smoothPlot$Signal)
#'  lines(rough$Time, rough$Signal, col="red")
#'  plot(smoothPlot$Time, smoothPlot$Noise)
#'  plot(smoothPlot$Time, smoothPlot$Signal+smoothPlot$Noise)
#'  lines(rough$Time, rough$Signal+rough$Noise, col='red')
#'
#' @export

StevesCoolRandomTS <- function(maxFlow=500*runif(1), maxNoise=500*runif(1), obs=10000, smoothed=TRUE, randomtimes = FALSE)
{

  end <- Sys.Date() - obs/24
  start <- runif(1, min=end-2*365, max=end)

  # make some random data
  randomMinutes <- rnorm(obs, mean = 60, sd = 30) # create 1000 random intervals
  randomMinutes <- randomMinutes[randomMinutes > 10] # exclude intervals less than 10 minutes
  randomMinutes <- as.POSIXct(cumsum(randomMinutes)*60, origin=as.Date(start, origin="1970-01-01"))
  # build a data frame of times and values
  df <- data.frame(Time = randomMinutes,
                   Signal = cumsum(rnorm(length(randomMinutes), mean = 0, sd = 1)) *
                     cumsum(rnorm(length(randomMinutes), mean = 0, sd = 1)) *
                     cumsum(rnorm(length(randomMinutes), mean = 0, sd = 1))^2
  )
  df$Signal <- df$Signal * c(0, df$Signal[-nrow(df)])  # autocorrelate
  df$Signal <- round(df$Signal, digits = 3)
  df$Signal <- df$Signal * c(0, df$Signal[-nrow(df)])


  if(smoothed)
  {
    hourly <- changeInterval(data.frame(df$Time, df$Signal), Interval = "Hourly", option="inst")
    godin <- godinFilter(hourly$Inst)

    if(!randomtimes){
      df <- data.frame(Time = hourly$Date, Signal = godin)
    }else{
      df <- data.frame(Time = df$Time, Signal = approx(hourly$Date, hourly$Inst, df$Time)$y)
    }
    df <- na.omit(df)

  }

  t <- as.numeric( df[,1] ) /60/60/24
  # Moon driven cycle
  d1 <- 5 * sin( 24/(12+(25/60)) * 2*pi * t) *
    1 * sin( ( 1/28 ) * 2*pi * t) *
    10 * sin( ( 1/364 ) * 2*pi * t+180)
  # Sun driven cycle
  d2 <- 2 * sin( 2*pi * t) *
    10 * sin( ( 1/365.25 ) * 2*pi * t+180)
  df$Noise <- d1+d2

  # Define parameters
  days <- 1:365         # Days of the year
  amplitude <- 1        # Amplitude of the sine wave
  phase_shift <- 0      # Phase shift (in radians)

  # Calculate the sine wave
  annual <- sin((2 * pi / 365) * t + 0) ^ 2
  df$Signal <- df$Signal + ( annual + abs(min(annual) ))

  df$Signal <- round( df$Signal / (max(abs(df$Signal)) / maxFlow) , 3 )
  df$Noise <- round( df$Noise / (max(df$Noise)/maxNoise), 3 )

  return(df)

}




