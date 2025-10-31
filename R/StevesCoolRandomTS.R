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
StevesCoolRandomTS <- function(maxFlow=500*runif(1), maxNoise=500*runif(1), obs=1000, smoothed=TRUE, randomtimes = FALSE, tideInteractions = TRUE, d1d2ratio = NULL, eventmagnitude = NULL)
{

  end <- Sys.Date()
  start <- end - obs/24

  # hourly
  hourlyseq <- seq(as.POSIXct(start), as.POSIXct(end), by = 60*60)

  posnumbers <- 0
  while(posnumbers < (obs/4)){
    # build a data frame of times and values
    df <- data.frame(Time = hourlyseq,
                     Signal = cumsum(rnorm(length(hourlyseq), mean = 0, sd = 1)) *
                       cumsum(rnorm(length(hourlyseq), mean = 0, sd = 1)) #*
                     #cumsum(rnorm(length(hourlyseq), mean = 0, sd = 1))^2
    )

    # zero any negative numbers
    choppedends <- df$Signal
    startend <- as.integer(obs*1/5)
    finishend <- as.integer(obs*4/5)

    choppedends[1:startend] <- 0
    choppedends[finishend:obs] <- 0

    posnumbers <- sum( choppedends > 0 )
  }

  df$Signal[df$Signal < 0] <- 0

  #df$Signal <- ( df$Signal / max(df$Signal) )

  df <- df %>% mutate(Signal2 = Signal + lead(Signal)) %>%
    mutate(Signal = Signal2 * lead(Signal,2)) %>%
    dplyr::select(-Signal2)

  #df$Signal <- df$Signal * c(0, df$Signal[-nrow(df)])  # autocorrelate
  #df$Signal <- round(df$Signal, digits = 2)
  #df$Signal <- df$Signal * c(0, df$Signal[-nrow(df)])



  if(smoothed)
  {
    #hourly <- changeInterval(data.frame(df$Time, df$Signal), Interval = "Hourly", option="inst")
    #godin <- godinFilter(hourly$Inst)

    # 24 hr smoothing filter
    n_hours <- 24
    filter_coeffs <- rep(1 / n_hours, n_hours)
    df$Signal <- stats::filter(df$Signal, filter = filter_coeffs, sides = 2)
    df <- na.omit(df)
  }
  if(randomtimes){

    teminseq <- seq(as.POSIXct(start), as.POSIXct(end), by = 10*60)
    randomtimestamps <- data.frame(times = teminseq,
                                   offset = rnorm(length(teminseq), mean = 5, sd = 5)) %>%
      dplyr::filter(offset > 1 & offset < 10) %>%
      mutate(times = times - 60*offset) %>% arrange(times) %>% pull(times)

    df <- data.frame(Time = randomtimestamps, Signal = approx(df$Time, df$Signal, randomtimestamps)$y)
    df <- na.omit(df)
  }

  t <- as.numeric( df[,1] ) /60/60/24
  # Sun driven cycle  #### S2
  d1 <- 2 * sin( 2*pi *  t + runif(1,0,1) * 365 )
  # Moon driven cycle ##### M2
  d2 <- 5 * sin( 24/(12+(25/60)) * 2*pi *  t + runif(1,0,1) * 28  ) # moon revolves every 12hours 25 minutes.

  # sun/moon interactions
  randomannualtime <- runif(1,0,1) * 365
  annualsun <- sin( ( 1/365.25 ) * 2 * pi *  t + randomannualtime ) + 1 # >0
  randomannualtime <- runif(1,0,1) * 365
  annualmoon <- 2 * sin( ( 1/364 ) * 2 * pi *  t + randomannualtime  )
  randomannualtime <- runif(1,0,1) * 365
  monthly <- sin( ( 1/29.53 ) * 2*pi * t + randomannualtime )  # monthly cycle is 29.5

  #if(is.null(d1d2ratio) | length(d1d2ratio) != 2)
  #{
  bias <- rnorm(2,1,0.16)
  #}else{
  #  bias <- d1d2ratio
  #}
  ratio <- bias[1]/bias[2]
  print(ratio)

  df$Noise <- bias[1]*d1*annualsun + bias[2]*d2*annualmoon*monthly #+
  # 0.1 * log((1/(bias[1]*d1*annualsun) * 1/(bias[2]*d2*annualmoon*monthly)) )# Add noise to tide

  # detrend noise with godin filtered data
  hourseq <- seq(df$Time[1], df$Time[nrow(df)], by = 60*60 )

  n_hours <- 72
  filter_coeffs <- rep(1 / n_hours, n_hours)
  trend <- approx(df$Time, df$Noise, hourseq)$y %>% stats::filter(filter = filter_coeffs, sides = 2)

  df <- df %>% mutate(Noise = Noise - approx(hourseq, trend, Time, rule = 2)$y)
  #trend <- approx(df$Time, df$Noise, hourseq)$y %>% stats::filter(1/24, sides = 2)
  #df <- df %>% mutate(Noise = Noise - trend)
  plot(df$Signal)

  df <- df %>% mutate(Noise = round( Noise * ( maxNoise / max(Noise)), 3 )) %>%
    mutate( Signal = round( Signal * ( maxFlow / max(Signal)) , 3 ))


  if(tideInteractions)
  {

    # dampen tide during events
    df <- df %>% mutate(signalratio = Signal / ( max(df$Signal) - min(df$Signal) ) ) #%>%
    #mutate(signalratio = signalratio - min(signalratio))

    range(df$signalratio)


    100 / (1000-500)

    if(is.null(eventmagnitude))
    {
      eventmagnitude <- rnorm(1, mean = 0.5, sd = 1/6)
      eventmagnitude <- min(max(eventmagnitude, 0), 1)
    }

    #df <- df %>% mutate(newnoise = Noise * (1-signalratio) ) %>% # it's the completement of the signal ratio
    #  mutate(TideInteraction = newnoise - Noise) %>% # it's a negative influence, so subtract the original noise %>%
    #  dplyr::select(-newnoise)

    df <- df %>% mutate(signalratio = signalratio * eventmagnitude)

    # add lag
    maxlag = 180 # 180 minutes

    laggedts <- df %>% mutate(lagminutes = signalratio ^ 0.5  * maxlag) %>%
      mutate(Time = Time + lagminutes * 60) %>%
      mutate(TideInteraction = Noise * ( 1 - signalratio ))

    # add the lagged interacted tide to original data
    # And subtract the original noise to get the mumeric offset
    df$TideInteraction <- approx(laggedts$Time, laggedts$TideInteraction, df$Time)$y - df$Noise

    df <- df %>% dplyr::select(-signalratio)

    df$eventmagnitude <- eventmagnitude


  }


  noise <- rnorm(length(t), mean = 0, sd = maxNoise/100)  # Standard Gaussian noise
  df$Noise <- df$Noise + noise
  df$d1d2ratio <- ratio

  # add some white noise
  #df$Signal <- df$Signal + rnorm(length(t), mean = 0.01, sd = maxFlow/10000)  # Standard Gaussian noise
  df$Signal[df$Signal < 0] <- 0

  plot(df$Signal)

  return(df)

}





