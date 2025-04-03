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
StevesCoolRandomTS <- function(maxFlow=1000*runif(1), maxNoise=500*runif(1), obs=1000, smoothed=TRUE, randomtimes = FALSE, tideInteractions = TRUE)
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
      df <- data.frame(Time = df$Time, Signal = approx(hourly$Date, godin, df$Time)$y)
    }
    df <- na.omit(df)

  }

  t <- as.numeric( df[,1] ) /60/60/24
  # Sun driven cycle
  d1 <- 2*sin( 2*pi * t)
  # Moon driven cycle
  d2 <- 5*sin( 24/(12+(25/60)) * 2*pi * t) # moon revolves every 12hours 25 minutes.

  # sun/moon interactions
  annualsun <- 5*sin( ( 1/365.25 ) * 2*pi * t+180)
  annualmoon <- 5*sin( ( 1/364 ) * 2*pi * t)
  monthly <- sin( ( 1/29.53 ) * 2*pi * t)  # monthly cycle is 29.5


  bias <- rnorm(2,1,0.16)
  ratio <- bias[1]/bias[2]
  print(ratio)

  df$Noise <- bias[1]*d1*annualsun + bias[2]*d2*annualmoon*monthly  # Add noise to tide

  # detrend noise with godin filtered data
  hourseq <- seq(df$Time[1], df$Time[nrow(df)], by = 60*60 )
  trend <- approx(df$Time, df$Noise, hourseq)$y %>% godinFilter
  df <- df %>% mutate(Noise = Noise - approx(hourseq, trend, Time, rule = 2)$y)

  noiseratio <- max(df$Noise) / maxNoise

  df$Signal <- df$Signal + ( annual + abs(min(annual) ))
  df$Signal <- round( df$Signal / (max(abs(df$Signal)) / maxFlow) , 3 )


  if(tideInteractions)
  {

    # dampen tide during events
    df <- df %>% mutate(signalratio = Signal / ( max(df$Signal) - min(df$Signal) ) ) %>%
      mutate(signalratio = signalratio - min(signalratio))

    df <- df %>% mutate(TideInteraction = ( Noise * ( 1 - signalratio ) ) - Noise ) # new noise minus noise

    # add lag
    maxlag = 120 # 120 minutes
    laggedts <- df %>% mutate(lagminutes = signalratio ^ 0.5  * maxlag) %>%
      mutate(Time = Time + lagminutes * 60)

    # lagged noise minus existing noise
    laggednoise <- approx(laggedts$Time, (laggedts$Noise+laggedts$TideInteraction), df$Time, rule = 2)$y
    laggdiff <- laggednoise - (df$Noise+df$TideInteraction)

    df$TideInteraction <- (df$TideInteraction + laggdiff) %>% round(3)

    df <- df %>% dplyr::select(-signalratio)
    df$TideInteraction <- df$TideInteraction / noiseratio

  }

  # round and normalise to desired values
  df <- df %>% mutate(Noise = round( Noise / noiseratio, 3 ))

  noise <- rnorm(length(t), mean = 0, sd = maxNoise/100)  # Standard Gaussian noise
  df$Noise <- df$Noise + noise
  df$d1d2ratio <- ratio

  # add some white noise
  df$Signal <- df$Signal + rnorm(length(t), mean = 0.01, sd = maxFlow/10000)  # Standard Gaussian noise
  df$Signal[df$Signal < 0] <- 0


  return(df)

}

df <- StevesCoolRandomTS()
head(df)
plot(df$Time, df$Noise)
xts(df$Noise+df$Signal+df$TideInteraction, df$Time) %>% dygraph
#xts(df$Signal, df$Time) %>% dygraph


library(Rlibeemd)

#imfs <- ceemdan(df$Noise)




# Sample time series (e.g., tidal or lunar cycle data)
#set.seed(123)
#t <- seq(0, 10, by = 0.1)  # Time vector
#y <- 5*sin(2 * pi * 1.93 * t) + 2*sin(2 * pi * 0.5 * t) + rnorm(length(t), sd = 0.5)  # Simulated data
y <- df$Noise


# Compute FFT
fft_result <- fft(y)
n <- length(y)
freqs <- seq(0, 1, length.out = n) * (1 / (t[2] - t[1]))  # Convert indices to frequencies

# Get amplitudes (magnitude of complex FFT values)
amplitudes <- Mod(fft_result) / n  # Normalize

# Only take the first half (positive frequencies)
half_n <- floor(n / 2)
freqs <- freqs[1:half_n]
amplitudes <- amplitudes[1:half_n]

# Plot spectrum
plot(freqs, amplitudes, type = "h", col = "blue", lwd = 2,
     xlab = "Frequency (cycles per unit time)", ylab = "Amplitude",
     main = "Frequency Spectrum of Time Series")

library(signal)  # Install if necessary: install.packages("signal")
library(pracma)
peaks <- findpeaks(amplitudes, nups = 1, ndowns = 1, npeaks = 5, sortstr = TRUE)
dominant_freqs <- freqs[peaks[,2]]
dominant_amplitudes <- peaks[,1]

# Print results
data.frame(Frequency = dominant_freqs, Amplitude = dominant_amplitudes)



