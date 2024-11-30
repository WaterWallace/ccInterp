#' Combined Cumulative Interpolation Filter
#'
#' Derivative of overlapping trapezoidally integrated datasets
#'
#' Useful for detiding or denoising of data
#' The overlapping averaging windows allow multiple recalculations of the same data
#' Which provides an estimate of the uncertainty in applying the filter
#'
#' @param data dataframe of posixct time (time in seconds) and an instantaneous value (per second i.e. cumecs)
#'
#' @return  dataframe with timestamp as posixct a trapezoidal integrated rate and quality code if included.
#'
#' @examples
#'
#' # set.seed(999)
#' # Generate some random data
#' df <- StevesCoolRandomTS()
#'
#' # must be hourly intervals
#' bwf <- butterworthFilter(data.frame(df$Time, df$Signal+df$Noise))
#' plot(Signal+Noise ~ Time, data = df, col="grey")
#' lines(Signal ~ Time, data = df, col="red")
#' lines(Filtered ~ Date, data = bwf, col = "blue")
#'
#' @export

butterworthFilter <- function(data)
{

  Rp <- 1#; % dB
  Rs <- 10#; % dB

  t <- data[,1]
  tr <- range(diff(t))
  if(tr[1] != tr [2]){
    message("Warning: non interval timeseries used")
    message(paste("Range:", tr[1], "to", tr[2], "Median" ,median(diff(t))))
  }
  #% From Rulh and Simpson 2005
  #% Use butterworth filter with
  #% a 30-hour stop period and
  #% a 40-hour pass period

  dt =  as.numeric( median((diff(t))) , units="secs" ) #;% sampling timestep (seconds)
  Fs = 1/(dt)#; % Sampling Frequency (Hz = samples/sec)
  Fn = Fs/2#; % Nyquist Frequency
  Wp = (1/(30*60*60))/Fn#; % Filter Passband (Normalised)
  Ws = (1/(40*60*60))/Fn#; % Filter Stopband (Normalised)

  btord <- signal::buttord(Wp,Ws,Rp,Rs)
  bw <- signal::butter(btord)

  BW <- signal::filtfilt(bw, data[,2])

  return(data.frame( Date = data[,1], Filtered = (data[,2] - BW)))

}
