
if(FALSE)
{
  library(ccInterp)
  library(dygraphs)
  library(xts)
  ################################
  # start/end of data
  #
  randomts <- StevesCoolRandomTS(maxFlow = 1000, obs = 300, maxNoise = 200)

  randomts <- randomts %>% mutate(TidedSignal = Signal + Noise)
  plot(randomts$Time, randomts$TidedSignal)

  bwf <- butterworthFilter(data.frame(randomts$Time,  randomts$TidedSignal))
  randomts$bwf <- bwf$Filtered
  cci <- ccInterpFilter(data.frame(randomts$Time,  randomts$TidedSignal), type = "spline")
  randomts$cci <- approx(cci$Date, cci$avg, randomts$Time)$y
  randomts$godin <- godinFilter(randomts$TidedSignal)


  plot(TidedSignal ~ Time, data = randomts, type = "l", ylab = "Discharge", col = "grey")
  points(cci ~ Time, data = randomts, col = "blue", pch = 4, cex = 1)
  points(bwf ~ Time, data = randomts, col = "orange", pch = 3, cex = 1)
  points(godin ~ Time, data = randomts, col = "darkgreen", cex = 1)
  lines(Signal ~ Time, data = randomts, col = "grey15", cex = 1)

  legend("topleft", legend = c("Unfiltered", "Input Signal", "cci","bwf","godin"),
         col = c("grey", "black", "blue", "orange", "darkgreen"), lty = c(1,1,NA,NA,NA), pch = c(NA, NA, 4,3,1))


  #####################################
  # Response functions
  #
  randomts <- StevesCoolRandomTS(maxFlow = 1000, obs = 10000, maxNoise = 200)

  xts(randomts$Signal, randomts$Time) %>% dygraph

  plot(randomts$Time, randomts$Signal)


  randomts


  par(mfrow=c(1,2))
  plot(randomts$Time, randomts$Noise)
  plot(randomts$Time, randomts$Signal)

  randomts$Time

  # Plot the input signal
  plot(randomts$Time, randomts$Signal, type = "l", col = "blue", main = "Input Signal",
       xlab = "Time (hours)", ylab = "Amplitude")

  randomts <- randomts %>% mutate(TidedSignal = Signal + Noise)
  plot(randomts$Time, randomts$TidedSignal)

  #######################
  # apply tidal filter

  # Define a simple moving average filter
  filter_length <- 25  # Length of the moving average
  filtered_signal <- stats::filter(randomts$TidedSignal, rep(1 / filter_length, filter_length), sides = 2)
  randomts$ma25Filter <- filtered_signal

  bwf <- butterworthFilter(data.frame(randomts$Time,  randomts$TidedSignal))
  randomts$bwf <- bwf$Filtered

  cci <- ccInterpFilter(data.frame(randomts$Time,  randomts$TidedSignal), type = "spline")
  randomts$cci <- approx(cci$Date, cci$avg, randomts$Time)$y

  # apply godin filter
  godin <- godinFilter(randomts$TidedSignal)
  randomts$Godin <- godin

  par(mfrow = c(2,1))

  # Plot the input signal
  plot(Signal ~ Time, data = randomts, type = "l", col = "grey10", main = "Input Signal",
       xlab = "Time (hours)", ylab = "Discharge")

  # Plot the filtered signal
  lines(bwf ~ Time, data = randomts, col = "orange", lwd = 2)
  lines(cci ~ Time, data = randomts, col = "blue", lwd = 2)
  lines(Godin ~ Time, data = randomts, col = "darkgreen", lwd = 2)
  legend("topleft", legend = c("Input Signal", "cci","bwf","godin"),
         col = c("black", "blue", "orange", "darkgreen"), lty = 1)

  maxloc <- which.max(randomts$TidedSignal)
  ts_subset <- randomts[(maxloc-100):(maxloc+100),]

  # Plot the input signal
  plot(Signal ~ Time, data = ts_subset, type = "l", col = "grey10", main = "Input Signal",
       xlab = "Time (hours)", ylab = "Discharge")

  # Plot the filtered signal
  lines(bwf ~ Time, data = ts_subset, col = "orange", lwd = 2)
  lines(cci ~ Time, data = ts_subset, col = "blue", lwd = 2)
  lines(Godin ~ Time, data = ts_subset, col = "darkgreen", lwd = 2)
  legend("topleft", legend = c("Input Signal", "cci","bwf","godin"),
         col = c("black", "blue", "orange", "darkgreen"), lty = 1)


  # remove NA's either side
  randomts <- na.omit(randomts)

  # Spectrum of input signal
  input_spectrum <- spectrum(randomts$TidedSignal , plot = FALSE)
  filtered_spectrum <- spectrum(randomts$bwf, plot = FALSE)

  dt <- as.numeric(median(diff(randomts$Time)))

  # Extract frequencies and amplitudes
  input_freq <- input_spectrum$freq / dt       # Adjust frequencies for sampling rate
  input_amp <- sqrt(input_spectrum$spec)       # Convert power to amplitude
  input_period <- 1/input_freq                 # convert frequency to period

  spectrum(randomts$cci)
  spectrum(randomts$bwf)
  spectrum(randomts$TidedSignal)

  filtered_spectrum_cci <- spectrum(randomts$cci, plot = FALSE) # spectral analysis into frequency and spectral densities
  filtered_freq_cci <- filtered_spectrum_cci$freq / dt          # divide frequency by time interval
  filtered_amp_cci <- sqrt(filtered_spectrum_cci$spec)          # calculate amplitude from spectral density
  filtered_period_cci <- 1/filtered_freq_cci                    # convert from frequency to period

  filtered_spectrum_bwf <- spectrum(randomts$bwf, plot = FALSE)
  filtered_freq_bwf <- filtered_spectrum_bwf$freq / dt
  filtered_amp_bwf <- sqrt(filtered_spectrum_bwf$spec)
  filtered_period_bwf <- 1/filtered_freq_bwf

  filtered_spectrum_gdn <- spectrum(randomts$Godin, plot = FALSE)
  filtered_freq_gdn <- filtered_spectrum_gdn$freq / dt
  filtered_amp_gdn <- sqrt(filtered_spectrum_gdn$spec)
  filtered_period_gdn <- 1/filtered_freq_gdn


  par(mfrow = c(2,2))
  plot(input_period, input_amp, log = "xy", main = "input amplitude")
  plot(filtered_period_bwf, filtered_amp_bwf, log = "xy", main = "filtered amplitude (butterworth)")
  plot(filtered_period_cci, filtered_amp_bwf, log = "xy", main = "filtered amplitude (cci)")
  plot(filtered_period_gdn, filtered_amp_gdn, log = "xy", main = "filtered amplitude (godin)")

  # Interpolate amplitudes to match frequencies and tehn
  # Divide the filtered amplitude by the input amplitude
  response_cci <- abs(approx(x = filtered_period_cci, y = filtered_amp_cci, xout = input_period)$y) / abs(input_amp)
  response_bwf <- abs(approx(x = filtered_period_bwf, y = filtered_amp_bwf, xout = input_period)$y) / abs(input_amp)
  response_gdn <- abs(approx(x = filtered_period_gdn, y = filtered_amp_gdn, xout = input_period)$y) / abs(input_amp)

  # save results for synthetic (for later comparison to real discharge)
  input_period_synthetic <- input_period
  response_bwf_synthetic <- response_bwf
  response_cci_synthetic <- response_cci

  par(mfrow = c(1,2))

  # Plot the response function
  plot(input_period, response_bwf, type = "l", col = "orange", lwd = 1,
       main = "Response Function (log y axis)", xlab = "Period (hours)",
       ylab = "Response (Amplitude Ratio)",
       xlim = c(0,60), log = "y")
  abline(h = 1, col = "gray", lty = 2)
  # Plot the response function
  lines(input_period, response_cci, type = "l", col = "blue", lwd = 1)
  abline(h = 1, col = "gray", lty = 2)
  lines(input_period, response_gdn, type = "l", col = "darkgreen", lwd = 1)
  abline(h = 1, col = "gray", lty = 2)
  legend("bottomright", legend = c("Input Signal", "cci","bwf","godin"),
         col = c("black", "blue", "orange", "darkgreen"), lty = 1)


  head(input_period)
  tail(input_period)

  # Plot the response function
  plot(input_period, response_bwf, type = "l", col = "orange", lwd = 1,
       main = "Response Function", xlab = "Period (hours)",
       ylab = "Response (Amplitude Ratio)",
       xlim = c(15,1000), ylim = c(0,1), log = "x")
  abline(h = 1, col = "gray", lty = 2)
  # Plot the response function
  lines(input_period, response_cci, type = "l", col = "blue", lwd = 1)
  abline(h = 1, col = "gray", lty = 2)
  lines(input_period, response_gdn, type = "l", col = "darkgreen", lwd = 1)
  abline(h = 1, col = "gray", lty = 2)
  legend("bottomright", legend = c("Input Signal", "cci","bwf","godin"),
         col = c("black", "blue", "orange", "darkgreen"), lty = 1)

  ####################################


  plot(randomts$Signal, randomts$bwf, col = "orange", log = "xy",
       xlab = "Synthetic Input",
       ylab = "Filtered",
       main = "Comparing to synthetic input (log)"
  )
  points(randomts$Signal, randomts$Godin, col = "darkgreen")
  points(randomts$Signal, randomts$cci, col = "blue")
  legend("bottomright", legend = c("1:1", "cci","bwf","godin"),
         col = c("black", "blue", "orange", "darkgreen"), lty = 1)
  abline(0, 1)


  plot(randomts$Signal, randomts$bwf, col = "orange",
       xlab = "Synthetic Input",
       ylab = "Filtered",
       main = "Comparing to synthetic input"
  )
  points(randomts$Signal, randomts$Godin, col = "darkgreen")
  points(randomts$Signal, randomts$cci, col = "blue")
  legend("bottomright", legend = c("1:1", "cci","bwf","godin"),
         col = c("black", "blue", "orange", "darkgreen"), lty = 1)

  abline(0, 1)

  library(data.table)
  library(ggplot2)
  melt(cci[3000:5000,], id.vars = "Date") %>% ggplot(aes(x = Date, y = value, colour = variable)) +
    geom_line()

  maxloc <- which.max(cci$avg)
  cci_subset <- cci[(maxloc-50):(maxloc+50),]

  melt(cci_subset, id.vars = "Date") %>% ggplot(aes(x = Date, y = value, colour = variable)) +
    geom_line()

  library(matrixStats)
  stdev <- rowSds(as.matrix(cci[2:(length(cci) - 1)]))

  #xts(randomts$Signal, randomts$Time) %>% dygraph

  cciConfi <- data.frame(Date = cci$Date, avg = cci$avg,
                         raw = approx(randomts$Time, randomts$TidedSignal, cci$Date)$y,
                         synthetic = approx(randomts$Time, randomts$Signal, cci$Date)$y
  )
  cciConfi <- cbind(cciConfi, upper_95 = cciConfi$avg + 2*stdev)
  cciConfi <- cbind(cciConfi, lower_95 = cciConfi$avg - 2*stdev)


  melt(cciConfi[3000:5000,], id.vars = "Date") %>% ggplot(aes(x = Date, y = value, colour = variable)) +
    geom_line()

  maxloc <- which.max(cciConfi$avg)
  cci_subset <- cciConfi[(maxloc-50):(maxloc+50),]

  melt(cci_subset, id.vars = "Date") %>% ggplot(aes(x = Date, y = value, colour = variable)) +
    geom_line()


  melt(cciConfi[3000:4000,], id.vars = "Date") %>% ggplot(aes(x = Date, y = value, colour = variable)) +
    geom_line()


  ######################
  ################################
  # real data

  mrdQ <- readRDS("data/MulgraveQ.rds")
  range(mrdQ$time)

  #randomts <- StevesCoolRandomTS()
  #randomts

  par(mfrow=c(1,2))
  plot(mrdQ$time, mrdQ$value_Discharge)


  #######################
  # apply tidal filter

  randomts <- changeInterval(data.frame(mrdQ$time, mrdQ$value_Discharge), Interval = "Hourly", option = "inst")
  names(randomts) <- c("Time", "TidedSignal")


  # Define a simple moving average filter
  filter_length <- 25  # Length of the moving average

  filtered_signal <- stats::filter(randomts$TidedSignal, rep(1 / filter_length, filter_length), sides = 2)
  randomts$ma25Filter <- filtered_signal

  bwf <- butterworthFilter(data.frame(randomts$Time,  randomts$TidedSignal))
  randomts$bwf <- bwf$Filtered


  cci <- ccInterpFilter(data.frame(mrdQ$time,  mrdQ$value_Discharge), type = "spline")
  head(cci)

  #lines(cci$Date, cci$avg, col = "blue")
  #plot(cci$Date, cci$avg)
  randomts$cci <- approx(cci$Date, cci$avg, randomts$Time)$y

  #View(randomts)

  # Plot the input signal
  #plot(randomts$Time, randomts$Signal, type = "l", col = "blue", main = "Input Signal",
  #     xlab = "Time (hours)", ylab = "Amplitude")


  # Plot the filtered signal
  lines(randomts$Time, randomts$ma25Filter, col = "red", lwd = 2)
  lines(randomts$Time, randomts$bwf, col = "orange", lwd = 2)
  lines(cci$Date, cci$avg, col = "blue")


  godin <- godinFilter(randomts$TidedSignal)
  randomts$Godin <- godin


  randomts <- na.omit(randomts)


  # Spectrum of input signal
  input_spectrum <- spectrum(randomts$TidedSignal , plot = FALSE)
  filtered_spectrum <- spectrum(randomts$bwf, plot = FALSE)

  # Extract frequencies and amplitudes
  input_freq <- input_spectrum$freq / dt       # Adjust frequencies for sampling rate
  input_amp <- sqrt(input_spectrum$spec)       # Convert power to amplitude
  input_period <- 1/input_freq

  filtered_freq <- filtered_spectrum$freq / dt
  filtered_amp <- sqrt(filtered_spectrum$spec)
  filtered_period <- 1/filtered_freq


  #Interpolate amplitudes to match frequencies
  response <- approx(x = filtered_period, y = filtered_amp, xout = input_period)$y / input_amp

  # Plot the response function
  plot(input_period, response, type = "l", col = "darkgreen", lwd = 2,
       main = "Response Function", xlab = "Frequency (cycles per hour)",
       ylab = "Response (Amplitude Ratio)",
       log = "x", ylim = c(0,1))
  abline(h = 1, col = "gray", lty = 2)  # Ideal response line

  ###################################

  # Spectrum of input signal
  #input_spectrum <- spectrum(randomts$TidedSignal , plot = FALSE)

  # Extract frequencies and amplitudes
  input_freq <- input_spectrum$freq / dt       # Adjust frequencies for sampling rate
  input_amp <- sqrt(input_spectrum$spec)       # Convert power to amplitude
  input_period <- 1/input_freq

  filtered_spectrum_cci <- spectrum(randomts$cci, plot = FALSE)
  filtered_freq_cci <- filtered_spectrum_cci$freq / dt
  filtered_amp_cci <- sqrt(filtered_spectrum_cci$spec)
  filtered_period_cci <- 1/filtered_freq_cci

  filtered_spectrum_bwf <- spectrum(randomts$bwf, plot = FALSE)
  filtered_freq_bwf <- filtered_spectrum_bwf$freq / dt
  filtered_amp_bwf <- sqrt(filtered_spectrum_bwf$spec)
  filtered_period_bwf <- 1/filtered_freq_bwf

  filtered_spectrum_gdn <- spectrum(randomts$Godin, plot = FALSE)
  filtered_freq_gdn <- filtered_spectrum_gdn$freq / dt
  filtered_amp_gdn <- sqrt(filtered_spectrum_gdn$spec)
  filtered_period_gdn <- 1/filtered_freq_gdn


  #Interpolate amplitudes to match frequencies
  response_cci <- approx(x = filtered_period_cci, y = filtered_amp_cci, xout = input_period)$y / input_amp
  response_bwf <- approx(x = filtered_period_bwf, y = filtered_amp_bwf, xout = input_period)$y / input_amp
  response_gdn <- approx(x = filtered_period_gdn, y = filtered_amp_gdn, xout = input_period)$y / input_amp


  input_period_real <- input_period
  response_bwf_real <- response_bwf
  response_cci_real <- response_cci


  # Plot the response function
  plot(input_period, response_bwf, type = "l", col = "orange", lwd = 1,
       main = "Response Function (Real Discharge)", xlab = "Period (hours)",
       ylab = "Response (Amplitude Ratio)",
       xlim = c(0,48), log = "y")

  abline(h = 1, col = "gray", lty = 2)  # Ideal response line
  # Plot the response function
  lines(input_period, response_cci, type = "l", col = "blue", lwd = 1,
        main = "Response Function", xlab = "Period (hours)",
        ylab = "Response (Amplitude Ratio)")
  abline(h = 1, col = "gray", lty = 2)  # Ideal response line
  lines(input_period, response_gdn, type = "l", col = "darkgreen", lwd = 1,
        main = "Response Function", xlab = "Period (hours)",
        ylab = "Response (Amplitude Ratio)")
  abline(h = 1, col = "gray", lty = 2)  # Ideal response line


  # Plot the response function
  plot(input_period, response_bwf, type = "l", col = "orange", lwd = 1,
       main = "Response Function (Real Discharge)", xlab = "Period (hours)",
       ylab = "Response (Amplitude Ratio)",
       xlim = c(15,1000), ylim = c(0,1), log = "x")
  abline(h = 1, col = "gray", lty = 2)  # Ideal response line
  # Plot the response function
  lines(input_period, response_cci, type = "l", col = "blue", lwd = 1,
        main = "Response Function", xlab = "Period (hours)",
        ylab = "Response (Amplitude Ratio)")
  abline(h = 1, col = "gray", lty = 2)  # Ideal response line
  lines(input_period, response_gdn, type = "l", col = "darkgreen", lwd = 1,
        main = "Response Function", xlab = "Period (hours)",
        ylab = "Response (Amplitude Ratio)")
  abline(h = 1, col = "gray", lty = 2)  # Ideal response line

  ####################################

  library(xts)
  library(dygraphs)

  cbind(
    xts(randomts$TidedSignal, randomts$Time),
    xts(randomts$cci, randomts$Time),
    xts(randomts$bwf, randomts$Time),
    xts(randomts$Godin, randomts$Time)
    #xts(randomts$ma25Filter, randomts$Time)
  ) %>% dygraph



  # Plot the response function
  plot(input_period_synthetic, response_bwf_synthetic, type = "l", col = "orange", lwd = 2,
       main = "Butterworth", xlab = "Period (hours)",
       ylab = "Response (Amplitude Ratio)",
       xlim = c(0,48), log = "y")



  # Plot the response function
  points(input_period_real, response_bwf_real, type = "l", col = "blue", lwd = 2,
         main = "Butterworth", xlab = "Period (hours)",
         ylab = "Response (Amplitude Ratio)")


  # Plot the response function
  plot(input_period_synthetic, response_cci_synthetic, type = "l", col = "orange", lwd = 2,
       main = "CCI", xlab = "Period (hours)",
       ylab = "Response (Amplitude Ratio)",
       xlim = c(0,48), log = "y")

  # Plot the response function
  points(input_period_real, response_cci_real, type = "l", col = "blue", lwd = 2,
         main = "CCI", xlab = "Period (hours)",
         ylab = "Response (Amplitude Ratio)")


  ############################################
  #
  # Plot the response function
  plot(input_period_synthetic, response_bwf_synthetic, type = "l", col = "orange", lwd = 2,
       main = "Butterworth", xlab = "Period (hours)",
       ylab = "Response (Amplitude Ratio)",
       xlim = c(15,1000), ylim = c(0,1), log = "x")
  abline(h = 1, col = "gray", lty = 2)  # Ideal response line

  # Plot the response function
  points(input_period_real, response_bwf_real, type = "l", col = "blue", lwd = 2,
         main = "Butterworth", xlab = "Period (hours)",
         ylab = "Response (Amplitude Ratio)")


  # Plot the response function
  plot(input_period_synthetic, response_cci_synthetic, type = "l", col = "orange", lwd = 2,
       main = "CCI", xlab = "Period (hours)",
       ylab = "Response (Amplitude Ratio)",
       xlim = c(15,1000), ylim = c(0,1), log = "x")
  abline(h = 1, col = "gray", lty = 2)  # Ideal response line

  # Plot the response function
  points(input_period_real, response_cci_real, type = "l", col = "blue", lwd = 2,
         main = "CCI", xlab = "Period (hours)",
         ylab = "Response (Amplitude Ratio)")



}
