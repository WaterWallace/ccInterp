
if(FALSE)
{
  library(ccInterp)
  library(dygraphs)
  library(xts)
  library(dplyr)
  ################################
  # start/end of data
  #

  # 9912 ringing around zero
  # 3339 just a good comparable example of bwf and cci being well matched
  #randomts <- StevesCoolRandomTS(maxFlow = 1000, obs = 600, maxNoise = 200, smoothed = TRUE, randomtimes = TRUE)


  randomSeed <- sample(1:10000, 1)
  set.seed(randomSeed)

  #set.seed(959)

  # 5698 rapid rises causes ringing at the beginning of an event.
  # 5905 rapid rises causes ringing at the beginning of an event. # good one
  randomts <- StevesCoolRandomTS(maxFlow = 500, obs = 200, maxNoise = 200, smoothed = FALSE, randomtimes = TRUE)
  #seq(randomts$Time[1], randomts$Time[1])

  randomts <- randomts %>% mutate(TidedSignal = Signal + Noise)
  plot(randomts$Time, randomts$TidedSignal)

  randomtsHrly <- randomts %>% dplyr::select(Time, TidedSignal) %>% changeInterval(Interval = "Hourly", option = "inst")
  points(randomtsHrly, col = "red")

  bwf <- butterworthFilter(randomtsHrly)
  randomtsHrly$bwf <- bwf$Filtered
  #head(randomtsHrly)
  #head(bwf)

  cci <- ccInterpFilter(data.frame(randomts$Time,  randomts$TidedSignal), type = "spinterp")

  randomtsHrly$cci <- approx(cci$Date, cci$avg, randomtsHrly$Date)$y
  randomtsHrly$godin <- godinFilter(randomtsHrly$Inst)

  plot(TidedSignal ~ Time, data = randomts, type = "l", ylab = "Discharge", col = "grey")
  points(cci ~ Date, data = randomtsHrly, col = "blue", pch = 4, cex = 1)
  points(bwf ~ Date, data = randomtsHrly, col = "orange", pch = 3, cex = 1)
  points(godin ~ Date, data = randomtsHrly, col = "darkgreen", cex = 1)
  lines(Signal ~ Time, data = randomts, col = "grey15", cex = 1)

  legend("topleft", legend = c("Unfiltered", "Input Signal", "cci","bwf","godin"),
         col = c("grey", "black", "blue", "orange", "darkgreen"), lty = c(1,1,NA,NA,NA), pch = c(NA, NA, 4,3,1))

  print(randomSeed)


  #####################################
  # Response functions
  #

  randomSeed <- sample(1:10000, 1)

  randomSeed <- 959
  set.seed(randomSeed)


  randomts <- StevesCoolRandomTS(maxFlow = 1000, obs = 10000, maxNoise = 200, smoothed = TRUE, randomtimes = TRUE)

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

  randomtsHrly <- randomts %>%
    dplyr::select(c(Time, TidedSignal)) %>%
    changeInterval(option = "inst", Interval = "Hourly")

  plot(TidedSignal ~ Time, data = randomts[150:200,])
  lines(randomtsHrly)


  #######################
  # apply tidal filter

  # Define a simple moving average filter
  filter_length <- 25  # Length of the moving average
  filtered_signal <- stats::filter(randomtsHrly$Inst, rep(1 / filter_length, filter_length), sides = 2)
  randomtsHrly$ma25Filter <- filtered_signal

  bwf <- butterworthFilter(randomtsHrly)
  randomtsHrly$bwf <- bwf$Filtered

  cci <- ccInterpFilter(data.frame(randomts$Time,  randomts$TidedSignal), type = "spinterp")
  randomtsHrly$cci <- approx(cci$Date, cci$avg, randomtsHrly$Date)$y

  # apply godin filter
  godin <- godinFilter(randomtsHrly$Inst)
  randomtsHrly$Godin <- godin

  par(mfrow = c(2,1))

  # Plot the input signal
  plot(Signal ~ Time, data = randomts, type = "l", col = "grey10", main = "Input Signal",
       xlab = "Time (hours)", ylab = "Discharge")

  # Plot the filtered signal
  lines(bwf ~ Date, data = randomtsHrly, col = "orange", lwd = 2)
  lines(cci ~ Date, data = randomtsHrly, col = "blue", lwd = 2)
  lines(Godin ~ Date, data = randomtsHrly, col = "darkgreen", lwd = 2)
  legend("topleft", legend = c("Input Signal", "cci","bwf","godin"),
         col = c("black", "blue", "orange", "darkgreen"), lty = 1)

  maxloc <- which.max(randomts$TidedSignal)
  ts_subset <- randomts[(maxloc-100):(maxloc+100),]

  # Plot the input signal
  plot(Signal ~ Time, data = ts_subset, type = "l", col = "grey10", main = "Input Signal",
       xlab = "Time (hours)", ylab = "Discharge")

  # Plot the filtered signal
  lines(bwf ~ Date, data = randomtsHrly, col = "orange", lwd = 2)
  lines(cci ~ Date, data = randomtsHrly, col = "blue", lwd = 2)
  lines(Godin ~ Date, data = randomtsHrly, col = "darkgreen", lwd = 2)
  legend("topleft", legend = c("Input Signal", "cci","bwf","godin"),
         col = c("black", "blue", "orange", "darkgreen"), lty = 1)

  # remove NA's either side
  randomtsHrly <- na.omit(randomtsHrly)

  # Spectrum of input signal
  input_spectrum <- spectrum(randomtsHrly$Inst , plot = FALSE)
  filtered_spectrum <- spectrum(randomtsHrly$bwf, plot = FALSE)

  dt <- as.numeric(median(diff(randomtsHrly$Date)))

  # Extract frequencies and amplitudes
  input_freq <- input_spectrum$freq / dt       # Adjust frequencies for sampling rate
  input_amp <- sqrt(input_spectrum$spec)       # Convert power to amplitude
  input_period <- 1/input_freq                 # convert frequency to period

  spectrum(randomtsHrly$cci)
  spectrum(randomtsHrly$bwf)
  spectrum(randomtsHrly$Inst)

  filtered_spectrum_cci <- spectrum(randomtsHrly$cci, plot = FALSE) # spectral analysis into frequency and spectral densities
  filtered_freq_cci <- filtered_spectrum_cci$freq / dt          # divide frequency by time interval
  filtered_amp_cci <- sqrt(filtered_spectrum_cci$spec)          # calculate amplitude from spectral density
  filtered_period_cci <- 1/filtered_freq_cci                    # convert from frequency to period

  filtered_spectrum_bwf <- spectrum(randomtsHrly$bwf, plot = FALSE)
  filtered_freq_bwf <- filtered_spectrum_bwf$freq / dt
  filtered_amp_bwf <- sqrt(filtered_spectrum_bwf$spec)
  filtered_period_bwf <- 1/filtered_freq_bwf

  filtered_spectrum_gdn <- spectrum(randomtsHrly$Godin, plot = FALSE)
  filtered_freq_gdn <- filtered_spectrum_gdn$freq / dt
  filtered_amp_gdn <- sqrt(filtered_spectrum_gdn$spec)
  filtered_period_gdn <- 1/filtered_freq_gdn


  par(mfrow = c(2,2))
  plot(input_period, input_amp, log = "xy", main = "input amplitude")
  plot(filtered_period_bwf, filtered_amp_bwf, log = "xy", main = "filtered amplitude (butterworth)")
  plot(filtered_period_cci, filtered_amp_cci, log = "xy", main = "filtered amplitude (cci)")
  plot(filtered_period_gdn, filtered_amp_gdn, log = "xy", main = "filtered amplitude (godin)")


  Periodograms <- bind_rows(
  data.frame (  Period = input_period, Amplitude = input_amp, Dataset = "Input" ),
  data.frame ( Period = filtered_period_bwf, Amplitude = filtered_amp_bwf, Dataset = "Butterworth" ),
  data.frame ( Period = filtered_period_cci, Amplitude = filtered_amp_cci, Dataset = "CCI" ),
  data.frame ( Period = filtered_period_gdn, Amplitude = filtered_amp_gdn, Dataset = "Godin" ))

  head(Periodograms)
  Periodograms %>% ggplot(aes(x = Period, y = Amplitude)) +
    geom_point()  +
    scale_x_continuous(trans = "log2", breaks = c(3, 6, 12, 24, 48, 168, 672, 8760),
                       minor_b) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    theme_minimal() +
    #annotation_logticks() +
    facet_wrap(vars(Dataset), scales = "free_y")


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
  plot(input_period, response_bwf^2, type = "l", col = "orange", lwd = 1,
       main = "Response Function", xlab = "Period (hours)",
       ylab = "Response (Amplitude Ratio)",
       xlim = c(15,1000), ylim = c(0,1), log = "x")
  abline(h = 1, col = "gray", lty = 2)
  # Plot the response function
  lines(input_period, response_cci^2, type = "l", col = "blue", lwd = 1)
  abline(h = 1, col = "gray", lty = 2)
  lines(input_period, response_gdn^2, type = "l", col = "darkgreen", lwd = 1)
  abline(h = 1, col = "gray", lty = 2)
  legend("bottomright", legend = c("Input Signal", "cci","bwf","godin"),
         col = c("black", "blue", "orange", "darkgreen"), lty = 1)

  ####################################


  plot(randomts$Signal, approx(randomtsHrly$Date, randomtsHrly$bwf, randomts$Time)$y, col = "orange", log = "xy",
       xlab = "Synthetic Input",
       ylab = "Filtered",
       main = "Comparing to synthetic input (log)"
  )
  points(randomts$Signal, approx(randomtsHrly$Date, randomtsHrly$Godin, randomts$Time)$y, col = "darkgreen")
  points(randomts$Signal, approx(randomtsHrly$Date, randomtsHrly$cci, randomts$Time)$y, col = "blue")
  legend("bottomright", legend = c("1:1", "cci","bwf","godin"),
         col = c("black", "blue", "orange", "darkgreen"), lty = 1)
  abline(0, 1)

  plot(randomts$Signal, approx(randomtsHrly$Date, randomtsHrly$bwf, randomts$Time)$y, col = "orange",
       xlab = "Synthetic Input",
       ylab = "Filtered",
       main = "Comparing to synthetic input"
  )
  points(randomts$Signal, approx(randomtsHrly$Date, randomtsHrly$Godin, randomts$Time)$y, col = "darkgreen")
  points(randomts$Signal, approx(randomtsHrly$Date, randomtsHrly$cci, randomts$Time)$y, col = "blue")
  legend("bottomright", legend = c("1:1", "cci","bwf","godin"),
         col = c("black", "blue", "orange", "darkgreen"), lty = 1)

  abline(0, 1)

  print(randomSeed)


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

  #JRIQ <- getTSDB(client, "1120053", var = "Discharge", bucket = "tsdata3",
  #                start = as.POSIXct("2021-01-01 00:00"),
  #               stop = as.POSIXct("2024-01-01 00:00"))
  #q2024 <- RREQ[[1]] %>% dplyr::select(time, value_Discharge, QC_Discharge)
  #qPre2024 <- RREQ[[2]] %>% dplyr::select(time, value_Discharge, QC_Discharge)
  #qMerged <- mergeTS(qPre2024, q2024)
  #saveRDS(qMerged, "RREQ.rds")

  realQ <- readRDS("data/JRIQ.rds")
  #realQ <- readRDS("data/MulgraveQ.rds")

  range(realQ$time)
  range(realQ$value_Discharge)


  #randomts <- StevesCoolRandomTS()
  #randomts

  par(mfrow=c(1,2))
  plot(realQ$time, realQ$value_Discharge,
       ylab = "Johnstone Discharge (cumecs)",
       xlab = "Date",
       type = "l")


  #######################
  # apply tidal filter

  randomts <- changeInterval(data.frame(realQ$time, realQ$value_Discharge), Interval = "Hourly", option = "inst")
  names(randomts) <- c("Time", "TidedSignal")


  # Define a simple moving average filter
  filter_length <- 25  # Length of the moving average

  filtered_signal <- stats::filter(randomts$TidedSignal, rep(1 / filter_length, filter_length), sides = 2)
  randomts$ma25Filter <- filtered_signal

  bwf <- butterworthFilter(data.frame(randomts$Time,  randomts$TidedSignal))
  randomts$bwf <- bwf$Filtered


  cci <- ccInterpFilter(data.frame(realQ$time,  realQ$value_Discharge), type = "spinterp")
  head(cci)

  #lines(cci$Date, cci$avg, col = "blue")
  #plot(cci$Date, cci$avg)
  randomts$cci <- approx(cci$Date, cci$avg, randomts$Time)$y

  #View(randomts)

  # Plot the input signal
  plot(randomts$Time, randomts$TidedSignal, type = "l", col = "blue", main = "Input Signal",
       xlab = "Time (hours)", ylab = "Amplitude")


  # Plot the filtered signal
  lines(randomts$Time, randomts$ma25Filter, col = "red", lwd = 2)
  lines(randomts$Time, randomts$bwf, col = "orange", lwd = 2)
  lines(cci$Date, cci$avg, col = "blue")


  godin <- godinFilter(randomts$TidedSignal)
  randomts$Godin <- godin
  randomts <- na.omit(randomts)

  # flood event peak
  randomts %>% dplyr::filter(Time > "2023-12-15" & Time < "2023-12-20") %>%
    melt(id.vars = "Time") %>%
    ggplot(aes(x= Time, y = value, colour = variable)) +
    ggtitle("Russell Event Peak") +
    geom_line()

  # Dip before rise
  randomts %>% dplyr::filter(Time > "2024-01-09" & Time < "2024-01-16") %>%
    melt(id.vars = "Time") %>%
    ggplot(aes(x= Time, y = value, colour = variable)) +
    ggtitle("Dip Before Rise") +
    geom_line()

  # Excessive ringing during low flows
  randomts %>%
    dplyr::select(-c(TidedSignal, ma25Filter)) %>%
    dplyr::filter(Time > "2022-10-02" & Time < "2022-12-16") %>%
    melt(id.vars = "Time") %>%
    ggplot(aes(x= Time, y = value, colour = variable)) +
    ggtitle("Excessive ringing during low flows") +
    #lims(y = c(0,50)) +
    geom_line()


  # Spectrum of input signal
  input_spectrum <- spectrum(randomts$TidedSignal , plot = TRUE)

  par(mfrow = c(1,1))
  plot(1/input_spectrum$freq, input_spectrum$spec, log = "xy", type = "l")



  filtered_spectrum <- spectrum(randomts$bwf, plot = TRUE)
  spectrum(randomts$cci, plot = TRUE)

  dt <- 1
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

  filtered_spectrum_ma25 <- spectrum(randomts$ma25Filter, plot = FALSE)
  filtered_freq_ma25 <- filtered_spectrum_ma25$freq / dt
  filtered_amp_ma25 <- sqrt(filtered_spectrum_ma25$spec)
  filtered_period_ma25 <- 1/filtered_freq_ma25


  #Interpolate amplitudes to match frequencies
  response_cci <- approx(x = filtered_period_cci, y = filtered_amp_cci, xout = input_period)$y / input_amp
  response_bwf <- approx(x = filtered_period_bwf, y = filtered_amp_bwf, xout = input_period)$y / input_amp
  response_gdn <- approx(x = filtered_period_gdn, y = filtered_amp_gdn, xout = input_period)$y / input_amp
  response_ma25 <- approx(x = filtered_period_ma25, y = filtered_amp_ma25, xout = input_period)$y / input_amp


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
  lines(input_period, response_ma25, type = "l", col = "darkgrey", lwd = 1,
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

  realVsSynthetic <- bind_rows(
  data.frame( Input = input_period_synthetic, Response = response_bwf_synthetic, Filter = "Butterworth", Dataset = "Synthetic" ),
  data.frame( Input = input_period_real, Response = response_bwf_real, Filter = "Butterworth", Dataset = "Real" ),
  data.frame( Input = input_period_synthetic, Response = response_cci_synthetic, Filter = "CCI", Dataset = "Synthetic" ),
  data.frame( Input = input_period_real, Response = response_cci_real, Filter = "CCI", Dataset = "Real" ))

  require(MASS) # to access Animals data sets
  require(scales) # to access break formatting functions
  #data(Animals) # load data
  realVsSynthetic %>% ggplot(aes(x = Input, y = Response, colour = Filter)) +
    geom_line()+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                   labels = trans_format("log10", math_format(10^.x)), limits = c(1,1e3)) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    theme_bw() +
    #annotation_logticks() +
    facet_grid(rows = vars(Filter), cols = vars(Dataset))
    #xlim(NA, 100)




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



  godinsum <- changeInterval(dplyr::select(randomts, c(Time, Godin)), option = "sum")
  bwfsum <-changeInterval(dplyr::select(randomts, c(Time, bwf)), option = "sum")
  ccisum <-changeInterval(dplyr::select(randomts, c(Time, cci)), option = "sum")
  rawsum <-changeInterval(dplyr::select(randomts, c(Time, TidedSignal )), option = "sum")

  plot(rawsum$Date, bwfsum$accum)
  points(rawsum$Date,   godinsum$accum, col ="blue")
  points(rawsum$Date,   ccisum$accum, col = "red")

  ( bwfsum[nrow(bwfsum),] - rawsum[nrow(rawsum),] ) /  rawsum[nrow(rawsum),] * 100
  ( ccisum[nrow(ccisum),] - rawsum[nrow(rawsum),]  ) /  rawsum[nrow(rawsum),] * 100
  ( godinsum[nrow(godinsum),] - rawsum[nrow(rawsum),] ) /  rawsum[nrow(rawsum),] * 100
  ( rawsum[nrow(rawsum),] - rawsum[nrow(rawsum),]  ) /  rawsum[nrow(rawsum),] * 100

  "cci = 0.013%"
  "bwf = 0.014%"
  "godin = 0.016%"


  # Impulse response as per Roberts & Roberts 1978
  # To examinet he transientr esponsea, n impulseo f magnitude
  # 50 at time 200 has been applied to each of the filters. The
  # results are shown in Figure 3.

  df <- data.frame(Time = seq(Sys.time() %>% round("hour"), by = 60*60, length.out = 512),
             Q = rep(0,512))
  df[200,2] <- 50
  plot(df, ylim = c(-1,4), type = "l")

  #lines(butterworthFilter(df), col = "orange")
  bwf <- butterworthFilter(df)
  df$bwf <- bwf$Filtered
  cci <- ccInterpFilter(df)
  df$cci <- approx(cci$Date, cci$avg, df$Time)$y

  df$godin <- godinFilter(df$Q)

  df <- melt(df, id.vars = "Time")

  # impulse response
  head(df)
  df %>% ggplot(aes(x = Time, y = value, color = variable)) +
    geom_line() +
    ylim (c(-1,5))


}
