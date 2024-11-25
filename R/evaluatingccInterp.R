if(FALSE)
{
library(ccInterp)

################################
# synthetic
randomts <- StevesCoolRandomTS(maxFlow = 1000, obs = 10000, maxNoise = 200)
#plot(randomts$Time, randomts$Signal)

par(mfrow=c(1,2))
plot(randomts$Time, randomts$Noise)
plot(randomts$Time, randomts$Signal)

randomts$Time

# Plot the input signal
plot(randomts$Time, randomts$Signal, type = "l", col = "blue", main = "Input Signal",
     xlab = "Time (hours)", ylab = "Amplitude")


randomts$TidedSignal <- randomts$Signal + randomts$Noise
randomts

plot(randomts$Time, randomts$TidedSignal)

#######################
# apply tidal filter

# Define a simple moving average filter
filter_length <- 25  # Length of the moving average
filtered_signal <- stats::filter(randomts$TidedSignal, rep(1 / filter_length, filter_length), sides = 2)
randomts$ma25Filter <- filtered_signal

bwf <- butterworthFilter(data.frame(randomts$Time,  randomts$TidedSignal))
randomts$bwf <- bwf$Filtered


cci <- ccInterpFilter(data.frame(randomts$Time,  randomts$TidedSignal))

lines(cci$Date, cci$avg, col = "blue")
#plot(cci$Date, cci$avg)

randomts$cci <- approx(cci$Date, cci$avg, randomts$Time)$y
#View(randomts)

# Plot the input signal
plot(randomts$Time, randomts$Signal, type = "l", col = "black", main = "Input Signal",
     xlab = "Time (hours)", ylab = "Amplitude")


# Plot the filtered signal
#lines(randomts$Time, randomts$ma25Filter, col = "red", lwd = 2)
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

input_period_synthetic <- input_period
response_bwf_synthetic <- response_bwf
response_cci_synthetic <- response_cci

# Plot the response function
plot(input_period, response_bwf, type = "l", col = "orange", lwd = 2,
     main = "Response Function", xlab = "Period (hours)",
     ylab = "Response (Amplitude Ratio)",
     xlim = c(0,60), log = "y")
abline(h = 1, col = "gray", lty = 2)  # Ideal response line
# Plot the response function
lines(input_period, response_cci, type = "l", col = "blue", lwd = 2,
     main = "Response Function", xlab = "Period (hours)",
     ylab = "Response (Amplitude Ratio)")
abline(h = 1, col = "gray", lty = 2)  # Ideal response line
lines(input_period, response_gdn, type = "l", col = "red", lwd = 2,
      main = "Response Function", xlab = "Period (hours)",
      ylab = "Response (Amplitude Ratio)")
abline(h = 1, col = "gray", lty = 2)  # Ideal response line


# Plot the response function
plot(input_period, response_bwf, type = "l", col = "orange", lwd = 2,
     main = "Response Function", xlab = "Period (hours)",
     ylab = "Response (Amplitude Ratio)",
     xlim = c(15,1000), ylim = c(0,1), log = "x")
abline(h = 1, col = "gray", lty = 2)  # Ideal response line
# Plot the response function
lines(input_period, response_cci, type = "l", col = "blue", lwd = 2,
     main = "Response Function", xlab = "Period (hours)",
     ylab = "Response (Amplitude Ratio)")
abline(h = 1, col = "gray", lty = 2)  # Ideal response line
lines(input_period, response_gdn, type = "l", col = "red", lwd = 2,
      main = "Response Function", xlab = "Period (hours)",
      ylab = "Response (Amplitude Ratio)")
abline(h = 1, col = "gray", lty = 2)  # Ideal response line

####################################


plot(randomts$Signal, randomts$bwf, col = "orange")
points(randomts$Signal, randomts$Godin, col = "red")
points(randomts$Signal, randomts$cci, col = "blue")
abline(0, 1)



######################
################################
# real data

mrdQ <- readRDS("data/MulgraveQ.rds")

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


cci <- ccInterpFilter(data.frame(mrdQ$time,  mrdQ$value_Discharge))
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
plot(input_period, response_bwf, type = "l", col = "orange", lwd = 2,
     main = "Response Function", xlab = "Period (hours)",
     ylab = "Response (Amplitude Ratio)",
     xlim = c(0,48), log = "y")

abline(h = 1, col = "gray", lty = 2)  # Ideal response line
# Plot the response function
lines(input_period, response_cci, type = "l", col = "blue", lwd = 2,
      main = "Response Function", xlab = "Period (hours)",
      ylab = "Response (Amplitude Ratio)")
abline(h = 1, col = "gray", lty = 2)  # Ideal response line
lines(input_period, response_gdn, type = "l", col = "red", lwd = 2,
      main = "Response Function", xlab = "Period (hours)",
      ylab = "Response (Amplitude Ratio)")
abline(h = 1, col = "gray", lty = 2)  # Ideal response line


# Plot the response function
plot(input_period, response_bwf, type = "l", col = "orange", lwd = 2,
     main = "Response Function", xlab = "Period (hours)",
     ylab = "Response (Amplitude Ratio)",
     xlim = c(15,1000), ylim = c(0,1), log = "x")
abline(h = 1, col = "gray", lty = 2)  # Ideal response line
# Plot the response function
lines(input_period, response_cci, type = "l", col = "blue", lwd = 2,
      main = "Response Function", xlab = "Period (hours)",
      ylab = "Response (Amplitude Ratio)")
abline(h = 1, col = "gray", lty = 2)  # Ideal response line
lines(input_period, response_gdn, type = "l", col = "red", lwd = 2,
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
