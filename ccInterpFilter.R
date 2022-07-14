
#WQI_timeseriesfunctions

# function for looking up the maximum y-value in a reference dataset, for a time at time of newts
# newts will be a forward looking mean, i.e. 10:00, 11:00, 12:00.. 10:00 will be the mean between 10:00 and 11:00
# And the function result will be the maximum in the reference period between 10:00 and 11:00

# This is useful for quality coded data that has been averaged.
# the maximum quality code can be extracted from the original data and applied to the averaged data.

# also useful for a daily maximum value
# or a daily minimum value

maxminfun <- function(refx, refy, newts, option="max", dt=2)
{ 
  #refx <- ts$Timestamp
  #refy <- ts$QC
  #newts <- spinterpData$Date
  #dt=1
  #option <- "max"
  
  originalnewts <- newts
  newts <- as.data.frame(newts)
  #newts[1]
  #newts <- spinterpData$Date
  #qcdata <- maxminfun(ts$Timestamp, ts$QC, , dt=1 )
  if (dt==1)
  {
    subtract <- (newts[2,1] -newts[1,1])/2
    #newts<- as.numeric(newts)
    #shift the target timeseries back by half an interval
    # function assumes default dt of 2, which is forward means.
    newts[,1] <- newts[,1] - (subtract )
  }
  
  #refx - reference x lookup
  #refy - reference y lookup i.e. a quality code
  #newts - a timeseries dataframe, with at least a column of timestamps
  #option - "max" or "min"
  
  maxfun <- function(t, y, x, option="max") #look forward and back to pick the max (or min)
  {

    # function looks to either side of point x and outputs the minimum
    f.qual0 <-  approxfun(t, y, method = "constant", ties=option,  f=0)
    f.qual1 <- approxfun(t, y, method = "constant", ties=option, f=1)
    
    if (option=="max") {
      return( pmax(f.qual0(x),  f.qual1(x)) )
    }else if(option=="min")
      return( pmin(f.qual0(x),  f.qual1(x)) )
  }
  
  oldts <- data.frame(refx, refy)
  #newts <- data.frame(Timestamp = newts ) # why did I add this?
  #colnames(newts)[1]
  
  # Function of hourly hours 
  f.hourlyHours <- approxfun(newts[,1], newts[,1], method="constant")
  # which hour does this belong
  oldts <- cbind (oldts, hour = f.hourlyHours(oldts[,1]) )
  # add max QC for period to df
  newts <- cbind(newts, QC = maxfun(oldts[,3], oldts[,2], newts[,1], option=option))

  #f.newts <- approxfun(newts, method="constant")
  #return(f.newts)
  
  if (dt==1)
  {
    newts[,1] <-  newts[,1]  + subtract
  }
  
  return(newts)
  
}


longts <- seq( 0, 500, by=1 )
# create a random timeseries
randts <- runif(length(longts), min = 0, max = 1000)

plot(longts, randts)
lines(longts, randts, col="grey")

intervalts <- seq( 0, 500, by=20 )
intervalMax <- maxminfun(longts, randts, intervalts)
f.max <- approxfun(intervalMax, method="constant")

lines(  f.max(longts), col="red" )

intervalMin <- maxminfun(longts, randts, intervalts, option="min")
f.min <- approxfun(intervalMin, method="constant")

lines(  f.min(longts), col="red" )




# spinterpConvert

# used for upsampling interval data
# desinged for outputing daily values as hourly.


library(dplyr)
library(zoo)
library(pracma)

# start is the start date of the averaging period, in days
# rate is the rate, ie megalitres per day
# outputInt is the new time interval required, i.e. hourly, being 1/24 days.
spinterpConvert <- function(start, rate, outputInt=(1/24), type="spinterp", dt=2)
{
  # inputs must be stepped averages (default, dt = 2)
  # If point data is used (not preferred as precision can be lost), set dt = 1
  
  # not appropriate for downsampling intervals from hourly back to daily
  # Use changeInterval() for that. 
  # (even though theoretically spinterpConvert will work fine for downsampling, and in some cases technically preferred)
  # Using "spinterp" Interpolation over long gaps could have unpredicatable results.
  # If necessary, smooth the input data first before using spinterp interpolation, 
  # or else the data noise could sway the result too much
  
  # start is the start date of the averaging period i.e. start of the day, start of hour.
  # time in days
  
  # spinterpConvert simply gets the derivative of the cumulative discharge
  
  
  if (dt==1)# if data is of point type, convert to an interval mean at from start point to end point
  {
    #start <- df$days
    #rate <- df$df...2.
    
    means <- data.frame(start, rate)
    means$end <- c(means$start[-1],0)
    means$rate2 <- c(means$rate[-1],0)
    means <- head(means,-1)
    tail(means)
    
    #get average
    means$avg <- (means$rate+means$rate2)/2
    start <- means$start
    rate <- means$avg
  }
  if (dt==3) # mean, average flow until this point
  {
    start <- start - (start[2] - start[1])
  }
  
  start <- as.numeric(start)
  nt <- length(start)
  t <- c(start, start[nt]+1)
  
  # duration between input data points # in days
  dur <- c(diff(t),0)   #t[2]-t[1]  #c(0,diff(t))
  
  # create a time sequence for output data
  xp <- seq(t[1], t[nt], by=outputInt)
  
  #length(dur)
  #length(rate )
  #length(t )
  #length(c(0,rate))
  #(length(cumdaily))
  negativeOffset <- min(rate)
  rate <- rate - negativeOffset
  
  # cumulative sum of discharge
  # from "cubic metres per second" into "cubic metres"
  DAYSTOSECONDS <- 60*60*24
  cumdaily <- cumsum( c(0, rate * head(dur, -1))) # * DAYSTOSECONDS)
  
  #t<- head(t,-1)
  #length(cumdaily)
  #plot(t, cumdaily)
  
  # output at new time interval
  if(type=="spinterp")
  {
    #print("spinterp")
    yp <- spinterp(t, cumdaily, xp)
  }else if(type=="spline")
  {
    #print("spline")
    f.yp <- splinefun(t, cumdaily, method="hyman")
    yp <- f.yp(xp)
  }else if(type=="schum")
  {
    #print("schum")
    SchumSpline = schumaker::Schumaker(t,cumdaily)
    yp <- SchumSpline$Spline(xp)
  }else if(type=="linear")
  {
    # Not useful to interpolate linearly
    f.yp <- approxfun(t, cumdaily)
    yp <- f.yp(xp)
  }
  
  length(rate)
  #plot(start, c(rate) )
  
  #plot(D)
  #plot(xp, yp)
  
  # derivative of spinterp
  f.spinterp <- splinefun(xp, yp)
  f.spinterp(xp, deriv=1)
  
  
  #f.dur <- approxfun(t, dur, method="constant", f=0)
  #df <- data.frame( Date=xp,  Data=   f.spinterp(xp, deriv=1), dur= round(f.dur(xp),4)  )
  df <- data.frame( Date=xp,  Data = f.spinterp(xp, deriv=1)+negativeOffset  )     #  /DAYSTOSECONDS  )
  #plot(xp, f.dur(xp))
  
  return(df)
  
  #dygraph(df)
  
  #plot(convertdata)
  #lines(df, col="red")
  
}



# changeInterval
# Stephen Wallace 17/03/2022

# changeInterval will change the timestep of any timeseries input to an evenly spaced timestep.
# inputs don't need to be evenly spaced, and it will use spline interpolation for gaps in data or NA values.

# ts input is a dataframe of posixct time (time in seconds) and an instantaneous value (per second i.e. cumecs)

# dt (datatype) of the input data, 1 = inst; 2 = fmean
# default is 1 for instantaenous, the inputs can be randomly spaced in time.
# 2 for fmean when you are converting evenly spaced forward means such as daily averages.

# Interval "Daily", "Hourly" or an integer (number of minutes)

# start is a timestamp written, that can be converted to posixct. 
# i.e. "2016-12-17 05:00" or numeric, seconds since "1970-01-01" i.e. 1481914800
# end is the same as start

# default to zero for start and end to juse use the same time as inputs

# offset is a positive or negative integer, is the output offset, in minutes
#Interval <- "Daily"
#start <- 0
#end <- 0

# option = "fmean", "inst", "sum"
# fmean will output a forward looking mean at that timestamp, 
# i.e. a value for 05/01/2010 will be the mean between 05/01/2010 to 06/01/2010.
# inst will output an instantaenous mean, i.e a midday reading on 05/01/2010 as a midpoint between.

# (sum needs testing, but theoretically should be right)
# sum will calculate the total change from the beginning of the period to the end of the period


changeInterval <- function(ts, dt=1, Interval="Daily", start=0, end=0, offset=0, option="fmean", rounded=TRUE)
{
 
  inputts <- ts
  
  #function omits na's, so will interpolate between
  ts <- na.omit(ts)
  ts <- ts[,1:2]

  
  if (start==0)
  {
    if (Interval == "Daily")
    {
      # round to full day
      start <- round(ts[1,1], units="days")
    }else if (Interval == "Hourly")
    {
      # round to full hour
      start <- round(ts[1,1], units="hours")
    }else{
      # start at raw value
      
      if (rounded == TRUE){
        start <- round(ts[1,1], units="hours")
      }else{
        start <- ts[1,1]
      }
    }
    
  }else
  {
    start <- as.POSIXct(start, origin="1970-01-01")
    
  }
  if ( end == 0)
  {
    end <- ts[nrow(ts),1]
  }else{
    end <- as.POSIXct(end, origin="1970-01-01")
  }
  
  
  
  #duration between this point and the previous point
  ts$dur <- c(0,diff(as.numeric(ts[,1] ) ) )  #duration in seconds
  
  
  # calculate a volume
  if (dt==2)
  {
    # means
    ts$megalitres <- ts$dur*ts[,2]/1000
  }
  if (dt==1)
  {
    # point data, or interval data (means, with value at centre point)
    v <- c(ts$dur,0) * c(ts[,2],0) + c(ts$dur,0) *  c(0,ts[,2])
    v <- v/2
    v <- head(v,-1)
    
    ts$megalitres <- v/1000
  }
  
  # cumulative sum
  ts$accum <- cumsum(ts$megalitres)
  
  # cumulative sum lookup function
  f.accum <- splinefun(ts[,1], ts$accum)
  
  
  # set time step based on interval
  if (Interval == "Daily"){ timestep = 24*60*60}else 
    if (Interval == "Hourly"){ timestep = 1*60*60}else{
    # assume an integer
      timestep = Interval*60}
  
  # new time sequence for new interval
  newintTS <- seq(start+offset*60, end, by= timestep )
  
  # delta y
  dy <-  c(diff(f.accum(newintTS) ),0 ) 
  # delta y / delta time
  if(option=="fmean"){
    df <- data.frame(Date=newintTS, FMean = round(dy/(timestep/1000),3))
  }else if(option=="sum")
  {
    # output total difference in cumulative 
    df <- data.frame(Date=newintTS, Sum = round(dy,3) )
  }else if(option=="inst")
  {
    #add half a timestep to report an instantaneous mean
    df <- data.frame(Date=newintTS+(timestep/2), Inst = round(dy/(timestep/1000),3))
  }
  
  if ( length(inputts) == 3 )
  {
    df <- maxminfun(inputts[,1], inputts[,3], df, option="max" )
  }
  
  df <- na.trim(df)
  return(df)
  
  
  #plot(df$Date, df$cumecs)
  
}



library(matrixStats)
library(zoo)

# Combined cumulative interpolation filter
# Stephen Wallace
# 18/03/2022

# designed for removing tidal signal for discharge data
# can be used as a standard low pass filter
# can remove noise 

# requires spinterpConvert and changeInterval

ccInterpFilter <- function(ts, hours = 24, discardbelowzero = FALSE, centred = FALSE, type="spinterp")
{
  # converts time series to daily, and cumulates, interpolates, and then gets the derivative.
  # does this for each hour, i.e. for 24 hour averaging it will do it 24 times offset an hour each time
  # the result is the average of all hours interpolated.
  
  # ts = dataframe of posixct time, and a value.
  
  #ccInterpFilter(tsdf, 24, centred=FALSE)
  #ts <- tsdf
  #hours <- 24
  #centred <- FALSE
  #discardbelowzero = FALSE
  #type="spinterp"
  
  spinterpData <- 0 #initialize
  
  #Loop through all daily interpolations
  for(i in 0:(hours-1))
  {
    #print(i)
    # for each 24 hour interpolation, increment the offset by an hour each time
    daily <- changeInterval(ts, Interval=hours*60, offset=i*60)
    offset <- 0
    if (discardbelowzero==FALSE){
      # spinterp can't do negative cumulative dishcharge, so add an offset to prevent any negative flows
      offset <- min(daily$FMean)
      daily$FMean <- daily$FMean - offset
    }
    
    daily$FMean[daily$FMean < 0] <- 0 
    daily[is.na(daily)] <- 0
    #plot(daily)
    # convert to numeric
    numdate <- as.numeric(as.POSIXct(daily$Date, format=format)) /60/60/24
    # create hourly spinterp data ( interpolation of cumulative daily discharge )
    spinterpSegment <- spinterpConvert(numdate, daily$FMean, type=type)
    # remove offset again
    spinterpSegment$Data <- spinterpSegment$Data+offset
    
    if (i>0){# after first row
      
      pad <- nrow(spinterpData) - length(c(rep(NA, i), head(spinterpSegment$Data, -i)) )
      if ( pad > 0)
      { # pad with na's
        spinterpData <- cbind(spinterpData, i=c(rep(NA, i), head(spinterpSegment$Data, -i), rep(NA, pad)   ))
      }else{
        spinterpData <- cbind(spinterpData, i=c(rep(NA, i), head(spinterpSegment$Data, -i)))
      }
      
    }else{# first row
      spinterpData <- spinterpSegment
    }
  }
  
  setcolnames <- paste("Col", 0:(hours-1), sep = "")
  setcolnames <- c("Date", setcolnames)
  
  colnames(spinterpData) <- setcolnames
  # trim na's
  spinterpData <- na.trim(spinterpData)
  
  # centred == TRUE = bias towards inner averages, tested but not useful for tidal filter
  if (centred == TRUE)
  {
    avgS <- round(hours/4, 0)
    avgF <- hours-avgS
    
    spinterpData <- cbind(spinterpData, avg = rowMeans(spinterpData[avgS+2:avgF+1]))
  }else{
    spinterpData <- cbind(spinterpData, avg = rowMeans(spinterpData[-1]))
  }
  
  
  spinterpData$Date <- as.POSIXct(spinterpData$Date*60*60*24, origin="1970-01-01")
  return(spinterpData)
}




###################################################################
# spinterpConvert example

# daily times and daily means
D <- list( x=c(16832.96,16833.96,16834.96,16835.96,16836.96,16837.96,16838.96,16839.96,16840.96,16841.96,16842.96,16843.96,16844.96,16845.96,16846.96,16847.96,16848.96,16849.96,16850.96,16851.96,16852.96,16853.96,16854.96,16855.96,16856.96,16857.96,16858.96,16859.96,16860.96,16861.96,16862.96,16863.96,16864.96,16865.96,16866.96,16867.96,16868.96,16869.96,16870.96,16871.96,16872.96,16873.96,16874.96,16875.96,16876.96,16877.96,16878.96,16879.96,16880.96,16881.96,16882.96,16883.96,16884.96,16885.96,16886.96,16887.96,16888.96,16889.96,16890.96,16891.96),
           y=c(  0,2.707,1.3795,211.932,22.565,4.409,2.536,1.162,0.66,0.358,0.291,0.1065,0.099,0.0265,0.0055,0,0,0,0,2.195,11.208,3.4125,1.3465,0.657,0.3255,0.2915,0.1235,0.0985,0.058,3.629,111.4225,148.32,144.8815,99.8055,45.8355,43.858,13.0915,5.876,3.6715,1.5845,0.8965,13.5315,2.902,1.608,0.7085,0.428,0.1965,0.2055,0.07,0.057,0.0235,0.007,0.006,0.001,0,0,0,0,0,0)
)

f.D <- approxfun(D, method="constant")

plot(D, ylim=c(0,450))
ts <- seq (first(D$x),last(D$x), by = (1/24 ) )
lines(ts, f.D(ts), ylim=c(0,450) )
# upsample to hourly 
# default = spinterp
DHourly <- spinterpConvert(D$x, D$y)
lines(DHourly, col="red")

# type = spline
DHourly <- spinterpConvert(D$x, D$y, type = "spline")
lines(DHourly, col="blue")

# type = linear
# DHourly <- spinterpConvert(D$x, D$y, type = "linear")
# lines(DHourly, col="green")



###################################################################
# changeInterval example


DDays <- data.frame(
  x = c(17151.61,17151.62,17151.63,17151.64,17151.65,17151.66,17151.67,17151.69,17151.7,17151.71,17151.72,17151.73,17151.74,17151.75,17151.76,17151.77,17152.07,17152.08,17152.09,17152.1,17152.11,17152.12,17152.13,17152.14,17152.15,17152.16,17152.17,17152.19,17152.2,17152.21,17152.22,17152.23,17152.24,17152.25,17152.26,17152.27,17152.28,17152.29,17152.3,17152.31,17152.32,17152.33,17152.34,17152.64,17152.65,17152.66,17152.67,17152.69,17152.7,17152.71,17152.72,17152.73,17152.74,17152.75,17152.76,17152.77,17152.78,17152.79,17152.8,17153.11,17153.12,17153.13,17153.14,17153.15,17153.16,17153.17,17153.19,17153.2,17153.21,17153.22,17153.23,17153.24,17153.25,17153.26,17153.27,17153.28,17153.29,17153.3,17153.31,17153.32,17153.33,17153.34,17153.35,17153.36,17153.69,17153.7,17153.71,17153.72,17153.73,17153.74,17153.75,17153.76,17153.77,17153.78,17153.79,17153.8,17153.81,17153.82,17153.83,17154.13,17154.14,17154.15,17154.16,17154.17,17154.19,17154.2,17154.21,17154.22,17154.23,17154.24,17154.25,17154.26,17154.27,17154.28,17154.29,17154.3,17154.31,17154.32,17154.33,17154.34,17154.35,17154.36,17154.37,17154.38,17154.39,17154.4,17154.72,17154.73,17154.74,17154.75,17154.76,17154.77,17154.78,17154.79,17154.8,17154.81,17154.82,17154.83,17154.84,17154.85,17154.86,17154.87,17154.88,17154.89,17155.17,17155.19,17155.2,17155.21,17155.22,17155.23,17155.24,17155.25,17155.26,17155.27,17155.28,17155.29,17155.3,17155.31,17155.32,17155.33,17155.34,17155.35,17155.36,17155.37,17155.38,17155.39,17155.4,17155.41,17155.42,17155.44,17155.77,17155.78,17155.79,17155.8,17155.81,17155.82,17155.83,17155.84,17155.85,17155.86,17155.87,17155.88,17155.89,17155.9,17155.91,17155.92,17155.94,17155.95,17156.24,17156.25,17156.26,17156.27,17156.28,17156.29,17156.3,17156.31,17156.32,17156.33,17156.34,17156.35,17156.36,17156.37,17156.38,17156.39,17156.4,17156.41,17156.42,17156.44,17156.45,17156.46,17156.47,17156.81,17156.82,17156.83,17156.84,17156.85,17156.86,17156.87,17156.88,17156.89,17156.9,17156.91,17156.92,17156.94,17156.95,17156.96,17156.97,17156.98,17156.99,17157,17157.29,17157.3,17157.31,17157.32,17157.33,17157.34,17157.35,17157.36,17157.37,17157.38,17157.39,17157.4,17157.41,17157.42,17157.44,17157.45,17157.46,17157.47,17157.48,17157.49,17157.5,17157.51,17157.84,17157.85,17157.86,17157.87,17157.88,17157.89,17157.9,17157.91,17157.92,17157.94,17157.95,17157.96,17157.97,17157.98,17157.99,17158,17158.01,17158.02,17158.03,17158.04,17158.05,17158.06,17158.33,17158.34,17158.35,17158.36,17158.37,17158.38,17158.39,17158.4,17158.41,17158.42,17158.44,17158.45,17158.46,17158.47,17158.48,17158.49,17158.5,17158.51,17158.52,17158.53,17158.54,17158.86,17158.87,17158.88,17158.89,17158.9,17158.91,17158.92,17158.94,17158.95,17158.96,17158.97,17158.98,17158.99,17159,17159.01,17159.02,17159.03,17159.04,17159.05,17159.06,17159.07,17159.08,17159.09,17159.36,17159.37,17159.38,17159.39,17159.4,17159.41,17159.42,17159.44,17159.45,17159.46,17159.47,17159.48,17159.49,17159.5,17159.51,17159.52,17159.53,17159.54,17159.55,17159.56,17159.57,17159.88,17159.89,17159.9,17159.91,17159.92,17159.94,17159.95,17159.96,17159.97,17159.98,17159.99,17160,17160.01,17160.02,17160.03,17160.04,17160.05,17160.06,17160.07,17160.08,17160.09,17160.1,17160.11,17160.12,17160.4,17160.41,17160.42,17160.44,17160.45,17160.46,17160.47,17160.48,17160.49,17160.5,17160.51,17160.52,17160.53,17160.54,17160.55,17160.56,17160.57,17160.58,17160.59,17160.9,17160.91,17160.92,17160.94,17160.95,17160.96,17160.97,17160.98,17160.99,17161,17161.01,17161.02,17161.03,17161.04,17161.05,17161.06,17161.07,17161.08,17161.09,17161.1,17161.11,17161.12,17161.13,17161.14,17161.15,17161.42,17161.44,17161.45,17161.46,17161.47,17161.48,17161.49,17161.5,17161.51,17161.52,17161.53,17161.54,17161.55,17161.56,17161.57,17161.58,17161.59,17161.6,17161.61,17161.62,17161.92,17161.94,17161.95,17161.96,17161.97,17161.98,17161.99,17162,17162.01,17162.02,17162.03,17162.04,17162.05,17162.06,17162.07,17162.08,17162.09,17162.1,17162.11,17162.12,17162.13,17162.14,17162.15,17162.16,17162.17,17162.19,17162.45,17162.46,17162.47,17162.48),
  y= c(-18.653,-22.085,-24.206,-23.637,-28.192,-26.725,-22.768,-12.975,8.319,35.432,47.791,44.786,41.467,35.165,30.382,25.992,-46.246,-54.979,-61.189,-59.834,-56.147,-53.484,-60.475,-65.256,-66.094,-69.688,-56.973,-55.308,-45.003,-36.956,-4.272,38.82,69.643,73.893,76.769,74.809,60.562,52.592,53.264,43.279,40.532,37.519,34.601,-18.422,-22.107,-25.234,-21.748,-24.197,-28.48,-24.147,-19.258,-12.939,8.822,30.679,36.583,33.512,32.056,27.805,24.003,-40.06,-45.823,-50.406,-48.847,-44.919,-45.859,-44.851,-44.193,-55.121,-49.648,-40.006,-31.659,-17.017,11.996,38.762,52.052,56.798,50.211,56.764,50.937,41.486,40.168,34.349,28.823,26.661,-19.918,-23.127,-22.321,-22.068,-25.69,-27.084,-27.054,-19.449,-6.605,15.779,31.56,32.972,32.734,28.59,24.8,-34.918,-39.973,-40.649,-39.271,-37.607,-34.766,-34.863,-33.661,-37.288,-36.661,-34.074,-36.47,-24.59,-19.154,-4.755,12.045,38.719,43.05,44.297,43.105,39.773,38.471,32.6,29.635,27.276,25.083,23.04,-17.794,-20.639,-21.169,-22.428,-25.567,-27.086,-26.901,-24.971,-19.366,-11.661,7.522,28.536,34.032,35.441,31.242,24.028,21.176,18.528,-24.36,-27.324,-24.681,-26.942,-28.918,-30.315,-24.444,-28.297,-28.172,-28.771,-30.975,-27.921,-21.731,-11.926,-2.526,25.194,32.471,35.17,34.006,36.147,32.109,29.899,26.714,28.599,26.338,24.223,-19.922,-22.893,-25.238,-29.382,-27.249,-27.755,-22.954,-18.401,-13.114,-5.045,15.734,27.33,29.579,31.721,24.39,23.546,20.87,18.394,-24.945,-28.008,-28.297,-28.817,-28.139,-28.293,-24.707,-23.145,-22.493,-16.977,-21.233,-19.091,-7.921,9.781,24.137,31.313,35.027,34.345,33.964,29.236,26.307,23.864,21.58,-25.33,-28.539,-30.156,-28.955,-30.575,-33.141,-29.936,-27.099,-23.409,-8.369,11.395,28.159,36.01,36.517,29.526,26.168,22.7,19.974,17.506,-24.569,-27.928,-27.238,-28.711,-30.122,-31.377,-31.697,-24.081,-22.533,-16.341,-11.964,-0.933,9.693,17.118,24.347,25.299,29.882,26.483,23.418,22.459,19.974,17.633,-28.628,-32.029,-33.529,-41.895,-40.883,-41.28,-44.114,-41.396,-33.762,-29.614,-15.798,1.942,24.487,33.692,40.948,36.355,39.003,33.888,28.442,25.433,24.942,22.349,-24.229,-28.098,-32.408,-29.211,-25.911,-23.632,-20.709,-20.092,-18.368,-13.046,-16.152,-11.59,3.499,15.755,31.07,30.751,31.846,30.849,27.461,24.585,21.416,-35.649,-40.609,-39.35,-45.037,-45.136,-42.603,-50.455,-44.622,-43.647,-36.158,-28.342,-13.979,2.094,27.724,39.423,45.112,38.664,43.163,36.584,33.462,26.619,25.746,22.686,-30.14,-34.817,-31.33,-33.771,-33.053,-28.005,-29.163,-24.753,-24.347,-26.183,-18.007,-12.427,1.28,21.562,32.515,37.375,36.015,28.871,28.338,24.539,21.033,-38.049,-43.806,-44.321,-48.568,-49.435,-47.334,-44.669,-54.716,-49.784,-47.462,-44.683,-26.356,-11.207,19.272,36.543,48.833,44.999,44.683,43.047,36.999,34.196,29.485,26.606,23.937,-25.615,-29.839,-29.219,-28.334,-27.333,-26.002,-28.523,-26.55,-21.499,-14.828,-7.22,4.43,22.857,39.064,41.272,35.283,30.65,28.275,24.193,-45.726,-53.44,-59.432,-50.086,-49.383,-47.175,-48.456,-56.536,-47.769,-46.784,-45.151,-36.815,-20.215,-2.805,27.936,54.861,43.438,45.742,48.601,50.477,45.91,37.939,35.27,31.669,28.385,-29.957,-35.23,-31.174,-39.106,-31.89,-22.075,-13.965,-18.587,-17.09,-22.42,-22.194,1.313,25.126,32.812,38.521,33.45,32.8,31.583,27.376,23.504,-42.185,-48.698,-54.481,-53.145,-51.758,-55.786,-55.003,-55.742,-55.03,-55.094,-43.074,-38.821,-31.601,-5.848,22.325,53.364,62.151,72.12,62.2,54.882,47.956,51.141,39.74,34.343,31.073,28.095,-26.911,-31.507,-29.552,-25.613)
)

D <- DDays
D$x <- as.POSIXct(D$x*60*60*24, origin="1970-01-01")
plot(D[50:100,])

# daily average of forward means
daily <- changeInterval(D, Interval="Daily")
head(daily)
lines(daily)


#spinterpConvert example

# Convert a numeric "days", output int 1/24 days = "hours", as spinterp, with a data type of 1, "Instantaneous"
# dailySpinterp <- spinterpConvert(as.numeric(D$x)/24/60/60, D$y, outputInt = (1/24), type = "spinterp", dt=1)
# lines(dailySpinterp$Date*24*60*60, dailySpinterp$Data, col="green")
# dailySpinterp <- spinterpConvert(as.numeric(D$x)/24/60/60, D$y, outputInt = (1/24), type = "spline", dt=1)
# lines(dailySpinterp$Date*24*60*60, dailySpinterp$Data, col="green")


# the same interval, but without automatically adjusting the offset to output by using "Daily" as interval
twentyfourhour <- changeInterval(D, Interval = 24*60)
head(twentyfourhour)


# output an hourly forward mean (timestamp at beginning of period ), 
hourly <- changeInterval(D, Interval="Hourly", option="fmean")
f.hourly <- approxfun(hourly$Date, hourly$FMean, method = "constant")
ts <- seq(first(hourly$Date), last(hourly$Date), by=10*60)
lines(ts, f.hourly(ts), col="red")

# output an hourly instanaenous (timestamp in middle of averaging period ), 
# offset average window by 30 minutes to output datapoints on the hour.
# essentially it is the average of the period 30 mins before and 30 mins after the timestamp.
hourlyInst <- changeInterval(D, Interval="Hourly", offset=0, option="inst")
head(hourlyInst)
lines(hourlyInst, col="orange")

#hourlyInst <- changeInterval(D, Interval="Hourly", offset=-30, option="inst")
#head(hourlyInst)
#lines(hourlyInst, col="blue")




###################################################################
# ccInterpFilter example

plot(D)

#as a tidal filter
spinterpData <- ccInterpFilter(D, 24, centred=FALSE)
lines(spinterpData$Date, spinterpData$avg, col="red")

# 95% confidence
stdev <- rowSds(as.matrix(spinterpData[2:(length(spinterpData) -1 )]  ))
spinterpConfi <- data.frame(spinterpData$Date, spinterpData$avg)
spinterpConfi <- cbind(spinterpConfi, upper = spinterpData$avg + stdev )
spinterpConfi <- cbind(spinterpConfi, lower = spinterpData$avg - stdev )

lines(spinterpConfi$spinterpData.Date, spinterpConfi$upper)
lines(spinterpConfi$spinterpData.Date, spinterpConfi$lower)




#as a noise removing filter
#keep tide but smooth
spinterpData <- ccInterpFilter(D, 2, centred=FALSE)
plot(D[50:100,])
lines(spinterpData$Date, spinterpData$avg, col="red")


# 95% confidence
stdev <- rowSds(as.matrix(spinterpData[2:(length(spinterpData) -1 )]  ))
spinterpConfi <- data.frame(spinterpData$Date, spinterpData$avg)
spinterpConfi <- cbind(spinterpConfi, upper = spinterpData$avg + stdev )
spinterpConfi <- cbind(spinterpConfi, lower = spinterpData$avg - stdev )

lines(spinterpConfi$spinterpData.Date, spinterpConfi$upper)
lines(spinterpConfi$spinterpData.Date, spinterpConfi$lower)

#dygraph(spinterpConfi)


#############################################
# Validation
#
# rm(PagendamPercivalSyntheticData)
# http://dx.doi.org/10.4225/08/543B4B1A57C33
# 
##############################################
#
#
#

library(dygraphs)
load("C:/Users/wallacesw/Downloads/RData_file_used_in_statistical_analyses/SyntheticData.RData")

PPS <- data.frame(Time = PagendamPercivalSyntheticData[,1])
PPS$Synthetic <- PagendamPercivalSyntheticData[,2]
PPS$Fresh <- PagendamPercivalSyntheticData[,3]
PPS$Tide <- PPS$Synthetic - PPS$Fresh

dygraph(PPS)

# Times didn't make sense, 4 tide cycles made up 6.15 days instead of an estaimted 4.13 days for 4 cycles
# adjusted accordingly, but largely arbitrary
PPS$Time <- PPS$Time * ((4*((24+5/6)/24))/6.1466)
dygraph(PPS)



# Convert to posixctime from days
PPS$Time <- as.POSIXct(PPS$Time*24*60*60 + as.numeric(as.POSIXct("2008-03-14")) , origin="1970-01-01")

# Plot
plot(PPS$Time, PPS$Synthetic)
lines(PPS$Time, PPS$Tide, col="grey" )
lines(PPS$Time, PPS$Fresh, col="red")

# Set na's to zero
PPS$Synthetic[is.na(PPS$Synthetic)] <- 0

# CCInterp Filter, daily means
ccInterpData <- ccInterpFilter(data.frame(PPS$Time, PPS$Synthetic), hours = 24, type="spinterp")
lines(ccInterpData$Date, ccInterpData$avg, col="green")

# 95% confidence
stdev <- rowSds(as.matrix(ccInterpData[2:(length(ccInterpData) -1 )]  ))
spinterpConfi <- data.frame(ccInterpData$Date, ccInterpData$avg)
spinterpConfi <- cbind(spinterpConfi, upper = ccInterpData$avg + stdev )
spinterpConfi <- cbind(spinterpConfi, lower = ccInterpData$avg - stdev )

# Plot Upper and Lower bounds of 95% confidence levels
lines(spinterpConfi$ccInterpData.Date, spinterpConfi$upper)
lines(spinterpConfi$ccInterpData.Date, spinterpConfi$lower)


# Butterworth for comparison
Hourly <- changeInterval(data.frame(PPS$Time, PPS$Synthetic), Interval = "Hourly", option = "inst")
lines(Hourly$Date, butterworthFilter(Hourly$Inst), col="orange")



#fapproxfun(ccInterpData$Date, ccInterpData$avg)

plot(PPS$Fresh, ccInterpData$avg)

# Timestamps, every ten minutes
TenMinTS <- seq(PPS$Time[1], PPS$Time[nrow(PPS)], by = 10*60 )
SixtyMinTS <- seq(PPS$Time[1], PPS$Time[nrow(PPS)], by = 60*60 )

# lookup functions
f.Raw <- approxfun(PPS$Time, PPS$Synthetic)
f.Actual <- approxfun(PPS$Time, PPS$Fresh)
f.CCInterp <- approxfun(ccInterpData$Date, ccInterpData$avg)
f.Butter <- approxfun(Hourly$Date, butterworthFilter(Hourly$Inst) )
f.Godin <- approxfun(Hourly$Date, godinFilter(Hourly$Inst) )
f.Upper <- approxfun(spinterpConfi$ccInterpData.Date, spinterpConfi$upper )
f.Lower <- approxfun(spinterpConfi$ccInterpData.Date, spinterpConfi$lower )

par(mfrow= c(3,1))
plot(f.Actual(SixtyMinTS), f.CCInterp(SixtyMinTS)-f.tidebias(SixtyMinTS), ylim=c(-50, 450), type="l")
lines(f.Actual(SixtyMinTS), f.CCInterp(SixtyMinTS), ylim=c(-50, 450), col="grey")
abline(0,1,col="red")
plot(f.Actual(SixtyMinTS), f.Butter(SixtyMinTS)-f.tidebias(SixtyMinTS), ylim=c(-50, 450), type="l")
lines(f.Actual(SixtyMinTS), f.Butter(SixtyMinTS), ylim=c(-50, 450), col="grey")
abline(0,1,col="red")
plot(f.Actual(SixtyMinTS), f.Godin(SixtyMinTS)-f.tidebias(SixtyMinTS),ylim=c(-50, 450), type="l")
lines(f.Actual(SixtyMinTS), f.Godin(SixtyMinTS),ylim=c(-50, 450), col="grey")
abline(0,1,col="red")


model <- lm( (f.CCInterp(TenMinTS)-f.tidebias(TenMinTS)) ~ (f.Actual(TenMinTS)) )
summary(model)
sqrt(mean(model$residuals^2))

model <- lm( (f.Butter(TenMinTS)-f.tidebias(TenMinTS)) ~ f.Actual(TenMinTS) )
summary(model)

model <- lm( (f.Godin(TenMinTS)-f.tidebias(TenMinTS)) ~ f.Actual(TenMinTS) )
summary(model)


library(xts)
# create a dygraph xts object for visualising
dyPlot <- xts(x = f.Raw(TenMinTS), order.by = TenMinTS)
dyPlot <- cbind(Raw = dyPlot, CCInterp = xts(x = f.CCInterp(TenMinTS), order.by = TenMinTS) )
dyPlot <- cbind(dyPlot, Butter = xts(x = f.Butter(TenMinTS), order.by = TenMinTS) )
dyPlot <- cbind(dyPlot, Godin = xts(x = f.Godin(TenMinTS), order.by = TenMinTS) )
dyPlot <- cbind(dyPlot, Actual = xts(x = f.Actual(TenMinTS), order.by = TenMinTS) )
dyPlot <- cbind(dyPlot, UpperBound = xts(x = f.Upper(TenMinTS), order.by = TenMinTS) )
dyPlot <- cbind(dyPlot, LowerBound = xts(x = f.Lower(TenMinTS), order.by = TenMinTS) )

#dygraph(dyPlot)
dygraph(dyPlot[,-1])


#####################################################
#
# without tide tide bias removed

# subtract daily tide bias
f.CCInterp <- approxfun(ccInterpData$Date, ccInterpData$avg)
plot(TenMinTS, f.CCInterp( TenMinTS )  )
lines(TenMinTS,  f.Upper( TenMinTS ))
lines(TenMinTS, f.Lower((TenMinTS) ))

# data frame
rawDischarge <- data.frame( Date = TenMinTS,  CCIFiltered = f.CCInterp( TenMinTS )  )
rawDischarge <- cbind(rawDischarge,  Upper =  f.Upper( TenMinTS ))
rawDischarge <- cbind(rawDischarge,  Lower = f.Lower( (TenMinTS) ))
rawDischarge <- cbind(rawDischarge,  Actual = f.Actual( (TenMinTS) ) )
# trim missing edges
rawDischarge <- na.trim(rawDischarge)

TimeAboveUpper <- (rawDischarge$Actual - rawDischarge$Upper )
TimeAboveUpper[TimeAboveUpper <= 0] <- 0
plot(TimeAboveUpper)


TimeBelowLower <- (rawDischarge$Actual - rawDischarge$Lower )
TimeBelowLower[TimeBelowLower >= 0] <- 0
plot(TimeBelowLower)

#TimeAboveUpper[ is.na(TimeAboveUpper)  ] <- 0
#TimeBelowLower[ is.na(TimeBelowLower)  ] <- 0
plot( cumsum(TimeAboveUpper) )
lines( abs(cumsum(TimeBelowLower) ))

# Plot total volume (black), adding total volume above (blue), subtracting total volume below (red)
plot(rawDischarge$Date, cumsum(rawDischarge$Actual), type="l")
lines(rawDischarge$Date, cumsum(rawDischarge$Actual)+cumsum(TimeAboveUpper), col="blue" )
lines(rawDischarge$Date, cumsum(rawDischarge$Actual)-abs(cumsum(TimeBelowLower) ), col="red")


# Total cumulative absolute discharge, of periods above and below the filtered discharge
(sum(TimeAboveUpper)  + abs(sum(TimeBelowLower))) / sum(rawDischarge$Actual)


message(paste( (sum(TimeAboveUpper)  + abs(sum(TimeBelowLower))) / sum(rawDischarge$Actual) * 100, "%", sep="" ))
# 8.631436%
#####################################################
# However, The tide data is biased, due to collection issues
# Theoretically, an unbiased tide will have a daily average of zero

# tide bias
plot (TenMinTS,  f.Raw(TenMinTS) - f.Actual(TenMinTS) )

# convert ten minute "tide only" to a daily value
# this is the extent the bias applies that day
DailyTideBias <- ( changeInterval(data.frame( as.POSIXct(TenMinTS), f.Raw(TenMinTS) - f.Actual(TenMinTS) )) )
lines(DailyTideBias$Date, DailyTideBias$FMean, col="blue")

# plot of cumulative daily bias
plot( cumsum(DailyTideBias$FMean ) )

# lookup function
?approxfun
head(TenMinTS)
f.tidebias <- approxfun(DailyTideBias$Date, DailyTideBias$FMean, method="constant", f=0)

#####################################################
#
#f.CCInterp <- approxfun(ccInterpData$Date, ccInterpData$avg)

# Plot
plot(TenMinTS, f.CCInterp( TenMinTS ) - f.tidebias(TenMinTS)  )
lines(TenMinTS,  f.Upper( TenMinTS )- f.tidebias(TenMinTS))
lines(TenMinTS, f.Lower((TenMinTS) )- f.tidebias(TenMinTS))
lines(TenMinTS, f.Actual(TenMinTS), col="green")

# dataframe
biasRemoved <- data.frame( Date = TenMinTS,  CCIFiltered = f.CCInterp( TenMinTS ) - f.tidebias(TenMinTS)  )
biasRemoved <- cbind(biasRemoved,  Upper =  f.Upper( TenMinTS )- f.tidebias(TenMinTS))
biasRemoved <- cbind(biasRemoved,  Lower = f.Lower( (TenMinTS) )- f.tidebias(TenMinTS))
biasRemoved <- cbind(biasRemoved,  Actual = f.Actual( (TenMinTS) ) )

# trim edges
biasRemoved <- na.trim(biasRemoved)

# Plot before and after
plot(biasRemoved$Date, cumsum(biasRemoved$Actual), ylim = c(-150000, 1800000) )
lines(biasRemoved$Date, cumsum(biasRemoved$CCIFiltered + f.tidebias(biasRemoved$Date)), col="blue")
lines(biasRemoved$Date, cumsum(f.tidebias(biasRemoved$Date)), col="red")
lines(biasRemoved$Date, cumsum(biasRemoved$CCIFiltered), col="red")

# create timeseries of time actual flow above upper bounds
TimeAboveUpper <- (biasRemoved$Actual - biasRemoved$Upper )
TimeAboveUpper[TimeAboveUpper <= 0] <- 0
plot(TimeAboveUpper)

# create a timeseries of time actual flow below lower bounds
TimeBelowLower <- (biasRemoved$Actual - biasRemoved$Lower )
TimeBelowLower[TimeBelowLower >= 0] <- 0
lines(abs(TimeBelowLower), col="red" )

#
plot(biasRemoved$Date, cumsum(biasRemoved$Actual), type="l")
lines(biasRemoved$Date, cumsum(biasRemoved$Actual)+cumsum(TimeAboveUpper), col="blue" )
lines(biasRemoved$Date, cumsum(biasRemoved$Actual)-abs(cumsum(TimeBelowLower) ), col="red")

# Total cumulative absolute discharge, of periods above and below the filtered discharge
(sum(TimeAboveUpper)  + abs(sum(TimeBelowLower))) / sum(biasRemoved$Actual)

message(paste( (sum(TimeAboveUpper)  + abs(sum(TimeBelowLower))) / sum(biasRemoved$Actual) * 100, "%", sep="" ))
# 3.822%

sum(biasRemoved$Actual)

#sum(biasRemoved$Actual)



#WQI::EIO_Node(APIKEY = api)

