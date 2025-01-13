changeInterval <- function(ts, dt = 1, Interval = "Daily", start = 0,
                           end = 0, offset = 0, option = "fmean",
                           rounded = TRUE)
{

  stopifnot("duplicate timestamps" = length(ts[,1]) == length(unique(ts[,1])) )

  names(ts)[1] <- "Date"
  names(ts)[2] <- "value"

  ts$Date <- as.POSIXct(ts$Date)
  inputts <- ts

  if (start == 0)
  {
    if (Interval == "Daily")
    {
      start <- round.POSIXt(ts[1, 1], units = "days")
      Interval = 60*60*24
    }else if ( Interval == "Hourly")
    {
      start <- round.POSIXt(ts[1, 1], units = "hours")
      Interval = 60*60
    }else if (rounded == TRUE){
      start <- round.POSIXt(ts[1, 1], units = "hours")
      Interval = Interval*60
    }else{
      start <- ts[1, 1]
      Interval = Interval*60
    }
  }else{
    # start should already be posixct time
  }

  if (end == 0) {
    end <- ts[nrow(ts), 1]
  } else {
    end <- as.POSIXct(end, origin = "1970-01-01")
  }

  # set to numeric
  ts <- ts %>% mutate( numDate = as.numeric(Date))

  #range(ts$Date) %>% round("hours")

  # new time sequence
  newintTS <- seq(start+offset*60, end , by = Interval)

  if(option == "resample")
  {
    return( data.frame(Date = newintTS,
                       Inst = approx(ts$Date, ts$value, newintTS)$y
    )
    )
  }

  if(dt == 1){
    # merge new timestamps into original dataset with linear interpolation
    f.timeties <- approxfun(ts$Date, ts$value)
    linearts <- data.frame(Date = newintTS, value = f.timeties(newintTS))
    merged <- rbind(linearts, dplyr::select(ts, c(Date, value)))
    merged <- merged[order(merged$Date), ]
    ts <- na.omit(distinct(merged))

    # convert rate to cumulative volume
    #ts <- ts %>% mutate( accum = cumsum(value * c(0, diff(numDate))))
    ts <- ts %>% mutate( accum = cumtrapz(Date %>% as.numeric, value) ) # this works

    # linear lookup accum
    f.accum <- approxfun(ts$Date, ts$accum)

    # new dataframe for output timestep
    newts <- data.frame(Date = newintTS,
                        accum = f.accum(newintTS))

  }else{
    stopifnot("other than dt of 1 must be interval data, i.e. daily forward mean" = max(diff(ts$numDate)) == mean(diff(ts$numDate)))
    stopifnot("select dt of 1(inst), 2(forward mean) or 3(trailing mean)" = max(diff(ts$numDate)) == mean(diff(ts$numDate)))

    # accumulate first
    ts <- ts %>% mutate( accum = cumsum(value * c(0,diff(numDate))))
    #ts <- ts %>% mutate( accum = cumtrapz(numDate, value) )

    if(dt == 2)
    {
      f.accum <- approxfun(ts$numDate + 1 * mean(diff(ts$numDate)), ts$accum)
    }else{
      f.accum <- approxfun(ts$numDate - 0 * mean(diff(ts$numDate)), ts$accum)
    }

    # new dataframe for output timestep
    newts <- data.frame(Date = newintTS,
                        accum = f.accum(newintTS))

  }

  newts <- na.omit(newts)

  # add half a timestep for instantaneous
  if(option == "inst"){


    #f.spline <- splinefun(newts$Date + 0.5 * mean(diff(newts$Date)), newts$accum  )
    #f.spline <- splinefun(newts$Date, newts$accum )
    #newts$Inst <- f.spline(newintTS, deriv = 1)
    #newts <- newts %>% mutate(splineInst = f.spline(Date, deriv = 1))
    #newts <- newts %>% dplyr::select(c(Date, Inst))
    newts$Date <- as.numeric(newts$Date)

    if(dt == 1)
    {
      # Convert Date to numeric if not already
      newts$Date <- as.numeric(newts$Date)
      newts$Inst <- inversecumtrapz(newts$Date, newts$accum, ts$value[1])

    }else{
      f.spline <- splinefun(newts$Date, newts$accum)
      newts$Inst <- f.spline(newintTS, deriv = 1)
      newts <- newts %>% mutate(splineInst = f.spline(Date, deriv = 1))
    }

    newts <- newts %>% dplyr::select(c(Date, Inst))

    ##########################
    # Print the result
    print(newts)


  }else if(option == "fmean")
  {
    newts$FMean <- c(diff(newts$accum)/Interval,NA )
    newts <- newts %>% dplyr::select(c(Date, FMean))
  }else if(option == "tmean")
  {
    newts$TMean <- c(NA,diff(newts$accum)/Interval )
    newts <- newts %>% dplyr::select(c(Date, TMean))
  }


  newts$Date <- as.POSIXct(newts$Date)

  if (length(inputts) == 3) {
    newts <- maxminfun(inputts[, 1], inputts[, 3], newts, option = "max")
  }
  newts[,2] <- round(newts[,2], 3)
  newts <- na.omit(newts)

  return(newts)

}


ts <- StevesCoolRandomTS(smoothed =  FALSE)
ts <- data.frame(ts$Time, ts$Signal)
Interval <- "Hourly"
option <- "inst"
start <- as.POSIXct("2023-06-01 00:00")
end <- 0
offset <- 0
rounded <- TRUE

hourly <- changeInterval(data.frame(ts$Time, ts$Signal), Interval = "Hourly", option = "inst")




library(pracma)

# Example: cumulative integral data
x <- seq(0, 10, length.out = 100)
y <- sin(x)  # Original function
cum_integral <- cumtrapz(x, y)  # Cumulative integral

# Reverse: approximate derivative
reverse_derivative <- diff(cum_integral) / diff(x)

# Add NA to align lengths (optional)
reverse_derivative <- c(NA, reverse_derivative)

# Compare with the original y
plot(x, y, type = "l", col = "blue", lwd = 2, ylim = c(-1.5, 1.5), main = "Original vs Reverse of cumtrapz")
lines(x, reverse_derivative, col = "red", lty = 2, lwd = 2)
legend("topright", legend = c("Original y", "Reverse cumtrapz"), col = c("blue", "red"), lty = c(1, 2))

#####

reverse_cumtrapz <- function(x, cum_integral) {
  n <- length(cum_integral)
  y <- numeric(n)

  # Assume the first derivative value y[1] is zero (or another known value)
  y[1] <- 0

  # Reconstruct y[i+1] from cum_integral
  for (i in 1:(n - 1)) {
    delta_x <- x[i + 1] - x[i]
    y[i + 1] <- 2 * (cum_integral[i + 1] - cum_integral[i]) / delta_x - y[i]
  }

  return(y)
}

# Example usage
library(pracma)

ts <- StevesCoolRandomTS(smoothed = FALSE, obs = 100)

x <- ts$Time %>% as.numeric()
y_original <- ts$Noise  # Original function
cum_integral <- cumtrapz(x, y_original)  # Cumulative integral

# Reconstruct the original function
y_reconstructed <- reverse_cumtrapz(x, cum_integral)

# Plot for comparison
plot(x, y_original, type = "l", col = "blue", lwd = 2, main = "Original vs Exact Reverse of cumtrapz")
lines(x, y_reconstructed, col = "red", lty = 2, lwd = 2)
legend("topright", legend = c("Original y", "Reconstructed y"), col = c("blue", "red"), lty = c(1, 2))






x = c(1,3.5,4,5,6)
y = c(1,2,5.545,6,7.5)
y2 = c(y[-1])


plot(x,y)
# area under curve

x = c(1,3,4,5,6)
ydx = c(0, diff(x)   *    ( y[-1]  +  y[-length(y)]  ) / 2 )
accum <- cumsum(ydx)
print(accum)
cumtrapz(x, y)




# now reverse it

# decumulate
decum <- c(0, diff(accum))

plot(x, decum)

recreated_y <- y[1]

( ( recreated_y + next_y ) / 2 )  *  diff(x)[1] == decum[2]


#(next_y) is actulaly 2
#next_y <- ( decum[2] * 2  ) - recreated_y

#something <- 2

#( diff(x)[1] * ( something + recreated_y ) ) / 2 == decum[2]
#decum[2] * 2 =  diff(x)[1] * ( something + recreated_y )
#( decum[2] * 2 ) / diff(x)[1] == something + recreated_y


###########################################
# this code works, so don't delete

decum <- c(0, diff(accum))

recreated_y <- y[1] # make recreated

i <- 1
recreated_y <- ( ( decum[i+1] * 2 ) / diff(x)[i] ) - recreated_y

i <- 2
recreated_y <- ( ( decum[i+1] * 2 ) / diff(x)[i] ) - recreated_y

i <- 3
recreated_y <- ( ( decum[i+1] * 2 ) / diff(x)[i] ) - recreated_y

i <- 4
recreated_y <- ( ( decum[i+1] * 2 ) / diff(x)[i] ) - recreated_y

#########################################
# Now put it in a loop

decum <- c(0, diff(accum))
recreated_y <- rep(0, length(decum))
recreated_y[1] <- 1
for(i in 1:(length(decum)-1) )
{
  recreated_y[i+1] <- ( ( decum[i+1] * 2 ) / diff(x)[i] ) - recreated_y[i]
}
y

##################


#########################

inversecumtrapz <- function( x, accum, startpos)
{
  # accum is the trapezoidal method of integration
  decum <- c(0, diff(accum))
  recreated_y <- rep(0, length(decum))
  recreated_y[1] <- startpos
  for(i in 1:(length(decum)-1) )
  {
    recreated_y[i+1] <- ( ( decum[i+1] * 2 ) / diff(x)[i] ) - recreated_y[i]
  }
  return(recreated_y)
}


# random x
x <- rnorm(50, mean = 10, sd = 5)
x[x < 0] <- 1
x <- cumsum(x)

# sequential x
#x <- seq(1, 50)

# random y
y <- rnorm(50, mean = 0, sd = 1)
y <- cumsum(y)

plot(x,y)
# area under curve
#x = c(1,3,4,5,6)

intervalx <- seq(range(x)[1], range(x)[2], by = 10 )
#plot(x,y)
ydx = c(0, diff(x)   *    ( y[-1]  +  y[-length(y)]  ) / 2 )
accum <- cumsum(ydx)
##print(accum)

intervalx <- c(x, intervalx)
intervalx <- unique(intervalx)
intervaly <- approx(x, y, intervalx)$y


df <- data.frame(intervalx, intervaly)
df <- df[order(df$intervalx),]

plot(df)


accum <- cumtrapz(df$intervalx, df$intervaly)

#ydx = c(0, diff(df$intervalx)   *    ( df$intervaly[-1]  +  df$intervaly[-length(df$intervaly)]  ) / 2 )
#accum <- cumsum(ydx)

plot(df$intervalx, accum)

intervalx <- seq(range(x)[1], range(x)[2], by = 10 )
intervalcumy <- approx(df$intervalx, accum, intervalx)$y

plot(intervalx, intervalcumy)

recreatedinterval <- inversecumtrapz(intervalx, intervalcumy,  df$intervaly[1])


plot(x, y)
points(intervalx, recreatedinterval, col = "blue")


intervalx <- seq(range(x)[1], range(x)[2], by = 10 )
recreatedintervaly <- approx(df$intervalx, recreatedinterval, intervalx )$y

plot(intervalx, recreatedintervaly)
lines(intervalx, recreatedintervaly)
points(x, y, col = "blue")

intervalcumy <- approx(x, accum, intervalx)$y

plot(x, accum)
points(intervalx, intervalcumy, col = "red")


plot(intervalx, intervalcumy)
recreatedinterval <- inversecumtrapz(intervalx, intervalcumy, y[1])
plot(intervalx, recreatedinterval)

plot(x, y)
lines(x,y)
#points(intervalx, recreatedinterval, col = "blue")
#intervalx
points(intervalx, intervaly, col = "blue")
points(intervalx, recreatedinterval, col = "red")
#points(intervalx ,c(0, rollmean(recreatedinterval, 2)), col = "orange")

print(intervaly)
print(recreatedinterval)


f.accum <- splinefun(intervalx, intervalcumy, method = "hyman")

lines(intervalx+5, f.accum(intervalx, deriv = 1), col = "red")

intervalcumy


#####################


accum <- intervalcumy
x <- intervalx

plot(x, accum)

###########################################
# this code works, so don't delete

decum <- c(0, diff(accum))

recreated_y <- y[1] # make recreated

i <- 1
recreated_y <- ( ( decum[i+1] * 2 ) / diff(x)[i] ) - recreated_y

i <- 2
recreated_y <- ( ( decum[i+1] * 2 ) / diff(x)[i] ) - recreated_y

i <- 3
recreated_y <- ( ( decum[i+1] * 2 ) / diff(x)[i] ) - recreated_y

i <- 4
recreated_y <- ( ( decum[i+1] * 2 ) / diff(x)[i] ) - recreated_y

#########################################



#########################################
# Now put it in a function
inversecumtrapz <- function( x, accum, startpos)
{
  # accum is the trapezoidal method of integration
  decum <- c(0, diff(accum))
  recreated_y <- rep(0, length(decum))
  recreated_y[1] <- startpos
  for(i in 1:(length(decum)-1) )
  {
    recreated_y[i+1] <- ( ( decum[i+1] * 2 ) / diff(x)[i] ) - recreated_y[i]
  }
  return(recreated_y)
}

################################################


randomts <- StevesCoolRandomTS(smoothed = FALSE)

# insert new interval before accumulating

newintTS <- seq(ceiling(randomts$Time[1]))

accum <- cumtrapz(randomts$Time %>% as.numeric, randomts$Noise)

#plot(randomts$Time, accum)
plot(randomts$Time, randomts$Noise)



recreated <- inversecumtrapz(randomts$Time %>% as.numeric, accum, randomts$Noise[1] )
points(randomts$Time, recreated, col = "red")


############################################


randomts <- StevesCoolRandomTS(smoothed = FALSE)

accum <- cumtrapz(randomts$Time %>% as.numeric, randomts$Signal)
tail(accum)
x <- randomts$Time %>% as.numeric
y <- randomts$Signal
y <- round(y, 1)
#plot(y)

ydx = c(0, diff(x)   *    ( y[-1]  +  y[-length(y)]  ) / 2 )
accum <- cumsum(ydx)
tail(accum)

plot(randomts$Time, accum)
plot(randomts$Time, randomts$Signal)


#head(randomts)


hourlyts <- seq(from = lubridate::ceiling_date(randomts$Time[1], unit = "hours"),  to = randomts$Time[nrow(randomts)], by = 60*60 )
#hourly <- approx(randomts$Time, accum, hourlyts) %>% as.data.frame
f.accum <- splinefun(randomts$Time, accum)
hourly <- data.frame(x = hourlyts, y = f.accum(hourlyts))


plot(randomts$Time, accum)
points(hourly, col = "red")

head(hourly)
tail(hourly)

#start value
approx(randomts$Time, randomts$Signal, hourly$x[1] )$y

recreated <- inversecumtrapz(hourly$x %>% as.numeric,
                             hourly$y,
                             approx(randomts$Time, randomts$Signal, hourly$x[1] )$y
                             )

hourly[10045,]$y

head(recreated)

plot(randomts$Time, randomts$Signal)
lines(hourly$x, recreated)


library(xts)
library(dygraphs)
cbind(
  xts(randomts$Signal, randomts$Time),
  xts(recreated, hourly$x)
) %>% dygraph

########################################
# This looks promising...
#
# duplicating above, but interpolating on the hour before accumulating


randomts <- StevesCoolRandomTS(smoothed = FALSE)
randomts <- randomts %>% dplyr::select(c(Time, Signal))

oldway <- changeInterval(randomts, option = "inst", Interval = "Hourly")

hourlyts <- seq(from = lubridate::ceiling_date(randomts$Time[1], unit = "hours"),  to = randomts$Time[nrow(randomts)], by = 60*60 )+0.5*60*60



timeties <- approx(randomts, xout = hourlyts)$y
randomts <- data.frame(Date = c(randomts$Time, hourlyts),
                       value = c(randomts$Signal, timeties))

#plot(randomts)
randomts <- randomts %>% distinct(Date, .keep_all = TRUE)
randomts <- randomts[order(randomts$Date),]

accum <- cumtrapz(randomts$Date %>% as.numeric, randomts$value)
tail(accum)
x <- randomts$Date %>% as.numeric
y <- randomts$value


#plot(y)


#ydx = c(0, diff(x)   *    ( y[-1]  +  y[-length(y)]  ) / 2 )
#accum <- cumsum(ydx)
#tail(accum)

plot(randomts$Date, accum)
plot(randomts$Date, randomts$value)


#hourlyts <- seq(from = lubridate::ceiling_date(randomts$Time[1], unit = "hours"),  to = randomts$Time[nrow(randomts)], by = 60*60 )
#hourly <- approx(randomts$Time, accum, hourlyts) %>% as.data.frame
f.accum <- approxfun(randomts$Date, accum)
hourly <- data.frame(x = hourlyts, y = f.accum(hourlyts))

plot(randomts$Date, accum)
points(hourly, col = "red")

head(hourly)
tail(hourly)

#start value
approx(randomts$Date, randomts$value, hourly$x[1] )$y

recreated <- inversecumtrapz(hourly$x %>% as.numeric,
                             hourly$y,
                             approx(randomts$Date, randomts$value, hourly$x[1] )$y
)

#hourly[10045,]$y

head(recreated)
tail(recreated)

plot(randomts$Date, randomts$value)
lines(hourly$x[-1], rollmean(recreated, 2), col = 'red')


library(xts)
library(dygraphs)
cbind(
  xts(randomts$value, randomts$Date),
  xts(rollmean(recreated, 2), (hourly$x[-1] - 0.5*60*60)),
  xts(oldway$Inst, oldway$Date)
) %>% dygraph


###
# comparing to old way


######################################################

#View(recreated)
#tail(recreated)




nrow(hourly)
length(c(0,hourly$x))
length(recreated)

length( rollmean(recreated, 2) )

cbind(
xts(randomts$Signal, randomts$Time),
xts(rollmean(recreated,2), hourly$x[-1]-(0.5*60*60) )
) %>% dygraph

############
# So close #
############


#plot(randomts$Time, randomts$Signal),
#points(hourly$x[-1], , col = "red" ))



#######################################


randomts <- StevesCoolRandomTS(smoothed = TRUE)

accum <- cumtrapz(randomts$Time %>% as.numeric, randomts$Signal)

plot(randomts$Time, accum)
plot(randomts$Time, randomts$Signal)


head(randomts)

accum <- round(accum,3)

hourlyts <- seq(from = round(randomts$Time[1], "hours"),  to = randomts$Time[nrow(randomts)], by = 60*60 )
#hourly <- approx(randomts$Time, accum, hourlyts) %>% as.data.frame
f.accum <- splinefun(randomts$Time, accum)
hourly <- data.frame(x = hourlyts, y = f.accum(hourlyts))


plot(randomts$Time, accum)
points(hourly, col = "red")

head(hourly)
tail(hourly)

f.hourly <- splinefun(hourly)
hourly$inst <- f.hourly(hourly$x, deriv = 1)
hourly$inst <- round(hourly$inst,3)


plot(randomts$Time, randomts$Signal)
points(hourly$x, hourly$inst, col= "red")


cbind(
  xts(hourly$inst, hourly$x),
  xts(randomts$Signal, randomts$Time)) %>% dygraph


















###################################


library(xts)
library(dygraphs)

cbind(
xts(recreated, hourly$x),
xts(randomts$Signal, randomts$Time)) %>% dygraph


###################################################################


hourlyinputs <- StevesCoolRandomTS(smoothed = TRUE)


hourlyinputs$accum <- cumtrapz(hourlyinputs$Time %>% as.numeric, hourlyinputs$Signal)
recreated <- inversecumtrapz(hourlyinputs$Time %>% as.numeric,
                             hourlyinputs$accum,
                             approx(hourlyinputs$Time, hourlyinputs$accum, hourlyinputs$Time[1])$y
)



plot(hourlyinputs$Time, recreated)


plot(randomts$Time, randomts$Signal)
points(hourly$x, recreated, col = "red")



head(randomts)
head(recreated)







#############################################

( recreated_y + something ) / 2 = 3
diff(x)[1]





plot(x[-length(x)], cumsum(ydx))

cumtrapz(x,y)

##############################################
# 2024-12-02
# gave up and tried again

ts <- StevesCoolRandomTS(smoothed =  FALSE)
ts <- data.frame(ts$Time, ts$Noise)
Interval <- "Hourly"
option <- "inst"
start <- 0
end <- 0
offset <- 0
rounded <- TRUE

changeInterval <- function(ts, dt = 1, Interval = "Daily", start = 0,
                           end = 0, offset = 0, option = "fmean",
                           rounded = TRUE)
{

  stopifnot("duplicate timestamps" = length(ts[,1]) == length(unique(ts[,1])) )

  names(ts)[1] <- "Date"
  names(ts)[2] <- "value"

  ts$Date <- as.POSIXct(ts$Date)
  inputts <- ts

  negoffset <- ts$value %>% min
  ts$value <- ts$value - negoffset

  if (start == 0)
  {
    if (Interval == "Daily")
    {
      start <- round.POSIXt(ts[1, 1], units = "days")
      Interval = 60*60*24
    }else if ( Interval == "Hourly")
    {
      start <- round.POSIXt(ts[1, 1], units = "hours")
      Interval = 60*60
    }else if (rounded == TRUE){
      start <- round.POSIXt(ts[1, 1], units = "hours")
      Interval = Interval*60
    }else{
      start <- ts[1, 1]
      Interval = Interval*60
    }
  }else{
    # start should already be posixct time
  }

  if (end == 0) {
    end <- ts[nrow(ts), 1]
  } else {
    end <- as.POSIXct(end, origin = "1970-01-01")
  }

  # set to numeric
  ts <- ts %>% mutate( numDate = as.numeric(Date))

  #range(ts$Date) %>% round("hours")

  # new time sequence
  newintTS <- seq(start+offset*60, end , by = Interval)

  if(dt == 1){
    # merge new timestamps into original dataset with linear interpolation
    f.timeties <- approxfun(ts$Date, ts$value)
    linearts <- data.frame(Date = newintTS, value = f.timeties(newintTS))
    merged <- rbind(linearts, dplyr::select(ts, c(Date, value)))
    merged <- merged[order(merged$Date), ]
    ts <- na.omit(distinct(merged))


    # convert rate to cumulative volume
    #ts <- ts %>% mutate( accum = cumsum(value * c(0, diff(numDate))))
    ts <- ts %>% mutate( accum = cumtrapz(Date %>% as.numeric, value) )

    # linear lookup accum
    f.accum <- approxfun(ts$Date, ts$accum)

    # new dataframe for output timestep
    newts <- data.frame(Date = newintTS,
                        accum = f.accum(newintTS))


  }else{
    stopifnot("other than dt of 1 must be interval data, i.e. daily forward mean" = max(diff(ts$numDate)) == mean(diff(ts$numDate)))
    stopifnot("select dt of 1(inst), 2(forward mean) or 3(trailing mean)" = max(diff(ts$numDate)) == mean(diff(ts$numDate)))

    # accumulate first
    ts <- ts %>% mutate( accum = cumsum(value * c(0,diff(numDate))))
    #ts <- ts %>% mutate( accum = cumtrapz(numDate, value) )

    if(dt == 2)
    {
      f.accum <- approxfun(ts$numDate + 1 * mean(diff(ts$numDate)), ts$accum)
    }else{
      f.accum <- approxfun(ts$numDate - 0 * mean(diff(ts$numDate)), ts$accum)
    }

    # new dataframe for output timestep
    newts <- data.frame(Date = newintTS,
                        accum = f.accum(newintTS))

  }

  # add half a timestep for instantaneous
  if(option == "inst"){
    #f.spline <- splinefun(newts$Date + 0.5 * mean(diff(newts$Date)), newts$accum  )

    f.spline <- splinefun(as.numeric ( newts$Date ), newts$accum, method = "hyman" )
    #newts$Inst <- f.spline(newintTS, deriv = 1)
    newts <- newts %>% mutate(Inst = f.spline( as.numeric( Date ), deriv = 1))
    newts <- newts %>% dplyr::select(c(Date, Inst))
    plot(newts)

  }else if(option == "fmean")
  {
    newts$FMean <- c(diff(newts$accum)/Interval,NA )
    newts <- newts %>% dplyr::select(c(Date, FMean))
  }else if(option == "tmean")
  {
    newts$TMean <- c(NA,diff(newts$accum)/Interval )
    newts <- newts %>% dplyr::select(c(Date, TMean))
  }
  newts$Date <- as.POSIXct(newts$Date)

  if (length(inputts) == 3) {
    newts <- maxminfun(inputts[, 1], inputts[, 3], newts, option = "max")
  }

  newts[,2] <- newts[,2] + negoffset
  newts[,2] <- round(newts[,2], 3)
  newts <- na.omit(newts)

  return(newts)

}

ts <- StevesCoolRandomTS(smoothed =  TRUE)
hourly <- changeInterval(data.frame(ts$Time, ts$Signal), Interval = "Hourly", option = "inst"     )
oldhourly <- ccInterp::changeInterval(data.frame(ts$Time, ts$Signal), Interval = "Hourly", option = "inst"     )



hourlyfmean <- changeInterval(data.frame(ts$Time, ts$Signal), Interval = "Hourly", option = "fmean", offset = 30)
hourlyfmean$Date <- hourlyfmean$Date + 0.5*(median(diff(hourlyfmean$Date)))

library(dygraphs)
library(xts)
cbind(
  xts(ts$Signal, ts$Time),
  xts(hourly$Inst, hourly$Date),
  xts(hourlyfmean$FMean, hourlyfmean$Date)
)%>% dygraph


cbind(
  xts(ts$Signal, ts$Time),
  xts(oldhourly$Inst, oldhourly$Date),
  xts(hourly$Inst, hourly$Date)
)%>% dygraph



