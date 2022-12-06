#' Spinterp Convert
#'
#' Simply gets the derivative of the cumulative discharge
#' not appropriate for downsampling intervals from hourly back to daily
#' Use changeInterval() for that.
#' (even though theoretically spinterpConvert will work fine for downsampling, and in some cases technically preferred)
#' Using "spinterp" Interpolation over long gaps could have unpredicatable results.
#'
#' @param start timestamp at start of the forward mean, time in days
#' inputs must be stepped averages (default, dt = 2)
#' If point data is used (not preferred as precision can be lost), set dt = 1
#' @param rate rate of the parameter, i.e. distance/time
#' @param ouputInt output interval, time in days, i.e. hourly = 1/24
#' @param type "spinterp", "spline"
#' @param dt datatrans, default 2 = forward mean
#'
#' @return  dataframe with time in days and an interpolated rate.
#'
#' @examples
#'
#'
#'# daily times and daily means
#'D <- list( x=c(16832.96,16833.96,16834.96,16835.96,16836.96,16837.96,16838.96,16839.96,16840.96,16841.96,16842.96,16843.96,16844.96,16845.96,16846.96,16847.96,16848.96,16849.96,16850.96,16851.96,16852.96,16853.96,16854.96,16855.96,16856.96,16857.96,16858.96,16859.96,16860.96,16861.96,16862.96,16863.96,16864.96,16865.96,16866.96,16867.96,16868.96,16869.96,16870.96,16871.96,16872.96,16873.96,16874.96,16875.96,16876.96,16877.96,16878.96,16879.96,16880.96,16881.96,16882.96,16883.96,16884.96,16885.96,16886.96,16887.96,16888.96,16889.96,16890.96,16891.96),
#'           y=c(  0,2.707,1.3795,211.932,22.565,4.409,2.536,1.162,0.66,0.358,0.291,0.1065,0.099,0.0265,0.0055,0,0,0,0,2.195,11.208,3.4125,1.3465,0.657,0.3255,0.2915,0.1235,0.0985,0.058,3.629,111.4225,148.32,144.8815,99.8055,45.8355,43.858,13.0915,5.876,3.6715,1.5845,0.8965,13.5315,2.902,1.608,0.7085,0.428,0.1965,0.2055,0.07,0.057,0.0235,0.007,0.006,0.001,0,0,0,0,0,0)
#')
#'
#'f.D <- approxfun(D, method="constant")
#'
#'plot(D, ylim=c(0,450))
#'ts <- seq (first(D$x),last(D$x), by = (1/24 ) )
#'lines(ts, f.D(ts), ylim=c(0,450) )
#'# upsample to hourly
#'# default = spinterp
#'DHourly <- spinterpConvert(D$x, D$y)
#'lines(DHourly, col="blue")
#'
#'# type = spline
#'DHourly <- spinterpConvert(D$x, D$y, type = "spline")
#'lines(DHourly, col="blue")
#'
#'# type = linear
#'# DHourly <- spinterpConvert(D$x, D$y, type = "linear")
#'# lines(DHourly, col="green")
#' @export
spinterpConvert <- function(start, rate, outputInt=(1/24), type="spinterp", dt=2)
{

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



}
