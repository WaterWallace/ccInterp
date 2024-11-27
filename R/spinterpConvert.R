#' Spinterp Convert
#'
#' Simply gets the derivative of the cumulative discharge
#' not appropriate for downsampling intervals from hourly back to daily
#' Use changeInterval() for that.
#' (even though theoretically spinterpConvert will work fine for downsampling,
#'  and in some cases technically preferred)
#' Using "spinterp" Interpolation over long gaps could have unpredicatable results.
#'
#' @param start timestamp at start of the forward mean, time in days
#' inputs must be stepped averages (default, dt = 2)
#' If point data is used (not preferred as precision can be lost), set dt = 1
#' @param rate rate of the parameter, i.e. distance/time
#' @param outputInt output interval, time in days, i.e. hourly = 1/24
#' @param type "spinterp", "spline"
#' @param dt datatrans, default 2 = forward mean
#'
#' @return  dataframe with time in days and an interpolated rate.
#'
#' @examples
#'
#'
#' # daily times and daily means
#' data(xy2_data)
#' D <- xy2_data
#'
#' f.D <- approxfun(D, method = "constant")
#'
#' plot(D, ylim = c(0, 450))
#' ts <- seq(D$x[1], D$x[length(D$x)], by = (1 / 24))
#' lines(ts, f.D(ts), ylim = c(0, 450))
#' # upsample to hourly
#' # default = spinterp
#' DHourly <- spinterpConvert(D$x, D$y, type = "spinterp")
#' lines(DHourly, col = "red")
#'
#' # type = spline
#' DHourly <- spinterpConvert(D$x, D$y, type = "spline")
#' lines(DHourly, col = "blue")
#'
#' # type = linear
#' # DHourly <- spinterpConvert(D$x, D$y, type = "linear")
#' # lines(DHourly, col="green")
#' @export
spinterpConvert <- function(start, rate, outputInt = (1 / 24), type = "spinterp", dt = 2) {
  if (dt == 1) # if data is of point type, convert to an interval mean at from start point to end point
  {
    # start <- df$days
    # rate <- df$df...2.
    means <- data.frame(start, rate)
    means$end <- c(means$start[-1], 0)
    means$rate2 <- c(means$rate[-1], 0)
    means <- head(means, -1)
    tail(means)

    # get average
    means$avg <- (means$rate + means$rate2) / 2
    start <- means$start
    rate <- means$avg
  }
  if (dt == 3) # mean, average flow until this point
  {
    start <- start - (start[2] - start[1])
  }

  start <- as.numeric(start)
  nt <- length(start)

  t <- c(start, start[nt] +  median(diff(start)))

    # duration between input data points # in days
  dur <- c(diff(t), 0) # t[2]-t[1]  #c(0,diff(t))

  # create a time sequence for output data
  xp <- seq(t[1], t[nt+1], by = outputInt)

  negativeOffset <- min(rate)
  rate <- rate - negativeOffset

  # cumulative sum of discharge
  # from "cubic metres per second" into "cubic metres"
  DAYSTOSECONDS <- 60 * 60 * 24
  cumdaily <- cumsum(c(0, rate * dur[-length(dur)] )) # * DAYSTOSECONDS)

  # output at new time interval
  if (type == "spinterp") {
    # print("spinterp")
    #plot(t, cumdaily)
    tfactor <- median(diff(t)) # transform dates so not a huge gap
    #tfactor <- 1
    yp <- spinterp(t/tfactor, cumdaily, xp/tfactor)
    #plot(xp, yp)
  } else if (type == "spline") {
    # print("spline")
    f.yp <- splinefun(t, cumdaily, method = "hyman")
    yp <- f.yp(xp)
    #plot(yp)
  } else if (type == "schum") {
    # print("schum")
    SchumSpline <- schumaker::Schumaker(t, cumdaily)
    yp <- SchumSpline$Spline(xp)
  } else if (type == "linear") {
    # Not useful to interpolate linearly
    f.yp <- approxfun(t, cumdaily)
    yp <- f.yp(xp)
  }


  xpyp <- data.frame(xp, yp) %>% na.omit

  # derivative of spinterp
  f.spinterp <- splinefun(xpyp, method = "monoH.FC")
  df <- data.frame(Date = xpyp$xp, Data = f.spinterp(xpyp$xp, deriv = 1) + negativeOffset) #  /DAYSTOSECONDS  )


  return(df)
}

