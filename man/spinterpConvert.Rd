% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spinterpConvert.R
\name{spinterpConvert}
\alias{spinterpConvert}
\title{Spinterp Convert}
\usage{
spinterpConvert(start, rate, outputInt = (1/24), type = "spinterp", dt = 2)
}
\arguments{
\item{start}{timestamp at start of the forward mean, time in days
inputs must be stepped averages (default, dt = 2)
If point data is used (not preferred as precision can be lost), set dt = 1}

\item{rate}{rate of the parameter, i.e. distance/time}

\item{outputInt}{output interval, time in days, i.e. hourly = 1/24}

\item{type}{"spinterp", "spline"}

\item{dt}{datatrans, default 2 = forward mean}
}
\value{
dataframe with time in days and an interpolated rate.
}
\description{
Simply gets the derivative of the cumulative discharge
not appropriate for downsampling intervals from hourly back to daily
Use changeInterval() for that.
(even though theoretically spinterpConvert will work fine for downsampling,
 and in some cases technically preferred)
Using "spinterp" Interpolation over long gaps could have unpredicatable results.
}
\examples{


# daily times and daily means
data(xy2_data)
D <- xy2_data

f.D <- approxfun(D, method = "constant")

plot(D, ylim = c(0, 450))
ts <- seq(D$x[1], D$x[length(D$x)], by = (1 / 24))
lines(ts, f.D(ts), ylim = c(0, 450))
# upsample to hourly
# default = spinterp
DHourly <- spinterpConvert(D$x, D$y, type = "spinterp")
lines(DHourly, col = "red")

# type = spline
DHourly <- spinterpConvert(D$x, D$y, type = "spline")
lines(DHourly, col = "blue")

# type = linear
# DHourly <- spinterpConvert(D$x, D$y, type = "linear")
# lines(DHourly, col="green")
}
