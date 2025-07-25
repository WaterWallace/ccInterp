% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maxminfun.R
\name{maxminfun}
\alias{maxminfun}
\title{MaxMin Function}
\usage{
maxminfun(refx, refy, newts, option = "max", dt = 1)
}
\arguments{
\item{refx}{source x input.}

\item{refy}{source y input.}

\item{newts}{target x interval i.e. hourly/daily.}

\item{option}{"max" or "min" to return interval max or interval min.}

\item{dt}{data type, 1 = inst, 2 = forward mean, 6 = insttotal}
}
\value{
dataframe with newts xaxis, and a forward looking max or min of refy.
}
\description{
function for looking up the maximum y-value in a reference dataset, for a time at time of newts
}
\details{
newts will be a forward looking mean, i.e. 10:00, 11:00, 12:00.. 10:00 will be the mean between 10:00 and 11:00
And the function result will be the maximum in the reference period between 10:00 and 11:00
This is useful for quality coded data that has been averaged.
the maximum quality code can be extracted from the original data and applied to the averaged data.
also useful for a daily maximum value
or a daily minimum value
}
\examples{

longts <- seq(0, 500, by = 1)
# create a random timeseries
randts <- runif(length(longts), min = 0, max = 1000)

plot(longts, randts)
lines(longts, randts, col = "grey")

# new interval timeseries
intervalts <- seq(0, 500, by = 20)

intervalMax <- maxminfun(longts, randts, intervalts)
f.max <- approxfun(intervalMax, method = "constant")
lines(f.max(longts), col = "red")

intervalMin <- maxminfun(longts, randts, intervalts, option = "min")
f.min <- approxfun(intervalMin, method = "constant")
lines(f.min(longts), col = "red")

}
