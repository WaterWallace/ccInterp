#' Combined Cumulative Interpolation Filter
#'
#' Derivative of overlapping trapezoidally integrated datasets
#'
#' Useful for detiding or denoising of data
#' The overlapping averaging windows allow multiple recalculations of the same data
#' Which provides an estimate of the uncertainty in applying the filter
#'
#' @param x vector of equally spaced interval data
#'
#' @return  dataframe with timestamp as posixct a trapezoidal integrated rate and quality code if included.
#'
#' @examples
#' data <- cumsum(rnorm(1000))
#' data <- data * c(0, data[-length(data)])
#' plot(data)
#' lines(godinFilter(data), col="red")
#' @export
godinFilter <- function(x)
{
  filt1 <- stats::filter(x, rep(1/24,24), method="convolution", sides=2)
  filt1 <- c(NA, filt1)
  filt2 <- stats::filter(filt1, rep(1/24,24), method="convolution", sides=2)
  filt3 <- stats::filter(filt2, rep(1/25,25), method="convolution", sides=2)
  filtered <- filt3[-length(filt3)]

  return (filtered)
}
