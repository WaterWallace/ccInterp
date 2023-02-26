#' Combined Cumulative Interpolation Filter
#'
#' Filters a source timeseries to remove noise/tides/etc
#' Output data as Hourly.
#'
#'
#' @docType package
#'
#' @author Stephen Wallace \email{stephen.wallace@des.qld.gov.au}
#'
#' @name ccInterp
NULL

#' @import dplyr
#' @import zoo
#' @import pracma
#' @import schumaker
#' @import matrixStats
#' @import signal
#' @import stats
#' @import lubridate
NULL

#' Xy data
#'
#' An data set used in the examples
#'
#' A data set of 457 rows and 2 columns
#' \describe{
#'   \item{x}{A numeric vector}
#'   \item{y}{A numeric vector}
#' }
"xy_data"

#' Xy2 data
#'
#' An data set used in the examples
#'
#' A data set of 457 rows and 2 columns
#' \describe{
#'   \item{x}{A numeric vector}
#'   \item{y}{A numeric vector}
#' }
"xy2_data"
