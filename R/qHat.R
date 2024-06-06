
#' Estimated theta quantile based on a given series
#'
#' @description Computing the estimated quantile according to a given series for a given level of theta
#' @param xt A vector representing a series.
#' @param theta A value from 0 to 1. Default is 0.95.
#'
#' @return A numeric number
#' @export
#' @import stats
#'
#'
#'
#'
qHat <- function(xt, theta = 0.95){
  if (requireNamespace(c('stats'), quietly = TRUE)){
    if ((theta < 0) | (theta > 1)){
      message("Parameter `theta` should be between 0 and 1.")
      return()
    }
    F_x <- ecdf(xt)
    return( min( xt[which(F_x(xt) >= theta)]) )
  }
}
