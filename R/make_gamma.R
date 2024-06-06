
#' Aggregation of estimated error
#'
#' @description Computing the aggregated estimated error at some index for constructing asymptotic normality
#' @param E A matrix representing the estimated error matrix at some time t of dim (p,q).
#' @param type Character input, choice from one of 'mu', 'alpha' and 'beta'. Default is 'mu'.
#' @param ind integer denoting the index of interest, only used when type is 'alpha' or 'beta'. Default is 1.
#'
#' @return A numeric number
#' @export
#'
#'
#'
#'
make_gamma <- function(E, type = 'mu', ind = 1){
  p <- dim(E)[1]
  q <- dim(E)[2]

  if (type == 'mu'){
    return(sqrt(sum( E^2) / (p*q) ))
  } else if (type == 'alpha'){
    if ((ind < 1)|(ind > p)){
      message("Parameter `ind` out of range.")
      return()
    }
    return(sqrt( sum(E[ind,]^2)/q ))
  } else if (type == 'beta'){
    if ((ind < 1)|(ind > q)){
      message("Parameter `ind` out of range.")
      return()
    }
    return(sqrt( sum(E[,ind]^2)/p ))
  } else {
    message("Parameter `type` should be `mu`, `alpha` or `beta`.")
    return()
  }
}
