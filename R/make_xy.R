
#' Construction of series for testing MEFM
#'
#' @description Constructing x or y series for the MEFM testing
#' @param E An array representing the sequence of estimated error matrix of dim (T,p,q).
#' @param type Character input, either 'alpha' or 'beta'. Default is 'alpha'.
#'
#' @return A vector representing the constructed x or y series
#' @export
#'
#'
#'
#'
make_xy <- function(E, type = 'alpha'){
  e_series <- rep(0, dim(E)[1])

  for (t in 1:(dim(E)[1]) ){
    E_mat <- E[t,,]
    if (type == 'alpha'){
      e_series[t] <- max(diag( E_mat %*% t(E_mat) )) / ncol(E_mat)
    } else if (type == 'beta'){
      e_series[t] <- max(diag( t(E_mat) %*% E_mat )) / nrow(E_mat)
    } else {
      message("Parameter `type` should be `alpha` or `beta`.")
      return()
    }
  }

  return(e_series)
}
