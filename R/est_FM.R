
#' Estimation of factor models on matrix time series
#'
#' @description Estimate the FM structure on the given matrix time series
#' @param Yt demeaned matrix time series, written in an array with dimension 3 and the first dimension for time.
#' @param r Rank of core factors for the common component, written in a vector of length 2. First value as 0 is to denote unknown rank which would be automatically estimated using ratio-based estimators. Default is 0.
#' @param delta Non-negative number as the correction parameter for rank estimation. Default is 0.2.
#'
#'
#' @return A list containing the following:
#' r: a vector representing either the given rank or the estimated rank, with length 2;
#' A: a list of the estimated row and column factor loading matrices;
#' Ft: the estimated core factor series, as multi-dimensional array with dimension 3, where mode-1 is the time mode;
#' Ct: the estimated common component time series, as multi-dimensional array with dimension 3, where mode-1 is the time mode;
#' covMatrix: a list of the estimated row and column covariance matrices which are used to estimate loading matrices;
#'
#'
#' @export
#' @import tensorMiss
#'
#' @examples
#' TT = 40;
#' d = c(40,40);
#' r = c(2,2);
#' re = c(2,2);
#' eta = list(c(0,0), c(0,0));
#' coef_f = c(0.7, 0.3, -0.4, 0.2, -0.1);
#' coef_fe = c(-0.7, -0.3, -0.4, 0.2, 0.1);
#' coef_e = c(0.8, 0.4, -0.4, 0.2, -0.1);
#' param_mu = c(0,1);
#' param_alpha = c(0,1);
#' param_beta = c(0,1);
#' data_example = gen_MEFM(TT,d,r,re,eta, coef_f, coef_fe, coef_e, param_mu, param_alpha, param_beta);
#' est_FM(data_example$FM);
#'
#'
est_FM <- function(Yt, r=0, delta=0.2){

  if (requireNamespace('tensorMiss', quietly = TRUE)){

    TT <- dim(Yt)[1]
    p <- dim(Yt)[2]
    q <- dim(Yt)[3]

    cov_mat <- list()
    for (k in 2:3){
      cov_mat_k <- unfold(Yt, k) %*% t(unfold(Yt, k))
      cov_mat[[k-1]] <- cov_mat_k/TT
    }

    #rank estimation if necessary
    if ( (r[1])==0 ){
      xi_1 <- delta * p*q*((TT*q)^(-0.5) + p^(-0.5))
      e1_val <- svd( cov_mat[[1]] )$d[1:floor(p/2)] + xi_1
      r1 <- which.min( ((e1_val[-1])/(e1_val[-length(e1_val)])) )

      xi_2 <- delta * p*q*((TT*p)^(-0.5) + q^(-0.5))
      e2_val <- svd( cov_mat[[2]] )$d[1:floor(q/2)] + xi_2
      r2 <- which.min( ((e2_val[-1])/(e2_val[-length(e2_val)])) )
    } else {
      r1 <- r[1]
      r2 <- r[2]
    } # rank estimation end here

    #estimate A1 and A2
    hat.A1 <- matrix(svd( cov_mat[[1]] )$u[, (1:r1)], nrow=p, ncol=r1)
    hat.A2 <- matrix(svd( cov_mat[[2]] )$u[, (1:r2)], nrow=q, ncol=r2)

    #estimate Ft
    hat.Ft <- ttm(Yt, t(hat.A1), 2)
    hat.Ft <- ttm(hat.Ft, t(hat.A2), 3)

    #construct Ct
    hat.Ct <- ttm(hat.Ft, hat.A1, 2)
    hat.Ct <- ttm(hat.Ct, hat.A2, 3)

    return(list(r=c(r1,r2), A=list(hat.A1,hat.A2), Ft=hat.Ft, Ct=hat.Ct, covMatrix=cov_mat))
  }
}
