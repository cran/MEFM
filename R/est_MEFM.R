
#' Estimation of MEFM on matrix time series
#'
#' @description Estimate the MEFM structure on the given matrix time series
#' @param Yt demeaned matrix time series, written in an array with dimension 3 and the first dimension for time.
#' @param r Rank of core factors for the common component, written in a vector of length 2. First value as 0 is to denote unknown rank which would be automatically estimated using ratio-based estimators. Default is 0.
#' @param delta Non-negative number as the correction parameter for rank estimation. Default is 0.2.
#'
#'
#' @return A list containing the following:
#' r: a vector representing either the given rank or the estimated rank, with length 2;
#' mu: a vector representing the estimated time-varying grand mean series;
#' alpha: a matrix representing the estimated time-varying row effect series, where the row index denotes time index;
#' beta: a matrix representing the estimated time-varying column effect series, where the row index denotes time index;
#' A: a list of the estimated row and column factor loading matrices;
#' Ft: the estimated core factor series, as multi-dimensional array with dimension 3, where mode-1 is the time mode;
#' Ct: the estimated common component time series, as multi-dimensional array with dimension 3, where mode-1 is the time mode;
#' Yt: the estimated matrix time series, as multi-dimensional array with dimension 3, where mode-1 is the time mode;
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
#' est_MEFM(data_example$MEFM);
#'
#'
est_MEFM <- function(Yt, r=0, delta=0.2){

  if (requireNamespace('tensorMiss', quietly = TRUE)){
    TT <- dim(Yt)[1]
    p <- dim(Yt)[2]
    q <- dim(Yt)[3]

    #estimate mu_t
    hat.mu_t <- rep(0, TT)
    for (t in 1:TT){
      hat.mu_t[t] <- t(rep(1,p)) %*% Yt[t,,] %*% rep(1,q) /(p*q)
    }

    #estimate alpha_t
    hat.alpha_t <- array(0 , dim=c(TT,p))
    for (t in 1:TT){
      hat.alpha_t[t,] <- (Yt[t,,] %*% rep(1,q) /q) - (hat.mu_t[t] * rep(1,p))
    }

    #estimate beta_t
    hat.beta_t <- array(0 , dim=c(TT,q))
    for (t in 1:TT){
      hat.beta_t[t,] <- ( t(Yt[t,,]) %*% rep(1,p) /p) - (hat.mu_t[t] * rep(1,q))
    }

    #estimate L_t
    hat.Lt <- array(0, dim=c(TT,p,q))
    for (t in 1:TT){
      hat.Lt[t,,] <- Yt[t,,] - (hat.mu_t[t]*matrix(1, nrow=p,ncol=q)) -
        hat.alpha_t[t,] %*% t(rep(1,q)) -
        rep(1,p) %*% t(hat.beta_t[t,])
    }

    cov_mat <- list()
    for (k in 2:3){
      cov_mat_k <- unfold(hat.Lt, k) %*% t(unfold(hat.Lt, k))
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
    hat.Ft <- ttm(hat.Lt, t(hat.A1), 2)
    hat.Ft <- ttm(hat.Ft, t(hat.A2), 3)


    #construct Ct
    hat.Ct <- ttm(hat.Ft, hat.A1, 2)
    hat.Ct <- ttm(hat.Ct, hat.A2, 3)

    #construct Yt
    hat.Yt <- array(0, dim=c(TT,p,q))
    for (t in 1:TT){
      hat.Yt[t,,] <- (hat.mu_t[t]*matrix(1, nrow=p, ncol=q)) +
        hat.alpha_t[t,] %*% t(rep(1,q)) +
        rep(1,p) %*% t(hat.beta_t[t,]) + hat.Ct[t,,]
    }

    return(list(r=c(r1,r2), mu=hat.mu_t, alpha=hat.alpha_t, beta=hat.beta_t, A=list(hat.A1,hat.A2),
                Ft=hat.Ft, Ct=hat.Ct, Yt=hat.Yt, covMatrix=cov_mat))
  }
}
