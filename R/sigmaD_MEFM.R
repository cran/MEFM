
#' HAC covariance estimator for asymptotic normality on each row j of loading matrix estimator
#'
#' @description Computing the HAC covariance estimator for asymptotic normality on each row j of the row or column loading matrix estimator
#' @param k Integer to choose the mode of loading matrix, either 1 or 2.
#' @param D Eigenvalue matrix of sample covariance matrix, with dimension rk by rk.
#' @param Q Estimated row (k=1) or column (k=2) loading matrix, with dimension p (for k=1) or q (for k=2) by rk.
#' @param C Estimated common component series, written in an array with dimension (T,p,q) where the first dimension denotes time.
#' @param E Estimated error matrix time series, written in an array with the same dimension as C.
#' @param j Integer representing the row of loading matrix. Value should be integers from minimum 1 to maximum p (for k=1) or q (for k=2).
#' @param beta Lag parameter of the HAC type. Default is 0.
#'
#' @return A matrix of dimension rk by rk
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
#' est_result = est_MEFM(data_example$MEFM, r=r);
#' D2 <- diag(x=(svd(est_result$covMatrix[[2]])$d)[1:r[2]], nrow=r[2], ncol=r[2]);
#' sigmaD_MEFM(2, D2, est_result$A[[2]], est_result$Ct, data_example$MEFM - est_result$Yt, 1, 0);
#'
#'
sigmaD_MEFM <- function(k, D, Q, C, E, j, beta=0){

  if (requireNamespace('tensorMiss', quietly = TRUE)){

    TT <- dim(C)[1]
    rk <- dim(D)[1]
    p <- dim(C)[2]
    q <- dim(C)[3]

    if (length(dim(C)) != length(dim(E))){
      message('The number of dimensions for C and E are inconsistent.')
      return()
    } else if ( prod(dim(C)) != prod(dim(E)) ){
      message('The dimensions of C and E are inconsistent.')
      return()
    }

    if ( (k < 1)|(k > 2) ){
      message('The parameter is invalid with the given matrix series.')
      return()
    }
    if ( (beta < 0)|(beta > (TT-1)) ){
      message('Parameter beta is either negative or too large.')
      return()
    }

    #construct DQCC and CE first;
    if (k==1){
      CC <- unfold(C,2) %*% t(unfold(C,2))
      DQCC <- solve(D) %*% t(Q) %*% CC / TT
      CE_j <- matrix(0, nrow=p, ncol=TT)
      for (t in 1:TT){
        CE_j[,t] <- C[t,,] %*% E[t,j,]
      }
    } else {
      CC <- unfold(C,3) %*% t(unfold(C,3))
      DQCC <- solve(D) %*% t(Q) %*% CC / TT
      CE_j <- matrix(0, nrow=q, ncol=TT)
      for (t in 1:TT){
        CE_j[,t] <- t(C[t,,]) %*% E[t,,j]
      }
    }


    #D_0 -----------------------------------------
    D_nu <- matrix(0, nrow = rk, ncol = rk)
    #estimation block ----------------------------
    if (k==1){
      for (t in 1:TT){
        first_term <- matrix(0, nrow=rk, ncol=1)
        for (i in 1:p){
          first_term <- first_term + (CE_j[i,t] * DQCC[,i])
        }
        D_nu <- D_nu + (first_term %*% t(first_term))
      }
    } else {
      for (t in 1:TT){
        first_term <- matrix(0, nrow=rk, ncol=1)
        for (i in 1:q){
          first_term <- first_term + (CE_j[i,t] * DQCC[,i])
        }
        D_nu <- D_nu + (first_term %*% t(first_term))
      }
    }#block ends ---------------------------------
    Sigma_HAC <- t(D_nu) + D_nu


    #D_nu -----------------------------------------
    if (beta > 0){
      for (nu in 1:beta){
        D_nu <- matrix(0, nrow = rk, ncol = rk)
        #estimation block ----------------------------
        if (k==1){
          for (t in (1+nu):TT){
            first_term <- matrix(0, nrow=rk, ncol=1)
            for (i in 1:p){
              first_term <- first_term + (CE_j[i,t] * DQCC[,i])
            }
            second_term <- matrix(0, nrow=rk, ncol=1)
            for (i in 1:p){
              second_term <- second_term + (CE_j[i,t-nu] * DQCC[,i])
            }
            D_nu <- D_nu + (first_term %*% t(second_term))
          }
        } else {
          for (t in (1+nu):TT){
            first_term <- matrix(0, nrow=rk, ncol=1)
            for (i in 1:q){
              first_term <- first_term + (CE_j[i,t] * DQCC[,i])
            }
            second_term <- matrix(0, nrow=rk, ncol=1)
            for (i in 1:q){
              second_term <- second_term + (CE_j[i,t-nu] * DQCC[,i])
            }
            D_nu <- D_nu + (first_term %*% t(second_term))
          }
        }#block ends ---------------------------------
        D_nu <- (1- (nu/(1+beta))) * (t(D_nu) + D_nu)
        Sigma_HAC <- Sigma_HAC + D_nu
      }
    }
    return(Sigma_HAC)
  }
}
