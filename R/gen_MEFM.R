
#' Data generation of matrix time series with MEFM structure
#'
#' @description Generate a matrix time series with MEFM at each time t, with the first mode as the time mode, the second as the row mode and the third as the column mode
#' @param TT Length of time series.
#' @param d Dimensions of the matrix at time t, written in a vector of length 2 where the first number denotes the number of rows p and the second denoted the number of columns q.
#' @param r Rank of the core factors, written in a vector of length 2.
#' @param re re: Rank of the cross-sectional common error core factors, written in a vector of length 2.
#' @param eta Quantities controlling factor strengths in each factor loading matrix, written in a list of 2 vectors.
#' @param coef_f AR(5) coefficients for the factor series, written in a vector of length 5.
#' @param coef_fe AR(5) coefficients for the common component in error series, written in a vector of length 5.
#' @param coef_e AR(5) coefficients for the idiosyncratic component in error series, written in a vector of length 5.
#' @param param_mu If rademacher = TRUE, represent parameters of normal distribution to generate grand mean series mu_t, written in a vector of length 2 representing the mean and standard deviation. Otherwise written in a scalar multiplied by the generated Rademacher random variable.
#' @param param_alpha If rademacher = TRUE, represent parameters of normal distribution to generate row effect series alpha_t, written in a vector of length 2 representing the mean and standard deviation. Otherwise written in a scalar multiplied by the generated Rademacher random variable.
#' @param param_beta If rademacher = TRUE, represent parameters of normal distribution to generate column effect series beta_t, written in a vector of length 2 representing the mean and standard deviation. Otherwise written in a scalar multiplied by the generated Rademacher random variable.
#' @param heavy_tailed Whether to generate data from heavy-tailed distribution. If FALSE, generate from N(0,1); if TRUE, generate from t-distribution. Default is FALSE.
#' @param t_df The degree of freedom for t-distribution if heavy_tailed = TRUE. Default is 3.
#' @param rademacher Mechanism to generate mu, alpha and beta. If FALSE, generate from normal with param_mu, param_alpha and param_beta; if TRUE, generate from Rademacher distribution and scaled by param_mu, param_alpha and param_beta. Default is FALSE.
#' @param seed Random seed required for reproducibility. Default is 2024.
#'
#'
#' @return A list containing the following:
#' mu: the generated time-varying grand mean series, as a vector of length TT;
#' alpha: the generated time-varying row effect series, as a matrix of dimension (TT,p);
#' beta: the generated time-varying column effect series, as a matrix of dimension (TT,q);
#' A: a list of 2 factor loading matrices;
#' C: the generated common component time series, as multi-dimensional array with dimension 3, where mode-1 is the time mode, mode-2 is for rows and mode-3 is for columns;
#' Ft: the generated core factor series, as multi-dimensional array with dimension 3, where mode-1 is the time mode, mode-2 is for rows and mode-3 is for columns;
#' MEFM: the generated matrix time series with MEFM structure, as multi-dimensional array with dimension 3, where mode-1 is the time mode, mode-2 is for rows and mode-3 is for columns;
#' FM: the generated matrix time series with only traditional factor structure, as multi-dimensional array with dimension 3, where mode-1 is the time mode, mode-2 is for rows and mode-3 is for columns;
#' E: the generated error time series with factor structure, as multi-dimensional array with dimension 3, where mode-1 is the time mode, mode-2 is for rows and mode-3 is for columns;
#'
#'
#' @export
#' @import tensorMiss
#' @import stats
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
#' gen_MEFM(TT,d,r,re,eta, coef_f, coef_fe, coef_e, param_mu, param_alpha, param_beta);
#'
#'
#'
gen_MEFM <- function(TT,d,r,re,eta, coef_f, coef_fe, coef_e, param_mu, param_alpha, param_beta, heavy_tailed = FALSE, t_df = 3, rademacher = FALSE, seed = 2024){

  if (requireNamespace(c('tensorMiss', 'stats'), quietly = TRUE)){

    dt <- tensor_gen(2,TT,d,r,re,eta, coef_f, coef_fe, coef_e, heavy_tailed, t_df, seed)
    p <- d[1]
    q <- d[2]

    set.seed(seed)

    #generate mu_t
    if (rademacher){
      if (param_mu[1] == 0){ # faster computation
        mu_t <- rep(0, TT)
      } else {
        mu_t <- (2* rbinom(TT, 1, 0.5) - 1) * (param_mu[1])
      }
    } else {
      mu_t <- rnorm(TT, mean= param_mu[1], sd= param_mu[2])
    }

    #generate alpha_t
    if (rademacher){
      if (param_alpha[1] == 0){ # faster computation
        alpha_t <- array(0, dim=c(TT,p))
      } else {
        alpha_t <- array(((2* rbinom(p*TT, 1, 0.5) - 1) * (param_alpha[1])), dim=c(TT,p))
      }
    } else {
      alpha_t <- array(rnorm(p*TT, mean= param_alpha[1], sd= param_alpha[2]) , dim=c(TT,p))
    }
    for (t in 1:TT){
      alpha_t[t,] <- alpha_t[t,] - (sum(alpha_t[t,])/p)
    }

    #generate beta_t
    if (rademacher){
      if (param_alpha[1] == 0){ # faster computation
        beta_t <- array(0, dim=c(TT,q))
      } else {
        beta_t <- array(((2* rbinom(q*TT, 1, 0.5) - 1) * (param_beta[1])), dim=c(TT,q))
      }
    } else {
      beta_t <- array(rnorm(q*TT, mean= param_beta[1], sd= param_beta[2]) , dim=c(TT,q))
    }
    for (t in 1:TT){
      beta_t[t,] <- beta_t[t,] - (sum(beta_t[t,])/q)
    }

    #generate Ft
    Ft <- dt$Ft

    #generate A1, A2
    A1 <- dt$A[[1]]
    A2 <- dt$A[[2]]
    for (i in 1:(r[1])){
      A1[,i] <- A1[,i] - (sum(A1[,i])/p)
    }
    for (i in 1:(r[2])){
      A2[,i] <- A2[,i] - (sum(A2[,i])/q)
    }

    #generate common component Ct
    Ct <- ttm(Ft, A1, 2)
    Ct <- ttm(Ct, A2, 3)

    #generate Et
    Et <- dt$X - dt$C

    #generate FM_Yt
    FM_Yt <- array(0, dim=c(TT,p,q))
    for (t in 1:TT){
      FM_Yt[t,,] <- Ct[t,,] + Et[t,,]
    }

    #generate MEFM_Yt
    MEFM_Yt <- array(0, dim=c(TT,p,q))
    for (t in 1:TT){
      MEFM_Yt[t,,] <- (mu_t[t]*matrix(1, nrow=p,ncol=q)) +
        alpha_t[t,] %*% t(rep(1,q)) +
        rep(1,p) %*% t(beta_t[t,]) +
        FM_Yt[t,,]
    }

    return(list(mu = mu_t, alpha = alpha_t, beta = beta_t, A = list(A1, A2), C = Ct,
                Ft = Ft, MEFM = MEFM_Yt, FM = FM_Yt, E = Et))
  }
}
