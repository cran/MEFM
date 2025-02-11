% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_MEFM.R
\name{gen_MEFM}
\alias{gen_MEFM}
\title{Data generation of matrix time series with MEFM structure}
\usage{
gen_MEFM(
  TT,
  d,
  r,
  re,
  eta,
  coef_f,
  coef_fe,
  coef_e,
  param_mu,
  param_alpha,
  param_beta,
  heavy_tailed = FALSE,
  t_df = 3,
  rademacher = FALSE,
  seed = 2024
)
}
\arguments{
\item{TT}{Length of time series.}

\item{d}{Dimensions of the matrix at time t, written in a vector of length 2 where the first number denotes the number of rows p and the second denoted the number of columns q.}

\item{r}{Rank of the core factors, written in a vector of length 2.}

\item{re}{re: Rank of the cross-sectional common error core factors, written in a vector of length 2.}

\item{eta}{Quantities controlling factor strengths in each factor loading matrix, written in a list of 2 vectors.}

\item{coef_f}{AR(5) coefficients for the factor series, written in a vector of length 5.}

\item{coef_fe}{AR(5) coefficients for the common component in error series, written in a vector of length 5.}

\item{coef_e}{AR(5) coefficients for the idiosyncratic component in error series, written in a vector of length 5.}

\item{param_mu}{If rademacher = TRUE, represent parameters of normal distribution to generate grand mean series mu_t, written in a vector of length 2 representing the mean and standard deviation. Otherwise written in a scalar multiplied by the generated Rademacher random variable.}

\item{param_alpha}{If rademacher = TRUE, represent parameters of normal distribution to generate row effect series alpha_t, written in a vector of length 2 representing the mean and standard deviation. Otherwise written in a scalar multiplied by the generated Rademacher random variable.}

\item{param_beta}{If rademacher = TRUE, represent parameters of normal distribution to generate column effect series beta_t, written in a vector of length 2 representing the mean and standard deviation. Otherwise written in a scalar multiplied by the generated Rademacher random variable.}

\item{heavy_tailed}{Whether to generate data from heavy-tailed distribution. If FALSE, generate from N(0,1); if TRUE, generate from t-distribution. Default is FALSE.}

\item{t_df}{The degree of freedom for t-distribution if heavy_tailed = TRUE. Default is 3.}

\item{rademacher}{Mechanism to generate mu, alpha and beta. If FALSE, generate from normal with param_mu, param_alpha and param_beta; if TRUE, generate from Rademacher distribution and scaled by param_mu, param_alpha and param_beta. Default is FALSE.}

\item{seed}{Random seed required for reproducibility. Default is 2024.}
}
\value{
A list containing the following:
mu: the generated time-varying grand mean series, as a vector of length TT;
alpha: the generated time-varying row effect series, as a matrix of dimension (TT,p);
beta: the generated time-varying column effect series, as a matrix of dimension (TT,q);
A: a list of 2 factor loading matrices;
C: the generated common component time series, as multi-dimensional array with dimension 3, where mode-1 is the time mode, mode-2 is for rows and mode-3 is for columns;
Ft: the generated core factor series, as multi-dimensional array with dimension 3, where mode-1 is the time mode, mode-2 is for rows and mode-3 is for columns;
MEFM: the generated matrix time series with MEFM structure, as multi-dimensional array with dimension 3, where mode-1 is the time mode, mode-2 is for rows and mode-3 is for columns;
FM: the generated matrix time series with only traditional factor structure, as multi-dimensional array with dimension 3, where mode-1 is the time mode, mode-2 is for rows and mode-3 is for columns;
E: the generated error time series with factor structure, as multi-dimensional array with dimension 3, where mode-1 is the time mode, mode-2 is for rows and mode-3 is for columns;
}
\description{
Generate a matrix time series with MEFM at each time t, with the first mode as the time mode, the second as the row mode and the third as the column mode
}
\examples{
TT = 40;
d = c(40,40);
r = c(2,2);
re = c(2,2);
eta = list(c(0,0), c(0,0));
coef_f = c(0.7, 0.3, -0.4, 0.2, -0.1);
coef_fe = c(-0.7, -0.3, -0.4, 0.2, 0.1);
coef_e = c(0.8, 0.4, -0.4, 0.2, -0.1);
param_mu = c(0,1);
param_alpha = c(0,1);
param_beta = c(0,1);
gen_MEFM(TT,d,r,re,eta, coef_f, coef_fe, coef_e, param_mu, param_alpha, param_beta);



}
