## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- warning=F, message=FALSE, eval=TRUE-------------------------------------
library(MEFM)

## -----------------------------------------------------------------------------
TT <- 60 #time length
d <- c(60,60) #spatial dimensions
r <- c(2,2) #rank of core tensor
re <- c(2,2) #rank of common component in error
eta <- list(c(0,0), c(0,0)) #strong factors
coef_f <- c(0.7, 0.3, -0.4, 0.2, -0.1) #AR(5) coefficients for core factor
coef_fe <- c(-0.7, -0.3, -0.4, 0.2, 0.1) #AR(5) coefficients for common component in error
coef_e <- c(0.8, 0.4, -0.4, 0.2, -0.1) #AR(5) coefficients for idiosyncratic part in error
param_mu <- c(0,1) #Normal mean and variance for mu
param_alpha <- c(0,1) #Normal mean and variance for alpha
param_beta <- c(0,1) #Normal mean and variance for beta
data_example <- gen_MEFM(TT,d,r,re,eta, coef_f, coef_fe, coef_e, param_mu, param_alpha, param_beta)

## -----------------------------------------------------------------------------
est_result <- est_MEFM(data_example$MEFM)

paste0('Estimated number of factors: ', est_result$r[1], ', ', est_result$r[2])
paste0('Estimation error for mu: ', sum((est_result$mu - data_example$mu)^2)/TT)
paste0('Estimation error for alpha: ', sum((est_result$alpha - data_example$alpha)^2)/prod(dim(data_example$alpha)))
paste0('Estimation error for beta: ', sum((est_result$beta - data_example$beta)^2)/prod(dim(data_example$beta)))
paste0('Estimation error for row factor loading: ', tensorMiss::fle(est_result$A[[1]], data_example$A[[1]]))
paste0('Estimation error for column factor loading: ', tensorMiss::fle(est_result$A[[2]], data_example$A[[2]]))
paste0('Estimation error for Yt: ', sum((est_result$Yt - (data_example$MEFM - data_example$E) )^2)/prod(dim(est_result$Yt)))

## -----------------------------------------------------------------------------
time_select <- 10 #select a time point
hat.E <- (est_result$Yt - data_example$MEFM)[time_select,,]
gamma_mu <- make_gamma(hat.E, type = 'mu')
asymp_i <- prod(d)^0.5 * gamma_mu^(-1) * (est_result$mu[time_select] - (-1.15))

paste0('The true mu at time 10 is: ', data_example$mu[time_select])
paste0('The test statistic is: ', asymp_i )

## -----------------------------------------------------------------------------
gamma_alpha_inv <- matrix(0, nrow=3, ncol=3)
for (j in 1:3){ gamma_alpha_inv[j,j] <- (make_gamma(hat.E, type = 'alpha', j))^(-1) }
asymp_alpha_example <- (d[2])^0.5 * gamma_alpha_inv %*% (est_result$alpha[time_select, 1:3] - data_example$alpha[time_select, 1:3])

gamma_beta_inv <- matrix(0, nrow=3, ncol=3)
for (j in 1:3){ gamma_beta_inv[j,j] <- (make_gamma(hat.E, type = 'beta', j))^(-1) }
asymp_beta_example <- (d[1])^0.5 * gamma_beta_inv %*% (est_result$beta[time_select, 1:3] - data_example$beta[time_select, 1:3])

paste0('The test statistic using true alpha: ', asymp_alpha_example[1], ', ', asymp_alpha_example[2], ', ', asymp_alpha_example[3])
paste0('The test statistic using true beta: ', asymp_beta_example[1], ', ', asymp_beta_example[2], ', ', asymp_beta_example[3])

## -----------------------------------------------------------------------------
# computing the covariance matrix estimator
r2 <- r[2]
eta <- floor(0.2*((TT * prod(d))^0.25))
D2 <- diag(x=(svd(est_result$covMatrix[[2]])$d)[1:r2], nrow=r2, ncol=r2)
# HAC_cov: HAC-type covariance matrix estimator
HAC_cov <- sigmaD_MEFM(2, D2, est_result$A[[2]], est_result$Ct, data_example$MEFM - est_result$Yt, 1, eta)

# computing the rotation matrix
A2 <- data_example$A[[2]]
Amk_i <- data_example$A[[1]]
R_ast_i <- 0
for (tt in 1:TT){
  R_ast_i <- R_ast_i + tensorMiss::unfold(matrix(data_example$Ft[tt,,], ncol=r2),2) %*% t(Amk_i) %*% Amk_i %*% t(tensorMiss::unfold(matrix(data_example$Ft[tt,,], ncol=r2),2))
}
R_ast_i <- A2 %*% R_ast_i %*% t(A2)
R_ast_i <- R_ast_i/TT
Z2_i <- diag(x = diag(t(A2) %*% A2), nrow=r2, ncol=r2)
Q2_i <- A2 %*% diag(x=diag(solve(Z2_i))^0.5, nrow=r2, ncol=r2)
# H2_i: rotation matrix
H2_i <- solve(D2) %*% t(est_result$A[[2]]) %*% R_ast_i %*% Q2_i %*% solve(t(Q2_i)%*% Q2_i)

## -----------------------------------------------------------------------------
HAC_cov.eigen <- eigen(HAC_cov)
HAC_cov.sqrt <- HAC_cov.eigen$vectors %*% diag(sqrt(HAC_cov.eigen$values)) %*% solve(HAC_cov.eigen$vectors)
A2_1 <- TT * (solve(HAC_cov.sqrt) %*% D2) %*% (matrix(est_result$A[[2]], nrow=d[2], ncol=r2)[1,] - (H2_i %*% Q2_i[1,]))
A2_1

## -----------------------------------------------------------------------------
paste0('MSE of estimating MEFM on FM-generated data: ', sum((est_MEFM(data_example$FM)$Yt -  data_example$FM)^2) / sum(data_example$FM^2))
paste0('MSE of estimating FM on FM-generated data: ', sum((est_FM(data_example$FM)$Ct -  data_example$FM)^2) / sum(data_example$FM^2))

## -----------------------------------------------------------------------------
paste0('MSE of estimating MEFM on MEFM-generated data: ', sum((est_MEFM(data_example$MEFM)$Yt -  data_example$MEFM)^2) / sum(data_example$MEFM^2))
paste0('MSE of estimating FM on MEFM-generated data (using number of factors as (3,3)): ', sum((est_FM(data_example$MEFM)$Ct -  data_example$MEFM)^2) / sum(data_example$MEFM^2))
paste0('MSE of estimating FM on MEFM-generated data (using number of factors as (30,30)): ', sum((est_FM(data_example$MEFM, r=c(30,30))$Ct -  data_example$MEFM)^2) / sum(data_example$MEFM^2))

## -----------------------------------------------------------------------------
MEFM_on_MEFM <- est_MEFM(data_example$MEFM)
FM_on_MEFM <- est_FM(data_example$MEFM, r=c(MEFM_on_MEFM$r[1] +1, MEFM_on_MEFM$r[1] +1))

x_alpha <- make_xy(MEFM_on_MEFM$Yt - data_example$MEFM, type = 'alpha')
y_alpha <- make_xy(FM_on_MEFM$Ct - data_example$MEFM, type = 'alpha')
x_beta <- make_xy(MEFM_on_MEFM$Yt - data_example$MEFM, type = 'beta')
y_beta <- make_xy(FM_on_MEFM$Ct - data_example$MEFM, type = 'beta')

paste0('Computed alpha_reject: ', sum(y_alpha >= qHat(x_alpha, 0.95))/length(y_alpha))
paste0('Computed beta_reject: ', sum(y_beta >= qHat(x_beta, 0.95))/length(y_beta))

