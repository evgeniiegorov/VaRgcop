#'  @description Derives P(X+Y<quantile) using Gaussian copula, normal margins. If X and Y are position P/Ls,
#'  then Value-at-Risk is equal to minus quantile.
#'  @param quantile portfolio P/L qunatile, number
#'  @param mu1 mean of first marginal P/L distribution, number
#'  @param mu2 mean of second marginal P/L distribution, number
#'  @param sigma1 standard deviation of first marginal P/L distribution, number, greater than zero
#'  @param sigma2 standard deviation of second marginal P/L distribution, number, greater than zero
#'  @param rho measure of correlation between two assets, number, must be in $[0;1)$
#'  @param nsteps the number of steps used in the copula approximation. 
#'  Increase in the number of steps increases the accuracy, but also increases the computation time.
#'  @title CDF of portfolio
cdfsum.copula.gauss <- function(quantile, mu1, mu2, sigma1, sigma2, rho, nsteps){
  
  if (!((sigma1 > 0) & (sigma2 > 0))){
    stop('Standard deviation must be greater than zero.') 
  } 
  
  if (!((rho >= 0) & (rho < 1))){
    stop('Correlation coefficient rho must be in [0;1)')
  }
  
  if (!(nsteps > 0)){
    stop('Number of steps in copula approximation nsteps must be greater than zero')
  }
  
  
  w.min <- 0.001
  w.max <- 0.999
  dw <- 0.001
  w <- seq(w.min, w.max, by = dw )
  y <- pnorm(quantile-qnorm(w,mu1,sigma1),mu2,sigma2)
  
  arg <- cbind(w,y)
  first.copula <- apply(arg, 1, function(x){
    copula.gauss(x[1],x[2],rho,nsteps)
  })
  
  arg <- cbind(w-dw,y)
  second.copula <- apply(arg[-1,], 1, function(x){
    copula.gauss(x[1],x[2],rho,nsteps)
  })
  
  approx.diff.copula <- (first.copula[-1]-second.copula)/dw
  
  return(sum(approx.diff.copula)*dw)
}