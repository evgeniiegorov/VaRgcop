#'  @description Derives Value-at-Risk using bivariate gaussian copula, for normal margins
#'  @param quantile portfolio P/L qunatile, number
#'  @param mu1 mean of first marginal P/L distribution, number
#'  @param mu2 mean of second marginal P/L distribution, number
#'  @param sigma1 standard deviation of first marginal P/L distribution, number, greater than zero
#'  @param sigma2 standard deviation of second marginal P/L distribution, number, greater than zero
#'  @param rho measure of correlation between two assets, number, must be in $[0;1)$
#'  @param nsteps the number of steps used in the copula approximation. 
#'  Increase in the number of steps increases the accuracy, but also increases the computation time.
#'  @param cl Value-at-Risk confidence level
#'  @param bound number of portfolio sd from portfolio mean to search VaR by bisection, 3 by default
#'  @param tol tolerance level to stop search VaR by bisection, 0.001 by default
#'  @return list of two numbers: VaR estimated by copula and  VaR estimated by variance-covariance
#'  @title Value-at-Risk with guass copula
var.copula.gauss <- function(mu1, mu2, sigma1, sigma2, rho, nsteps, cl,
                             bound = 3, tol = 0.001){
  if (!((sigma1 > 0) & (sigma2 > 0))){
    stop('Standard deviation must be greater than zero.') 
  } 
  
  if (!((rho >= 0) & (rho < 1))){
    stop('Correlation coefficient rho must be in [0;1)')
  }
  
  if (!(nsteps > 0)){
    stop('Number of steps in copula approximation nsteps must be greater than zero')
  }
  
  if (!((cl > 0) & (cl < 1))){
    stop('Confidence level cl must be in (0;1)')
  }
  
  p <- 1-cl
  
  p.mu <- mu1+mu2
  cov.m <- matrix(c(sigma1^2, rho*sigma1*sigma2, rho*sigma1*sigma2, sigma2^2),
                  nrow = 2, ncol = 2)
  p.var <- as.numeric(c(1,1) %*% cov.m %*% c(1,1))
  VaR.varcovar <- -p.mu-sqrt(p.var)*qnorm(p, mean = 0, sd = 1)
  
  L <- -p.mu-bound*sqrt(p.var)
  fL <- cdfsum.copula.gauss(L, mu1, mu2, sigma1, sigma2, rho, nsteps) - p
  
  U <- -p.mu+bound*sqrt(p.var)
  fU <- cdfsum.copula.gauss(U, mu1, mu2, sigma1, sigma2, rho, nsteps) - p
  
  if (sign(fL) == sign(fU)){
    stop('Assumed bounds are not include answer')
  }

  while(U-L > tol){
    x <- (L+U)/2
    cum.prob <- cdfsum.copula.gauss(x, mu1, mu2, sigma1, sigma2, rho, nsteps)
    fx <- cum.prob-p
    if (sign(fx) == sign(fL)){
      L <- x
      fL <- fx
    } else {
      U <- x
      fU <- fx
    }
  }
  result <- list(copulaVaR = -x, variancecovarianceVaR = VaR.varcovar)
  return(result)
}