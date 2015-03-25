#'  @description Derives value of gaussian copula
#'  @param u value of frist marginal for random variable X: $P(X<x)=u$, number, $u in (0;1)$
#'  @param v value of second marginal for random variable Y: $P(Y<v)=u$, number, $v in (0;1)$
#'  @param rho correlation between X and Y, $rho in [0;1)$
#'  @param nsteps number of steps used for each x and y
#'  @return value of copula in point $(u,v)$, number
#'  @title Gauss copula 
copula.gauss <- function(u,v,rho,nsteps){
    
  if (!((u > 0) & (u < 1))){
    stop('First marginal Prob(X<x) u must be in (0;1)')
  }
  
  if (!((v > 0) & (v < 1))){
    stop('Second marginal Prob(Y<y) u must be in (0;1)')
  }
  
  if (!((rho >= 0) & (rho < 1))){
    stop('Correlation coefficient rho must be in [0;1)')
  }
  
  if (!(nsteps > 0)){
    stop('Number of steps in copula approximation nsteps must be greater than zero')
  }
  
  x.min <- -3
  x.max <- qnorm(u, mean = 0, sd = 1)
  dx <- (x.max-x.min)/nsteps
  x <- seq(x.min, x.max, by = dx)
  
  y.min <- -3
  y.max <- qnorm(v, mean = 0, sd = 1)
  dy <- (y.max-y.min)/nsteps
  y <- seq(y.min, y.max, by = dy)
 
  term <- matrix(, nrow = length(x), ncol = length(y))
  for(i in 1:length(x)){
    for(j in 1:length(y)){
      term[i,j] <- exp(-1*(x[i]^2-2*rho*x[i]*y[j]+y[j]^2)/(2*(1-rho^2)))
    }
  }
  triangle <- function(x,d){
    0.5*d*(2*sum(x)-x[1]-x[length(x)])
  }
  sum.part <- apply(term,1,triangle, d = dy)
  sum.total <- triangle(sum.part, d = dx)
  return(sum.total/(2*pi*sqrt(1-rho^2)))
}