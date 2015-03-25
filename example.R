#please install package from source 
#link to donwload 
#https://github.com/evgeniiegorov/VaRgcop/blob/master/VaRgcop_1.0.tar.gz
library('VaRgcop')

#cdf gaussian copula 
copula.gauss(u = 0.1, v = 0.2, rho = 0.3, nsteps = 10) 

#cdf P/L of portfolio returns with two assets 
cdfsum.copula.gauss(quantile = 2, mu1 = 0, mu2 = 1, sigma1 = 1, sigma2 = 3, rho = 0.9, nsteps = 10)

#Value-at-Risk of two assets
var.copula.gauss(mu1 = 0, mu2 = 0, sigma1 = 1, sigma2 = 1, rho = 0, nsteps = 10, cl = 0.95)

#manipulate with level of correlation and see how cdf gaussian copula will change
install.packages('manipulate')
library('manipulate')
cop.rho <- function(rho){
  u <- v <- seq(from = 0.01, to = 0.99, by = 0.01)
  cop.g <- function(x,y,rho){
    arg <- cbind(x,y)
    apply(arg, 1, function(x){
      copula.gauss(u = x[1], v = x[2], rho, nsteps = 10 )
    })  
  }
  z <- outer(u, v, cop.g, rho)
  persp(u, v, z, col="lightblue", expand = 0.5,shade = 0.2) 
}
manipulate(cop.rho(rhoScale), rhoScale = manipulate::slider(0,0.9,step = 0.1))
  
           

