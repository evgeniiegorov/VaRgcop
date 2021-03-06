% Generated by roxygen2 (4.0.2): do not edit by hand
\name{var.copula.gauss}
\alias{var.copula.gauss}
\title{Value-at-Risk with guass copula}
\usage{
var.copula.gauss(mu1, mu2, sigma1, sigma2, rho, nsteps, cl, bound = 3,
  tol = 0.001)
}
\arguments{
\item{mu1}{mean of first marginal P/L distribution, number}

\item{mu2}{mean of second marginal P/L distribution, number}

\item{sigma1}{standard deviation of first marginal P/L distribution, number, greater than zero}

\item{sigma2}{standard deviation of second marginal P/L distribution, number, greater than zero}

\item{rho}{measure of correlation between two assets, number, must be in $[0;1)$}

\item{nsteps}{the number of steps used in the copula approximation.
Increase in the number of steps increases the accuracy, but also increases the computation time.}

\item{cl}{Value-at-Risk confidence level}

\item{bound}{number of portfolio sd from portfolio mean to search VaR by bisection, 3 by default}

\item{tol}{tolerance level to stop search VaR by bisection, 0.001 by default}

\item{quantile}{portfolio P/L qunatile, number}
}
\value{
list of two numbers: VaR estimated by copula and  VaR estimated by variance-covariance
}
\description{
Derives Value-at-Risk using bivariate gaussian copula, for normal margins
}
\details{

}

