\name{el.CS.Hz}
\alias{el.CS.Hz}
\title{
Current Status Data Empirical Likellihood Test for a Parameter
Defined by Weighted Cumulative Hazard
}
\description{
Given n current status data, we may estimate the CDF F(t) by NPMLE (e.g. by isotNEW2 in this package).
By the 1 to 1 correspondence, \eqn{ \Lambda (t)= - \log [1-F(t)]}, 
we can also calculate the NPMLE of the cumulative hazard function \eqn{ \Lambda(t)}.
We are interested in a parameter defined by weighted average of cumulative hazard function
(or integrated cumulative hazard) defined below.
Empirical likelihood ratio test is performed and return a Wilks statistics.
The Wilks statistics, -2 log likelihood ratio, under H0 is approximately chi square DF=1 distributed.
}
\usage{ 
el.CS.Hz(ti, di, Pfun, thetaMU, error=1e-11, maxit=25) 
}
\arguments{
\item{ti}{The inspection times, a vector of length n.}
\item{di}{Either 0 or 1. I[yi <= ti]. length n.}
\item{Pfun}{The function used in calculate the weighted cumulative hazard. See an example below.}
\item{thetaMU}{The null hypothesis value of the weighted cumulative hazard to be tested.}
\item{error}{Default to 1e-11. Control the SQP iteration.}
\item{maxit}{Default to 25. Control the SQP iteration.}
}
\details{
This function tests the null hypothesis that the parameter of weighted (or integrated)
cumulative hazard function equal to the given value (thetaMU).
We assume the data given is current status censored data.

The (unknown) true value of the parameter of interest is defined by
\deqn{
    \theta(\Lambda_0) = \int_0^M \Lambda_0(s) d \Psi(s) 
}
where \eqn{\Lambda_0(s)} is the true
cumulative hazard function. 
And \eqn{\Psi(s)} is the weight function given(Pfun, used to generate weights).
The function \eqn{\Psi (s)} is assumed to be continuous and at least piecewise differtiable.

The NPMLE of cumulative hazard, \eqn{\hat \Lambda_n(t)}, is obtained via the NPMLE of CDF \eqn{\hat F_n(t)}:
\deqn{
\hat \Lambda_n(t) = - \log [ 1- \hat F_n(t)] ~.
}

It goes without saying that we assume the parameter is well defined (not equal to infinity)
and its NPMLE has finite asymptotic variance.
}
\value{
It returns a list containing
\item{"-2LLR"}{The Wilks statistics of the EL test, has approximate chi SQ DF=1 distribution under null hypothesis.}
\item{location}{The locations of jumps of the NPMLE of cumulative hazard function (also of CDF).}
\item{Haz}{The constrained NPMLE of cumulative hazard function.}
\item{iter}{Number of iterations done.}
\item{error}{The error at final iteration, = sum(abs(difference))}
\item{Loglik1}{The log lik value achieved by the constrained NPMLE.}
\item{Check}{The weighted cumulative hazard by the constrained NPMLE. 
             This should equal to thetaMU (of input). A check of convergence.}
\item{thetaMLE}{The NPMLE of the (non-constrained) weighted cumulative hazard. 
                The input value thetaMU should not deviate from this value by too much.}
}
\references{
Zhou, M. (2026). \emph{Empirical Likelihood Method in Survival Analysis 2nd Edition} Chapman & Hall/CRC
}
\author{
Mai Zhou <maizhou@gmail.com>.
}
\examples{
#### An example of the Pfun (or Psi(s) or weight function) and calculation ####
mydgfun2 <- function(t){0.3*t*(10-t)*as.numeric(t<=10)}

N <- 3000
set.seed(12345)
itime <- sort(rexp(N, rate=0.1))       #### inspection times      
Stime <- rexp(N, rate=0.1)             #### survival times
delta <- as.numeric(Stime <= itime)    ####  current status censoring

el.CS.Hz(ti=itime, di=delta, Pfun=mydgfun2, thetaMU= -5) ## -5 is the true value of parameter.

#### You should get 
## $`-2LLR`
## [1] 1.04782          #### and more.
}
\keyword{nonparametric}
