\name{el.CS.prob}
\alias{el.CS.prob}
\title{
Current Status Data Empirical Likellihood Test for the Probability \eqn{F(t_0)}.
}
\description{
Given n current status data, we may estimate the CDF F(t) by NPMLE (e.g. by \code{isotNEW2()} function in this package).
This function, \code{el.CS.prob}, uses empirical likelihood to test the hypothesis that F(t) at a given location (t0)
equal to a given value (Ft0): i.e. H0: F(t0) = Ft0.

Empirical likelihood ratio test returns the Wilks statistics, -2LLR.
According to our working conjecture, 
the -2 log likelihood ratio times (5/3) under H0 is approximately chi square DF=1 distributed.
See references below.
}
\usage{ 
el.CS.prob(ti, di, t0=0.5, Ft0=0.5) 
}
\arguments{
\item{ti}{The inspection times, a vector of length n.}
\item{di}{Either 0 or 1. I[yi <= ti]. length n.}
\item{t0}{The given time where F() value is tested.}
\item{Ft0}{The hypothesized value of F(t0). Must be within (0, 1).}
}
\details{
This function tests the null hypothesis that F(t0) = Ft0 versus not equal.
We assume the data given is current status censored data.

We require t0 equals to one of the inspection times. If not, you have
to do something by the right continuity of the NPMLE (change t0 to the closest ti on the left).

The NPMLE F(t) is convergent at cubic root speed and the -2LLR times (5/3) has chi square DF=1
null distribution, according to our working conjecture.

It goes without saying that we assume the NPMLE has finite asymptotic variance (when normalized by cubic root n).
}
\value{
It returns a list containing
\item{"-2LLR"}{The Wilks statistics of the EL test, when multiply by (5/3) has approximate
                              chi SQ DF=1 distribution under null hypothesis.}
\item{LogLik0}{The log lik value achieved by the un-constrained NPMLE.}
\item{LogLik1}{The log lik value achieved by the constrained NPMLE.}
}
\references{
Banerjee, M., and Wellner, J. (2001). Likelihood ratio tests for monotone
functions. Annals of Statistics 29, 1699-1731.

Banerjee, M., and Wellner, J. (2005). Score statistics for current sta-
tus data: Comparisons with likelihood ratio and Wald statistics. The
International Journal of Biostatistics 1, 1-29.

Zhou, M. (2026). \emph{Empirical Likelihood Method in Survival Analysis 2nd Edition} Chapman & Hall/CRC

Sun, J. (2006). \emph{The Statistical Analysis of Interval-Censored Failure Time Data} Springer, New York.
}
\author{
Mai Zhou <maizhou@gmail.com>.
}
\examples{
N <- 300
set.seed(12345)
itime <- sort(c(rexp(N-1), 0.5) )       #### inspection times      
Stime <- rexp(N)             #### survival times
delta <- as.numeric(Stime <= itime)    ####  current status censoring

el.CS.prob( ti=itime, di=delta, t0=0.5, Ft0=pexp(0.5) ) 

#### You should get 
## $`-2LLR`
## [1] 1.867655     #### and more.
}
\keyword{nonparametric}
