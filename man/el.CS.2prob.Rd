\name{el.CS.2prob}
\alias{el.CS.2prob}
\title{
Current Status Data Empirical Likellihood Test for Two Parameters: \eqn{F(t_1)} and \eqn{F(t_2)} Jointly
}
\description{
Given n current status data, we may estimate the CDF F(t) by NPMLE
(e.g. by isotNEW2() function in this package).
This function, \code{el.CS.2prob}, uses empirical likelihood to test the hypothesis
that F(t) at two given locations(t01 and t02)
equal to two given values (Ft01 and Ft02): i.e. H0: F(t01) = Ft01 and F(t02) = Ft02 jointly.

Empirical likelihood ratio test returns the Wilks statistics, -2LLR.
By our working conjecture, 
the -2 log likelihood ratio times (5/3) under H0 is approximately chi square DF=2 distributed.
See reference below.
}
\usage{ 
el.CS.2prob(ti, di, t01, Ft01, t02, Ft02) 
}
\arguments{
\item{ti}{The inspection times, a vector of length n.}
\item{di}{Either 0 or 1. I[yi <= ti]. length n.}
\item{t01}{The given time where F() value is tested.}
\item{Ft01}{The hypothesized value of F(t01). Must be within (0, 1).}
\item{t02}{The given time where F() value is tested.}
\item{Ft02}{The hypothesized value of F(t02). Must be within (0, 1).}
}
\details{
This function tests the null hypothesis that F(t01) = Ft01 and F(t02) = Ft02 versus at least one not equal.
We assume the data given is current status censored data.

We require t01 (and also t02) be equal to one of the inspection times. If not, you have
to do something using the right continuity of the NPMLE (change t01 to the closest ti on the left, etc.).

The NPMLE F(t) is convergent at cubic root n speed and the -2LLR times (5/3) has chi square DF=2
null distribution, according to the working conjecture.

It goes without saying that we assume the NPMLE has finite asymptotic variance (when normalized by cubic root n).
}
\value{
It returns a list containing
\item{"-2LLR"}{The Wilks statistics of the EL test, after multiply by (5/3) 
                has approximate chi SQ DF=2 distribution under null hypothesis.}
\item{LogLik0}{The log lik value achieved by the un-constrained NPMLE.}
\item{LogLik1}{The log lik value achieved by the constrained NPMLE.}
}
\references{
Banerjee, M., and Wellner, J. (2001). Likelihood ratio tests for monotone
functions. Annals of Statistics 29, 1699-1731.

Banerjee, M., and Wellner, J. (2005). Score statistics for current sta-
tus data: Comparisons with likelihood ratio and Wald statistics. The
International Journal of Biostatistics 1, 1-29.

Zhou, M. (2026). \emph{Empirical Likelihood Method in Survival Analysis 2nd Edition} Chapman & Hall/CRC.

Sun, J. (2006). \emph{The Statistical Analysis of Interval-Censored Failure Time Data} Springer, New York.
}
\author{
Mai Zhou <maizhou@gmail.com>.
}
\examples{
N <- 300
set.seed(12345)
itime <- sort(c(rexp(N-2), 0.3, 0.6) )       #### inspection times      
Stime <- rexp(N)                         #### survival times
delta <- as.numeric(Stime <= itime)      ####  current status censoring
}
\keyword{nonparametric}
