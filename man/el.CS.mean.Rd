\name{el.CS.mean}
\alias{el.CS.mean}
\title{
Current Status Data Empirical Likellihood Test for the Parameter of Mean: mu(F)
}
\description{
Given n current status data, we may estimate the CDF F(t) by NPMLE
(e.g. by the \code{isotNEW2()} function in this package). 
Based on the NPMLE \eqn{ \hat F_n(t)} we can estimate the mean mu(F).

This function, \code{el.CS.mean}, uses empirical likelihood to test the hypothesis that mu(F)
equal to a given value(mu): i.e. H0: mu(F) = mu.

Empirical likelihood ratio test returns the Wilks statistics, -2LLR.
The -2 log likelihood ratio under H0 is asymptotically chi square DF=1 distributed.
See reference below.
}
\usage{ 
el.CS.mean(mu, Itime, delta, Pfun) 
}
\arguments{
\item{mu}{The hypothesized mean value.}
\item{Itime}{The inspection times, a vector of length n.}
\item{delta}{Either 0 or 1. I[yi <= ti]. length n.}
\item{Pfun}{A given function, \eqn{\Psi(s)}, used to define the (weighted) mean.}
}
\details{
This function tests the null hypothesis that mu(F) = mu versus not equal.
We assume the data given is current status censored data.

The definition of the mean, mu(F) is 
\deqn{
 \mu(F) = \int_0^M [1- F(t)] d \Psi(t) 
}
and its estimator based on \eqn{ (\delta_i, t_i) }
or \eqn{ \hat F_n } is (assume \eqn{ \min(t_i) =0} or \eqn{t_{(1)} =0})
\deqn{
  {\mu(\hat F_n )} = \sum_{i=1}^n [1-\hat F_n(t_{(i)})] \Delta \Psi(t_{(i)})~,
}
where \eqn{ \Psi (t) } is a given function and
\eqn{ \Delta \Psi(t_{(i)})= \Psi (t_{(i+1)}) - \Psi(t_{(i)}) }.
If \eqn{ \Psi(t) =t } in the above, then this
is the ordinary mean (assuming F(t) has support (0, M) ).

The NPMLE \eqn{\hat F_n(t)} is convergent at cubic root n speed, but the mean estimator is 
convergent at ordinary root n speed. The -2LLR has chi square DF=1 null distribution asymptotically.

It goes without saying that we assume the NPMLE mu(hat F)
has finite asymptotic variance (when normalized by root n).
}
\value{
It returns a list containing
\item{"-2LLR"}{The Wilks statistics of the EL test, has approximate chi SQ DF=1 distribution under null hypothesis.}
}
\references{
Zhou, M. (2026). \emph{Empirical Likelihood Method in Survival Analysis 2nd Edition} Chapman & Hall/CRC

Huang, J. and Wellner, J. (1995). \emph{Asymptotic normality of the NPMLE
of linear functionals for interval censored data, case 1}
Statistica Neerlandica \bold{49}, 2, 153--163.

Sun, J. (2006). \emph{The Statistical Analysis of Interval-Censored Failure Time Data} Springer, New York.
}
\author{
Mai Zhou <maizhou@gmail.com>.
}
\examples{
N <- 300
set.seed(12345)
itime <- sort(c(rexp(N-1), 0.5) )       #### inspection times      
Stime <- rexp(N)                       #### survival times
delta <- as.numeric(Stime <= itime)    ####  current status censoring
}
\keyword{nonparametric}
