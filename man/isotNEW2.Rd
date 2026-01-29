\name{isotNEW2}
\alias{isotNEW2}
\title{
Given Current Status Data, Calculate the NPMLE of CDF F(t)
by Calling the monotone Function from monotone Package
}
\description{
Using improved PAVA algorithm to calculate the NPMLE of CDF F(t). Input inspection times (x) can
have ties. We require two extra inputs a and b to make sure the output is a proper CDF:
F(t(1)-a) =0, and F(t(n)+b) =1 and the rest are increase from 0 increase to 1.
}
\usage{
isotNEW2(x, y, a, b, LONG=TRUE) 
}
\arguments{
\item{x}{Inspection time of the current status data.}
\item{y}{Equivalent to delta in the current status data. Either 0 or 1 as in I[survTi <= xi]. length N.}
\item{a}{To make sure the output F(.) is a proper CDF: F(x[1]-a) = 0 always.}
\item{b}{To make sure the output F(.) is a proper CDF: F(x[n]+b) = 1 always.}
\item{LONG}{Should the output in the LONG format or not? Default is TRUE.}
}
\details{
Usually, the current status data are stored in either long format or short format.
The short format is often used when there are many tied inspection times.
This function, \code{isotNEW2}, takes in the current status data in the long form:
\code{x=(t1, t2, ... tN)} and \code{y=(0, 1, ... 1)}. The only values of the \code{y} are 0 or 1.

For more details please refer to monotone package.

May be we should put a=1, b=1 as default. 

Since the NPMLE of F(t) has very few number of jumps
(for sample size N, the number of positive jumps are about cubic root N), 
a short format output can same space.
For current status data of sample size 1000, usually the NPMLE F(t) has about 10 jumps.
So the short format has length about 10 and the long format has length 1000.
}
\value{
It returns a list in either the long format or short format containing
\item{x}{The ordered inspection times, including x[1]-a, and x[n]+b.}
\item{y}{The NPMLE of F(x) at the x time, F(x[1]-a)=0 always; F(x[n]+b)=1 always.}
}
\references{
Busing, F. (2022). \emph{Monotone regression: A simple and fast O(n) PAVA
implementation.} Journal of Statistical Software,
Code Snippets 102, 1 (2022), 1--25. \doi{https://doi.org/10.18637/jss.v102.c01}

Zhou, M. (2026). \emph{Empirical Likelihood Method in Survival Analysis 2nd Edition} Chapman & Hall/CRC
}
\author{
Mai Zhou <maizhou@gmail.com>.
}
\examples{ 
itime <- c(1,  2,  3,  4,  5,  6,  7,  8)
delta <- c(0,  1,  0,  1,  1,  1,  0,  1)

isotNEW2(x=itime, y=delta, a=0.5, b=0.5)
## $x
##  [1] 0.5 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 8.5
## 
## $y
##  [1] 0.00 0.00 0.50 0.50 0.75 0.75 0.75 0.75 1.00 1.00
####  the correct answer is  F(t) =  0,  .5, .5, .75, .75, .75, .75,  1
####  at the ordered itime, augumented by a and b.
isotNEW2(x=itime, y=delta, a=0.5, b=0.5, LONG=FALSE)
## $x
## [1] 2 4 8
## 
## $y
## [1] 0.50 0.75 1.00  #### for time t < 2, F(t) = 0 by right cont.
}
\keyword{nonparametric}
