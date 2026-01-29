\name{CSdataclean}
\alias{CSdataclean}
\title{
Order the Given Current Status Data and Output with Some Additional Information
}
\description{
Given N current status data, order according to inspection times, and
find the positions (in the ordered data) of the first 1 and last 0 of the delta.
}
\usage{
CSdataclean(itime, delta) 
}
\arguments{
\item{itime}{The inspection times. Length N.}
\item{delta}{Either 0 or 1. I[yi <= itimei]. Length N.}
}
\details{
When calculate NPMLE of the CDF F(t) from current status data, it is obvious that 
F(t) = 0  before the location (itime) of first delta=1 occur in the delta list, 
and similarly, F(t) = 1 (one position) after the last delta=0 occur. So in the calculation of NPMLE
we need only to consider F(.) when time t is from itime[first] to itime[last]. 

We take the definition of NPMLE \eqn{ \hat F_n(t) } as right continuous.

Usually, the current status data are stored in either long format or short format.
The short format is often used when there are many tied inspection times.
This function, \code{CSdataclean}, takes input the current status data in the long form:
\code{itime=(t1, t2, ... tN)} and \code{delta=(0,1,... 1)}. The only values in the delta are 0 or 1.
}
\value{
It returns a list containing
\item{itime}{The ordered inspection times.}
\item{delta}{The delta, ordered according to itime.}
\item{Istart}{delta[Istart] is the first 1 in the ordered delta output.}
\item{Iend}{delta[Iend] is the last 0 in the ordered delta output.}
}
\references{
Zhou, M. (2026). \emph{Empirical Likelihood Method in Survival Analysis 2nd Edition} Chapman & Hall/CRC
}
\author{
Mai Zhou <maizhou@gmail.com>.
}
\examples{ 
y <- c(10, 209, 273, 279, 324, 391, 566, 785)
x <- c(21, 38, 39, 51, 77, 185, 240, 289, 524)
}
\keyword{nonparametric}
