\name{L2S}
\alias{L2S}
\title{
Convert the Current Status Data from Long Format to Short Format
}
\description{
Given current status data (two vectors of length n): one vector of inspection times,
the other vector = I[Yi <= ti]. 
We assume the input is ordered according to itime. Convert the same data
into short format. Useful when there are a lot of tied inspection times. (saves space).
}
\usage{
L2S(itime, delta) 
}
\arguments{
\item{itime}{The inspection times. Length n.}
\item{delta}{Either 0 or 1. I[yi <= itimei]. Length n.}
}
\details{
Same data set of 850 cases Hepatitis A tests from Bulgaria recorded
in long format \code{hepABulg} in the \code{R} package \code{csci} compare to
the short format 
\code{hepatitisA} in the \code{R} package \code{curstatCI}.

This function (and the sister function \code{S2L})
can be used to convert between the short and long format.
}
\value{
It returns a list containing
\item{Sitime}{The ordered inspection times. Short format. Tied itimes removed.}
\item{fi}{Number of delta's=1, with identical itime ti.}
\item{ni}{Number of itimes tied at ti.}
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