\name{S2L}
\alias{S2L}
\title{
Convert the Current Status Data from Short Format to Long Format
}
\description{
Given current status data in short format (three vectors),
convert it to the long format (one vector of inspection times,
another vector=I[Yi <= ti]). 
We assume the input is ordered according to inspection times itime.
}
\usage{
S2L(itime, fi, ni, OneFirst=TRUE) 
}
\arguments{
\item{itime}{The unique inspection times. Tie removed. Length k.}
\item{fi}{Frequency of 1's in the delta.  Definition is I[yi <= itimei]. Length k.}
\item{ni}{Frequency of inspection times ti tied at this value. Length k.}
\item{OneFirst}{In the output of long format data, within the same inspection times, delta=1's comes before 0's.
                Default to TRUE.}
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
\item{Litime}{The ordered inspection times. Long format.}
\item{delta}{This is I[Yi <= ti].}
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