\name{CSbj}
\alias{CSbj}
\title{
Current Status Data Buckley-James Estimator For Linear Regression Models
}
\description{
In an AFT regression model, when the reponses are current status censored (observe yi either > ti or <= ti),
we may still estimate the regression coefficients by the Buckley-James (extension from right censored case).
We assume the inspection time ti have a larger support to cover the support of error epsilon,
which is assumed iid.
}
\usage{
CSbj(x, delta, Itime, maxiter = 99, error = 0.0001) 
}
\arguments{
\item{x}{Design matrix N row p col. }
\item{delta}{Either 0 or 1. I[yi <= Itimei]. Length N. yi = beta xi + ei}
\item{Itime}{The inspection times. Length N.}
\item{maxiter}{Default to 99. Control the iteration.}
\item{error}{Default to 0.0001. Control the iteration.}
}
\details{
This function is an implementation of the Buckley-James estimator for the
regression parameter beta in the 
AFT regression model when the observed responses are current status censored. Similar
to the Binary Choice model in econometrics where all the inspection times are fixed at zero.
I wrote an S-plus function for the binary choice model (name bibj).
It is easily adapted to the current status situation, and this is the function.
The AFT model we considered here has a intercept term. But we try to estimate the regression
parameter beta, without intercept term.
The estimator of intercept can be obtained as the mean of the iid error term
after we got the estimator of the slope terms.

Depends on the functions \code{monotone} from package monotone
and \code{lsfit} from the basic stats package.

At this point, we do not have a good estimate for the variance for 
the Buckley-James estimators.
Bootstrap is one method one can try.
}
\value{
It returns a list containing
\item{est}{The Buckley-James estimator of the regression coefficients.}
\item{iterN}{Number of iterations done.}
\item{distFx}{Locations of the jumps of final estimator of the error distribution.}
\item{distFy}{Probabilities of the final estimator of the error distribution at jump locations.
      Mean of this error distribution is the intercept term est. of the regression model.}
}
\references{
Zhou, M. (2026). \emph{Empirical Likelihood Method in Survival Analysis 2nd Edition} Chapman & Hall/CRC

Wang, W., and Zhou, M. (1995). \emph{Iterative least squares estimator of binary
choice models: a semi-parametric approach.} Tech. Report, University of
Kentucky.  \url{https://www.ms.uky.edu/~mai/research/eco5wz.pdf}

Buckley, J. J., and James, I. R. (1979). \emph{Linear regression with censored
data.} Biometrika \bold{66}, 429--436.
}
\author{
Mai Zhou <maizhou@gmail.com>.
}
\examples{ 
y <- c(10, 209, 273, 279, 324, 391, 566, 785)
x <- c(21, 38, 39, 51, 77, 185, 240, 289, 524)
}
\keyword{nonparametric}
