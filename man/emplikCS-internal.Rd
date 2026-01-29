\name{emplikCS-internal}
\alias{CSQP}
\alias{CSQPprob}
\alias{CSQP2prob}
\alias{CSQPmean}
\alias{CSQPmean2}
\title{Internal emplikCS Functions}
\description{
These are internal functions called by various other functions in the package \code{emplikCS}.
They are R codes of Sequential Quadratic Programing to compute the (various) constrained NPMLE
of \eqn{\hat F_n(t)} with current status data. They all call \code{solve.QP( )} from quadprog package.
They are not intended to be called directly by ordinary users.
}
\usage{
CSQP(Itime, d, w0, error=1e-11, maxit=25)
CSQPprob(Itime, d, w0, t0=0.50, Ft0=0.5, error=1e-11, maxit=25) 
CSQP2prob(Itime, d, w0, t01=0.4, Ft01=0.4, t02=0.6, Ft02=0.6, error=1e-11, maxit=25)
CSQPmean(Itime,d,w0=rep(.5,length(d)),MU,dp=rep(1,length(d)),error=1e-11,maxit=25)
CSQPmean2(Itime,d,w0=rep(.5,length(d)),MU,dp=rep(1,length(d)),error=1e-11,maxit=25) 
}
\details{
\code{CSQP( )} calculate the NPMLE of \eqn{ \hat F_n(t) } by using the SQP method to maximize the
empirical likelihood function subject to monotonicity.  This should give the same result as PAVA algorithm
provided by \code{ isotNEW2( ) } in this package. Should also be the same as \code{ ComputeMLE( ) }
from csci package (from slope of GCM algorithm?). We did some tests and all are equal, but
may be more tests are due.

This function is slower than \code{isotNEW2( )}, therefore not used often, 
other than to proof the SQP algorithm works.

\code{CSQPprob( )} similar to function \code{CSQP( )}, but with an extra constraint of F(t0) = Ft0.
The SQP can easily add one equality constraint in addition to the n-1 monotone inequality constraints.
The PAVA and slope of GCM algorithm may also add a constraint, but SQP seem the easiest to do.


\code{CSQP2prob( )} similar to \code{CSQPprob( )} but with two constraints. We may write one for three
constraints similarly.

\code{CSQPmean( )} similar to \code{CSQPprob} but the constraint is in terms of the mean
         (or a linear functional) of the CDF F(t). The PAVA and slope of GCM algorithm seems not usable here. 

\code{CSQPmean2( )} similar to \code{CSQPmean( )}.
}
\keyword{internal}
