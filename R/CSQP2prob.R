CSQP2prob <- function(Itime, d, w0, t01=0.4, Ft01=0.4, t02=0.6, Ft02=0.6, error=1e-11, maxit=25) 
{
    ####  Here we try to use SQP to find the NPMLE of CDF with 
    ####  current status data, AND with TWO extra constraint:
    ####   F(t01)=theta1=Ft01,  F(t02)=theta2=Ft02 .
  
    ####  w0 should be the initial prob, preferable from monotone() or CSQP(), (that
    ####  is the NPMLE without constraint) but should also work, but less well, 
    ####  if w0 == rep(0.5, n). 
    ####  Itime is inspection times, will be sorted later, assume no ties.
    ####  d is the indicator, I(T <= Itime) as usual. 

    ####  We assume t01 (and also t02) is actually equal to one and only one
    ####  of the inspection times. (Too restrictive??)

    ####  In general, smallest Itime(s) with d=0 needs to be deleted, since here w=0.
    ####  similarly, largest Itime(s) with d=1 needs to be deleted, since here w = 1.
    ####  And we have to compute log(w), log(1-w), so w=0/w=1 will be problematic.
    ####  In the LogLik contributions, these terms will be zero (Lik will be x1).
    ####  The logLik is sum d log(w) + (1-d) log(1-w) .

    Ivec <- as.vector(Itime)
    n <- length(Ivec)
    if (length(d) != n) 
        stop("length of Itime and d must agree")
    if (length(w0) != n)
        stop("length of w0 must be same as length of d")
    if (any((d != 0) & (d != 1))) 
        stop("d must be 0(<itime) or 1(>itime)")
    if (!is.numeric(Ivec)) 
        stop("Itime must be numeric")
          
    Iorder <- order(Ivec)
    sortedIvec <- Ivec[Iorder]
    sortedd <- d[Iorder]
    sortedw0 <- w0[Iorder]
 
    if(sortedd[1] == 0) stop("the smallest Itime has d=0")
    if(sortedd[n] == 1) stop("the largest Itime has d=1")
    LogLik0 <- sum(sortedd*log(sortedw0)+(1-sortedd)*log(1-sortedw0)) ##starting loglik

    #### May be the best is to sort everything (according to Itime)
    #### before entering the fun. and collaps the tie, with weight. Using Wdataclean()??
    #### But first, assume no tie. and assume input already ordered according to Itime.
    #### print(sum(w0))

    Dmat0 <- d*w0 + (1-d)*(1-w0)         ## second derivative

    dvec0 <- d/w0 - (1-d)/(1-w0)         ## first derivative

#### Assume t01 is equal to one (and only one) of the inspection times, sortedIvec.  
#### This time, we want to add 2 extra constraints in the form of  F(t0) = Ft0.
#### That is assume t0 == one of the itime. 

#### find the position kk1/2 in sortedIvec that equal to t01/2. 
    kk1 <- which(sortedIvec == t01)
    if(length(kk1) != 1) stop("no equal or more than one equal t01?")

    kk2 <- which(sortedIvec == t02)
    if(length(kk2) != 1) stop("no equal or more than one equal t02?")

    Amat0 <- diag(rep(-1, n))[,-n] + diag(rep(1, n))[, -1]

    bvec0 <- - (t(Amat0)%*%w0)

    Exconst1 <- rep(0, n)
    Exconst1[kk1] <- 1

    Exconst2 <- rep(0, n)
    Exconst2[kk2] <- 1
 
    Amat <- cbind(Exconst1, Exconst2, Amat0)             ### 2 constraints for w[kk1], w[kk2]
    bvec0 <- c(Ft01-w0[kk1], Ft02-w0[kk2], bvec0)         ### 2 constraints, w[kk12] = Ft012

    value0 <- solve.QP(diag(Dmat0), dvec0, Amat, bvec0, meq=2, factorized = TRUE)
    w <- w0 + value0$solution
 
  ### if (any(w <= 0)) 
  ###     stop("There is no probability satisfying the constraints")

  diff <- 10
  m <- 0
  while ((diff > error) & (m < maxit)) {
       dvec <- d/w - (1-d)/(1-w)
       Dmat <- d*w + (1-d)*(1-w)
       bvec <- - (t(Amat0)%*%w)
       bvec <- c(0, 0, bvec)
       value0 <- solve.QP(diag(Dmat), dvec, Amat, bvec, meq=2, factorized = TRUE)
       w <- w + value0$solution
       diff <- sum(abs(value0$solution))
       m <- m + 1
  }
   LogLik1 <- sum(d*log(w)+(1-d)*log(1-w))
  list(prob=w, iter=m, error=diff, LogLik11=LogLik1, ObsTime=sortedIvec, index=c(kk1,kk2),
                   checktime=c(sortedIvec[kk1], sortedIvec[kk2]), Check=c(w[kk1], w[kk2]) )  
}

