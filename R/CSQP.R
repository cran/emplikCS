CSQP <- function(Itime, d, w0, error = 1e-11, maxit = 25) 
{
    ####  Here we try to use SQP to find the NPMLE of CDF with 
    ####  current status data, (later with an extra (mean/prob) constraint).
    ####  But FIRST, let us do it without extra constraint and 
    ####  compare to monotone() to make sure we done SQP right.
    ####  w0 should be the initial prob, preferable from monotone(), but
    ####  also should work if w0 == 0.5 etc. 
    ####  Itime is inspection times, will be sorted later, assume no ties.
    ####  d is the indicator, I[T <= Itime] as usual.
    ####  F(t) for T is the parameter of interest.

    ####  Test data: Itime: 1,   2,  3,   4,   5,   6,   7,  8.
    ####                d : 0,   1,  0,   1,   1,   1,   0,  1.
    ####  the answer is     0,  .5, .5, .75, .75, .75, .75,  1.
    ####  This is checked out by isotNEW() and isotNEW2(); however we
    ####  need to get rid of the first observation here (because it has d=0)
    ####  Similarly, we need to get rid of the last observation (because it has d=1)
    ####  In general, smallest Itime(s) with d=0 needs to be deleted, since w here will be 0.
    ####  similarly, largest Itime(s) with d=1 needs to be deleted, since here w = 1.
    ####  And we have to compute log(w), log(1-w), so w=0/w=1 will be problematic.
    ####  In the LogLik contributions, these terms will be zero (Lik will be x 1).
    ####  The logLik is sum d log(w) + (1-d) log(1-w) .
    ####  We need to CSdataclean() input data before call this function.
    Ivec <- as.vector(Itime)
    n <- length(Ivec)
    if (length(d) != n) 
        stop("length of Itime and d must agree")
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
    LogLik0 <- sum(sortedd*log(sortedw0)+(1-sortedd)*log(1-sortedw0))

    #### May be the best is to sort everything (according to Itime)
    #### before entering the fun. and collaps the tie, with weight. Using Wdataclean()??
    #### But first, assume no tie.
    ###############################################

    Dmat0 <- d*w0 + (1-d)*(1-w0) 

    dvec0 <- d/w0 - (1-d)/(1-w0)

    Amat <- diag(rep(-1, n))[,-n] + diag(rep(1, n))[, -1]
    #### Amat <- cbind( c(1, rep(0, n-1)), Amat)
    bvec0 <- - (t(Amat)%*%w0)
    
    value0 <- solve.QP(diag(Dmat0), dvec0, Amat, bvec0, meq=0, factorized = TRUE)
    w <- w0 + value0$solution

  ### if (any(w <= 0)) 
  ###     stop("There is no probability satisfying the constraints")

  diff <- 10
  m <- 0
  while ((diff > error) & (m < maxit)) {
       dvec <- d/w - (1-d)/(1-w)
       Dmat <- d*w + (1-d)*(1-w)
       bvec <- - (t(Amat)%*%w)
       value0 <- solve.QP(diag(Dmat), dvec, Amat, bvec, meq = 0, factorized = TRUE)
       w <- w + value0$solution
       diff <- sum(abs(value0$solution))
       m <- m + 1
  }
  LogLik1 <- sum(d*log(w)+(1-d)*log(1-w))
  list(prob=w, times=sortedIvec, D=sortedd, LogLik00=LogLik0, LogLik11=LogLik1, iteration=m, error=diff)  
}

#### Output do not include prob=0 or 1. 