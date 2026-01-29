CSQPprob <- function(Itime, d, w0, t0=0.50, Ft0=0.5, error = 1e-11, maxit = 25) 
{
    ####  Here we try to use SQP to find the NPMLE of CDF with 
    ####  current status data, AND with an extra constraint of: F(t0)=theta=Ft0, with t0 given too.
    ####  When there is no constraint, we do it with CSQP() OR use monotone() OR use isotNEW()/isotNEW2().
    ####  The one extra constraint is F(t0) = theta=Ft0. 
    ####  w0 should be the initial prob, preferable from monotone() or CSQP(), (that
    ####  is the NPMLE without constraint) but should also work if w0 == rep(0.5, n). 
    ####  Itime is inspection times, will be sorted later, assume no ties.
    ####  d is the indicator, I(T <= Itime) as usual. Dist. for T is parameter.

    ####  Test data: Itime: 1,   2,  3,   4,   5,   6,   7,  8.
    ####                d : 0,   1,  0,   1,   1,   1,   0,  1.
    ####  the answer is     0,  .5, .5, .75, .75, .75, .75,  1.
    ####  This is checked out by isotNEW() and isotNEW2(); however we
    ####  need to get rid of the first observation (because it has d=0)
    ####  Similarly, we need to get rid of the last observation (because it has d=1)
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

#### Assume t0 is equal to one (and only one) of the obs. times sortedIvec.  
#### This time, we want to add one extra constraint in the form of  F(t0) = Ft0.
#### That is assume t0 == one of the itime. 

#### find the position kk in sortedIvec that equal to t0. 
    kk <- which(sortedIvec == t0)
    if(length(kk) != 1) stop("no equal or more than one equal?")

    Amat0 <- diag(rep(-1, n))[,-n] + diag(rep(1, n))[, -1]

    bvec0 <- - (t(Amat0)%*%w0)

    Exconst <- rep(0, n)
    Exconst[kk] <- 1
 
    Amat <- cbind(Exconst, Amat0)             ### one constraint for w[kk]
    bvec0 <- c(Ft0-w0[kk], bvec0)             ### one constraint, w[kk] = Ft0

    value0 <- solve.QP(diag(Dmat0), dvec0, Amat, bvec0, meq=1, factorized = TRUE)
    w <- w0 + value0$solution
 
  ### if (any(w <= 0)) 
  ###     stop("There is no probability satisfying the constraints")

  diff <- 10
  m <- 0
  while ((diff > error) & (m < maxit)) {
       dvec <- d/w - (1-d)/(1-w)
       Dmat <- d*w + (1-d)*(1-w)
       bvec <- - (t(Amat0)%*%w)
       bvec <- c(0, bvec)
       value0 <- solve.QP(diag(Dmat), dvec, Amat, bvec, meq=1, factorized = TRUE)
       w <- w + value0$solution
       diff <- sum(abs(value0$solution))
       m <- m + 1
  }
   LogLik1 <- sum(d*log(w)+(1-d)*log(1-w))
  list(prob=w, iter=m, error=diff, LogLik11=LogLik1, ObsTime=sortedIvec, index=kk, Check=w[kk])  

  ###  lik00 <- sum(ww * log(dvec00))
  ###  tval <- 2 * (lik00 - sum(ww * log(w)))
  ###  list(llik00=lik00, llik11=sum(ww*log(w)), `-2LLR` = tval, 
  ###       Pval = 1 - pchisq(tval, df = kk), prob1 = w[dd == 1], 
  ###       xtime = xx[dd == 1], iteration = m, error = diff)
}
