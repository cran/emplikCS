CSQPmean2 <- function(Itime, d, w0=rep(0.5, length(d)), MU, dp=rep(1, length(d)), error = 1e-11, maxit = 25) 
{
    ####  Newer version with one more input dp = Delta \Psi(t_i).
    ####  Here we try to use SQP to find the (constrained) NPMLE of CDF with 
    ####  current status data, AND with an extra (mean) constraint defined by
    ####  sum[1-F(t_i)]dp_i = MU. ( ==> or  sum(dp_i) - MU = sum[F(t_i)dp_i] ). 
    ####
    ####  If we have equal-spaced inspection times 0,1,2,3, ..., n/(n+1); then Dt_i=1/(n+1).
    ####  we may let dp_i=1 and MU=true mean*(n+1). ( this ==> sum[1-F(t_i)]=mean*(n+1),
    ####  i.e.  sum[1-F(t_i)]/(n+1) = mean ).
    ####   One sticky point: the SQP only works for F(ti) that are 0<F(t0)<1 strictly. 
    ####  For F(ti)=0 or F(ti)=1 we need to manually adjust the MU value to account for that.
    ####  w0 should be the initial prob, preferable from monotone() or CSQP(),
    ####  (that is the NPMLE without constraint) but should also work (less well)
    ####  if w0 == rep(0.5, n). w0 do not need to satisfy constraint.
    ####  Itime is inspection times, will be sorted later, assume no ties.
    ####  d is the censor indicator I(lifetime <= Itime) as usual. 

    ####  Test data: Itime: 1,   2,  3,   4,   5,   6,   7,  8.
    ####                d : 0,   1,  0,   1,   1,   1,   0,  1.
    ####  the answer is     0,  .5, .5, .75, .75, .75, .75,  1.
    ####  This is checked out by isotNEW() and isotNEW2(); however we
    ####  need to get rid of the first observation (because it has d=0)
    ####  Similarly, we need to get rid of the last observation (because it has d=1)
    ####  In general, smallest Itime(s) with d=0 needs to be deleted, since here w=0.
    ####  similarly, largest Itime(s) with d=1 needs to be deleted, since here w = 1.
    ####  These observations do not contribute to loglik, but they do contribute to mean value.
    ####  And we have to compute log(w), log(1-w), so w=0/w=1 will be problematic.
    ####  In the LogLik contributions, these terms will be zero (Lik will be *1).
    ####  The logLik is sum[d*log(w) + (1-d)*log(1-w)] .


    Ivec <- as.vector(Itime)
    n <- length(Ivec)
    if (length(d) != n) 
        stop("length of Itime and d must agree")
    if (any((d != 0) & (d != 1))) 
        stop("d must be 0(<itime) or 1(>itime)")
    if (!is.numeric(Ivec)) 
        stop("Itime must be numeric")
    if (length(dp) != n) 
        stop("length of d and dp must agree")
     if (length(w0) != n) 
        stop("length of d and w0 must agree")
 
    Iorder <- order(Ivec)
    sortedIvec <- Ivec[Iorder]
    sortedd <- d[Iorder]
    sortedw0 <- w0[Iorder]
    sorteddp <- dp[Iorder]
 
    if(sortedd[1] == 0) stop("the smallest Itime has d=0")
    if(sortedd[n] == 1) stop("the largest Itime has d=1")
    LogLik0 <- sum(sortedd*log(sortedw0)+(1-sortedd)*log(1-sortedw0))

    #### May be the best is to sort everything (according to Itime)
    #### before entering the fun. and collaps the tie, with weight. Using Wdataclean()??
    #### But first, assume no tie.
    #### print(sum(w0))

    Dmat0 <- d*w0 + (1-d)*(1-w0) 
    dvec0 <- d/w0 - (1-d)/(1-w0)

#### This time, we want to add one more constraint in the form of sum(1-w)dp = MU.
#### However, we need to rewrite it as:  sum(w*dp) = sum(dp) - MU. 

    Amat0 <- diag(rep(-1, n))[,-n] + diag(rep(1, n))[, -1]

    bvec0 <- - (t(Amat0)%*%w0)
 

    Amat <- cbind(sorteddp, Amat0)         ### one more constraint, sum(w*dp)=sum(dp)-MU
    MU0 <- sum((1-sortedw0)*sorteddp)      ### initial mean estimator
    bvec0 <- c(MU0-MU, bvec0)              ### one more constraint value.

    value0 <- solve.QP(diag(Dmat0), dvec0, Amat, bvec0, meq=1, factorized = TRUE)
    w <- w0 + value0$solution
    #### print(sum(w)) 
 
  diff <- 10
  m <- 0
  while ((diff > error) & (m < maxit)) {
       dvec <- d/w - (1-d)/(1-w)
       Dmat <- d*w + (1-d)*(1-w)
       bvec <- - (t(Amat0)%*%w)
       bvec <- c(0, bvec)             #### after one iteration, already satisfy H0.
       value0 <- solve.QP(diag(Dmat), dvec, Amat, bvec, meq=1, factorized = TRUE)
       w <- w + value0$solution
       diff <- sum(abs(value0$solution))
       m <- m + 1
  }
  LogLik1 <- sum(d*log(w)+(1-d)*log(1-w))
  list(prob=w, iter=m, error=diff, LogLik11=LogLik1, MU0=MU0, Check=sum((1-w)*sorteddp))  

  ###  lik00 <- sum(ww * log(dvec00))
  ###  tval <- 2 * (lik00 - sum(ww * log(w)))
  ###  list(llik00=lik00, llik11=sum(ww*log(w)), `-2LLR` = tval, 
  ###       Pval = 1 - pchisq(tval, df = kk), prob1 = w[dd == 1], 
  ###       xtime = xx[dd == 1], MU0=MU0, iteration = m, error = diff)
}
