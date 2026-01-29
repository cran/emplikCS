el.CS.mean <- function(mu, Itime, delta, Pfun){
#### Itime include all inspection times and sorted. No tie. [no need to sort]
#### This function differs with el.CS.mean2 only by adding the 
#### CSDataclean() in the beginning, and eliminate the input start/end.
#### Also, no need to have CSDataclean()  in the DATAGenerat1/2().

temp0 <- CSdataclean(itime=Itime, delta=delta)   #### from input
Itime <- temp0$itime
delta <- temp0$delta
Istart <- temp0$Istart
Iend <- temp0$Iend
#### Istart=position of first 1's; Iend=position of last 0's in
#### the sorted delta vector. For ti below Itime[Istart], F(ti)=0. 
#### For ti over Itime[Iend], F(ti)=1.
 
if(min(Itime) < 0) 
     stop("we assume all inspection times >= 0. Supp(F)=(0, M)")
if(Istart > Iend) stop("data: (0's followed by 1's)")

Itime2 <- Itime[Istart:Iend]
delta2 <- delta[Istart:Iend]
 
Itime1 <- Itime[1:Istart]
if(min(Itime1) > 0) {Itime1 <- c(0, Itime1)}
mu1 <- 0 
if(length(Itime1) > 1) { 
      temp <- Pfun(Itime1)
      mu1 <- sum( diff(temp) )
      }
if(length(Itime) == Iend) {Itime3 <- c(Itime2, Itime[Iend]+1)}
#### May add other options other than +1 ?
if(length(Itime) > Iend) {Itime3 <- Itime[Istart:(Iend+1)]}
FunItime3 <- Pfun(Itime3)
dpvec <- diff(FunItime3)

temp00 <- CSQP(Itime=Itime2, d=delta2, w0=rep(0.5, length(delta2)))
LogLikH0 <- temp00$LogLik11

LogLikHA <- CSQPmean2(Itime=Itime2, d=delta2, w0=temp00$prob,
                                     MU=mu-mu1, dp=dpvec)$LogLik11            
list( "-2LLR"= 2*(LogLikH0 - LogLikHA) )
}

