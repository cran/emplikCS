el.CS.prob <- function(ti, di, t0=0.5, Ft0=0.5) {
## ti is the inspection time vector [do we require ordering? No, CSdataclean will do ordering]
## di is the associated indicator I[X_i < ti]


temp <- CSdataclean(itime=ti, delta=di)
itime2 <- temp$itime[temp$Istart:temp$Iend]
delta2 <- temp$delta[temp$Istart:temp$Iend]

tem <- CSQP(Itime=itime2, d=delta2, w0=rep(0.5, length(delta2)))
tem2 <- CSQPprob(Itime=itime2, d=delta2, w0=tem$prob, t0=t0, Ft0=Ft0)
Wtheta <- 2*(tem$LogLik11 - tem2$LogLik11)

list( "-2LLR"=Wtheta, LogLik0=tem$LogLik11, LogLik1=tem2$LogLik11)
}


