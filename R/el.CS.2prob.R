el.CS.2prob <- function(ti, di, t01, Ft01, t02, Ft02){
## ti is the inspection time vector [do we require ordering? No, CSdataclean will do ordering]
## di is the associated indicator I[X_i < ti]

temp <- CSdataclean(itime=ti, delta=di)
itime2 <- temp$itime[temp$Istart:temp$Iend]
delta2 <- temp$delta[temp$Istart:temp$Iend]

tem0 <- CSQP(Itime=itime2, d=delta2, w0=rep(0.5, length(delta2)))
tem2 <- CSQP2prob(Itime=itime2, d=delta2, w0=tem0$prob, t01=t01, Ft01=Ft01, t02=t02, Ft02=Ft02)
Wtheta <- 2*(tem0$LogLik11 - tem2$LogLik11)

list( "-2LLR"=Wtheta, LogLik0=tem0$LogLik11, LogLik1=tem2$LogLik11 )
}


