S2L <- function(itime, fi, ni, OneFirst=TRUE){
## Convert current status data from short format to long format.
## We assume the input itime is ordered. fi is the frequency
## of 1's (I[Y < ti]) at ti. ni is the total obs. at ti.
## OneFirst is logical. Affect output format.
## In output delta's, with identical itime, delta=1 comes before =0?

Litime <- rep(itime, ni)

k <- length(ni)

if( OneFirst ) {
OneZero <- rep(c(1,0), k)
OneZeroF <- as.vector( rbind(fi, (ni-fi)) )
delta <- rep(OneZero, OneZeroF)
}
if( !OneFirst ) {
ZeroOne <- rep(c(0,1), k)
ZeroOneF <- as.vector( rbind((ni-fi), fi) )
delta <- rep(ZeroOne, ZeroOneF)
}
list(Litime=Litime, delta=delta)
}
