L2S <- function(itime, delta){
## Convert current status data from long format to short format.
## Input itime should be ordered. delta=1 or 0, that is I[Yi <= ti]

temp <- rle( itime )
fi <- ni <- temp$lengths
indx <- cumsum(ni)

fi[1] <- sum( delta[1:indx[1]] )
for(i in 2:length(fi)){
   fi[i] <- sum( delta[(1+indx[i-1]):indx[i]])
}
list(Sitime=temp$values, fi=fi, ni=ni)
}