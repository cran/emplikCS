isotNEW2<- function(x,y,a,b, LONG=TRUE) {
### This version require the R package monotone, to do the heavy lifting.  
### the parameter a, b are to argument output dist. F
### so that F(x[1]-a)=0 and F(x[n]+b) =1; always (a proper dist.)
### y are 0's and 1's (in our application, but may not in other cases).
if ( length(x) != length(y) ) stop("check length of x and/or y")
if (any((y != 0) & (y != 1))) stop("y values, I[Si <= xi], must be either 0 or 1")

    xorder <- order(x)                          # 
    xsort <- x[xorder]                          # this block updated
    yxorder <- y[xorder]                        # March 29, 1995
    temp <- rle(xsort)                          # (replace a for loop)
    uniquexsort <- temp$values
    nsort <- temp$lengths
    ysort <- cumsum(yxorder)[cumsum(nsort)]
    ysort <- diff(c(0,ysort))
    m <- length(nsort)

p <- ysort/nsort
temp <- monotone(p, w=nsort)    
## The monotone() function from R package monotone. Avoid for() loop.
## output temp is of the same length as p.  Mai Z. 2025/6/10

longOUT <- list(x=c(uniquexsort[1]-a,uniquexsort,uniquexsort[m]+b), y=c(0, temp, 1))

if(LONG) return(longOUT)
        else {
   inde <- which(c(0, diff(longOUT$y)) > 0)    #### location of (y--jumps), [jump>0] 
   NPMLEp <- longOUT$y[inde]                   #### the values do not have tie   
   dival <- longOUT$x[inde]                
   shortOUT <- list(x=dival, y=NPMLEp)
   return(shortOUT) }
## In our applications (CDF from current status data), usually the
## jumps in y values are few and far between. For example, if sample size=1000,
## the number of jumps are around 10. So sometime we want to
## make the output in a short format. only out put locations where is a jump >0.
}

## Suppose longOUT is the output from isotNEW2() [long format]
## inde <- which(c(0, diff(longOUT$y)) > 0)    #### location of (y--jumps), [jump>0] 
## NPMLEp <- longOUT$y[inde]                   #### the values do not have tie   
## dival <- longOUT$x[inde]                    ##
## The short format, shortOUT, is:  (locations=dival,  jumps=NPMLEp) We assume right continuous
## of the resulting hat F(t), so there is no need to add a location of 0

