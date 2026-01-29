CSbj <- function(x, delta, Itime, maxiter = 99, error = 0.0001) {  
### x is a design matrix of N rows p cols. (We require each col sum to 0.
### that is x  <- (x - bar x). This also exclude an intercept term in x.) NOT? 
### Itime is a vector of length N. (inspection time) ti.
### delta is a vector of length N. (value +1 and -1)
### delta =1 mean the actual (response) obs (yi) is > ti.  delta = -1 means yi < ti. 
### This definition may be different from current status (which has +1, 0).
### Output: beta is a vector of length = no. of column(x) =p
### The model before censoring is yi = alpha + xi*beta + epsiloni.
### We are estimating beta only, not alpha. alpha+epsilon = baseline parameter.

dimofX <- dim(x)
if( dimofX[1] != length(delta) ) stop("check dim of x and/or delta")
x <- x - matrix(colMeans(x), ncol=ncol(x), nrow=nrow(x), byrow=TRUE)

      newtemp <- matrix(NA, ncol=dimofX[2], nrow=4)

    #### Get an initial beta est. if not given ####
        FakeY <- delta 
        posD1 <- which(delta == 1)
        posDn1 <- which(delta == -1)
       FakeY[posD1] <- Itime[posD1] + 1    ### or + some const?
       FakeY[posDn1] <- Itime[posDn1] - 1

      newtemp[1,] <- lsfit(x, FakeY, intercept=FALSE)$coef  
    ###  newtemp[1,]<-newtemp[1,]/newtemp[1,1]   ### seems do not need to normalize??
    ################################################### initial est#

      for (i in 2:4) { 
      newtemp[i,] <- CSiter2(x, delta, insp=Itime, beta=newtemp[i-1,])
    ####  newtemp[i,]<-newtemp[i,]/newtemp[i,1]
      } 

      num <- 3 
      while(num <= maxiter && error < max(abs(newtemp[3,]-newtemp[4,])) ) {
             newtemp[3,] <- newtemp[4,]
             newtemp[4,] <- CSiter2(x, delta, insp=Itime, beta=newtemp[3,])
      #######   newtemp[4,]<-newtemp[4,]/newtemp[4,1]
             num <- num + 1 
      } 
      delta[delta < 0.5] <- 0
      ut <- Itime - x %*% newtemp[4,]
      tempF <- isotNEW2(x=ut, y=1-delta, a=1, b=1)  ### Some other a, b?
      list( est=newtemp[4,], iterN=num,  distFx=tempF$x, distFy=tempF$y )
}

CSiter2 <- function(x, delta, insp, beta) {  
### x is a matrix of N rows, delta is a vector of length N.
### insp is a vector of length N.
### delta =1 mean the actual obs (y) is > insp(i). -1 otherwise.
### beta is a vector of length = no. of column(x)

N <- length(delta)
k <- length(as.vector(beta))

if (dim(x)[1] != N) stop("check dim of x")
if (dim(x)[2] != k) stop("check length of beta and dim of x") 
u <- x %*% beta 
ystar <- u      # should I just let  ystar <- delta ?
uu <- insp - u 

dd <- delta
dd[delta<0.5] <- 0 
temp <- isotNEW2(uu,1-dd,2,2)
tx <- temp$x
ty <- temp$y
TN <- length(tx)

   for(i in 1:N) {
              m <- length(tx[tx <= uu[i]])
              if(delta[i] > 0.5 && ty[m] < 1 ) {
                    ystar[i] <- u[i] + 0.5*sum((tx[(m+1):TN]+tx[m:(TN-1)])*diff(ty[m:TN])) / (1-ty[m]) 
                                }
              else if( ty[m] > 0 ) {
                    ystar[i] <- u[i] + 0.5*sum((tx[2:m]+tx[1:(m-1)])*diff(ty[1:m])) / ty[m]
                    }
    }
return( lsfit(x,ystar,intercept=FALSE)$coef ) 
}







