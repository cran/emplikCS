el.CS.Hz <- function(ti, di, Pfun, thetaMU, error=1e-11, maxit=25) {
###  thetaMU: the constrained value of sum(Lam*dgi) =\int Lam dgi 
###  Pfun:  used to calculate dgi = Pfun(t[i+1]) - Pfun(t[i])
###  First, do a NPMLE of F(t) by PAVA using monotone() via R function
###  isotNEW2() [package emplikCS]

temp <- isotNEW2(x=ti, y=di, a=0.1, b=0.1)   #### a, b values can change.

### From the output temp, we find the values in ti that F(*) have a jump >0.

inde <- which(c(0, diff(temp$y)) > 0)    #### location of (y--jump) 
NPMLEp <- temp$y[inde]                   #### the values do not hace tie
NPMLEp <- NPMLEp[ NPMLEp < 1 ]      #### this only include 0 < F() < 1
Lam0 <- - log( 1- NPMLEp )          #### NPMLE of cumulative hazard Lam0(t)     
dival <- temp$x[inde]    ### at those locations, \hat F or \hat Lambda has jump>0. jump size in NPMLEp or Lam0

J <- length(dival) -1
ni <- rr <- rep(NA, J)
for(k in 1:J) {
     rr[k] <- sum( di[(dival[k] <= ti) & (ti < dival[k+1])] )
  ni[k] <- length( ti[(dival[k] <= ti) & (ti < dival[k+1])] )
}
  
dgi <- diff( Pfun(dival) )   #### constrain dti vector. such as function(t) = t*exp(-t)*as.numeric(t<=1.3)

##### Could take other function to replace exp(-9**). Should have a stop at t=0.8 etc.
##### Make sure the integration is finite. sum( Lam0*dgi )

theta <- sum( Lam0*dgi )    #### The constrain value at Lam0, MLE. Also the initial value. real one = theta +- 0.02 etc

###  Now data is reduced to 3 vector:  dival, rr, ni. (of length J+1; J; J)
###  (the CDF NPMLEp should be rr/ni in the interval [dival[k], dival[k+1])  )
###  We are ready to start hazard calculation. The logLik is
###  logLik = sum(over k=1:J)  rr[k]*log[1-exp(-Lam[k])] - (ni[k]-rr[k])*(Lam[k])
###  We need first and second derivative wrt Lam[k]
###  First derivative: (there are J of the partial derivatives ) = d logLik/d Lam[k]

dvec0 <- rr/(exp(Lam0) -1 ) - (ni - rr)

###  Second partial derivatives of logLik:  k=1 : J (diag matrix)
###  rr[k] * (exp(Lam[k]))/( exp(Lam[k]) -1)^2
###  Write the second derivative matrix as R'R 
###  R = diag(  vecR= sqrt(rr[k]) * (exp(Lam[k]/2))/(exp(Lam[k]) -1) )
###  Finally  the R^{-1} is just 
###  R{-1} =  diag( 1/(vecR) ) = diag(  (exp(Lam[k]) -1)/(sqrt(rr[k])*exp(Lam[k]/2)  )

Dmat0 <- (exp(Lam0) -1)/(sqrt(rr)*exp(Lam0/2) )    #### make it a diag later in the call to solve.QP

#### Constraint: first 0 < Lam[k] < Lam[k+1]  < infinity  (last one no need)
####  Constrain is formulated as: t(A) %*% Lam >= bvec
#### So the t(A) matrix is  of JxJ.  Later, with 1 more constraint, t(A) will be (J+1) by J.

m <- diag( rep(-1, (J-1)) )
k <- cbind( m, rep(0, (J-1)) )
mk <- rbind(rep(0,J), k)
Amat0 <- t( diag( rep(1, J) ) + mk )

bvec0 <- - (t(Amat0)%*%Lam0)

#######################this part for one more constraint ###############
Amat1 <- cbind(dgi,  Amat0)
thetaD = thetaMU - theta
bvec1 <- c(thetaD, bvec0)        #### thetaD = thetaMU - theta = desired(theta) - initial(theta.Lam0)
#######################################################################
####   What is the initial lam?   ==> - log( 1- F(tk) ) ??? YES  =Haz=Lam0

value0 <- solve.QP(diag(Dmat0), dvec0, Amat1, bvec1, meq=1, factorized = TRUE)
Lam <- Lam0 + value0$solution

#### This new Lam  will satisfy the constrain of sum(dgi*Lam) = thetaMU, thetaMU given above.

#### Now do some iteration with solve.QP(), each time with derivative computed at new Lam ###
diff <- 10                                       ##### No change to Amat
m <- 0
  while ((diff > error) & (m < maxit)) {
        dvec <- rr/(exp(Lam) -1 ) - (ni - rr)
        Dmat <- (exp(Lam) -1)/(sqrt(rr)*exp(Lam/2) ) 
        bvec <- - (t(Amat0)%*%Lam)
        bvec <- c(0, bvec)         #### already satisfy constrain, so no change in this constrain
        value0 <- solve.QP(diag(Dmat), dvec, Amat1, bvec, meq=1, factorized = TRUE)
        Lam <- Lam + value0$solution
        diff <- sum(abs(value0$solution))
        m <- m + 1
  }
 #####  Lam output is final solution
Loglik1 <- sum( rr*log(1-exp(-Lam)) - (ni-rr)*(Lam) )
Loglik0 <- sum( rr*log(1-exp(-Lam0)) - (ni-rr)*(Lam0) )
Wilks <- 2*(Loglik0 - Loglik1)

list("-2LLR"=Wilks, location=dival, Haz=Lam, iter=m, error=diff,
              LogLik11=Loglik1, Check=sum(Lam*dgi), thetaMLE=theta) 
}

######  Simulation to show chi square null dist. in a qq plot for el.CS.Hz( ) ###############
#
#set.seed(123)
#### Assume the true Hazard is 0.1t (i.e. rexp(n, rate=0.1) , coupled with 2 function below
#mydgfun <- function(t){0.3*t*(12-t)*as.numeric(t<=12) }   #### for this fun true theta = -8.64
#mydgfun2 <- function(t){0.3*t*(10-t)*as.numeric(t<=10) }   #### for this fun has true theta = -5
#
#SimuHaz <- function(N=600, Pfun=mydgfun2, thetaMU= -5) {
#           itime <- sort( rexp(N, rate=0.1) )            
#           Stime <- rexp(N, rate=0.1)               #### Haz(t) = 0.1t 
#           yy <- as.numeric(Stime <= itime)
#
#           temp <- el.CS.Hz(ti=itime, di=yy, Pfun=Pfun,  thetaMU=thetaMU)
#           return(temp$"-2LLR")
#}
#
