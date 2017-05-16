foi_renewal<-function(I,N,Theta,rSI,T_prior)
{
  I0 <- matrix(0,N,T_prior)
  I0[,1] <- Theta[(N+1):(2*N)]
  
  SI_trunc <- length(rSI)-1
  for (i in 2:T_prior){
    f <- max(c(1,(i-SI_trunc)))
    I0[,i] <- Theta[1:N]*(I0[,f:i]%*%rSI[((SI_trunc+1)-(i-f)):(SI_trunc+1)])
  }
  I_tot<- cbind(I0,I)
  
  lambda=matrix(NA,N,ncol(I))
  for (k in (T_prior+1):ncol(I_tot)){
    f=k-SI_trunc
    lambda[,k-T_prior]=I_tot[,f:k]%*%rSI[((SI_trunc+1)-(k-f)):(SI_trunc+1)]
  }
  return(lambda)
}