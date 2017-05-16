Likelihood_renewal<-function(foi,I,R0)
{
  # L <- sum(rowSums(-R0*foi+I*log(R0*foi),na.rm=TRUE),na.rm=TRUE)
  L <- sum(dpois(I,R0*foi,log=T),na.rm=T)
  return(L)
}