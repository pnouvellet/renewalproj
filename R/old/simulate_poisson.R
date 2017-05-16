simulate_poisson <- function(N_Week_Est,N_location,R,Mean_SI,CV_SI,SI_trunc,I_start)
{
  t_sim <- 7*N_Week_Est+1
  SI_Distr <- vector()
  for(i in 0:SI_trunc) SI_Distr[i+1] <- DiscrSI(i, Mean_SI, CV_SI*Mean_SI)
  SI_Distr<-SI_Distr/sum(SI_Distr)
  # plot(0:SI_trunc,SI_Distr)
  ws <- rev(SI_Distr) # reversed to use as past infectivity
  
  # simulate incidence
  I <- matrix(0,N_location,t_sim) # declare simulated incidence
  I[,1:ncol(I_start)] <- I_start
  for (i in (ncol(I_start)+1):t_sim){
    f=max(c(1,(i-SI_trunc))) # truncated the incidence according to SI truncation
    lambda <- I[,f:i]%*%ws[((SI_trunc+1)-(i-f)):(SI_trunc+1)]
    I[,i]=rpois(N_location,R*lambda)
  }
  # plot(I[1,])
  # lines(I[2,],type = "p",col='red')
  return(list(Incidence=I,SI_Distr=SI_Distr))
}