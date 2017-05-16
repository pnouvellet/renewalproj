#' simulate renewal process
#'
#'
#'  simulate poisson renewal based process in independent location
#'
#'
#' @author Pierre Nouvellet (\email{p.nouvellet@imperial.ac.uk})
#'
#' @export
#'
#' @param N_Week_Est (a single integer) number of week of simulation. the simulation output is daily numbers
#' This can be any positive number.
#'
#' @param N_location (a single integer) number of of independent location.
#'
#' @param R: (a vector of real positive(s) of size N_location) the reproduction number at each location.
#'
#' @param Mean_SI (a single real positive) mean of the serial interval.
#'
#' @param CV_SI (a single real positive) coefficient of variation of the serial interval.
#'
#' @param SI_trunc (a single real positive) threshold/cut-off for the serial interval calculation  
#'
#' @param I_start: (a vector of integer(s) of size N_location) initial size of the outbreak(s). 
#'Tthe output excludes those and return simulated incidence from the following day
#'
#'
#' @return
#'  The function returns simulated time series. the list return include $Incidence: the daily incidence, 
#'  and $SI_Distr  daily value of the serial interval from the same day of symptoms onset up to 
#'  SI_trunc+1 after the day of symptoms onset
#'
#'
#' @examples
#'
#' 
#'
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