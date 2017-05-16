#' force of infection for renewal process
#'
#' ~internal function
#'
#' @author Pierre Nouvellet (\email{p.nouvellet@imperial.ac.uk})
#'
#' @export
#'
#' @param I (a matrix of integer) daily time-series at the different location
#' This can be any positive number.
#'
#' @param N (a single integer) number of of independent location.
#'
#' @param Theta: (a vector of real positive(s) of size 2 x N), from 1:N values of the reproduction number,
#' and from (N+1):(2N) values of the back calculated incdence at time T_prior- first date of time window
#'
#' @param rSI (a vector of real positive) reversed serial interval distribution.
#'
#' @param T_prior (a single real positive) number of backward step for which incidence is backcalculated.
#'
#' @return
#'  The function returns sum_{s=(t-SI_trunc)}^{t}(I_{t-s}w_s) for each day with the time window of interest.
#'
#'
#' @examples
#'
#' 
#'
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