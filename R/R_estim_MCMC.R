#' Estimate of R's (and initial states) by MCMC sampling
#'
#'  gh
#'
#' @author Pierre Nouvellet (\email{p.nouvellet@imperial.ac.uk})
#'
#' @export
#'
#' @param rep (a matrix of integer) daily time-series at the different location
#' This can be any positive number.
#'
#' @param I (a single integer) number of of independent location.
#'
#' @param Mean_SI (a single real positive) mean of the serial interval.
#'
#' @param CV_SI (a single real positive) coefficient of variation of the serial interval.
#'
#' @param SI_trunc (a single real positive) threshold/cut-off for the serial interval calculation  
#'
#' @param sigma_Prop: (a vector of real positive(s) of size 2 x N), from 1:N values of the reproduction number,
#' and from (N+1):(2N) values of the back calculated incdence at time T_prior- first date of time window
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
R_estim_MCMC<-function(rep,I,Mean_SI,CV_SI,SI_trunc,sigma_Prop,T_prior)
{
  
  N_sim <- nrow(I)
  SI_Distr <- vector()
  for(i in 0:SI_trunc) SI_Distr[i+1] <- DiscrSI(i, Mean_SI, CV_SI*Mean_SI)
  SI_Distr<-SI_Distr/sum(SI_Distr)
  # plot(0:SI_trunc,SI_Distr)
  ws <- rev(SI_Distr) # reversed to use as past infectivity
  
  ###########################################################
  ###########################################################
  # Estimation of R0 based on renewal equation with MCMC sampling
  ##################################################
  ############################################################
  
  # useful parameters and results
  
  Acc <- matrix(NA,N_sim,2)        # will contain the acceptance rate for MCMC
  thetas <- matrix(NA,rep,N_sim*2)  # parameters posterior, Rt+I0
  L <- thetas                      # likelihoods

  #############################################################################
  #############################################################################
  # initiallisation
  init<-function(mu){
    I0=matrix(0,1,SI_trunc+50)
    I0[,1]=mu
    for (i in 2:ncol(I0)){
      f=max(c(1,(i-SI_trunc)))
      I0[,i]=(I0[,f:i]%*%ws[((SI_trunc+1)-(i-f)):(SI_trunc+1)])
    }
    return(I0[1,ncol(I0)])
  }
  I0 <- rep(NA,N_sim)
  for (i in 1:N_sim){
    par_opt<-optim(7,function(x) (init(x)-mean(I[i,]))^2,method="Brent",lower=0,upper=1e5)
    I0[i]<-par_opt$par
  }

  
  theta0 <- c(rep(1,N_sim),I0)#rep(1,N.country),rep(30,N.country))
  lambda <- foi_renewal(I=I,N=N_sim,Theta=theta0,rSI=ws,T_prior=T_prior) 
  L1 <- Likelihood_renewal(lambda,I,theta0[1:N_sim])
  L[1,] <- L1
  thetas[1,] <- theta0       
  #############################################################################
  # sampling
  for (i in 2:rep){   
    print(i)
    #print(i)
    for (j in 1:(2*N_sim)){
      Ts <- theta0
      lambdaT <- lambda
      Ts[j] <- Ts[j]*exp(sigma_Prop[j]*rnorm(1,0,1))
      lambdaT <- foi_renewal(I=I,N=N_sim,Theta=Ts,rSI=ws,T_prior=1e2)
      Lint <- Likelihood_renewal(lambdaT,I,Ts[1:N_sim])
      if (j<3){
        r <- exp(Lint-L1)*Ts[j]/theta0[j]
      }else{
        Lint <- Likelihood_renewal(lambdaT,I,Ts[1:N_sim])
        r <- exp(Lint-L1)*Ts[j]/theta0[j]*exp(-Ts[j]/I0[j-N_sim])
      }
      if (runif(1,0,1)<= r){
        theta0[j] <- Ts[j]
        lambda <- lambdaT
        L1 <- Lint
      }  
      L[i,j] <- L1    
    }
    thetas[i,] <- theta0
  }
  Acc <- colSums(diff(thetas)!=0)/rep  
  
  res <- list(post_theta = thetas, post_logLik = L, Acc_rate = Acc, SI_Distr=S1$SI_Distr)
  
  return(res)
}