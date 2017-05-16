simulate_forward <- function(SI_Distr,thetaf,N_simul,T_prior,N_sim,N_Week_Forward,N_Week_Est)
{
  ###########################################################
  ############         Simulate           ###################
  ###########################################################
  
  I_predict_day <- array(0,c((N_Week_Forward+N_Week_Est)*7,N_sim,N_simul))  # days of simulation+nWeekEst week of estimation
  I_predict_week <- array(0,c((N_Week_Forward+N_Week_Est),N_sim,N_simul))
  
  ws <- rev(SI_Distr)
  SI_trunc <- length(ws)-1
  
  for (k in 1:N_simul){
    fR <- round(runif(1,1,nrow(thetaf)))
    R0 <- thetaf[fR,1:N_sim]
    I0 <- matrix(0,N_sim,T_prior)
    I0[,1] <- thetaf[fR,(N_sim+1):(N_sim*2)]
    
    # deal with pre-inference period
    for (i in 2:T_prior){
      f <- max(c(1,(i-SI_trunc)))
      I0[,i] <- R0*(I0[,f:i]%*%ws[((SI_trunc+1)-(i-f)):(SI_trunc+1)])
    }
    # deal with inference and post-inference period
    I_new <- cbind(I0,matrix(0,N_sim,(N_Week_Forward+N_Week_Est)*7))
    for (i in ((T_prior+1):ncol(I_new))){
      f <- max(c(1,(i-SI_trunc)))
      lambda <- I_new[,f:i]%*%ws[((SI_trunc+1)-(i-f)):(SI_trunc+1)]
      I_new[,i] <- rpois(N_sim,R0*lambda)
    }
    # store simulations
    I_predict_day[,,k] <- t(I_new[,(T_prior+1):ncol(I_new)])
    
    for (i in c(1:N_sim))  I_predict_week[,i,k]=colSums(matrix(I_predict_day[,i,k],7,N_Week_Forward+N_Week_Est))
    
  }
  res <- list(I_predict_day=I_predict_day,I_predict_week=I_predict_week)
  return(res)
}
