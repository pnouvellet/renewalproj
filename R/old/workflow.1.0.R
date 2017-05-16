rm(list=ls(all.names=TRUE))

library(EpiEstim)
source('simulate_poisson.R')

###############################
### simulate incidence  ######
##############################
# # parameters
# N_Week_Est <- 30
# t_sim <- 7*N_Week_Est+1 # time window for R estimation
# N_sim <- 2  # number of simulation, i.e. could be iteration of same place/different R/incidence, or different places altogheter
# # trnamsissibility
# R <-matrix(c(1.2,1.1),N_sim,1) # reproduction number
# Mean_SI <- 15
# CV_SI <- .8
# SI_trunc <- 30
# #initial conditions
# I_start <- matrix(c(10,5),N_sim,1)  #
# 
# # simulate incidence
# S1 <-simulate_poisson(N_Week_Est,N_sim,R,Mean_SI,CV_SI,SI_trunc,I_start)
# 
# #remove weird initial points
# S1$Incidence <- S1$Incidence[,-1]
# 
# #weekly
# S1$WeeklyIncidence <- matrix(NA,N_sim,N_Week_Est)
# for (i in 1:N_sim){
#   S1$WeeklyIncidence[i,] <- colSums(matrix(S1$Incidence[i,],7,N_Week_Est))
# }
# ### plot
# 
# 
# plot(colSums(matrix(S1$Incidence[1,],7,N_Week_Est)),ylim=c(0,max(S1$WeeklyIncidence)))
# lines(colSums(matrix(S1$Incidence[2,],7,N_Week_Est)),type = "p",col='red')
# 
# save.image(file='D.RData')
###########################################################
#######   MCMC #########
###########################################################

load('D.RData')
#  parameters
rep <- 1e4
I <- S1$Incidence
N_sim <- nrow(I)
Mean_SI <- 15
CV_SI <- .8
SI_trunc <- 30
# sigma_Prop <- rep(1.5e-1,N_sim*2)     # proposal variance
sigma_Prop <- c(1.5e-1,1.5e-1,5e-1,5e-1)     # proposal variances
T_prior <-1e2
#window of estimation
Wind_range <- c(5,28)
N_Week_Est <- diff(Wind_range)+1
I_window <- I[,(1+(Wind_range[1]-1)*7):((Wind_range[2])*7)]
I_window_week <- S1$WeeklyIncidence[,(Wind_range[1]):(Wind_range[2])]

source('foi_renewal.R')
source('Likelihood_renewal.R')
source('R_estim_MCMC.R')


res <- R_estim_MCMC(rep=rep, I=I_window, 
                    N_sim=N_sim, Mean_SI=Mean_SI, 
                    CV_SI=CV_SI, SI_trunc=SI_trunc, 
                    sigma_Prop=sigma_Prop,T_prior=T_prior)
  


res$Acc_rate
#########################
# check
burn <- round(rep/3)
Lf <- res$post_logLik[burn:rep,]
thetaf <- res$post_theta[burn:rep,]
# save(thetaf,Acc,Lf,file <- paste0("xxx.RData"))
plot(Lf[,1],ylab='log_likelihood')
for (i in 1:(2*N_sim)) {
  if (i>N_sim){
    plot(thetaf[,i],xlab='iteration',ylab='I0')
    title(paste0('location ',i-2))
  }else{
    plot(thetaf[,i],xlab='iteration',ylab='R')
    title(paste0('location ',i))
  }
    
  if (i<3) lines(c(0,rep),c(R[i],R[i]),col='red')
}


###########################################################
############   projecting cases forward       ###################
###########################################################
N_Week_Forward<-4
N_simul <- 1e3

source('simulate_forward.R')
I_predict <- simulate_forward(SI_Distr=res$SI_Distr, thetaf,
                              N_simul=N_simul, T_prior=T_prior,
                              N_sim=N_sim,N_Week_Forward=N_Week_Forward,
                              N_Week_Est=N_Week_Est)


##################################


# ###########################################################
# ############       Pull results         ###################
# ###########################################################


ResProj<-array(NA,c(3,N_Week_Forward+N_Week_Est,N_sim))
for (i in 1:N_sim){
  A <- I_predict$I_predict_week[,i,]
  ResProj[,,i] <- apply(A,1,function(x) quantile(x,c(.25,.5,.75),na.rm=T))
}


# ###########################################################
# ############       Plot             ###################
# ###########################################################


col.scenarios <- c("Sc1"="#4f81bd","Sc2"="#9bbb59",
                   "Sc3"="#c0504d","Sc3"="#e67e00")
col.scenarios.transp <- paste(col.scenarios,"55",sep="")
col.scenarios.transp2 <- paste(col.scenarios,"35",sep="")

# png("figs/TEST.png",width=800,height=800,res=200)

for (i in 1:N_sim){
  
  
  x <- Wind_range[1]:Wind_range[2]
  
  plot(S1$WeeklyIncidence[i,],
       pch=1,lwd=2,
       col='grey',
       xlim=c(0,length(S1$WeeklyIncidence[i,])+N_Week_Forward+1),
       ylim=c(0,max(S1$WeeklyIncidence)),
       xlab='',ylab='',
       main=paste('Site ',i,sep=''),
       xaxt  = 'n',yaxt='n',bty='n')
  
  axis(side = 1, at = seq(0,length(S1$WeeklyIncidence[i,])+N_Week_Forward+1,by=5))
  axis(side = 2, at = seq(0,max(S1$WeeklyIncidence),length=5))
  mtext("Incidence (weekly)",2,3,outer=FALSE,las=0,adj=-1.5)
  mtext("Weeks", 1, 4, outer=FALSE,padj=-2,adj=.25) 
  
  # time window for estimation
  points(x,I_window_week[i,],
         pch=16,lwd=2,
         col=col.scenarios[i])
  
  
  x2<-seq(Wind_range[1],Wind_range[2]+N_Week_Forward)
  lines(x2,ResProj[2,,i],
        type='l',lwd=2, col=col.scenarios[i])
  
  #95%CrI
  polygon(c(x2,rev(x2)),
          c(ResProj[1,,i],rev(ResProj[3,,i])),
          col = col.scenarios.transp[i],border=NA)
    
    
  
  # legend("topleft",c("Data","Fitted","IQR"),
  #            col=c(col.scenarios[i],col.scenarios[i],
  #                  col.scenarios.transp[i]),
  #        lty=c(-1,1,-1),pch=c(19,-1,15),pt.cex=c(1,1,2),
  #        lwd=c(2,2,1),cex=0.8,bty="n")
  legend("topleft",c("Data","Fitted","95%CrI"),
         col=c(col.scenarios[i],col.scenarios[i],
               col.scenarios.transp[i],
               col.scenarios.transp2[i]),
         lty=c(-1,1,-1,-1),pch=c(19,-1,15,15),pt.cex=c(1,1,2,2),
         lwd=c(2,2,1,1),cex=0.8,bty="n")
  
  
}

apply(thetaf,2,function(x) quantile(x,c(.25,.5,.75),na.rm=T))
