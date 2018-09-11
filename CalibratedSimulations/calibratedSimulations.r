#This script implements calibrated simulations, sampling outcomes from the empirical distribution of an actual experiment.
library(parallel)
library(tidyverse)
library(forcats)
library(xtable)

source("../CoreFunctions/welfareplotsFunctions.R")
source("../CoreFunctions/welfareplotsGraphics.R")
source("../CoreFunctions/SimulatedFunctions.R")

set.seed(12231983)
# prior precision
alpha0=1

# # Step 1: load experimental data.
# # for now: simulating data for sake of getting started.
# # include covariates later
# 
# loadPseudoData=function(n=100, #number of observations
#                         k=3 #number of treatments
#                         )
# {
#     theta=runif(k) #true average potential outcomes
#     
#     D=sample(1:k, n, replace=TRUE) #random treatment assignment
#     Y=simulatedSample(D,theta) #bernoulli draws with probability theta(D)
#     thetaSample=tapply(Y, D,mean) #sample mean outcomes - this is our pseudo-true reference parameter
#     
#     list(D=D, Y=Y, thetaSample=thetaSample)
# }



# Simulate2WaveDesign=function(N1,N2,C,theta, method="optimal"){
#   k=length(theta) #number of treatment arms
#   
#   D1=EqualAssignment(N1,k)
#   Y1=simulatedSample(D1,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution
#   
#   posterior=betaposterior(D1,Y1)
# 
#   D2=D2choice(posterior$A,posterior$B,C, N2, method)
# 
#   Y2=simulatedSample(D2,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution
#   
#   Regret(c(D1,D2),c(Y1,Y2),C,theta)
# }



SimulateTWaveDesign=function(NN,C,theta, method="modifiedthompson"){
  k=length(theta) #number of treatment arms
  theta=sample(theta) #randomly permute so as not to privilege any options in policy choice
  
  T=length(NN) #number of waves
  MM=cumsum(NN)
  
  D=rep(0,MM[T])
  Y=rep(0,MM[T])
  
  Dt=EqualAssignment(NN[1],k)
  D[1:NN[1]]=Dt
  Y[1:NN[1]]=simulatedSample(Dt,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution
  
  for (t in 2:T) {
    posterior=betaposterior(D[1:MM[t-1]],Y[1:MM[t-1]])
    Dt=D2choice(posterior$A,posterior$B,C, NN[t], method)
    D[(MM[t-1]+1):MM[t]]=Dt
    Y[(MM[t-1]+1):MM[t]]=simulatedSample(Dt,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution
  }
  Regret(D,Y,C,theta)
}




SimulateConventionalDesign=function(N,C,theta){
  k=length(theta) #number of treatment arms
  theta=sample(theta) #randomly permute so as not to privilege any options in policy choice
  
  D=EqualAssignment(N,k)
  
  Y=simulatedSample(D,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution

  Regret(D,Y,C,theta)
}



ExpectedRegret=function(NN,C,theta,R, filename=NULL){
  k=length(theta) #number of treatment arms

      #parallelize simulations
      no_cores = detectCores()
      clust = makeCluster(no_cores, type="FORK")  #forking requires mac or Linux OS!
    
  regretConventional=parSapply(clust, 1:R, function(i) SimulateConventionalDesign(sum(NN),C,theta))
  regretFrame=data_frame(Statistic="regret, conventional design",
                         Value=mean(regretConventional[1,]))
  shareoptFrame=data_frame(Statistic="share optimal, conventional design",
                           Value=mean(regretConventional[1,]==0))
  csvStats=c(0,
             mean(regretConventional[1,]),
             colMeans(outer(regretConventional[2,],(1:k), "==")),
             theta)
  
    
  twoWaveMethods=c("optimalhandpicked",
                   "optimalrandom",
                   "optimalhat",
                   "optimal",
                   "besthalf",
                   "thompson",
                   "modifiedthompson")   
  twoWaveNames=c("handpicked",
                 "random",
                 "estimated U",
                 "optimized",
                 "rule of thumb",
                 "Thompson",
                 "modified Thompson")
  for (i in c(5,6,7,1)) { #pick here which methods to simulate
    #browser()
    regret2Wave=parSapply(clust, 1:R, function(j) SimulateTWaveDesign(NN,C,theta, twoWaveMethods[i]))
    regretFrame=rbind(regretFrame,
                      data_frame(Statistic=paste("regret, T wave ", twoWaveNames[i], sep=""),
                      Value=mean(regret2Wave[1,])))
    shareoptFrame=rbind(shareoptFrame,
                        data_frame(Statistic=paste("share optimal, T wave ", twoWaveNames[i], sep=""),
                        Value=mean(regret2Wave[1,]==0)))
    
    csvStats=rbind(csvStats, c(i,
                               mean(regret2Wave[1,]),
                               colMeans(outer(regret2Wave[2,],(1:k), "==")),
                               theta))
#can move these outside loop; but useful for storing intermediate results:    

    # if (!is.null(filename)){
    #   write_csv(as.data.frame(csvStats), paste("../../Figures/Simulations/", filename, "Regret.csv", sep=""))
    # }
  }
  
 
  stopCluster(clust)
  
  rbind(regretFrame, shareoptFrame)
}



DesignTable=function(NN,ThetaList,R=100,columnames=NULL, filename=NULL) {
  for (i in 1:length(ThetaList)){
    theta=ThetaList[[i]]
    C=rep(0, length(theta))
    #this is where the action happens:
    simResults=ExpectedRegret(NN,C,theta,R, paste(filename,i,sep=""))
    
    if (is.null(columnames)){
      colnames(simResults)[2]=paste("$\\theta = (", toString(theta,sep=", "), ")$")
    } else {
      colnames(simResults)[2]=columnames[i]
    }
    if (i==1){
      RegretTable=simResults
    } else {
      RegretTable= RegretTable %>%
        bind_cols(simResults[2])
    }
  }
  
  if (!is.null(filename)){
    #write_csv(RegretTable, paste("../../Figures/", filename, "RegretTable.csv", sep=""))
    labtext=paste(filename,"_", paste(NN, collapse="_"),"_","RegretTable", sep="")
    rows=dim(RegretTable)[1]
    print(xtable(RegretTable, type = "latex",
                 caption=paste("Performance of alternative experimental designs, ",
                               R, " replications, wave sizes ", paste(NN, collapse = ", "), ".", sep=""),
                 label=paste("tab:", labtext, sep=""),
                 digits=3),
          hline.after = c(-1,0,rows/2,rows), #horizontal lines
          file = paste("../../Figures/Simulations/", labtext,".tex", sep=""),
          caption.placement = "top",
          latex.environments = "center", #centering the table and caption
          include.rownames=FALSE, #to omit row numbering
          sanitize.text.function=function(x){x}) #to maintain tex strings as intended
  }
}



#List of Thetas to consider
# N1=16
# N2=16
 R= 4000
# ThetaList=list(c(.2,.5,.8),
#                c(.6,.5,.4),
#                c(.2,.4,.6,.8),
#                c(.2,.4,.5,.6,.7,.8))
# DesignTable(N1,N2,ThetaList,R,filename="Test")

#another list
# NN=rep(36,2)
# ThetaList=list(seq(.2,length=8, by=.1),
#                seq(.4,length=8, by=.05),
#                seq(.4,length=8, by=.025),
#                seq(.4,length=8, by=.01))
# columnames=c("space .1", "space .05", "space .025", "space .01")
# DesignTable(NN,ThetaList,R,columnames,filename="8options")


#another list
NN=rep(36,2)
ThetaList=list(seq(.5,length=6, by=.05),
               c(.6, .59, .36, .35))
columnames=c("like Ashraf", "like Bryan")
DesignTable(NN,ThetaList,R,columnames,filename="almostEmpirical")
