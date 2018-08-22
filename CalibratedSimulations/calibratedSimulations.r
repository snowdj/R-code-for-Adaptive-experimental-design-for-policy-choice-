#This script implements calibrated simulations, sampling outcomes from the empirical distribution of an actual experiment.
library(parallel)
source("../CoreFunctions/welfareplotsFunctions.R")
source("../CoreFunctions/welfareplotsGraphics.R")

set.seed(12231983)

# Step 1: load experimental data.
# for now: simulating data for sake of getting started.
# include covariates later

loadPseudoData=function(n=100, #number of observations
                        k=3 #number of treatments
                        )
{
    theta=runif(k) #true average potential outcomes
    
    D=sample(1:k, n, replace=TRUE) #random treatment assignment
    Y=simulatedSample(D,theta) #bernoulli draws with probability theta(D)
    thetaSample=tapply(Y, D,mean) #sample mean outcomes - this is our pseudo-true reference parameter
    
    list(D=D, Y=Y, thetaSample=thetaSample)
}
    
# function for drawing outcomes from sampling distribution
simulatedSample=function(D, #treatment vector
                         theta #success probabilities
                         )
{rbinom(length(D), 1, theta[D]) #bernoulli draws with probability theta(D)
}

PolicyChoice=function(A,B,C){
    # policy choice maximizing expected average outcome of optimal treatment
    # A,B vectors with a, b parameters of beta distribution
    # C vector with cost of treatments
    which.max(A /(A+B) - C)
}

Regret=function(D, #vector of treatments across all waves
                Y, #vector of outcomes across all waves
                C, #treatment cost vector
                theta #true parameter vector
                )
{
    k=max(D) #number of treatment arms, assuming we start at 1
    A=tapply(Y, D,sum) + 1 #posterior parameters, starting with Beta(1,1) (uniform) prior
    B=tapply(1-Y, D,sum) + 1
    
    d=PolicyChoice(A,B,C) # policy chosen given experimental data
    
    regret=max(theta)-theta[d] # regret - difference in welfare between chosen and optimal option
    
    c(regret, d)
}



Simulate2WaveDesign=function(N1,N2,C,theta){
  k=length(theta) #number of treatment arms
  
  n1=rep(floor(N1/k),k)
  n1[seq_len(N1-sum(n1))]=n1[1:(N1-sum(n1))]+1 #equitable treatment assignment in period 1
  D1=as.vector(unlist(sapply(1:k, function(i) rep(i,n1[i])))) #this seems odd way to construct D2 - revisit?
  
  
  Y1=simulatedSample(D1,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution
  
  A=tapply(Y1, D1,sum) + 1 #posterior parameters, starting with Beta(1,1) (uniform) prior
  B=tapply(1-Y1, D1,sum) + 1
  
  USimplex=UoverSimplex(A,B,C,N2, Ufunction=U)
  noptimal=USimplex[which.max(USimplex$U), 1:k] #optimal treatment assignment in wave 2
  
  D2=unlist(sapply(1:k, function(i) rep(i,noptimal[i]))) #this seems odd way to construct D2 - revisit?
  Y2=simulatedSample(D2,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution
  
  D=c(D1,D2)
  Y=c(Y1,Y2)
  Regret(D,Y,C,theta)
}



SimulateConventionalDesign=function(N,C,theta){
  k=length(theta) #number of treatment arms
  
  n=rep(floor(N/k),k)
  n[1:(N-sum(n))]=n[1:(N-sum(n))]+1 #equitable treatment assignment in period 1
  D=unlist(sapply(1:k, function(i) rep(i,n[i]))) #this seems odd way to construct D2 - revisit?
  
  Y=simulatedSample(D,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution

  Regret(D,Y,C,theta)
}



ExpectedRegret=function(N1,N2,C,theta,R, filename=NULL){
  #parallelize simulations
    no_cores = detectCores()
    clust = makeCluster(no_cores, type="FORK")  #forking requires mac or Linux OS!

  regret2Wave=parSapply(clust, 1:R, function(i) Simulate2WaveDesign(N1,N2,C,theta))
  regretConventional=parSapply(clust, 1:R, function(i) SimulateConventionalDesign(N1+N2,C,theta))
    stopCluster(clust)
    
  RegretStats=data_frame(Statistic=c("regret, 2 Wave design", 
               "share optimal, 2 Wave design",
               "regret, conventional design",
               "share optimal, conventional design"),
             Value=c(mean(regret2Wave[1,]),
               mean(regret2Wave[1,]==0),
               mean(regretConventional[1,]),
               mean(regretConventional[1,]==0))
             )
  if (!is.null(filename)){
    write_csv(RegretStats, paste("../../Figures/Applications/", filename, "Regret.csv", sep=""))
  }
  
  RegretStats
}



DesignTable=function(N1,N2,ThetaList,R=100,columnames=NULL, filename=NULL) {
  for (i in 1:length(ThetaList)){
    theta=ThetaList[[i]]
    C=rep(0, length(theta))
    simResults=ExpectedRegret(N1,N2,C,theta,R)
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
    #write_csv(RegretTable, paste("../../Figures/Applications/", filename, "RegretTable.csv", sep=""))
    print(xtable(RegretTable, type = "latex"), 
          file = paste("../../Figures/Applications/", filename,"_", N1,"_",N2,"_","RegretTable.tex", sep=""),
          include.rownames=FALSE, #to omit row numbering
          sanitize.text.function=function(x){x}) #to maintain tex strings as intended
  }
}



#List of Thetas to consider
# ThetaList=list(c(.2,.5,.8),
#                c(.6,.5,.4),
#                c(.2,.4,.6,.8))
# DesignTable(N1,N2,ThetaList,R=100,filename="Test")
