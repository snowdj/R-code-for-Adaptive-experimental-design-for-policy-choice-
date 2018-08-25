#This script implements calibrated simulations, sampling outcomes from the empirical distribution of an actual experiment.
library(parallel)
library(tidyverse)
library(forcats)
library(xtable)

source("../CoreFunctions/welfareplotsFunctions.R")
source("../CoreFunctions/welfareplotsGraphics.R")

set.seed(12231983)
# prior precision
alpha0=3
# how many random subpoints of simplex to consider in optimization:
RR=400 #drop this line if full optimization is desired

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
    A=tapply(Y, D,sum) + alpha0 #posterior parameters, starting with Beta(alpha0,alpha0) (uniform) prior
    B=tapply(1-Y, D,sum) + alpha0
    
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
  

  A=tapply(Y1, D1,sum) + alpha0 #posterior parameters, starting with Beta(alpha0,alpha0) prior
  B=tapply(1-Y1, D1,sum) + alpha0

  USimplex=UoverSimplex(A,B,C,N2, Ufunction=U,RR)
  noptimal=USimplex[which.max(USimplex$U), 1:k] #optimal treatment assignment in wave 2
  
  D2=as.vector(unlist(sapply(1:k, function(i) rep(i,noptimal[i])))) #this seems odd way to construct D2 - revisit?
  Y2=simulatedSample(D2,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution
  
  D=c(D1,D2)
  Y=c(Y1,Y2)
  Regret(D,Y,C,theta)
}


SimulateRuleOfThumb=function(N1,N2,C,theta){
  k=length(theta) #number of treatment arms
  
  n1=rep(floor(N1/k),k)
  n1[seq_len(N1-sum(n1))]=n1[1:(N1-sum(n1))]+1 #equitable treatment assignment in period 1
  D1=as.vector(unlist(sapply(1:k, function(i) rep(i,n1[i])))) #this seems odd way to construct D2 - revisit?
  
  Y1=simulatedSample(D1,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution
  
  # MLE of theta
  thetahat=tapply(Y1, D1,mean)
  #get ordered list of treatments, starting with the highest MLE
  bestoptions=order(thetahat,decreasing=TRUE)
  k2=ceiling(k/2) #consider half the options in second round. maybe do something more sophisticated here?
  
  #start with n2 vector as if options were ordered, then reorder afterwards
  n2=c(rep(floor(N2/k2),k2), rep(0,k-k2))
  n2[seq_len(N2-sum(n2))]=n2[1:(N2-sum(n2))]+1
  n2[bestoptions]=n2
  
  D2=as.vector(unlist(sapply(1:k, function(i) rep(i,n2[i])))) #this seems odd way to construct D2 - revisit?
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
  regretRuleOfThumb=parSapply(clust, 1:R, function(i) SimulateRuleOfThumb(N1,N2,C,theta))
  regretConventional=parSapply(clust, 1:R, function(i) SimulateConventionalDesign(N1+N2,C,theta))
    stopCluster(clust)
    
  RegretStats=data_frame(Statistic=c("regret, 2 wave optimized",
                                     "regret, 2 wave rule of thumb",
                                     "regret, conventional design",
                                     "share optimal, 2 wave optimized",
                                     "share optimal, 2 wave rule of thumb",
                                     "share optimal, conventional design"),
             Value=c(mean(regret2Wave[1,]),
                     mean(regretRuleOfThumb[1,]),
                     mean(regretConventional[1,]),
                     mean(regret2Wave[1,]==0),
                     mean(regretRuleOfThumb[1,]==0),
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
    #this is where the action happens:
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
    labtext=paste(filename,"_", N1,"_",N2,"_","RegretTable", sep="")
    print(xtable(RegretTable, type = "latex",
                 caption=paste("Performance of alternative experimental designs, ",
                               "$N_1=", N1,
                               "$, $N_2=", N2, "$, ",
                               R, " replications"),
                 label=paste("tab:", labtext, sep=""),
                 digits=3),
          hline.after = c(-1,0,3,6), #horizontal lines
          file = paste("../../Figures/Applications/", labtext,".tex", sep=""),
          caption.placement = "top",
          latex.environments = "center", #centering the table and caption
          include.rownames=FALSE, #to omit row numbering
          sanitize.text.function=function(x){x}) #to maintain tex strings as intended
  }
}



#List of Thetas to consider
# N1=12
# N2=9
# ThetaList=list(c(.2,.5,.8),
#                c(.6,.5,.4),
#                c(.2,.4,.6,.8))
# DesignTable(N1,N2,ThetaList,R=400,filename="Test")
