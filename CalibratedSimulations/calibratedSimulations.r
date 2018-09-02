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
alpha0=3

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


#method may be "optimal", "besthalf", or thompson
Simulate2WaveDesign=function(N1,N2,C,theta, method="optimal"){
  k=length(theta) #number of treatment arms
  
  D1=EqualAssignment(N1,k)
  Y1=simulatedSample(D1,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution
  
  posterior=betaposterior(D1,Y1)

  D2=D2choice(posterior$A,posterior$B,C, N2, method)

  Y2=simulatedSample(D2,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution
  
  Regret(c(D1,D2),c(Y1,Y2),C,theta)
}





SimulateConventionalDesign=function(N,C,theta){
  k=length(theta) #number of treatment arms
  D=EqualAssignment(N,k)
  
  Y=simulatedSample(D,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution

  Regret(D,Y,C,theta)
}



ExpectedRegret=function(N1,N2,C,theta,R, filename=NULL){

    
    # # chunk for debugging:
    #   browser() #invoke command line for debugging
    # 
    #   Rprof()
    #   R=1
    #   sapply(1:R, function(i) Simulate2WaveDesign(N1,N2,C,theta, "optimal"))
    #   summaryRprof()
    #   Rprof()
    #   sapply(1:R, function(i) Simulate2WaveDesign(N1,N2,C,theta, "optimalhat"))
    #   summaryRprof()
    #   Rprof()
    #   sapply(1:R, function(i) Simulate2WaveDesign(N1,N2,C,theta, "besthalf"))
    #   summaryRprof()
    #   Rprof()
    #   sapply(1:R, function(i) Simulate2WaveDesign(N1,N2,C,theta, "thompson"))
    #   summaryRprof()
      
# notes on timing:
      # for optimal, most time spent in betabinomial
      # for optimalhat, most time in Vfunction
      # the latter scales better for larger V2
      # but both take very long
  
      # time for optimalhandpicked, optimalrandom, betterhalf, Thompson, is negligible
      # conundrum: why is calculation much faster in app?     
      
      #parallelize simulations
      no_cores = detectCores()
      clust = makeCluster(no_cores, type="FORK")  #forking requires mac or Linux OS!
    
  regretConventional=parSapply(clust, 1:R, function(i) SimulateConventionalDesign(N1+N2,C,theta))
      stopCluster(clust)      
      
      
  regret2Wave=parSapply(clust, 1:R, function(i) Simulate2WaveDesign(N1,N2,C,theta, "optimal"))
  regret2WaveHandpicked=parSapply(clust, 1:R, function(i) Simulate2WaveDesign(N1,N2,C,theta, "optimalhandpicked"))
  regret2WaveRandom=parSapply(clust, 1:R, function(i) Simulate2WaveDesign(N1,N2,C,theta, "optimalrandom"))
  regret2WaveHat=parSapply(clust, 1:R, function(i) Simulate2WaveDesign(N1,N2,C,theta, "optimalhat"))
  
  regretRuleOfThumb=parSapply(clust, 1:R, function(i) Simulate2WaveDesign(N1,N2,C,theta, "besthalf"))
  regretThompson=parSapply(clust, 1:R, function(i) Simulate2WaveDesign(N1,N2,C,theta, "thompson"))

    
  RegretStats=data_frame(Statistic=c("regret, 2 wave optimized",
                                     "regret, 2 wave handpicked",
                                     "regret, 2 wave random",
                                     "regret, 2 wave estimated U",
                                     
                                     "regret, 2 wave rule of thumb",
                                     "regret, 2 wave, thompson",
                                     "regret, conventional design",
                                     
                                     "share optimal, 2 wave optimized",
                                     "share optimal, 2 wave handpicked",
                                     "share optimal, 2 wave random",
                                     "share optimal, 2 wave estimated U",
                                     
                                     "share optimal, 2 wave rule of thumb",
                                     "share optimal, 2 wave, thompson",
                                     "share optimal, conventional design"),
             Value=c(mean(regret2Wave[1,]),
                     mean(regret2WaveHandpicked[1,]),
                     mean(regret2WaveRandom[1,]),
                     mean(regret2WaveHat[1,]),
                     
                     mean(regretRuleOfThumb[1,]),
                     mean(regretThompson[1,]),
                     mean(regretConventional[1,]),
                     
                     
                     mean(regret2Wave[1,]==0),
                     mean(regret2WaveHandpicked[1,]==0),
                     mean(regret2WaveRandom[1,]==0),
                     mean(regret2WaveHat[1,]==0),
                     
                     mean(regretRuleOfThumb[1,]==0),
                     mean(regretThompson[1,]==0),
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
    rows=dim(RegretTable)[1]
    print(xtable(RegretTable, type = "latex",
                 caption=paste("Performance of alternative experimental designs, ",
                               "$N_1=", N1,
                               "$, $N_2=", N2, "$, ",
                               R, " replications"),
                 label=paste("tab:", labtext, sep=""),
                 digits=3),
          hline.after = c(-1,0,rows/2,rows), #horizontal lines
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
