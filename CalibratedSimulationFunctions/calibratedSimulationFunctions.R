#This script implements calibrated simulations, sampling outcomes from the empirical distribution of an actual experiment.
library(parallel)
library(tidyverse)
library(forcats)
library(xtable)

set.seed(12231983)
# prior precision
alpha0=1



SimulateTWaveDesign=function(NN,C,theta, method="modifiedthompson"){
  k=length(theta) #number of treatment arms
  theta=sample(theta) #randomly permute so as not to privilege any options in policy choice
  
  if (method=="conventional") {
      D=EqualAssignment(sum(NN),k)
      Y=simulatedSample(D,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution
  } else {
      T=length(NN) #number of waves
      MM=cumsum(NN)
      
      D=rep(0,MM[T])
      Y=rep(0,MM[T])
      
      Dt=EqualAssignment(NN[1],k)
      D[1:NN[1]]=Dt
      Y[1:NN[1]]=simulatedSample(Dt,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution
  
      for (t in 2:T) {
        posterior=betaposterior(D[1:MM[t-1]],Y[1:MM[t-1]])
        Dt=Dtchoice(posterior$A,posterior$B,C, NN[t], method)
        D[(MM[t-1]+1):MM[t]]=Dt
        Y[(MM[t-1]+1):MM[t]]=simulatedSample(Dt,theta) #bernoulli draws with probability theta(D) - drawing from sample distribution
      }
  }
  Regret(D,Y,C,theta)
}


ExpectedRegret=function(NN,C,theta,methods,R){
  k=length(theta) #number of treatment arms

  #parallelize simulations
  no_cores = detectCores()
  clust = makeCluster(no_cores, type="FORK")  #forking requires mac or Linux OS!

  twoWaveMethods=c("conventional",
                   "optimalhandpicked",
                   "optimalrandom",
                   "optimalhat",
                   "optimal",
                   "besthalf",
                   "thompson",
                   "modifiedthompson")   
  twoWaveNames=c("conventional design",
                 "handpicked",
                 "random",
                 "estimated U",
                 "optimized",
                 "rule of thumb",
                 "Thompson",
                 "modified Thompson")

  for (i in methods) { #pick here which methods to simulate
    regret2Wave=parSapply(clust, 1:R, function(j) SimulateTWaveDesign(NN,C,theta, twoWaveMethods[i]))
    regretTable=rbind(get0("regretTable"),
                      tibble(Statistic=paste("regret, T wave ", twoWaveNames[i], sep=""),
                                Value=mean(regret2Wave[1,])))
    shareoptTable=rbind(get0("shareoptTable"),
                        tibble(Statistic=paste("share optimal, T wave ", twoWaveNames[i], sep=""),
                                Value=mean(regret2Wave[1,]==0)))
  }
  
 
  stopCluster(clust)
  
  rbind(regretTable, shareoptTable)

}



DesignTable=function(NN,DataList,methods,R=100,columnames=NULL,filename=NULL) {

  # Run Expected Regret simulations for each value of theta    
  RegretTableTemp=map(DataList, function(data) 
                    ExpectedRegret(NN,C=rep(0, length(data$theta)),data$theta,methods,R))      
    
  # Combine tables
  RegretTable=bind_cols(RegretTableTemp[[1]]["Statistic"], #row labels
                        map(RegretTableTemp, 2)) #keep only simulation values

  #change column names
  if (!is.null(columnames)) {
      colnames(RegretTable)=c("Statistic", columnames)
  } else {
      colnames(RegretTable)=c("Statistic", map(DataList, function(data) 
                          paste("$\\theta = (", toString(data$theta,sep=", "), ")$")))
  }
  
  #write to tex file
  if (!is.null(filename)) PrintRegretTable(RegretTable,NN,filename)
}



PrintRegretTable = function(RegretTable,NN,filename){
    #write_csv(RegretTable, paste("../Figures/", filename, "RegretTable.csv", sep=""))
    labtext=paste(filename,"_", paste(NN, collapse="_"),"_","RegretTable", sep="")
    rows=dim(RegretTable)[1]
    print(xtable(RegretTable, type = "latex",
                 caption=paste("Performance of alternative experimental designs, ",
                               R, " replications, wave sizes ", paste(NN, collapse = ", "), ".", sep=""),
                 label=paste("tab:", labtext, sep=""),
                 digits=3),
          hline.after = c(-1,0,rows/2,rows), #horizontal lines
          file = paste("../Figures/Simulations/", labtext,".tex", sep=""),
          caption.placement = "top",
          latex.environments = "center", #centering the table and caption
          include.rownames=FALSE, #to omit row numbering
          sanitize.text.function=function(x){x}) #to maintain tex strings as intended
}








