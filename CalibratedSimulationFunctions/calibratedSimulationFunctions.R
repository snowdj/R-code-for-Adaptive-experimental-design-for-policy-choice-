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

  Methods=c("conventional",
                   "optimalhandpicked",
                   "optimalrandom",
                   "optimalhat",
                   "optimalhatrandom",
                   "optimal",
                   "besthalf",
                   "thompson",
                   "modifiedthompson")   
  MethodNames=c("non-adaptive",
                 "handpicked",
                 "random",
                 "estimated U",
                 "estimated U, random set",
                 "optimized",
                 "best half",
                 "Thompson",
                 "modified Thompson")
  shareTreatments=list() #empty list to store vectors of shares assigned to each treatment
      
  for (i in methods) { #pick here which methods to simulate
      sink("status_ExpectedRegret.txt") #status file to track computations
      cat("waves ", NN, "\n", 
          "theta", theta, "\n",
          "method", Methods[i])
      sink()  
    #if (i==9) browser()  
      
    regret2Wave=parSapply(clust, 1:R, function(j) SimulateTWaveDesign(NN,C,theta, Methods[i]))
    regretTable=rbind(get0("regretTable"),
                      tibble(Statistic=paste("Regret, ", MethodNames[i], sep=""),
                                Value=mean(regret2Wave[1,])))
    shareoptTable=rbind(get0("shareoptTable"),
                        tibble(Statistic=paste("Share optimal, ", MethodNames[i], sep=""),
                                Value=mean(regret2Wave[1,]==0)))
    
    shareTreatments[[Methods[i]]]=table(factor(regret2Wave[1,], levels=max(theta)-theta)) #store shares assigned to each treatment for method i
  }
  
 
  stopCluster(clust)
  
  lastrows=tibble(Statistic=c("Units per wave", "Number of treatments"),
                 Value=c(NN[1], k))
  tablecolumns=rbind(regretTable, shareoptTable, lastrows)

  list(tablecolumns = tablecolumns, shareTreatments=shareTreatments)
}



DesignTable=function(DataList,methods,MC_replicates=100,columnames=NULL,filename=NULL) {

  # Run Expected Regret simulations for each value of theta    
  ResultsTemp=map(DataList, function(data) 
                    ExpectedRegret(data$NtN,C=rep(0, length(data$theta)),data$theta,methods,MC_replicates))      
    
  RegretTableTemp=map(ResultsTemp, "tablecolumns") #extract table columns
  shareTreatmentsList=map(ResultsTemp, "shareTreatments") #extract shares assigned to treatments
  
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
  waves=length(DataList[[1]]$NtN)
  if (!is.null(filename)) {
      PrintRegretTable(RegretTable,filename, MC_replicates, waves)
      PrintRegretHistogram(shareTreatmentsList,filename, MC_replicates, waves, columnames)
  }      
}



PrintRegretTable = function(RegretTable,filename, MC_replicates, waves){
    #write_csv(RegretTable, paste("../Figures/", filename, "RegretTable.csv", sep=""))
    labtext=paste(filename,"_", 
                  #paste(NN, collapse="_"),"_",
                  "RegretTable", sep="")
    rows=dim(RegretTable)[1]-2
    cols=dim(RegretTable)[2]
    digs=matrix(c(rep(3, rows*(cols+1)), rep(0,2*(cols+1))), nrow=rows+2, ncol=cols+1, byrow=T) #controling digits
    
    print(xtable(RegretTable, type = "latex",
                 caption=paste(MC_replicates, " replications, ",
                               waves, " waves.", sep=""),
                 label=paste("tab:", labtext, sep=""),
                 digits=digs),
          hline.after = c(-1,0,rows/2,rows,rows+2), #horizontal lines
          file = paste("../Figures/Simulations/", labtext,".tex", sep=""),
          caption.placement = "top",
          latex.environments = "center", #centering the table and caption
          include.rownames=FALSE, #to omit row numbering
          sanitize.text.function=function(x){x}) #to maintain tex strings as intended
}



PrintRegretHistogram=function(shareTreatmentsList,filename, MC_replicates, waves, columnames) {
    for (i in 1:length(columnames)) {
        pathname=paste("../Figures/Simulations/Histograms/",
                       filename,"_", columnames[i],
                       "_RegretHistogram.pdf", sep="")
        
        shareTreatments=shareTreatmentsList[[i]]
        
        histTibble=tibble(regret=as.double(names(shareTreatments[[1]])),
                           non_adaptive= as.double(shareTreatments[[1]] / MC_replicates),
                           modifiedThompson= as.double(shareTreatments[[4]] / MC_replicates)) #careful about coordinates being right when you change methods!

        ggplot(histTibble, aes(x=regret, y=modifiedThompson)) +
            geom_point(color=fillcolor, size=1) +
            geom_segment(aes(x=regret,xend=regret, y=non_adaptive,  yend=modifiedThompson), color=fillcolor, size=.5) +
            scale_y_continuous(limits=c(0,1)) +
            #coord_cartesian(ylim=c(0,1)) +
            theme_light() +
            theme(#panel.background = element_rect(fill = backcolor, colour = NA),
                #panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.ticks.x=element_blank(),
                axis.ticks.y=element_blank(),
                plot.caption=element_text(size=7)) +
            labs(title=paste(columnames[[i]], ", ", waves, " waves", sep=""),
                 x="Regret", y="Share of simulations")
        

        
        ggsave(pathname, width = 4, height = 3)
    }
}




