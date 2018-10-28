#This script implements calibrated simulations, sampling outcomes from the empirical distribution of an actual experiment.
library(parallel)
library(tidyverse)
library(forcats)
library(xtable)

set.seed(12231983)
# prior precision
alpha0=1


#TBC here!!!
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
    Dt=Dtchoice(posterior$A,posterior$B,C, NN[t], method)
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
  for (i in c(5,6,7)) { #pick here which methods to simulate
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
    #   write_csv(as.data.frame(csvStats), paste("../Figures/Simulations/", filename, "Regret.csv", sep=""))
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
}







DataToTheta=function(filename, dataname, dbar, strataVars, printFigures=FALSE){
    Data=read_csv(paste("../Datasets/Cleaned/", filename, ".csv", sep=""))
    head(Data)
    #check for missings?
    
    treatDummies=paste("treatment",1:dbar, sep="") #names of treatment variables
    
    Data=Data %>%
        mutate(treatment=factor(as.matrix(Data[treatDummies])%*%(1:dbar))) %>%
        mutate(Strata=(interaction(select(.,strataVars))))
    
    #recoding the levels in a reproducible way    
    oldlevels=levels(Data$Strata)  
    key=data_frame(Strata=oldlevels, strata=factor(1:length(oldlevels)))
    Data=Data %>%
        left_join(., key, by = "Strata") %>%
        select(-Strata)
    
    #average outcomes by treatment and stratum
    sumstats=Data %>% 
        group_by(treatment, strata) %>%
        summarize(meanout=mean(outcome), obs=n()) 
    
    #average outcomes by treatment alone
    theta=Data %>% 
        group_by(treatment) %>%
        summarize(meanout=mean(outcome), obs=n())
    
    stratasizes = sumstats %>%
        group_by(strata) %>%
        summarize(n=sum(obs))
    
    if (printFigures) {
        #careful: we are dropping "missing" strata from figures!
        ggplot(drop_na(stratasizes),aes(x=factor(strata, levels = rev(levels(strata))), y=n)) +
            geom_col(fill=fillcolor, width=.5) + 
            theme_light() +
            theme(#panel.background = element_rect(fill = backcolor, colour = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) +
            coord_flip() +
            labs(title="Strata size",
                 subtitle=dataname,
                 x="stratum", y="number of observations")
        
        ggsave(paste("../Figures/Applications/", filename, "strata.pdf", sep=""), width = 4, height =1+ .3*length(levels(sumstats$strata)))
        
        xmax=1
        if (max(sumstats$meanout) < .5) xmax= max(sumstats$meanout)+.02
        ggplot(drop_na(sumstats), aes(x=meanout, y=factor(treatment,levels=c(dbar:1)))) +
            geom_point(color=fillcolor) +
            geom_segment(aes(x=0, xend=meanout, yend=factor(treatment,levels=c(dbar:1))), color=fillcolor, size=.5) +
            #scale_y_discrete(labels=paste("treatment", 1:dbar)) + 
            facet_grid(strata~.) +
            scale_x_continuous(limits=c(0,xmax)) +
            theme_light() +
            theme(#panel.background = element_rect(fill = backcolor, colour = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.ticks.x=element_blank()) +
            labs(title="Average outcomes by treatment and stratum",
                 subtitle=dataname,
                 x="mean outcome", y="treatment")
        
        ggsave(paste("../Figures/Applications/", filename, ".pdf", sep=""), width = 4, height = 0.15*length(levels(sumstats$strata))*dbar+1.5)
        
        write_csv(sumstats, paste("../Figures/Applications/", filename, "Sumstats.csv", sep=""))
    }
    
    list(theta=theta$meanout, sumstats=sumstats, stratasizes =stratasizes, key=key)
}


RunCalibratedSimulations=function(){
    ThetaList=list()
    columnames=list()
    
    for (application in 1:3){
        #parameters for each simulation
        if (application==1){
            filename="Ashraf"
            dataname="Ashraf, Berry, and Shapiro (2010)" #,\n  \"Can Higher Prices Stimulate Product Use?  Evidence from a Field Experiment in Zambia.\""
            dbar=6 #number of treatment values
            strataVars=c("covar1") #, "covar2") #variables to stratify on. we might want to do be careful about keeping track of strata meaning, though
        } else if (application==2) {
            filename="Bryan"
            dataname="Bryan, Chowdhury, and Mobarak (2014)" #,\n \"Underinvestment in a Profitable Technology: The Case of Seasonal Migration in Bangladesh\""
            dbar=4
            strataVars=c("covar1") #, "covar2")
        } else if (application==3) {
            filename="KarlanList2"
            dataname="Karlan and List (2007)" #,\n \"Underinvestment in a Profitable Technology: The Case of Seasonal Migration in Bangladesh\""
            dbar=13
            strataVars=c("covar1")
        }
        
        #produce figures and get Thetas
        DataSummary=DataToTheta(filename, dataname, dbar, strataVars, printFigures=T)
        
        ThetaList[[application]]=DataSummary$theta
        columnames[[application]]=filename
    }
    
    #parameters of hypothetical experiment
    NN=rep(24,12)
    R=4000
    #DesignTable(NN,ThetaList,R,columnames,"Applications/CalibratedSimulations")
    
}



#List of Thetas to consider
# N1=16
# N2=16
# R= 4000
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
# NN=rep(12,10)
# ThetaList=list(seq(.5,length=6, by=.05),
#                c(.6, .59, .36, .35))
# columnames=c("like Ashraf", "like Bryan")
# DesignTable(NN,ThetaList,R,columnames,filename="almostEmpirical")
