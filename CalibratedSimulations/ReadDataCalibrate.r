rm(list = ls())
library(tidyverse)
library(forcats)
source("calibratedSimulations.r")

backcolor="azure2" #background color for plots
fillcolor="skyblue4"


DataToTheta=function(filename, dataname, dbar, strataVars){
  Data=read_csv(paste("../../Datasets/Cleaned/", filename, ".csv", sep=""))
  head(Data)
  #check for missings?
  
  treatDummies=paste("treatment",1:dbar, sep="") #names of treatment variables
  
  Data=Data %>%
    mutate(treatment=factor(as.matrix(Data[treatDummies])%*%(1:dbar))) %>%
    mutate(Strata=(interaction(select(.,strataVars)))) #this is useful for renaming levels but introduces randomness...
  
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
  
  ggplot(stratasizes, aes(x=factor(strata, levels = rev(levels(strata))), y=n)) +
    geom_col(fill=fillcolor, width=.5) + 
    theme(panel.background = element_rect(fill = backcolor, colour = NA)) +
    coord_flip() +
    labs(title="Strata size",
         subtitle=dataname,
         x="stratum", y="number of observations")
  
  ggsave(paste("../../Figures/Applications/", filename, "strata.pdf", sep=""), width = 5, height =1.5+ .2*length(levels(sumstats$strata)))
  
  
  ggplot(sumstats, aes(x=meanout, y=factor(treatment,levels=c(dbar:1)))) +
    geom_point(color=fillcolor) +
    #scale_y_discrete(labels=paste("treatment", 1:dbar)) + 
    facet_grid(strata~.) +
    scale_x_continuous(limits=c(0,1)) +
    theme(panel.background = element_rect(fill = backcolor, colour = NA)) +
    labs(title="Average outcomes by treatment and stratum",
         subtitle=dataname,
         x="outcome", y="treatment")
  
  ggsave(paste("../../Figures/Applications/", filename, ".pdf", sep=""), width = 5, height = 0.15*length(levels(sumstats$strata))*dbar+1.5)

  write_csv(sumstats, paste("../../Figures/Applications/", filename, "Sumstats.csv", sep=""))
  
  list(theta=theta$meanout, sumstats=sumstats, stratasizes =stratasizes, key=key)
}


RunCalibratedSimulations=function(){
  ThetaList=list()
  columnames=list()
  
  for (application in 1:2){
    #parameters for each simulation
    if (application==1){
      filename="Ashraf"
      dataname="Ashraf, Berry, and Shapiro (2010)" #,\n  \"Can Higher Prices Stimulate Product Use?  Evidence from a Field Experiment in Zambia.\""
      dbar=6 #number of treatment values
      strataVars=c("covar1", "covar2") #variables to stratify on. we might want to do be careful about keeping track of strata meaning, though
    } else if (application==2) {
      filename="Bryan"
      dataname="Bryan, Chowdhury, and Mobarak (2014)" #,\n \"Underinvestment in a Profitable Technology: The Case of Seasonal Migration in Bangladesh\""
      dbar=4
      strataVars=c("covar1", "covar2")
    }
    
    #produce figures and get Thetas
    DataSummary=DataToTheta(filename, dataname, dbar, strataVars)
    
    ThetaList[[application]]=DataSummary$theta
    columnames[[application]]=filename
  }
  
  #parameters of hypothetical experiment
  N1=12
  N2=12
  R=20
  DesignTable(N1,N2,ThetaList,R,columnames,"CalibratedSimulations")
  
}

   
RunCalibratedSimulations()

