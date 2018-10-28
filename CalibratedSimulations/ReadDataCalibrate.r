rm(list = ls())

source("calibratedSimulations.r")

backcolor="azure2" #background color for plots
gridcolor="azure1"
fillcolor="skyblue4"


DataToTheta=function(filename, dataname, dbar, strataVars, printFigures=FALSE){
  Data=read_csv(paste("../../Datasets/Cleaned/", filename, ".csv", sep=""))
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
    
    ggsave(paste("../../Figures/Applications/", filename, "strata.pdf", sep=""), width = 4, height =1+ .3*length(levels(sumstats$strata)))
    
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
    
    ggsave(paste("../../Figures/Applications/", filename, ".pdf", sep=""), width = 4, height = 0.15*length(levels(sumstats$strata))*dbar+1.5)
  
    write_csv(sumstats, paste("../../Figures/Applications/", filename, "Sumstats.csv", sep=""))
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

   
RunCalibratedSimulations()

