#.rs.restartR()
library(tidyverse)
source("ThompsonHierarchical.R")

DataToThetaCovariates=function(filename, dataname, dbar, strataVars){
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
        summarize(meanout=mean(outcome), sumout=sum(outcome), obs=n()) 
    
    #average outcomes by treatment alone
    theta=Data %>% 
        group_by(treatment) %>%
        summarize(meanout=mean(outcome), obs=n())
    
    stratasizes = sumstats %>%
        group_by(strata) %>%
        summarize(n=sum(obs))
 
    # converting to wide matrices of successes and trials, dropping NA strata
    nstrata=nrow(key)
    
    SS= sumstats  %>% #number of successes, treatments by row, strata by column
        select(treatment, strata, sumout) %>%
        spread(strata, sumout) %>%
        ungroup %>% 
        select(paste(1:nstrata)) %>%
        data.matrix()
    
    NN= sumstats  %>% #number of trials, treatments by row, strata by column
        select(treatment, strata, obs) %>%
        spread(strata, obs) %>%
        ungroup %>% 
        select(paste(1:nstrata)) %>%
        data.matrix()
    
    NX= stratasizes %>% #number of units per stratum
        slice(1:nstrata) %>%
        select(n) %>%
        data.matrix()
    
    PX=NX/sum(NX)
    
    list(SS=SS, NN=NN, PX=PX)
}






ReadAllDataThompson=function(){
    DataList=list()
    columnames=list()
    
    for (application in 1:2){
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
        }
        
        #produce figures and get Thetas
        DataList[[application]]=DataToThetaCovariates(filename, dataname, dbar, strataVars)
        columnames[[application]]=filename
    }
    
    # #parameters of hypothetical experiment
    # NN=rep(24,12)
    # R=4000
    # DesignTable(NN,ThetaList,R,columnames,"Applications/CalibratedSimulations")
    DataList
}





#TBC here!!!
SimulateTWaveDesignThompson=function(NN,C,theta, PX){
    k=length(C) #number of treatment arms
    
    T=length(NN) #number of waves
    MM=cumsum(NN)
    
    X=SimulateX(PX,MM[T])
    D=rep(0,MM[T])
    Y=rep(0,MM[T])
    
    Dt=EqualAssignment(NN[1],k) #Do we really want to start like this? better to stratify!
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




DataList=ReadAllDataThompson()
print(DataList)






