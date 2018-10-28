SimulateTWaveDesignThompson=function(Nt,C,theta, PX){
    k=length(C) #number of treatment arms
    nx = length(PX)
    
    T=length(Nt) #number of waves
    MM=cumsum(Nt)
    
    #X=SimulateX(PX,MM[T])  #random sampling from distribution PX
    X=1:MM[T] %% nx + 1 #making life easier by "stratified sampling" from X
    
    D=rep(0,MM[T])
    Y=rep(0,MM[T])
    
    Dt=StratifiedAssignment(X[1:Nt[1]], k, nx)
    Xt=X[1:Nt[1]]
    D[1:Nt[1]]=Dt
    Y[1:Nt[1]]=SimulateY(theta, Dt, Xt)
    
    for (t in seq(2, length=max(0,T-1))) {
        previous=1:MM[t-1]
        current=(MM[t-1]+1):MM[t]
        Dt=DtchoiceThompsonHierarchicalAveraged(Y[previous], D[previous], X[previous],
                                        k,nx, X[current], RR=5)
        D[current]=Dt
        Y[current]=SimulateY(theta, Dt, X[current]) 
    }
    
    thetahat=hierarchicalPosteriorMean(Y,D,X)
    
    #careful here - giving priority to lower indices...
    Dstar=apply(thetahat,2, which.max)
    
    regretX=apply(theta,2, max) - theta[cbind(Dstar, 1:nx)]
    
    list(X=X, D=D, Y=Y, thetahat=thetahat, Dstar=Dstar, regretX=regretX)
}




DataToThetaCovariates=function(filename, dataname, k, strataVars){
    Data=read_csv(paste("../Datasets/Cleaned/", filename, ".csv", sep=""))
    head(Data)
    #check for missings?
    
    treatDummies=paste("treatment",1:k, sep="") #names of treatment variables
    
    Data=Data %>%
        mutate(treatment=factor(as.matrix(Data[treatDummies])%*%(1:k))) %>%
        mutate(Strata=(interaction(select(.,strataVars)))) 
    
    #recoding the levels in a reproducible way    
    oldlevels=levels(Data$Strata)  
    key=data_frame(Strata=oldlevels, strata=factor(1:length(oldlevels)))
    Data=Data %>%
        left_join(., key, by = "Strata") %>%
        select(-Strata)
    browser()
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
    
    list(SS=SS, NN=NN, PX=PX, k=k)
}






ReadAllDataThompson=function(){
    DataList=list()
    columnames=list()
    
    for (application in 1:3){
        #parameters for each simulation
        if (application==1){
            filename="Ashraf"
            dataname="Ashraf, Berry, and Shapiro (2010)" #,\n  \"Can Higher Prices Stimulate Product Use?  Evidence from a Field Experiment in Zambia.\""
            k=6 #number of treatment values
            strataVars=c("covar1") #, "covar2") #variables to stratify on. we might want to do be careful about keeping track of strata meaning, though
        } else if (application==2) {
            filename="Bryan"
            dataname="Bryan, Chowdhury, and Mobarak (2014)" #,\n \"Underinvestment in a Profitable Technology: The Case of Seasonal Migration in Bangladesh\""
            k=4
            strataVars=c("covar1") #, "covar2")
        }  else if (application==3) {
          filename="KarlanList1"
          dataname="Karlan and List (2007)" #,\n \"Underinvestment in a Profitable Technology: The Case of Seasonal Migration in Bangladesh\""
          k=4
        }
        
        #produce figures and get Thetas
        DataList[[application]]=DataToThetaCovariates(filename, dataname, k, strataVars)
        columnames[[application]]=filename
    }
    

    DataList
}






RunAllSimulationsThompson=function(T = 4, #number of waves
                                   nt = 36, #units per wave
                                   RR = 1000){ #number of replications for simulation
    DataList=ReadAllDataThompson()
    print(DataList)
    
    simRegret=matrix(0,2,2) #columns for applications, rows for methods
    
    for (application in 1:2){
        theta=DataList[[application]]$SS / DataList[[application]]$NN
        PX=DataList[[application]]$PX
        #note: need enough units to observe each treatment/covariate combo in first round, for MLE

        C=rep(0,DataList[[application]]$k)
        
        for (r in 1:RR) {
            #note 2: permute theta rows for fair comparisons
            
            # Thompson design
            simDesign=SimulateTWaveDesignThompson(Nt = rep(nt,T),C,theta, PX)
            # 1 wave comparison with the same number of units. stratified assignment.
            simDesign1wave=SimulateTWaveDesignThompson(Nt = nt * T,C,theta, PX)
            
            simRegret[1,application] = simRegret[1,application] +  simDesign$regretX %*% PX
            simRegret[2,application] = simRegret[2,application] +  simDesign1wave$regretX %*% PX
            
        }    
    }
    simRegret=simRegret / RR
    print(simRegret)
    ## TBD: formatting, design choices, further comparisons?
}



#modify this to get repeated draws, and adjust X distribution / weighting of regret
# include stratified equalsplit as comparison, and fully random assignment.

