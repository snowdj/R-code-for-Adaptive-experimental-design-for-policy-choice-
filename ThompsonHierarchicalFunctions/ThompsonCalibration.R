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





RunAllSimulationsThompson=function(DataList,
                                   T = 4, #number of waves
                                   nt = 36, #units per wave
                                   RR = 1000){ #number of replications for simulation

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

