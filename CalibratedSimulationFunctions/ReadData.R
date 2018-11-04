DataToTheta=function(filename, dataname, k, strataVars, OutcomeName,TreatmentName, CovariatesNames, printFigures=FALSE){
    Data=read_csv(paste("../Datasets/Cleaned/", filename, ".csv", sep=""))
    N=nrow(Data)
    
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
    
    
    # Figures and tables
    if (printFigures) {
        FooterText=paste("Average of ", OutcomeName,
                         "\nfor each value of the treatment (", TreatmentName,
                         "),\nand each value of ", CovariatesNames, ".", sep="")
        PrintDataFigures(stratasizes,sumstats,filename,dataname,FooterText, k)
    }
    
    
    
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
    
    #think about unifying output of this function across old simulations / Thompson!
    list(filename=filename, dataname=dataname, N=N,
         theta=theta$meanout, sumstats=sumstats, stratasizes =stratasizes, key=key, #output for old simulations (no covariates)
         SS=SS, NN=NN, PX=PX, k=k) #output for thompson simulations (using covariates)
}



#produce figures and tables of calibrated parameter values
PrintDataFigures=function(stratasizes, sumstats, filename, dataname,FooterText, k){
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
    

    
    ggplot(drop_na(sumstats), aes(x=meanout, y=factor(treatment,levels=c(k:1)))) +
        geom_point(color=fillcolor) +
        geom_segment(aes(x=0, xend=meanout, yend=factor(treatment,levels=c(k:1))), color=fillcolor, size=.5) +
        #scale_y_discrete(labels=paste("treatment", 1:k)) +
        facet_grid(strata~.) +
        scale_x_continuous(limits=c(0,xmax)) +
        theme_light() +
        theme(#panel.background = element_rect(fill = backcolor, colour = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.x=element_blank(),
            plot.caption=element_text(size=7)) +
        labs(title="Average outcomes by treatment and stratum",
             subtitle=dataname,
             x="mean outcome", y="treatment",
             caption=FooterText)
    
    # alternative version, with sorting first by strata, then by treatment
    
    # ggplot(drop_na(sumstats), aes(x=meanout, y=factor(strata))) +
    #     geom_point(color=fillcolor) +
    #     geom_segment(aes(x=0, xend=meanout, yend=factor(strata)), color=fillcolor, size=.5) +
    #     #scale_y_discrete(labels=paste("treatment", 1:k)) + 
    #     facet_grid(factor(treatment, levels=(1:k)) ~.) +
    #     scale_x_continuous(limits=c(0,xmax)) +
    #     theme_light() +
    #     theme(#panel.background = element_rect(fill = backcolor, colour = NA),
    #         panel.grid.major = element_blank(),
    #         panel.grid.minor = element_blank(),
    #         axis.ticks.x=element_blank()) +
    #     labs(title="Average outcomes by treatment and stratum",
    #          subtitle=dataname,
    #          x="mean outcome", y="stratum", facet="treatment")
    
    ggsave(paste("../Figures/Applications/", filename, ".pdf", sep=""), width = 4, height = 0.15*length(levels(sumstats$strata))*k+1.5)
    
    write_csv(sumstats, paste("../Figures/Applications/", filename, "Sumstats.csv", sep=""))
    
}





ReadAllData=function(printFigures=F){
    #read table of all applications from csv file
    ApplicationTable=read_csv("../Datasets/ApplicationTable.csv")
  
    #read in data for each row of Application table, store output in each element of DataList
    DataList=pmap(ApplicationTable, function(filename, dataname, k, stratavars, OutcomeName,TreatmentName, CovariatesNames) #make sure ApplicationTable has these column names
                  DataToTheta(filename, dataname, k, stratavars, OutcomeName, TreatmentName, CovariatesNames, printFigures=printFigures))
        

    DataList
}






#parameters of hypothetical experiment
# NN=rep(24,12)
# R=4000
#DesignTable(NN,ThetaList,R,columnames,"Applications/CalibratedSimulations")


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