library(tidyverse)
library(magrittr)
#filepath="../../Datasets/TestDataApp/Ashraf.csv"

ReadDataApp=function(filepath,
                               key=NULL){ #key to consistently create strata
  Data=read_csv(filepath)
  varnames=names(Data)
  
  #extract number of treatments, number of covariates
  k=0
  while (paste("treatment",k+1, sep="") %in% varnames ) {
    k=k+1 
  }
  ncovs=0
  while (paste("covar",ncovs+1, sep="") %in% varnames ) {
    ncovs=ncovs+1 
  }
  
  treatDummies=paste("treatment",1:k, sep="") #names of treatment variables
  strataVars=paste("covar",1:ncovs, sep="") #names of control variables
  
  Data %<>% drop_na %>%  #drop rows with missing entries
            mutate(Strata=(interaction(select(.,strataVars)))) #create strata
  
  #recoding the strata levels in a reproducible way across datasets 
  oldlevels=levels(Data$Strata) 
  if (is.null(key)) {
    key=data_frame(Strata=oldlevels, strata=factor(1:length(oldlevels)))
  }
  Data %<>% left_join(., key, by = "Strata") %>%
            select(-Strata)
  
  #create factor variable from treatment
  if (k>0) {
    Data %<>% mutate(treatment=factor(as.matrix(Data[treatDummies])%*%(1:k)))
  } else if ("treatment" %in% varnames) { #if treatment is provided as categorical variable, rather than dummies
    Data %<>% mutate(treatment=factor(treatment))
  }    

  if ("outcome" %in% varnames) {
    retlist=list(Y=Data$outcome, 
                 D=Data$treatment, 
                 X=Data$strata,
                 k=k, nx=length(key$strata),
                 key=key)
  } else {
    retlist=Data[strataVars]
    retlist$Xt=Data$strata
  }
  
  
  retlist
}

