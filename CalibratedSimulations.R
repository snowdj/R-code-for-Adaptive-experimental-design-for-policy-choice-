source("OptimalAssignmentFunctions/WelfareFunctions.R")
source("OptimalAssignmentFunctions/welfareplotsGraphics.R")
source("OptimalAssignmentFunctions/SimulatedWelfareFunctions.R")

source("CalibratedSimulationFunctions/ReadData.R")
source("CalibratedSimulationFunctions/calibratedSimulationFunctions.R")

source("ThompsonHierarchicalFunctions/ThompsonHierarchical.R")
source("ThompsonHierarchicalFunctions/calibratedSimulationFunctionsCovariates.R")

#backcolor="azure2" #background color for plots
#gridcolor="azure1"
fillcolor="skyblue4"

DoNoCovariates=F
DoCovariates=T


DataList=ReadAllData(printFigures=T)
#for each elemant of DataList (each application), extract the filename sub-element
columnames=map(DataList, "filename")   

#number of replications
MC_replicates=5000

if (DoNoCovariates) {
    #which methods of treatment assignment to simulate
    methods=c(1,7,8,9)
    #Simulations without covariates
    for (Mult in c(.5,1)) { # multiples of original sample size
      for (waves in c(2,4,10)) { # number of waves
          for (i in 1:length(DataList)) {
              DataList[[i]]$wavesizes = rep( floor(DataList[[i]]$N * Mult / waves) ,waves) 
          }
          DesignTable(DataList,methods,MC_replicates,columnames,filename=paste("test_CalibratedSimulations_T_", waves, "_" ,Mult, sep=""))
      }
    }
}




if (DoCovariates) {
    #Simulations with covariates
    methods=c(1,2,3,4)
    
    for (Mult in c(.5,1)) { # multiples of original sample size
      for (waves in c(2,4,10)) { # number of waves
        for (i in 1:length(DataList)) {
          DataList[[i]]$wavesizes = rep( floor(DataList[[i]]$N * Mult / waves) ,waves) 
        }
        DesignTableCovariates(DataList,methods,MC_replicates,columnames,filename=paste("test_Covariates_CalibratedSimulations_T_", waves, "_" ,Mult, sep=""))
      }
    }
}







#Old List of Thetas to consider
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

