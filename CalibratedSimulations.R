source("OptimalAssignmentFunctions/WelfareFunctions.R")
source("OptimalAssignmentFunctions/welfareplotsGraphics.R")
source("OptimalAssignmentFunctions/SimulatedWelfareFunctions.R")

source("CalibratedSimulationFunctions/ReadData.R")
source("CalibratedSimulationFunctions/calibratedSimulationFunctions.R")


backcolor="azure2" #background color for plots
gridcolor="azure1"
fillcolor="skyblue4"


DataList=ReadAllData(printFigures=F)

#for each elemant of DataList (each application), extract the theta sub-element, and filename sub-element
#ThetaList=map(DataList, "theta")
columnames=map(DataList, "filename")   

#size of waves, and number of replications
NN=c(24,24)
R=1000

#which methods of treatment assignment to simulate
#indices corresponding to:
# c("conventional",
#   "optimalhandpicked",
#   "optimalrandom",
#   "optimalhat",
#   "optimal",
#   "besthalf",
#   "thompson",
#   "modifiedthompson")   
methods=c(1,2,7,8)

DesignTable(NN,DataList,methods,R,columnames,filename="CalibratedSimulations")



