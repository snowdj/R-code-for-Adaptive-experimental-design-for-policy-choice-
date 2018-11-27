#.rs.restartR()
library(tidyverse)

source("CalibratedSimulationFunctions/ReadData.R")
source("ThompsonHierarchicalFunctions/ThompsonHierarchical.R")
source("ThompsonHierarchicalFunctions/ThompsonCalibration.R")



DataList=ReadAllData(printFigures=F)

RunAllSimulationsThompson(DataList, RR=10)