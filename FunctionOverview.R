#DirectoryList=list.dirs(recursive = F)
library(purrr)


listFunctions <- function(filename) {
  temp.env <- new.env()
  sys.source(filename, envir = temp.env)
  functions <- lsf.str(envir=temp.env)
  rm(temp.env)
  return(functions)
}

printfunctions=function(script) {
  cat(script)
  cat(": \n\n")
  print(listFunctions(script))
  cat("\n\n\n")
}

printfunctions_insubfolders = function(){
  dirlist= list.dirs(recursive = F)
  for (directory in dirlist[-1]) {
    setwd(directory)
    print(directory)
    ScriptList=list.files(pattern="*.R", recursive = F)
    walk(ScriptList, printfunctions)
    setwd("../")
  }
}

ScriptList=list.files(pattern="*.R", recursive = T)
#dropscripts=c(4, 7, 9, 14)
corescripts=c(3,5,6)
thompsonscripts=c(10,11,13)

#how to deal with scripts invoking relative folders?
#calibrationscripts=c(1,2)

#webapp1wave just uses copy of core scripts
#webappscripts=c(14)

allscripts=c(corescripts, thompsonscripts)

sink("functionlist.txt")
walk(ScriptList[corescripts], printfunctions)
#walk(ScriptList[webappscripts], printfunctions)
walk(ScriptList[thompsonscripts], printfunctions)
#walk(ScriptList[calibrationscripts], printfunctions)

cat("missing scripts: \n")
print(ScriptList[-allscripts])

sink()






  
