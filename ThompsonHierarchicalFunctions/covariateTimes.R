sink("covariateTimes.txt")

for (i in 2:4) {
    Rprof("covariateTimes.Rprof")
    SimulateTWaveDesignCovariates(wavesizes,C,theta,PX, Methods[i])
    Rprof()
    
    cat("\n", Methods[i], "\n")
    times=summaryRprof("covariateTimes.Rprof")
    print(times)
}

sink()
