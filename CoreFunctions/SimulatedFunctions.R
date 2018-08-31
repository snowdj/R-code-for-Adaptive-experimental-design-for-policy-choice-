library(tidyverse)


#number of replication draws from simplex
#replace by NULL if full optimization is desired
#RR=400
RR=NULL

# number of replication draws for Uhat
RU=1000
# maximal number of units in a wave
Nmax=20
# number of treatment arms
k=3

# generating random draws once and for all
# using global assignment operator <<- 
Seed = function(A,B){
  set.seed(12231983)
  
  #draws from prior distribution of theta - dimensions (RU, k)
  theta <<- sapply(1:k, function(d)  rbeta(RU, A[d], B[d]))
  #draws of potential outcomes - dimensions (Nmax, RU, k)
  Yd <<- sapply(theta, function(y) rbernoulli(Nmax, p = y))
  Yd <<- array(Yd, c(Nmax, RU,k))
}


# a simulated approximation of U
Uhat=function(A,B,C,n, Vfunction=SWF){
  # number of successes for design n, for each of r=1..RU replicates
  SR=sapply(1:k, function(d) if (n[d]>1) {colSums(Yd[1:n[d],,d],1)}
                                else if (n[d]==1) {Yd[1,,d]}
                                else {rep(0,RU)} )
  # posterior expected social welfare for each realization of SR
  SW=sapply(1:RU, function(r) Vfunction(A+SR[r,] ,B+n-SR[r,] ,C))
  mean(SW)
}




#creating nmatrix that has all elements of simplex
simplex = function(N,k){
  Nplus1=N+1
  Nplus1tok=Nplus1^k
  if (is.null(RR)){ #computationally costly version of getting full simplex
    nmatrix=matrix(0, Nplus1tok, k)
    i=1:(Nplus1tok)
    for (j in k:1) {
      nmatrix[,j]=i%%Nplus1
      i=i%/%Nplus1
    }
    #keep only rows with correct rowsum
    nmatrix[rowSums(nmatrix)==N,,drop=FALSE]
  } else 
  { #getting random subsample of simplex
    nmatrix=matrix(runif(RR*k), RR, k) #sample from contiuous hypercube
    nmatrix=t(apply(nmatrix, 1, function(x) N*x/sum(x))) #row normalize. this gives more points at center!
    nmatrix=floor(nmatrix) #round down
    nmatrix+  #add missing numbers to get right rowsum again
      t(sapply(N-rowSums(floor(nmatrix)), function(dN) sample(c(rep(1, dN), rep(0, k-dN))))) 
  }
  
}




####
#for testing
# A=c(1,1,1)
# B=c(1,1,1)
# C=c(0,0,0)
# 
# Seed(A,B)
# Uhat(A,B,C,c(0,0,10))
# Uhat(A,B,C,c(10,0,0))
# Uhat(A,B,C,c(5,5,5))
# Uhat(A,B,C,c(10,10,10))
# Uhat(A,B,C,c(20,20,20))
