betabinomial = function(n,s,a,b) {
  #probability mass function of the beta-binomial distribution
  #n trials, s successes, beta parameters a,b
  #in logs for numerical reasons
  exp(lgamma(n + 1)-lgamma(s + 1)-lgamma(n - s + 1) +
        lbeta((a + s),(b + n - s)) - lbeta(a,b)) 
}

betabinomialvector = function(N,S,A,B) {
  #wrapper for betabinomial:
  #calculating betabinomial probability, product across components
  k=length(N)
  prod(sapply(1:k, function(j) betabinomial(N[j],S[j],A[j],B[j])))
}


SWF= function(A, B, C){
  #social welfare function for expected average outcome of optimal treatment
  # A,B vectors with a, b parameters of beta distribution
  # C vector with cost of treatments
  max(A /(A+B) - C)
}

#creating nmatrix that has all elements of simplex
simplex = function(N,k){
  Nplus1=N+1
  nmatrix=matrix(0, Nplus1^k, k)
  i=1:(Nplus1^k)
  for (j in k:1) {
    nmatrix[,j]=i%%Nplus1
    i=i%/%Nplus1
  }
  #keep only rows with correct rowsum
  nmatrix[rowSums(nmatrix)==N,,drop=FALSE]
}



# beginning of period value function U
# expected welfare as a fuction of design n
# for prior A, B, , cost C, end of period value function Vfunction
U=function(A,B,C,n, Vfunction=SWF){
  k=length(A) #number of treatment arms
  #possible success vectors
  nplus1=n+1;
  N=prod(nplus1)
  SM=matrix(0, N, k)
  i=0:(N-1)
  for (j in k:1) { #mapping each element of i into corresponding vector of successes
    SM[,j]=i%%nplus1[j]
    i=i%/%nplus1[j]
  }
  
  #probability of each combination of successes
  p=sapply(1:N, function(ii) betabinomialvector(n,SM[ii,],A,B))
  #posterior expected social welfare for each combination of successes
  SW=sapply(1:N, function(ii) Vfunction(A+SM[ii,] ,B+n-SM[ii,] ,C))
  
  #expected SW
  p %*% SW
}


# for each design calculate U, given sample size N
# maximum over these will give value function V
UoverSimplex=function(A,B,C,N, Ufunction){
  k=length(A)
  nmatrix=simplex(N,k) #number of units assigned to each of k treatments, summing to N
  
  USimplex=as.data.frame(nmatrix)
  Nassignments=nrow(nmatrix) #number of different treatment assignments
  names(USimplex)=paste("n", 1:k, sep="")
  
  #calculate Ufunction for each row of the nmatrix
  USimplex$U=sapply(1:Nassignments, function(i) Ufunction(A,B,C, nmatrix[i,]))

  USimplex
}


# beginning of period value function V
# maximizing over possible designs n
# for prior A, B, , cost C
# sample sizes NN for coming waves
V=function(A,B,C,NN) {
  if (length(NN)>1) {
    Ufunction=function(A,B,C,n) U(A,B,C,n, Vfunction=function(A,B,C) V(A,B,C,NN[-1]))
  } else {
    Ufunction=U
  }
  
  USimplex=UoverSimplex(A,B,C,NN[1], Ufunction)
  max(USimplex$U) 
}  
