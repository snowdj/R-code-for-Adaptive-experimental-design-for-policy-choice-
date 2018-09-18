betabinomialMLE = function(NN,SS) { #matrix of trials and successes, across treatments and strata
    k=dim(NN)[1]
    nx=dim(NN)[2]
    
    AB=matrix(0,k,2)
    LLH=rep(0,k)
    
    # for each treatment arm, find MLE of hyperparameters
    for (d in 1:k) {
        minusLLH=function(ab) { #negative of betabinomial log likelihood, up to constant
            -sum(mapply(lbeta, SS[d,]+ab[1], NN[d,]-SS[d,] + ab[2])) +  nx* lbeta(ab[1], ab[2])
        }
        findmin<-nlminb(objective=minusLLH, 
                        start=c(1,1),
                        lower=c(0,0),
                        upper=rep(500, 2),
                        control = list(eval.max = 10^5, 
                                       iter.max = 10^5, 
                                       abs.tol = 10^(-8)))
        AB[d,] = findmin$par
        LLH[d] = -findmin$objective
    }
    
    list(AB=AB, LLH=LLH)
}

hierarchicalPosteriorDraw = function(NN,SS, 
                                     LLH) { #vector of maximized log-likelihoods, up to constant, for each treatment arm
    k=dim(NN)[1]
    nx=dim(NN)[2]
    AB=matrix(0,k,2)
    f=rep(0,k)
    
    for (d in 1:k) { #loop over treatment values
        crit=-1
        while (crit<0) { #rejection sampling from the hyper-posterior
            nn=rchisq(1,3) #chi squared prior for precision
            mn=runif(1) #uniform prior for mean
            ab=c(mn*nn, (1-mn)*nn)
            f[d]=sum(mapply(lbeta, SS[d,]+ab[1], NN[d,]-SS[d,] + ab[2])) -  k* lbeta(ab[1], ab[2])  #negative of betabinomial log likelihood, up to constant
            crit= f[d]-LLH[d]- log(runif(1))
        }
        AB[d,]=ab
    }

    list(theta=matrix(mapply(rbeta, 1,  SS+AB[,1], NN-SS + AB[,2]),k,nx), #draws from conditional beta posterior
         A=AB[,1], B=AB[,2], f=f)
}



DtchoiceThompsonHierarchical=function(Y,D,X, #outcomes, treatments, and covariates thus far
                                      k,nx, #number of treatments and number of strata
                                      Xt){ # covariates for period t
  
    SS=tapply(Y,list(D,X),sum) #matrix of successes
    NN=tapply(Y,list(D,X),length) #matrix of trials
    
    MLE=betabinomialMLE(NN,SS)
    
    Nt=length(Xt)
    Dt=rep(0,Nt)
    previousD=rep(-Inf, nx) # auxiliary vector to avoid repeat assignments of same D
    
    for (i  in 1:Nt) {
        thetadraw=hierarchicalPosteriorDraw(NN,SS,MLE$LLH)$theta[,Xt[i]] #draw from posterior for covariate value Xt[i]
        Dt[i]=which.max(thetadraw)
        if (Dt[i] == previousD[Xt[i]]) {
            thetadraw[Dt[i]] = -Inf
            Dt[i]=which.max(thetadraw)
        }
        previousD[Xt[i]] = Dt[i]
    }
    
  Dt
}




SimulateY=function(theta, D, X, N){
    thetaDX=mapply(function(d,x) theta[d,x], D,X)
    Y=runif(N)<thetaDX    
}

SimulateX=function(PX,N){
    apply(rmultinom(N,1,PX),2, which.max) #probably not the most efficient way to do this...
}

#dividing N units equally (up to rounding) between k treatments  
EqualAssignment=function(N,k){
    floor((1:N-1)*k/N)+1
}




debug = function(){
  # simulate Y, D, X based on theta matrix
  N=1000
  k=3
  nx=5
  theta=matrix(runif(k*nx),k,nx)
  print("theta:")
  print(theta)

  D=sample(k,N,replace=TRUE)
  X=sample(nx,N,replace=TRUE)
  Y=SimulateY(theta, D, X, N)

  # test whether MLE works as intended
  SS=tapply(Y,list(D,X),sum) #matrix of successes
  NN=tapply(Y,list(D,X),length) #matrix of trials
  print("share successes:")
  print(SS/NN)

  MLE=betabinomialMLE(NN,SS)
  print("hyper parameter MLEs: ")
  print(MLE$AB)

  # test whether thetadraw works as intended
  thetadraw=hierarchicalPosteriorDraw(NN,SS,MLE$LLH)$theta
  print("posterior draw estimation error: ")
  print(theta-thetadraw)

  # test whether DtchoiceThompsonHierarchical works as intended
  Dt=DtchoiceThompsonHierarchical(Y,D,X,
                               k,nx, #number of treatments and number of strata
                               SimulateX(rep(1/nx,nx), 24))
  print("Thompson assignment:")
  print(Dt)

}

debug()














#OLD VERSION WITH ARGUMENT Nt: 
# DtchoiceThompsonHierarchical=function(Y,D,X, #outcomes, treatments, and covariates thus far
#                                       k,nx, #number of treatments and number of strata
#                                       Nt){ #vector of length nx of available units for each stratum
#     
#     SS=tapply(Y,list(D,X),sum) #matrix of successes
#     NN=tapply(Y,list(D,X),length) #matrix of trials
#     
#     MLE=betabinomialMLE(NN,SS)
#     
#     #treat each treatment arm separate, thus the outer for loop
#     #but set up hierarchical model across strata
#     
#     Dt=vector("list", nx)
#     for (x in 1:nx) {
#         Dt[[x]]=rep(0,Nt[x])
#         
#         for (i in 1:Nt[x]) {
#             thetadraw=hierarchicalPosteriorDraw(NN,SS,MLE$LLH)$theta[,x]
#             Dt[[x]][i]=which.max(thetadraw)
#             if (i>1) {  #this step is modification of basic Thompson method
#                 if (Dt[[x]][i] ==Dt[[x]][i-1]) {
#                     thetadraw[Dt[[x]][i]] = -Inf
#                     Dt[[x]][i]=which.max(thetadraw)
#                 }
#             }
#             
#         }
#     }
#     
#     Dt
# }

