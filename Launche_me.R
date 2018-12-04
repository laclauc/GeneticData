rm(list=ls())

### Loading of all the functions of NFLBM
source("Prog/Chargement.R")

### Parameters
g=5 # Nb cluster of row
m=3 # Nb cluster of column + noise cluster
a=4 # Prior on clusters
b=1 # Prior on alpha_kl
theta=list() 
theta$pi_k=1/g*rep(1,g) # Probability of row cluster
theta$tau_l=1/m*rep(1,m) # Probability of column cluster
theta$phi=0.8 # Probability of part of LBM
eps=0.2 # Difficulty
theta$alpha_kl=matrix(c(1-eps,eps  ,eps  ,
                        eps  ,1-eps,eps  ,
                        eps  ,1-eps,1-eps,
                        1-eps,1-eps,eps  ,
                        eps  ,eps  ,eps   ),ncol=m,byrow=TRUE)

n=500 # Number of row
d= 300 # Number of column
theta$lambda_j=runif(d) # Parameters for potential noise cluster

#### Simulation of data
data=BinBlocRnd(n,d,theta)
BinBlocVisu(data,z = data$xrow,data$xcol)

#### Vbayes algorithm
ResVBayes=BinBlocVBayes(data,g,m,verbatim=TRUE)
BinBlocVisu(data,ResVBayes$z,ResVBayes$w)
Crit_Vbayes=BinBlocICL(data,ResVBayes$z,ResVBayes$w)

#### Gibbs Sampler
ResGibbs=BinBlocGibbs(data,g,m,niter = 10000,verbatim=TRUE)
init=list(s_ik=!(ResGibbs$z%*%t(rep(1,g))-rep(1,n)%*%t(1:g)),
          t_jl=!(ResGibbs$w%*%t(rep(1,m+1))-rep(1,d)%*%t(1:(m+1))))
ResGibbsVBayes=BinBlocVBayes(data,g,m,ini=init,verbatim = TRUE)
BinBlocVisu(data,ResGibbsVBayes$z,ResGibbsVBayes$w)
Crit_Vbayes=BinBlocICL(data,ResGibbsVBayes$z,ResGibbsVBayes$w)

