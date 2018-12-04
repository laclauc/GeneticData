# -------------
# Copyright (c) 2018 V. Brault, C. Laclau
# -------------
#
# -------------
# License
# The code can be used for research purposes only.

BinBlocGVB=function(data,g,m,a,b,niter,ini=list(ale=FALSE),classement=TRUE,eps=10^(-16),niterVB=100000,epsi=10^(-16),epsi_int=0.1,niter_int=1){
  res=BinBlocGibbs(data=data,g=g,m=m,a=a,b=b,niter=niter,ini=ini,classement=classement,eps=eps)
  BinBlocVBayes(data=data,g=g,m=m,a=a,b=b,ini=list(theta=res$theta,t_jl=res$Part$w),niter=niterVB,eps=eps,epsi=epsi,epsi_int=epsi_int,niter_int=niter_int)
}