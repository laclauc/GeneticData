# -------------
# Copyright (c) 2018 V. Brault, C. Laclau
# -------------
#
# -------------
# License
# The code can be used for research purposes only.

MultBlocRnd =function(n,d,theta){
  
  Z=PartRnd(n,theta$pi_k);
  W=PartRnd(d,c(1-theta$phi,theta$phi*theta$tau_l));
  x_ij=matrix(0,n,d)
  w0=sum(W$z_ik[,1])
  W$z=W$z-1
  x_ij[,W$z!=0]=sapply(which(!W$z_ik[,1]),function(j){
    sapply(1:n,function(i){
      sample(x=1:r,size=1,prob=theta$alpha_kl[Z$z[i],W$z[j],],replace=TRUE)}
    )})
  x_ij[,W$z_ik[,1]]=sapply(which(W$z_ik[,1]),function(j){sample(x = 1:r,size = n,prob = theta$lambda[,j],replace = TRUE)})
  list(x=x_ij,xrow=Z$z,xcol=W$z+1)
}
