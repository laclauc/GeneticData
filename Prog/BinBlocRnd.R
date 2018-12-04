BinBlocRnd =function(n,d,theta){

Z=PartRnd(n,theta$pi_k);
if (theta$phi<1){
  W=PartRnd(d,c(1-theta$phi,theta$phi*theta$tau_l));
  x_ij=matrix(0,n,d)
  w0=sum(W$z_ik[,1])
  x_ij[,W$z!=1]=matrix(as.numeric(matrix(runif(n*(d-w0)),nrow=n)<Z$z_ik%*%theta$alpha_kl%*%t(W$z_ik[!W$z_ik[,1],-1])),nrow=n);
  x_ij[,W$z_ik[,1]]=matrix(as.numeric(matrix(runif(n*w0),nrow=n)<matrix(1,n,1)%*%theta$lambda_j[W$z_ik[,1]]),nrow=n)
}else{
  W=PartRnd(d,theta$tau_l);
  x_ij=matrix(0,n,d)
  x_ij=matrix(as.numeric(matrix(runif(n*d),nrow=n)<Z$z_ik%*%theta$alpha_kl%*%t(W$z_ik)),nrow=n);
}
list(x=x_ij,xrow=Z$z,xcol=W$z)
}
