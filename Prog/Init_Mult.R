Init_Mult=function(g,m,x_ij,ale=TRUE){
  eps=10^(-16)
  n=dim(x_ij)[1];
  d=dim(x_ij)[2]
  r=max(x_ij)
  # Initialisation with random sampling of centers
  
  theta=list()
  if (ale){
    theta$pi_k=runif(g)
    theta$pi_k=theta$pi_k/sum(theta$pi_k)
    theta$tau_l=runif(m+1)
    theta$tau_l=theta$tau_l/sum(theta$tau_l)
  }else{
    theta$pi_k=1/g *rep(1,g);
    theta$tau_l=1/m *rep(1,m);
  }
  Z=PartRnd(n,theta$pi_k);
  W=PartRnd(d,theta$tau_l);
  Norm=colSums(Z$z_ik)%*%t(colSums(W$z_ik))
  theta$alpha_kl=array(unlist(lapply(1:r,function(h){
    t(Z$z_ik)%*%(data$x==h)%*%W$z_ik/Norm
    })),c(g,m,r))
  theta$lambda=matrix(runif(d*r),ncol=d)
  theta$lambda=theta$lambda/matrix(rep(1,r),ncol=1)%*%apply(theta$lambda,MARGIN = 2,sum)
  theta$phi=1-mean(W$z_ik[,1])
  theta$tau_l=apply(W$z_ik[,-1],MARGIN = 2,mean)
  theta$tau_l=theta$tau_l/sum(theta$tau_l)
  list(t_jl=W$z_ik,theta=theta,w=W$z)
}