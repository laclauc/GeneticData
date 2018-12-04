Init=function(g,m,x_ij,ale=TRUE){
  eps=10^(-16)
  n=dim(x_ij)[1];
  d=dim(x_ij)[2]
  
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
  res=PartRnd(d,theta$tau_l);
  alpha_kl=matrix(runif(g*m),ncol=m)
  alpha_kl[alpha_kl<eps]=eps;
  alpha_kl[alpha_kl>1-eps]=1-eps;
  theta$alpha_kl=alpha_kl
  theta$lambda=runif(d)
  theta$phi=1-mean(res$z_ik[,1])
  theta$tau_l=apply(res$z_ik[,-1],MARGIN = 2,mean)
  theta$tau_l=theta$tau_l/sum(theta$tau_l)
  list(t_jl=res$z_ik,theta=theta,w=res$z)
}