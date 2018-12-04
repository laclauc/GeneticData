# -------------
# Copyright (c) 2018 V. Brault, C. Laclau
# -------------
#
# -------------
# License
# The code can be used for research purposes only.

MultBlocXemAlpha_kl =function(v_ijh,theta,t_jl,a,b,C,e,niter=100000,eps=10^(-16),epsi=10^(-16),epsi_int=0.1,niter_int=1){
  alpha_kl=theta$alpha_kl;
  pi_k=theta$pi_k;
  tau_l=theta$tau_l;
  phi=theta$phi
  
  g=dim(alpha_kl)[1]
  m=dim(alpha_kl)[2]
  r=dim(alpha_kl)[3]
  n=dim(v_ijh)[1]
  d=dim(v_ijh)[2]
  t_l=t(colSums(t_jl));
  T_jl=matrix(0,d,m+1)
  xplusj=t(sapply(1:r,function(h){colSums(v_ijh[,,h])}))
  lambda=xplusj+t(sapply(1:r,function(h){rep(e[h]-1,d)}))
  lambda=lambda/matrix(rep(1,r),ncol=1)%*%apply(lambda,MARGIN = 2,sum)
  log_lambda=log(lambda*(lambda>0)+(lambda<=0)*eps)
  
  # Loop over the algorithm iterations
  pi_old=pi_k
  tau_old=tau_l
  phi_old=phi
  alpha_old=alpha_kl
  for (iter in 1:niter){
    log_alpha_kl=log((alpha_kl*(alpha_kl>0)+(alpha_kl<=0)*eps));
    log_pi_k=log(pi_k*(pi_k>0)+(pi_k<=0)*eps);
    log_tau_l=log(tau_l*(tau_l>0)+(tau_l<=0)*eps);
    log_phi=log(phi*(phi>0)+(phi<=0)*eps)
    log_moinsphi=log((1-phi)*(phi<1)+(phi>=1)*eps)
    
    W2=-Inf;
    for (iter_int in 1:niter_int){
      # Computation of sik
      S_ik=matrix(1,n,1)%*%t(log_pi_k)+
        apply(array(unlist(lapply(1:r,function(h){
          (v_ijh[,,h]%*%t_jl[,-1])%*%t(log_alpha_kl[,,h])
        })),c(n,g,r)),MARGIN = c(1,2),sum)
      S_ikmax=apply(S_ik, 1,max)
      Sp_ik = exp(S_ik-S_ikmax);
      s_ik=Sp_ik/(rowSums(Sp_ik))
      s_k=t(colSums(s_ik));
      
      # Computation of t_jl
      T_jl[,-1]=matrix(rep(1,d),ncol=1)%*%log_tau_l+log_phi+apply(array(unlist(lapply(1:r,function(h){
        (t(v_ijh[,,h])%*%s_ik)%*%log_alpha_kl[,,h]
      })),c(d,m,r)),MARGIN = c(1,2),sum)
      T_jl[,1]=log_moinsphi+colSums(xplusj*log_lambda)
      T_jlmax=apply(T_jl, 1,max)
      Tp_jl = exp(T_jl-T_jlmax);
      t_jl=Tp_jl/(rowSums(Tp_jl));  
      t_l=t(colSums(t_jl));
    }
    # Computation of parameters
    pi_k=t((s_k+a-1)/(n+g*(a-1)));
    tau_l=t((t_l[-1]+a-1)/(sum(t_l[-1])+m*(a-1)));
    tau_l=t((t_l[-1]+a-1)/((sum(t_l[-1])+m*(a-1))*((sum(t_l[-1])+m*(a-1))>0)+
                             ((sum(t_l[-1])+m*(a-1))<=0)*eps));
    phi=1-(t_l[1]+C[1]-1)/(d+sum(C)-length(C))
    u_kl=array(unlist(lapply(1:r,function(h){t(s_ik)%*%v_ijh[,,h]%*%t_jl[,-1]+b-1})),c(g,m,r));
    n_kl=array(rep(t(s_k)%*%t_l[-1]+r*(b-1),r),c(g,m,r))
    u_kl=u_kl*(u_kl>0);
    n_kl=n_kl*(n_kl>0);
    alpha_kl=u_kl/n_kl
    if (any(is.na(alpha_kl))){
      alpha_kl=(u_kl*(u_kl>0)+(u_kl<=0)*eps)/(n_kl*(n_kl>0)+(n_kl<=0)*eps)
    }
    #Stopping criterion
    if (max(c(abs(pi_old-pi_k),abs(tau_old-tau_l),abs(alpha_old-alpha_kl),
              abs(phi-phi_old)))< epsi){
      break
    }else{
      pi_old=pi_k
      tau_old=tau_l
      phi_old=phi
      alpha_old=alpha_kl
    }
    
  }
  # W=sum(c(s_k)*log(pi_k))+t_l[1]*log(1-phi)+(d-t_l[1])*log(phi)+sum(t_l[-1]*log(tau_l))+
  #   sum((xplusj*log(lambda))%*%t_jl[,1])+
  #   sum(u_kl*log(alpha_kl))
  # if (is.na(W)){
    W=sum(c(s_k)*log(pi_k*(pi_k>0)+(pi_k<=0)*eps))+
      t_l[1]*log((1-phi)*(phi<1)+(phi>=1)*eps)+
      (d-t_l[1])*log(phi*(phi>0)+(phi<=0)*eps)+
      sum(t_l[-1]*log(tau_l*(tau_l>0)+(tau_l<=0)*eps))+
      sum((xplusj*log(lambda*(lambda>0)+(lambda<=0)*eps))%*%t_jl[,1])+
      sum(u_kl*log(alpha_kl*(alpha_kl>0)+(alpha_kl<=0)*eps))
  # }
  theta$alpha_kl=alpha_kl;
  theta$pi_k=pi_k;
  theta$tau_l=tau_l;
  theta$phi=phi
  theta$lambda=lambda
  z=apply(S_ik, 1, which.max)
  w=apply(T_jl, 1, which.max)
  empty_cluster=((length(unique(z))!=g)||(length(unique(w))!=(m+1)))
  list(W=W,theta=theta,s_ik=s_ik,t_jl=t_jl,iter=iter,empty_cluster=empty_cluster,z=z,w=w)
}