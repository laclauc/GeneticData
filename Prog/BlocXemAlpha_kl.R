# -------------
# Copyright (c) 2018 V. Brault, C. Laclau
# -------------
#
# -------------
# License
# The code can be used for research purposes only.

BlocXemAlpha_kl =function(x_ij,theta,t_jl,a,b,C,e,niter=100000,eps=10^(-16),epsi=10^(-16),epsi_int=0.1,niter_int=1){
  alpha_kl=theta$alpha_kl;
  pi_k=theta$pi_k;
  tau_l=theta$tau_l;
  phi=theta$phi
  
  g=dim(alpha_kl)[1]
  m=dim(alpha_kl)[2]
  n=dim(x_ij)[1]
  d=dim(x_ij)[2]
  t_l=t(colSums(t_jl));
  T_jl=matrix(0,d,m+1)
  xplusj=colSums(x_ij)
  lambda=(xplusj+e[1]-1)/(n+sum(e)-2)
  A_lambda=log(((lambda*(lambda>0)+(lambda<=0)*eps))/((1-lambda)*(lambda<1)+(lambda>=1)*eps))
  B_lambda=log((1-lambda)*(lambda<1)+(lambda>=1)*eps)
  # Loop over the algorithm iterations
  pi_old=pi_k
  tau_old=tau_l
  phi_old=phi
  alpha_old=alpha_kl
  for (iter in 1:niter){
    log_phi=log(phi*(phi>0)+(phi<=0)*eps)
    log_moinsphi=log((1-phi)*(phi<1)+(phi>=1)*eps)
    A_kl=log((alpha_kl*(alpha_kl>0)+(alpha_kl<=0)*eps)/((1-alpha_kl)*((1-alpha_kl)>0)+((1-alpha_kl)<=0)*eps));
    B_kl=log((1-alpha_kl)*(alpha_kl<1)+(alpha_kl>=1)*eps);
    log_pi_k=log(pi_k+(pi_k<=0)*eps);
    log_tau_l=log(tau_l+(tau_l<=0)*eps);
    
    W2=-Inf;
    for (iter_int in 1:niter_int){
      # Computation of sik
      S_ik=(x_ij%*%t_jl[,-1])%*%t(A_kl)+rep(1,n)%*%t(t_l[-1])%*%(t(B_kl))+matrix(1,n,1)%*%t(log_pi_k);
      S_ikmax=apply(S_ik, 1,max)
      Sp_ik = exp(S_ik-S_ikmax);
      s_ik=Sp_ik/(rowSums(Sp_ik))
      s_k=t(colSums(s_ik));
      
      # Computation of t_jl
      T_jl[,-1]=(t(x_ij)%*%s_ik)%*%A_kl+rep(1,d)%*%s_k%*%B_kl+matrix(rep(1,d),ncol=1)%*%log_tau_l+log_phi;
      T_jl[,1]=log_moinsphi+xplusj*A_lambda+n*B_lambda
      T_jlmax=apply(T_jl, 1,max)
      Tp_jl = exp(T_jl-T_jlmax);
      t_jl=Tp_jl/(rowSums(Tp_jl));  
      t_l=t(colSums(t_jl));
    }
    # Computation of parameters
    pi_k=t((s_k+a-1)/(n+g*(a-1)));
    tau_l=t((t_l[-1]+a-1)/(sum(t_l[-1])+m*(a-1)));
    if (any(is.nan(tau_l))){
      tau_l=t((t_l[-1]+a-1)/((sum(t_l[-1])+m*(a-1))*((sum(t_l[-1])+m*(a-1))>0)+
                               eps*((sum(t_l[-1])+m*(a-1))<=0)));
    }
    phi=1-(t_l[1]+C[1]-1)/(d+sum(C)-2)
    u_kl=t(s_ik)%*%x_ij%*%t_jl[,-1];
    n_kl=t(s_k)%*%t_l[-1];
    u_kl=u_kl*(u_kl>0);
    n_kl=n_kl*(n_kl>0);
    alpha_kl=(u_kl+b-1)/(n_kl+2*(b-1));
    if (any(is.nan(alpha_kl))){
      alpha_kl=((u_kl+b-1)*((u_kl+b-1)>0))/((n_kl+2*(b-1))*((n_kl+2*(b-1))>0)+
                                              eps*((n_kl+2*(b-1))<=0));
    }
    #Stopping criterion
    if ((max(c(abs(pi_old-pi_k),abs(tau_old-tau_l),abs(alpha_old-alpha_kl),
              abs(phi-phi_old)))< epsi)){
      break
    }else{
      pi_old=pi_k
      tau_old=tau_l
      phi_old=phi
      alpha_old=alpha_kl
    }
    
  }
  # W=sum(c(s_k)*log(pi_k))+t_l[1]*log(1-phi)+(d-t_l[1])*log(phi)+sum(t_l[-1]*log(tau_l))+
  #   sum(t_jl[,1]*(xplusj*log(lambda)+(n-xplusj)*log(1-lambda)))+
  #   sum(u_kl*log(alpha_kl)+(n_kl-u_kl)*log(1-alpha_kl))
  # if (is.na(W)){
    W=sum(c(s_k)*log(pi_k*(pi_k>0)+eps*(pi_k<=0)))+
      t_l[1]*log((1-phi)*(phi<1)+eps*(phi>=1))+(d-t_l[1])*log(phi*(phi>0)+eps*(phi<=0))+
      sum(t_l[-1]*log(tau_l*(tau_l>0)+eps*(tau_l<=0)))+
      sum(t_jl[,1]*(xplusj*log(lambda*(lambda>0)+eps*(lambda<=0))+
                      (n-xplusj)*log((1-lambda)*(lambda<1)+eps*(lambda>=1))))+
      sum(u_kl*log(alpha_kl*(alpha_kl>0)+eps*(alpha_kl<=0))+
            (n_kl-u_kl)*log((1-alpha_kl)*(alpha_kl<1)+eps*(alpha_kl>=1)))
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