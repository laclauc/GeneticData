# -------------
# Copyright (c) 2018 V. Brault, C. Laclau
# -------------
#
# -------------
# License
# The code can be used for research purposes only.

BlocGibbs_LBM=function(x_ij,theta,a,b,t_jl,niter,eps=10^(-16),verbatim=TRUE){
  g=dim(theta$alpha_kl)[1]
  m=dim(theta$alpha_kl)[2]
  n=dim(x_ij)[1]
  d=dim(x_ij)[2]

  ## Save all parameters
  z=matrix(0,n,1);w=matrix(0,d,1);
  alpha_kl_path=matrix(0,g,m);
  pi_k_path=rep(0,g);
  tau_l_path=rep(0,m);
  Part=list()
  Part$z=matrix(0,nrow=n,ncol=g);
  Part$w=matrix(0,nrow=d,ncol=m);

  alpha_kl=theta$alpha_kl;
  pi_k=theta$pi_k;
  tau_l=theta$tau_l;
  t_l=colSums(t_jl);
  T_jl=matrix(0,d,m)
  if (verbatim){
    pb=txtProgressBar(min=1,max=niter,style=3)
  }

  for (iter in 1:niter){
    A_kl=log((alpha_kl)/(1-alpha_kl));
    B_kl=log(1-alpha_kl);
    log_pi_k=log(pi_k);
    log_tau_l=log(tau_l);

    # Computation of z
    S_ik=(x_ij%*%t_jl)%*%t(A_kl)+rep(1,n)%*%t(t_l)%*%(t(B_kl))+matrix(1,n,1)%*%t(log_pi_k);
    S_ikmax=apply(S_ik, 1,max)
    Sp_ik = exp(S_ik-S_ikmax);
    s_ik=Sp_ik/(rowSums(Sp_ik))
    for (iter_i in 1:n){
      z[iter_i]=which(as.logical(rmultinom(1,size=1,prob=s_ik[iter_i,])));
    }
    s_ik=!(z%*%matrix(1,1,g)-matrix(1,n,1)%*%(1:g));
    s_k=t(colSums(s_ik));

    # Computation of w
    T_jl=(t(x_ij)%*%s_ik)%*%A_kl+rep(1,d)%*%s_k%*%B_kl+matrix(1,d,1)%*%t(log_tau_l);
    T_jlmax=apply(T_jl, 1,max)
    Tp_jl = exp(T_jl-T_jlmax);
    t_jl=Tp_jl/(rowSums(Tp_jl));
    for (iter_j in 1:d){
      w[iter_j]=which(as.logical(rmultinom(1,size=1,prob=t_jl[iter_j,])));
    }
    t_jl=!(w%*%matrix(1,1,m)-matrix(1,d,1)%*%(1:(m)));
    t_l=colSums(t_jl);   

    # Computation of parameters
    # pi_k
    nz=s_k+a;
    pi_k=rgamma(g,nz,1)
    pi_k=pi_k/(sum(pi_k));
    #tau_l
    dw=t_l+a;
    tau_l=rgamma(m,dw,1)
    tau_l=tau_l/(sum(tau_l));
    # alpha_kl
    u_kl=c(t(s_ik)%*%x_ij%*%t_jl+b);
    n_kl=c(t(s_k)%*%t_l-u_kl+2*b);
    alpha_kl=matrix(rbeta(g*m,u_kl,n_kl),ncol=m)
    
    # trajectories
    resorg=BinBlocClassement_LBM(theta=list(
      pi_k=pi_k,
      tau_l=tau_l,
      alpha_kl=alpha_kl),
      Part=list(
        z=s_ik,
        w=t_jl)
    )
    alpha_kl_path=alpha_kl_path+resorg$theta$alpha_kl;
    pi_k_path=pi_k_path+resorg$theta$pi_k;
    tau_l_path=tau_l_path+resorg$theta$tau_l;
    Part$z=Part$z+resorg$z;
    Part$w=Part$w+resorg$w;
    if (verbatim){
      setTxtProgressBar(pb,iter)
    }
  }


# output arguments
  theta=list()
  theta$alpha_kl=alpha_kl_path/niter;
  theta$pi_k=pi_k_path/niter;
  theta$tau_l=tau_l_path/niter;

  list(theta=theta,Part=Part)
}
  