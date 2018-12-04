BlocGibbsMult=function(v_ijh,theta,a,b,C,e,t_jl,niter,eps=10^(-16)){
  g=dim(theta$alpha_kl)[1]
  m=dim(theta$alpha_kl)[2]
  r=dim(theta$alpha_kl)[3]
  n=dim(v_ijh)[1]
  d=dim(v_ijh)[2]

  ## Save all parameters
  z=matrix(0,n,1);w=matrix(0,d,1);
  alpha_kl_path=array(0,c(g,m,r));
  pi_k_path=rep(0,g);
  tau_l_path=rep(0,m);
  phi_path=0;
  lambda_path=matrix(0,r,d)
  Part=list()
  Part$z=matrix(0,nrow=n,ncol=g);
  Part$w=matrix(0,nrow=d,ncol=m+1);

  alpha_kl=theta$alpha_kl;
  pi_k=theta$pi_k;
  tau_l=theta$tau_l;
  phi=theta$phi
  lambda=theta$lambda
  t_l=colSums(t_jl);
  T_jl=matrix(0,d,m+1)
  xplusj=t(sapply(1:r,function(h){colSums(v_ijh[,,h])}))
  xplusje=xplusj+e%*%t(rep(1,d))
  pb=txtProgressBar(min=1,max=niter,style=3)

  for (iter in 1:niter){
    log_alpha_kl=log((alpha_kl*(alpha_kl>0)+(alpha_kl<=0)*eps));
    log_pi_k=log(pi_k*(pi_k>0)+(pi_k<=0)*eps);
    log_tau_l=log(tau_l*(tau_l>0)+(tau_l<=0)*eps);
    log_lambda=log(lambda*(lambda>0)+(lambda<=0)*eps)

    # Computation of z
    S_ik=matrix(1,n,1)%*%t(log_pi_k)+
      apply(array(unlist(lapply(1:r,function(h){
        (v_ijh[,,h]%*%t_jl[,-1])%*%t(log_alpha_kl[,,h])
      })),c(n,g,r)),MARGIN = c(1,2),sum)
    S_ikmax=apply(S_ik, 1,max)
    Sp_ik = exp(S_ik-S_ikmax);
    s_ik=Sp_ik/(rowSums(Sp_ik))
    for (iter_i in 1:n){
      z[iter_i]=which(as.logical(rmultinom(1,size=1,prob=s_ik[iter_i,])));
    }
    s_ik=!(z%*%matrix(1,1,g)-matrix(1,n,1)%*%(1:g));
    s_k=t(colSums(s_ik));

    # Computation of w
    T_jl[,-1]=matrix(rep(1,d),ncol=1)%*%log_tau_l+log(phi)+apply(array(unlist(lapply(1:r,function(h){
      (t(v_ijh[,,h])%*%s_ik)%*%log_alpha_kl[,,h]
    })),c(d,m,r)),MARGIN = c(1,2),sum)
    T_jl[,1]=log(1-phi)+colSums(xplusj*log_lambda)
    T_jlmax=apply(T_jl, 1,max)
    Tp_jl = exp(T_jl-T_jlmax);
    t_jl=Tp_jl/(rowSums(Tp_jl));
    for (iter_j in 1:d){
      w[iter_j]=which(as.logical(rmultinom(1,size=1,prob=t_jl[iter_j,])));
    }
    t_jl=!(w%*%matrix(1,1,m+1)-matrix(1,d,1)%*%(1:(m+1)));
    t_l=colSums(t_jl);   

    # Computation of parameters
    # pi_k
    nz=s_k+a;
    pi_k=rgamma(g,nz,1)
    pi_k=pi_k/(sum(pi_k));
    #tau_l
    dw=t_l[-1]+a;
    tau_l=rgamma(m,dw,1)
    tau_l=tau_l/(sum(tau_l));
    # alpha_kl
    u_kl=array(unlist(lapply(1:r,function(h){t(s_ik)%*%v_ijh[,,h]%*%t_jl[,-1]+b})),c(g,m,r));
    alpha_kl=array(rgamma(g*m*r,u_kl,1),c(g,m,r))
    Norm=array(rep(sapply(1:m,function(l){sapply(1:g,function(k){sum(alpha_kl[k,l,])})}),r),
               c(g,m,r))
    alpha_kl=alpha_kl/Norm
    # phi
    phi=rbeta(1,d-t_l[1]+C[1],t_l[1]+C[2])
    # lambda
    lambda=matrix(rgamma(d*r,xplusje,1),r)
    lambda=lambda/(rep(1,r)%*%t(colSums(lambda)))
    
    # trajectories
    resorg=MultBlocClassement(theta=list(
      pi_k=pi_k,
      tau_l=tau_l,
      alpha_kl=alpha_kl,
      phi=phi,
      lambda=lambda),
      Part=list(
        z=s_ik,
        w=t_jl)
    )
    alpha_kl_path=alpha_kl_path+resorg$theta$alpha_kl;
    pi_k_path=pi_k_path+resorg$theta$pi_k;
    tau_l_path=tau_l_path+resorg$theta$tau_l;
    phi_path=phi_path+resorg$theta$phi;
    lambda_path=lambda_path+resorg$theta$lambda;
    Part$z=Part$z+resorg$z;
    Part$w=Part$w+resorg$w;
    setTxtProgressBar(pb,iter)
  }


# output arguments
  theta=list()
  theta$alpha_kl=alpha_kl_path/niter;
  theta$pi_k=pi_k_path/niter;
  theta$tau_l=tau_l_path/niter;
  theta$phi=phi_path/niter;
  theta$lambda=lambda_path/niter;

  list(theta=theta,Part=Part)
}
  