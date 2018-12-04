MultBlocVBayes=function(data,g,m,a=4,b=1,C=rep(1,2),e=rep(1,max(data$x)),ini=list(n_init=20,ale=TRUE),niter=1000,eps=10^(-16),epsi=10^(-16),epsi_int=0.1,niter_int=1){
  r=max(data$x)
  v_ijh=array(unlist(lapply(1:r,function(h){data$x==h})),c(n,d,r))
  if (any(names(ini)=="n_init")){
    if (!is.numeric(ini$n_init)){
      stop("For a random initialisation, n_init must be a numeric value")
    }
    if (ini$n_init<=0){
      stop("n_init must be positive")
    }
    if (!is.numeric(m)){
      stop("m must be a numeric value")
    }
    if (m<2){
      stop("m must be greater than 1")
    }
    W=-Inf
    Res=list()
    Empty=TRUE
    n=nrow(data$x)
    d=ncol(data$x)
    pb=txtProgressBar(min=1,max=ini$n_init,style=3)
    for (i in 1:ini$n_init){
      init=Init_Mult(g,m,data$x,ale=ini$ale)
      ResVBayesAle=MultBlocXemAlpha_kl(v_ijh,init$theta,init$t_jl,a=a,b=b,C=C,e=e,niter=niter,eps=eps,epsi=epsi,epsi_int=epsi_int,niter_int=niter_int)
      if (Empty){
        if (!ResVBayesAle$empty_cluster){
          Res=ResVBayesAle
          W=ResVBayesAle$W
          Empty=FALSE
        }else{
          if (ResVBayesAle$W>W){
            Res=ResVBayesAle
            W=ResVBayesAle$W
          }
        }
      }else{
        if (!ResVBayesAle$empty_cluster){
          if (ResVBayesAle$W>W){
            Res=ResVBayesAle
            W=ResVBayesAle$W
          }
        }
      }
      setTxtProgressBar(pb,i)
    }
  }else{
    if ((any(names(ini)=="t_jl"))&&(any(names(ini)=="theta"))){
      Res=MultBlocXemAlpha_kl(v_ijh,ini$theta,ini$t_jl,a=a,b=b,C=C,e=e,niter=niter,eps=eps,epsi=epsi,epsi_int=epsi_int,niter_int=niter_int)
    }else{
      if ((any(names(ini)=="t_jl"))&&(any(names(ini)=="s_ik"))){ 
        t_l=t(colSums(ini$t_jl));
        s_k=t(colSums(ini$s_ik));
        n=nrow(data$x)
        d=ncol(data$x)
        theta=list()
        # Computation of parameters
        theta$pi_k=t((s_k)/(n));
        theta$tau_l=t((t_l[-1])/(sum(t_l[-1])));
        theta$tau_l=t((t_l[-1])/((sum(t_l[-1]))*(sum(t_l[-1])>0)+
                                   (sum(t_l[-1])<=0)*eps));
        theta$phi=1-(t_l[1])/d
        u_kl=array(unlist(lapply(1:r,function(h){t(ini$s_ik)%*%v_ijh[,,h]%*%ini$t_jl[,-1]})),c(g,m,r));
        n_kl=array(rep(t(s_k)%*%t_l[-1],r),c(g,m,r))
        u_kl=u_kl*(u_kl>0);
        n_kl=n_kl*(n_kl>0);
        theta$alpha_kl=u_kl/n_kl;
        if (any(is.nan(theta$alpha_kl))){
          theta$alpha_kl=u_kl/(n_kl*(n_kl>0)+eps*(n_kl<=0));
        }
        xplusj=t(sapply(1:r,function(h){colSums(v_ijh[,,h])}))
        lambda=xplusj+t(sapply(1:r,function(h){rep(e[h]-1,d)}))
        theta$lambda=lambda/matrix(rep(1,r),ncol=1)%*%apply(lambda,MARGIN = 2,sum)
        Res=MultBlocXemAlpha_kl(v_ijh,theta,ini$t_jl,a=a,b=b,C=C,e=e,niter=niter,eps=eps,epsi=epsi,epsi_int=epsi_int,niter_int=niter_int)
      }
    }
  }
  return(Res)
}