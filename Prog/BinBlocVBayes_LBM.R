# -------------
# Copyright (c) 2018 V. Brault, C. Laclau
# -------------
#
# -------------
# License
# The code can be used for research purposes only.

BinBlocVBayes_LBM=function(data,g,m,a=4,b=1,ini=list(n_init=20,ale=TRUE),niter=1000,eps=10^(-16),epsi=10^(-16),epsi_int=0.1,niter_int=1,verbatim=TRUE){
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
    if (verbatim){
      pb=txtProgressBar(min=1,max=ini$n_init,style=3)
    }
    for (i in 1:ini$n_init){
      init=Init_LBM(g,m,data$x,ale=ini$ale)
      ResVBayesAle=BlocXemAlpha_kl_LBM(data$x,init$theta,init$t_jl,a=a,b=b,niter=niter,eps=eps,epsi=epsi,epsi_int=epsi_int,niter_int=niter_int)
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
      if (verbatim){
       setTxtProgressBar(pb,i)
      }
    }
  }else{
    if ((any(names(ini)=="t_jl"))&&(any(names(ini)=="theta"))){
      Res=BlocXemAlpha_kl_LBM(data$x,ini$theta,ini$t_jl,a=a,b=b,niter=niter,eps=eps,epsi=epsi,epsi_int=epsi_int,niter_int=niter_int)
    }else{
      if ((any(names(ini)=="t_jl"))&&(any(names(ini)=="s_ik"))){ 
        t_l=t(colSums(ini$t_jl));
        s_k=t(colSums(ini$s_ik));
        n=nrow(data$x)
        d=ncol(data$x)
        theta=list()
        # Computation of parameters
        theta$pi_k=t((s_k)/(n));
        theta$tau_l=t((t_l)/(sum(t_l)));
        u_kl=t(ini$s_ik)%*%data$x%*%ini$t_jl;
        n_kl=t(s_k)%*%t_l;
        u_kl=u_kl*(u_kl>0);
        n_kl=n_kl*(n_kl>0);
        theta$alpha_kl=u_kl/n_kl;
        if (any(is.nan(theta$alpha_kl))){
          theta$alpha_kl=u_kl/(n_kl*(n_kl>0)+eps*(n_kl<=0));
        }
        theta$lambda=colSums(data$x)/n
        Res=BlocXemAlpha_kl_LBM(data$x,theta,ini$t_jl,a=a,b=b,niter=niter,eps=eps,epsi=epsi,epsi_int=epsi_int,niter_int=niter_int)
      }
    }
  }
  return(Res)
}