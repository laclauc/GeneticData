MultBlocGibbs=function(data,g,m,a=4,b=1,C=rep(1,2),e=rep(1,max(data$x)),niter=1000,ini=list(ale=TRUE),eps=10^(-16)){
  r=max(data$x)
  v_ijh=array(unlist(lapply(1:r,function(h){data$x==h})),c(n,d,r))
  if (is.null(ini$theta)){
    init=Init_Mult(g,m,data$x,ale=ini$ale)
  }else{
    init=ini
  }
  res=BlocGibbsMult(v_ijh,init$theta,a=a,b=b,C=C,e=e,init$t_jl,niter=niter,eps=eps)
  # g=dim(res$theta_path$pi_k)[1]
  # m=dim(res$theta_path$rho_l)[1]
  # theta=list(pi_k=matrix(0,g,1),rho_l=matrix(0,m,1),alpha_kl=matrix(0,g,m))
  # if (classement){
  #   for (i in 1:niter){
  #     thetabis=list(alpha_kl=res$theta_path$alpha_kl[[i]],pi_k=res$theta_path$pi_k[,i],rho_l=res$theta_path$rho_l[,i])
  #     thetabisorg=BinBlocClassement(thetabis)
  #     theta$pi_k=theta$pi_k+thetabisorg$pi_k
  #     theta$rho_l=theta$rho_l+thetabisorg$rho_l
  #     theta$alpha_kl=theta$alpha_kl+thetabisorg$alpha_kl
  #   }
  # }else{
  #   for (i in 1:niter){
  #     theta$alpha_kl=theta$alpha_kl+res$theta_path$alpha_kl[[i]]
  #   }
  #   theta$pi_k=rowSums(res$theta_path$pi_k)
  #   theta$rho_l=rowSums(res$theta_path$rho_l)
  # }
  # theta$pi_k=theta$pi_k/niter
  # theta$rho_l=theta$rho_l/niter
  # theta$alpha_kl=theta$alpha_kl/niter
  z=apply(res$Part$z, 1, which.max)
  w=apply(res$Part$w, 1, which.max)
  list(theta=res$theta,z=z,w=w)
}