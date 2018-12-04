MultBlocClassement=function(theta,Part=list()){
  thetares=theta
  ind_i=order(theta$alpha_kl[,,1]%*%theta$tau_l)
  ind_j=order(t(theta$pi_k)%*%theta$alpha_kl[,,1])
  thetares$pi_k=theta$pi_k[ind_i]
  thetares$tau_l=theta$tau_l[ind_j]
  thetares$alpha_kl=theta$alpha_kl[ind_i,ind_j,]
  if (is.null(Part$z)||is.null(Part$w)){
    return(thetares)
  }else{
    if (is.vector(Part$z)){
      z=matrix(0,length(Part$z),1)
      maxz=max(Part$z)
      for (i in 1:maxz){
        z[Part$z==i]=ind_i[i]
      }
      w=Part$w
      ind_j=ind_j+1
      maxw=max(Part$w)-1
      for (j in 1:maxw){
        w[Part$w==(j+1)]=ind_j[j]
      }
    }else{
      z=matrix(0,nrow(Part$z),ncol(Part$z))
      for (i in 1:ncol(Part$z)){
        z[,i]=Part$z[,ind_i[i]]
      }
      w=Part$w
      ind_j=ind_j+1
      for (j in 2:ncol(Part$w)){
        w[,j]=Part$w[,ind_j[j-1]]
      }
    }
    list(theta=thetares,z=z,w=w)
  }
}