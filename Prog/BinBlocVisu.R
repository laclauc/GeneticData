# -------------
# Copyright (c) 2018 V. Brault, C. Laclau
# -------------
#
# -------------
# License
# The code can be used for research purposes only.

BinBlocVisu=function(data,z,w,lwd=2){
  ## Initialisation
  require(ade4)
  g=max(z);
  m=max(w);
  n=dim(data$x)[1]
  d=dim(data$x)[2]
  ii=c()
  jj=c()
  for (i in 1:g){
    ii=c(ii,which(z==i))
  }
  for (j in 1:m){
    jj=c(jj,which(w==j))
  }
  z_ik=((z%*%t(rep(1,g)))==(rep(1,n)%*%t(1:g)));
  w_jl=((w%*%t(rep(1,m)))==(rep(1,d)%*%t(1:m)));
  
  ## Display
  n_k=n-sort(cumsum(colSums(z_ik))-0.5,decreasing=TRUE);
  n_l=cumsum(colSums(w_jl))+0.5;
  
  # table.paint(data$x[ii,jj],clegend=0)
  truc=data$x[ii,jj]
  par(mar=c(0,0,0,0))
  image(t(truc[n-(0:(n-1)),]),col = gray(seq(0,1, length=256)),yaxt="n",xaxt="n",
        x = 1:d,y=1:n,xlab="",ylab="")
  
  for (i in n_k[2:g]){
    lines(c(n_l[1],(d+1)),c(i,i),col='blue',lwd=lwd)
  }
  for (i in n_l[2:(m-1)]){
    lines(c(i,i),c(0.5,(n+1)),col='blue',lwd=lwd)
  }
  lines(c(n_l[1],n_l[1]),c(0.5,(n+1)),col='red',lwd=lwd)
  
}