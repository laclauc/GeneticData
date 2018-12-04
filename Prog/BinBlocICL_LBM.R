# -------------
# Copyright (c) 2018 V. Brault, C. Laclau
# -------------
#
# -------------
# License
# The code can be used for research purposes only.

BinBlocICL_LBM =function (data,z1,w1,a=4,b=1){
  
  x=data$x;
  n=dim(x)[1]
  d=dim(x)[2]
  if (!is.matrix(z1)){
    g=max(z1);
    z=!(z1%*%t(rep(1,g))-rep(1,n)%*%t(1:g));
  }else{
    z=z1;
  }
  z=as.matrix(z[,colSums(z)!=0])
  g=ncol(z)
  if (!is.matrix(w1)){
    m=max(w1);
    w=!(w1%*%t(rep(1,m))-rep(1,d)%*%t(1:m));
  }else{
    w=w1;
  }
  m=ncol(w)
  
  nk=colSums(z);
  nl=colSums(w);
  Nkl=nk%*%t(nl);
  nkl=t(z)%*%x%*%w;
  xplusj=colSums(x)
  critere=lgamma(g*a)+lgamma(m*a)-(m+g)*lgamma(a)+m*g*(lgamma(2*b)-2*lgamma(b))-
    lgamma(n+g*a)-lgamma(d+m*a)+
    sum(lgamma(nk+a))+sum(lgamma(nl+a))+
    sum(lgamma(nkl+b)+lgamma(Nkl-nkl+b)-lgamma(Nkl+2*b))
  critere
}