# -------------
# Copyright (c) 2018 V. Brault, C. Laclau
# -------------
#
# -------------
# License
# The code can be used for research purposes only.

BinBlocICL =function (data,z1,w1,a=4,b=1,C=c(1,1),e=c(1,1)){
  
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
    if (m==1){
      warning("Only trash class")
      return(-Inf)
    }
    w=!(w1%*%t(rep(1,m))-rep(1,d)%*%t(1:m));
  }else{
    w=w1;
  }
  if (ncol(w)>2){
    w=as.matrix(w[,c(TRUE,colSums(w[,-1])!=0)])
  }
  m=ncol(w)-1
  
  nk=colSums(z);
  nl=colSums(w);
  Nkl=nk%*%t(nl[-1]);
  nkl=t(z)%*%x%*%w[,-1];
  xplusj=colSums(x)
  critere=lgamma(g*a)-g*lgamma(a)-lgamma(n+g*a)+sum(lgamma(nk+a))+
    lgamma(sum(C))-sum(lgamma(C))+lgamma(m*a)-m*lgamma(a)+
    lgamma(d-nl[1]+C[1])+lgamma(nl[1]+C[2])-lgamma(d+sum(C))+
    sum(lgamma(nl[-1]+a))-lgamma(d-nl[1]+m*a)+
    nl[1]*(lgamma(sum(e))-sum(lgamma(e))-lgamma(n+sum(e)))+
    sum(w[,1]*(lgamma(xplusj+e[1])+lgamma(n-xplusj+e[2])))+
    g*m*(lgamma(2*b)-2*lgamma(b))+
    sum(lgamma(nkl+b)+lgamma(Nkl-nkl+b)-lgamma(Nkl+2*b))
  critere
}