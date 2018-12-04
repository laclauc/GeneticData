# -------------
# Copyright (c) 2018 V. Brault, C. Laclau
# -------------
#
# -------------
# License
# The code can be used for research purposes only.

BinBlocVisuResum=function(data,z,w){
  ## Initialisation
  require(ade4)
  g=max(z);
  m=max(w);
  n=dim(data$x)[1]
  d=dim(data$x)[2]
  z_ik=((z%*%t(rep(1,g)))==(rep(1,n)%*%t(1:g)));
  w_jl=((w%*%t(rep(1,m)))==(rep(1,d)%*%t(1:m)));
  
  ## Display
  nk=colSums(z_ik);
  nl=colSums(w_jl);
  Nkl=nk%*%t(nl);
  nkl=t(z_ik)%*%data$x%*%w_jl;
  
  table.paint(nkl/Nkl)
  
}