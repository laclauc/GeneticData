# -------------
# Copyright (c) 2018 V. Brault, C. Laclau
# -------------
#
# -------------
# License
# The code can be used for research purposes only.

ARI=function(z,zprime){
  L_z=length(z)
  if (L_z!=length(zprime)){
    print('Both partitions must contain the same number of points.')
  } 
  N_z=max(z)  
  N_zprime=max(zprime)
  nz=matrix(0,N_z,N_zprime)
  for (i in 1:N_z){
    for (j in 1:N_zprime){
      G1 = which(z==i)
      G2 = which(zprime==j)
      nz[i,j] = length(intersect(G1,G2))
    }
  }
  
  ssm = 0
  sm1 = 0
  sm2 = 0
  for (i in 1:N_z){
    for (j in 1:N_zprime){
      ssm = ssm + choose(nz[i,j],2)
    }
  }
  temp = rowSums(nz)
  for (i in 1:N_z){
    sm1 = sm1 + choose(temp[i],2)
  }
  temp = colSums(nz)
  for (i in 1:N_zprime){
    sm2 = sm2 + choose(temp[i],2)
  }
  NN = ssm - (sm1*sm2)/choose(L_z,2)
  DD = (sm1 + sm2)/2 - (sm1*sm2)/choose(L_z,2)
  ari = NN/DD
  return(list(ari=ari, nz=nz))
  
}


## function defining the CARI
CARI=function(z,w,zprime,wprime){
  L_z=length(z)
  L_w=length(w)
  if (L_z!=length(zprime)){
    print('Both partitions must contain the same number of points.')
  }
  if (L_w!=length(wprime)){
    print('Both partitions must contain the same number of points.')
  }
  
  N_z=max(z)
  N_w=max(w)
  N_zprime=max(zprime)
  N_wprime=max(wprime)
  ari_z=ARI(z,zprime) 
  ari_w=ARI(w,wprime)
  nz=ari_z[[2]]
  nw=ari_w[[2]]
  nzw=kronecker(nz,nw)
  ssm = 0
  sm1 = 0
  sm2 = 0
  for (i in 1:(N_z*N_w)){
    for (j in 1:(N_zprime*N_wprime)){
      ssm = ssm + choose(nzw[i,j],2)
    }
  }
  temp1=rowSums(nz)
  temp2=rowSums(nw)
  temp=kronecker(temp1,temp2)
  for (i in 1:(N_z*N_w)){
    sm1 = sm1 + choose(temp[i],2)
  }
  temp1=colSums(nz)
  temp2=colSums(nw)
  temp = kronecker(temp1,temp2)
  for (i in 1:(N_zprime*N_wprime)){
    sm2 = sm2 + choose(temp[i],2)
  }
  NN = ssm - sm1*sm2/choose(L_z*L_w,2)
  DD = (sm1 + sm2)/2 - sm1*sm2/choose(L_z*L_w,2)
  cari = NN/DD
  return(list(cari=cari, nzw=nzw))
}