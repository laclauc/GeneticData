# -------------
# Copyright (c) 2018 V. Brault, C. Laclau
# -------------
#
# -------------
# License
# The code can be used for research purposes only.

library(permute)

SimpleError=function(zs,z){
  g=max(length(unique(zs)),length(unique(z)))
  n=length(z)
  if (max(z)>g){
    val=sort(unique(z),decreasing = FALSE)
    for (it in 1:length(val)){
      if (val[it]!=it){
        z[z==val[it]]=it
      }
    }
  }
  zsik=!(zs%*%t(rep(1,g))-rep(1,n)%*%t(1:g));
  zik=!(z%*%t(rep(1,g))-rep(1,n)%*%t(1:g));
  if (g<7){
    P<-permute::allPerms(g)
  }else{
    P1<-permute::allPerms(6)
    if (g>=7){
      P=cbind(rep(7,nrow(P1)),P1)
      for (it in 2:7){
        P2=matrix(7,nrow=nrow(P1),ncol=7)
        P2[,-it]=P1
        P<-rbind(P,P2)
      }
      if (g>=8){
        P1=P
        P=cbind(rep(8,nrow(P1)),P1)
        for (it in 2:8){
          P2=matrix(8,nrow=nrow(P1),ncol=8)
          P2[,-it]=P1
          P<-rbind(P,P2)
        }
        if (g==9){
          P1=P
          P=cbind(rep(9,nrow(P1)),P1)
          for (it in 2:9){
            P2=matrix(9,nrow=nrow(P1),ncol=9)
            P2[,-it]=P1
            P<-rbind(P,P2)
          }
        }
      }
    }
  }
  Crit=Inf
  for (iter in 1:nrow(P)){
    Val=sum(abs(zik-zsik[,P[iter,]]))/(2*n)
    if (Val<Crit){
      Crit=Val
      if (Val==0){
        break
      }
    }
  }
  Crit
}


## function defining the CARI
CoError=function(z,w,zprime,wprime){
  ez=SimpleError(z,zprime)
  ew=SimpleError(w,wprime)
  ez+ew-ez*ew
}