BinBlocInfo=function(res){
  z=res$z
  w=res$w
  g=max(z)
  m=max(w)
  resumz=matrix(0,g,1)
  resumw=matrix(0,m,1)
  for (i in 1:g){
    resumz[i]=sum(res$z==i)
  }
  for (j in 1:m){
    resumw[j]=sum(res$w==j)
  }
  
  list(Repz=resumz,Repw=resumw)
}