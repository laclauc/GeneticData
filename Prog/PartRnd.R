PartRnd = function(n,proba){

g=length(proba)

i=sample(1:n,n, replace = FALSE);
i1=i[1:g];
i2=i[(g+1):n];
z=rep(0,n);
z[i1]=1:g;

if (n > g) {
  z[i2]=(g+1)-colSums(rep(1,g)%*%t(runif(n-g)) < cumsum(proba)%*%t(rep(1,n-g)));
  z[z==(g+1)]=g
} 

z_ik=!(z%*%t(rep(1,g))-rep(1,n)%*%t(1:g));

list(z=z,z_ik=z_ik)
}