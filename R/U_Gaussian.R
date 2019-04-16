K2=function(u)
{
  v=abs(u)
  9/16*(v<2)*((1-v^2)*(2-v)-v*((1-v)^2-1)+1/3*(v^2-2)*((1-v)^3+1)+1/2*v*((1-v)^4-1)+1/5*((1-v)^5+1))
}

#plot(seq(-1,1,by=0.01),apply(matrix(seq(-1,1,by=0.01)),1,K2))

Gaussvec=function(m, n, h, mu)
{
  t=seq(0.01,1,by=0.01)
  p=length(t)
  x=matrix(0,p,n)
  lam=rep(0,m)
  lam[1:(m/2+1)]=(5/3)*apply((1/(p*h))*matrix((0:(m/2))),1,K2)
  lam[(m/2+2):m]=(5/3)*apply((1/(p*h))*matrix(((m/2-1):1)),1,K2)
  lam0=fft(lam);
  set.seed(1989+as.integer(100*h))
  urv=matrix(rnorm(m*n,0,1),m,n)
  aseq=matrix(0,m,n)
  aseq[1,]=sqrt(lam0[1])*urv[1,]/sqrt(m);
  aseq[m/2+1,]=sqrt(lam0[m/2+1])*urv[m/2+1,]/sqrt(m);
  aseq[2:(m/2),]=sqrt(lam0[2:(m/2)])*(urv[2:(m/2),]+1i*urv[(2+m/2):m,])/sqrt(2*m);
  aseq[m:(m/2+2),]=sqrt(lam0[2:(m/2)])*(urv[2:(m/2),]-1i*urv[(2+m/2):m,])/sqrt(2*m)
  aseq0=mvfft(aseq);
  x=Re(aseq0)[1:p,]+mu            
  return(x);
}
## Here n is the number of realizations of 100-dim random vector
## Gaussvec(2^8,5,0.3,0) will produce 5 columns, which corresponds to 5 realizations 
## of 100-dim random vectors.

#U=Gaussvec(2^10,100000,0.7,0)
#hist(U[1,],breaks=100,probability=TRUE)
#plot((U[,3])^2)
#T_standard=2*apply(U^2,2,mean)
#hist(T_standard,breaks=100,probability=TRUE)
#TQ=quantile(T_standard,c(0.99,0.98,0.97,0.96,0.95,0.9,0.85,0.8))
#write.table(T_standard, file = "T_standard_0.7.csv", sep = ",", col.names = NA)

#dnor=function(x)
#{
#  (1/(sqrt(2*pi)))*exp(-x^2/(2))
#}
#dchi=function(x,k=2)
#{
#  x^(k/2-1)*exp(-x/2)/(2^(k/2)*gamma(k/2))
#}
#curve(dnor,add=TRUE)
#curve(dchi,add=TRUE)