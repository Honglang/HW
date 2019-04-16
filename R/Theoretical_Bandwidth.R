################## Theoretical Bandwidth selection (6/22/2014)###################################
#    Since it's theoretical, i.e. we assume that we know the true model, we need to
#    modify the whole script when we have different simulation setup. And it will 
#    be never used in the real data analysis.
#    Here the true parameters are:
#    Y(t)=0.5*sin(t)*(1+2*exp(t)+N(0,1))+0.5*sin(t)*(3-4*t^2+N(0,1))+N(0,sigma2*rho^(100*t))
######################################################################################


#    Bias of the local linear kernel estimation of the coefficient functions betahat(t0)
bias=function(t0,h,X,t,mv)
{
  D=cbind(X,(t-t0)*X)
  K=(1/mv)*ker(t-t0,h)
  l=(t-t0)^2*X
  L=matrix(0,ncol=2*p,nrow=2*p)
  R=matrix(0,ncol=p,nrow=2*p)
  for(i in 1:(2*p))
  {
    for( j in 1:(2*p))
    {
      L[i,j]=sum(K*D[,i]*D[,j])
    }
  }
  
  for(i in 1:(2*p))
  {
    for( j in 1:p)
    {
      R[i,j]=sum(K*(D[,i]*l[,j]))
    }
  }
  #  Second deravative of the truth coefficient functions beta(t), which need to be estimated
  #  in estimated bandwidth selection procedure.
  beta.2=matrix(c(-0.5*sin(t0),-0.5*sin(t0)),ncol=1)
  betahat0=solve(L)%*%R%*%beta.2
  return(betahat0[1:p])
}

# biast=t(apply(as.matrix(t),1,function(x) bias(x,h, X, t, mv)))
# plot(t,biast[,1])
# plot(t,biast[,2])



#    The following function is used to calculate the covariance of local linear kernel
#    estimator betahat(t0).
subject.terms=function(id,V,t,t0,h,m)
#    id: the covariance matrix of betahat(t0) (p*p) is from the first corder of a (2p*2p)
#    covariance matrix of betahat and the first derivative of betahat. When we calculate
#    this (2p*2p) matrix, we calculate 4 corners separately. And id indicates which
#    corner we are calculating.
#    V: the covariance matrix of error for one individial times the X[i,j]%*%t(X[i,k]),
#    where i indicates the individual thus it's fixed here, X[i,j] is the j-th covariate
#    vector for the individual;  
{
  K=(1/m)*ker(t-t0,h)
  KK=(matrix(rep(K,length(K)),nrow=length(K)))*t(matrix(rep(K,length(K)),nrow=length(K)))
  if(id==1)
  {
    term=sum(V*KK)
  }
  if(id==2)
  {
    term=sum(V*KK*t(matrix(rep((t-t0),length(t)),nrow=length(t))))
  }
  if(id==3)
  {
    term=sum(V*KK*matrix(rep((t-t0),length(t)),nrow=length(t)))
  }
  if(id==4)
  {
    tt0=(matrix(rep((t-t0),length(t)),nrow=length(t)))*t(matrix(rep((t-t0),length(t)),nrow=length(t)))
    term=sum(V*KK*tt0)
  }
  return(term)
}

#    The covariance matrix of betahat(t0) (p*p)
cov=function(t0,h,X,t,n,p,mvinfo,mv)
{
  D=cbind(X,(t-t0)*X)
  K=(1/mv)*ker(t-t0,h)
  L=matrix(0,ncol=2*p,nrow=2*p)
  R=matrix(0,ncol=2*p,nrow=2*p)
  for(i in 1:(2*p))
  {
    for( j in 1:(2*p))
    {
      R[i,j]=sum(K*D[,i]*D[,j])
    }
  }
  
  #  The calculation of the four corners of the (2p*2p) covariance matrix of 
  #  betahat and the first derivative of betahat.
  for(i in 1:p)
  {
    for( j in 1:p)
    {
      count=0
      for(k in 1:n)
      {
        Xi=X[(count+1):(count+mvinfo[k]),i]
        Xj=X[(count+1):(count+mvinfo[k]),j]
        tk=t[(count+1):(count+mvinfo[k])]
        XX=(matrix(rep(Xi,mvinfo[k]),nrow=mvinfo[k]))*t(matrix(rep(Xj,mvinfo[k]),nrow=mvinfo[k]))
        tt=(matrix(rep(tk,mvinfo[k]),nrow=mvinfo[k]))-t(matrix(rep(tk,mvinfo[k]),nrow=mvinfo[k]))
        #SD=sqrt((sin(2*pi*t))^2)
        #SDSD=matrix(rep(SD,length(t)),nrow=length(t))*t(matrix(rep(SD,length(t)),nrow=length(t)))
        #V=XX*(rho^(abs(tt)))*SDSD
        V=XX*(rho^(abs(100*tt)))
        L[i,j]=L[i,j]+subject.terms(1,V,tk,t0,h,mvinfo[k])
        L[i,j+p]=L[i,j+p]+subject.terms(2,V,tk,t0,h,mvinfo[k])
        L[i+p,j]=L[i+p,j]+subject.terms(3,V,tk,t0,h,mvinfo[k])
        L[i+p,j+p]=L[i+p,j+p]+subject.terms(4,V,tk,t0,h,mvinfo[k])
        count=count+mvinfo[k]
      }
    }
  }
  
  #  Get the covariance of betahat(t0)
  Kr=cbind(diag(rep(1,p)),diag(rep(0,p)))
  cov=Kr%*%solve(R)%*%L%*%solve(R)%*%t(Kr)
  return(cov)
}


#Covariance=t(apply(as.matrix(t),1,function(x) cov(x,h,X,t,n,p,mvinfo,mv)))
#plot(t,Covariance[,1])

#     Mean Square Error for betahat(t0)
mset0=function(t0,h,X,t,n,p,mvinfo,mv)
{
  b=bias(t0,h,X,t,mv)
  c=cov(t0,h,X,t,n,p,mvinfo,mv)
  mse=t(b)%*%b+sum(diag(c))
  return(mse)
}

#    Integrated Mean Square Error for betahat
mse=function(h,X,t,n,p,mvinfo,mv)
{
  integrand=apply(as.matrix(seq(0.01,1,0.01)),1,function(x) mset0(x,h,X,t,n,p,mvinfo,mv))
  0.01*sum(integrand)
}

#    Theoretical bandwidth selection
h_theoretical=function(X,t,n,p,mvinfo,mv)
{
  h_temp=as.matrix(seq(0.2,1.0,by=0.3))
  MSE=apply(h_temp,1,function(x) mse(x,X,t,n,p,mvinfo,mv))
  plot(h_temp,MSE)
  h=h_temp[which(MSE==min(MSE))] 
  #h_the=function(h)
  #{
  #  mse(h,X,t,n,p,mvinfo,mv)
  #}
  #h=nlm(h_the,0.3)$estimate
  return(h)
}

##################The followings will not be used here, just standby###################
#     This is second moment function of our random covariates X(t).
omiga=function(t0)
{
  omiga=matrix(c(1+(1+2*exp(t0))^2,(1+2*exp(t0))*(3-4*t0^2),(1+2*exp(t0))*(3-4*t0^2),1+(3-4*t0^2)^2),nrow=2)
  return(omiga)
}

#     Omigat=t(apply(as.matrix(t),1,function(x) omiga(x)))

#     Mean Square Error for betahat(t0) but with X
mset0_with_X=function(t0,h,X,t,m)
{
  b=bias(t0,h,X,t,mv)
  o=omiga(t0)
  c=cov(t0,h,X,t,n,p,mvinfo,mv)
  mse=t(b)%*%o%*%b+sum(diag(o%*%c))
  return(mse)
}
