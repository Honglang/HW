#########################Bandwidth Estimation (6/22/2014)#####################################
# There are 4 pre-determined bandwidths: h_beta in the residuals, h_beta_2 in the bias
# h_rho in the correlation estimation, h_sig in the variance estimation
####################################################################################


#   Local Linear Kernel Estimator vector of beta(t0) (p).
betahat=function(t0,X,y,t,h,mv)
  #   X is a covariates matrix (N*p), y is the response vector (N), t is the time vector (N)
  #   For one variable, the longitudinals are vetorized together.
  #   h is the seleted bandwidth.
  #   mv is a vector (N) telling us how many repeated measurements for each individual
{
  D=cbind(X,as.vector(((t-t0)/h))*X)
  K=(1/mv)*ker(t-t0,h)
  L=matrix(0,ncol=2*p,nrow=2*p)
  R=matrix(0,ncol=1,nrow=2*p)
  for(i in 1:(2*p))
  {
    for( j in 1:(2*p))
    {
      L[i,j]=sum(K*D[,i]*D[,j])
    }
    R[i,1]=sum(K*(D[,i]*y))
  }
  
  betahat0=solve(L)%*%R
  return(betahat0[1:p])
}

#    Estimation of the second derivative of beta(t0)
beta.2.hat=function(t0,h,X,t,y,mv)
#    h here is a precribed bandwidth for the estimation of second derivative of beta(t0)  
{
  D=cbind(X,(t-t0)*X,((t-t0)^2)*X)
  K=(1/mv)*ker(t-t0,h)
  L=matrix(0,ncol=3*p,nrow=3*p)
  R=matrix(0,ncol=1,nrow=3*p)
  for(i in 1:(3*p))
  {
    for( j in 1:(3*p))
    {
      L[i,j]=sum(K*D[,i]*D[,j])
    }
    R[i,1]=sum(K*(D[,i]*y))
  }
  
  alphahat0=solve(L)%*%R
  return(alphahat0[(2*p+1):(3*p)])
}

#    Bias of the local linear kernel estimation of the coefficient functions betahat(t0)
bias.hat=function(t0,h,X,t,y,mv,h_beta_2)
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
  #  Estimated second deravative of beta(t).
  betapphat=matrix(beta.2.hat(t0,h_beta_2,X,t,y,mv),ncol=1)
  betahat0=solve(L)%*%R%*%betapphat
  return(betahat0[1:p])
}

#    Kernel Estimation of error variance
sigma.t0.hat=function(D,t0,t,h)
#    D here is just all of the estimated residuals, no need to differentiate among different individuals  
#    h: the second intitial bandwidth to choose for error variance
{  
  K=ker(t-t0,h)
  sig2.hat=sum(D^2*K)/sum(K)
  return(sig2.hat)
}

#    Polynomial Approximation of the Variance
##########Constant Approximation#################
sigma.t0.hat.constant.par=function(r.hat)
{
  ttt=matrix(rep(1,n*m),ncol=1)
  psi=solve(t(ttt)%*%ttt)%*%t(ttt)%*%r.hat^2
  return(psi)
}


sigma.t0.hat.constant=function(t0,psi)
{
  v=matrix(c(1),ncol=1)
  sig2.hat=t(psi)%*%v
  return(sig2.hat)
}

###########Quadratic Approximation#####################
sigma.t0.hat.quadratic.par=function(t,r.hat)
{
  ttt=cbind(rep(1,n*m),t,t^2)
  psi=solve(t(ttt)%*%ttt)%*%t(ttt)%*%r.hat^2
  return(psi)
}


sigma.t0.hat.quadratic=function(t0,psi)
{
  v=matrix(c(1,t0,t0^2),ncol=1)
  sig2.hat=t(psi)%*%v
  return(sig2.hat)
}



#    The following function is used for error correlation function estimation
rhohat_terms=function(id,Res_info,t0,h)
#    id: indicates which term we are calculating in the estimated correlation function
#    Res_info: includes two parts, the first part is the residual and the second part is the time information  
{
  m=length(Res_info)/2
  X=Res_info[1:m]
  t=Res_info[(1+m):(2*m)]
  #Xbar=mean(X)
  #XX=(matrix(rep(X-Xbar,length(X)),nrow=length(X)))*t(matrix(rep(X-Xbar,length(X)),nrow=length(X)))
  XX=(matrix(rep(X,length(X)),nrow=length(X)))*t(matrix(rep(X,length(X)),nrow=length(X)))
  tau=abs((matrix(rep(t,m),nrow=m))-t(matrix(rep(t,m),nrow=m)))
  K=matrix(ker(as.vector(tau)-t0,h),nrow=m)
  # This is the numerator component for estimation of the correlation.
  if(id==1)
  {
    #term=sum(XX*K)
    term=sum(XX*K)-sum(diag(XX*K))
  }
  #This is the denominator component for estimation of the correlation.
  if(id==0)
  {
    #term=sum(K)
    term=sum(K)-sum(diag(K))
  }
  return(term)
}

#    Kernel Estimation of the correlation function at t0.
rho.t0.hat=function(ehat,t,n,t0,h,mvinfo)
#    ehat is the whole vector of standardized residuals
#    h is the third initial bandwidth to choose for correlation estimation
{
  count=0
  up=0
  down=0
  for(i in 1:n)
  {
    up=up+rhohat_terms(1,c(ehat[(count+1):(count+mvinfo[i])],t[(count+1):(count+mvinfo[i])]),t0,h)
    down=down+rhohat_terms(0,c(ehat[(count+1):(count+mvinfo[i])],t[(count+1):(count+mvinfo[i])]),t0,h)
    count=count+mvinfo[i]
  }
  return(up/down)
}




#    The following two functions are used to calculate the exponential approxiation of the 
#    correlation funciton. We use Newton's method, so we have to have the following 
#    first and second derivatives.
rho_e_fun_p=function(theta,tau,ee)
{
  tautheta=theta[1]*(exp(theta[2]*tau))
  re1=(ee-tautheta)*(exp(theta[2]*tau))
  re2=re1*tau
  return(c(sum(re1),sum(re2)))
}

rho_e_fun_pp=function(theta,tau,ee)
{
  re11=-sum(exp(2*theta[2]*tau))
  re12=sum((ee*(exp(theta[2]*tau))-2*theta[1]*(exp(2*theta[2]*tau)))*tau)
  re21=-sum((exp(2*theta[2]*tau))*tau)
  re22=sum((ee*(exp(theta[2]*tau))-2*theta[1]*(exp(2*theta[2]*tau)))*tau^2)
  re=matrix(c(re11,re21,re12,re22),nrow=2)
  return(re)
}

#     The estimation of the two parameters in the exponential approximation, since it's parametric method
theta_haat=function(e.hat,t,n,mvinfo)
{
  count=0
  ehatmatrix=matrix(rep(e.hat[(count+1):(count+mvinfo[1])],mvinfo[1]),nrow=mvinfo[1])
  ee=matrix(as.vector(ehatmatrix*t(ehatmatrix)),ncol=1)
  tmatrix=matrix(rep(t[(count+1):(count+mvinfo[1])],mvinfo[1]),nrow=mvinfo[1])
  tau=matrix(as.vector(abs(tmatrix-t(tmatrix))),ncol=1)
  for(i in 2:n){
    ehatmatrix=matrix(rep(e.hat[(count+1):(count+mvinfo[i])],mvinfo[i]),nrow=mvinfo[i])
    ee_temp=matrix(as.vector(ehatmatrix*t(ehatmatrix)),ncol=1)
    tmatrix=matrix(rep(t[(count+1):(count+mvinfo[i])],mvinfo[i]),nrow=mvinfo[i])
    tau_temp=matrix(as.vector(abs(tmatrix-t(tmatrix))),ncol=1)
    ee=rbind(ee,ee_temp)
    tau=rbind(tau,tau_temp)
    count=count+mvinfo[i]
  }
  theta0=c(1,0)
  F0=rho_e_fun_p(theta0,tau,ee)
  F0=sqrt(sum(F0*F0))
  while( abs(F0) > 10^(-5) )
  {
    theta1=theta0-solve(rho_e_fun_pp(theta0,tau,ee))%*%rho_e_fun_p(theta0,tau,ee)
    theta0=theta1
    F0=rho_e_fun_p(theta0,tau,ee)
    F0=sqrt(sum(F0*F0))
  }
  
  return(theta0)
}

#    Exponential Estimation of the correlation function at t0.
rho.t0.exp.hat=function(t0,theta0)
{
  theta0[1]*exp(theta0[2]*t0)
}




#    Polynomial Approximation of the correlation function
rho.t0.poly.hat.par=function(e.hat,t,n,mvinfo)
{
  count=0
  ehatmatrix=matrix(rep(e.hat[(count+1):(count+mvinfo[1])],mvinfo[1]),nrow=mvinfo[1])
  ee=matrix(as.vector(ehatmatrix*t(ehatmatrix)),ncol=1)
  tmatrix=matrix(rep(t[(count+1):(count+mvinfo[1])],mvinfo[1]),nrow=mvinfo[1])
  tau=matrix(as.vector(abs(tmatrix-t(tmatrix))),ncol=1)
  for(i in 2:n){
    ehatmatrix=matrix(rep(e.hat[(count+1):(count+mvinfo[i])],mvinfo[i]),nrow=mvinfo[i])
    ee_temp=matrix(as.vector(ehatmatrix*t(ehatmatrix)),ncol=1)
    tmatrix=matrix(rep(t[(count+1):(count+mvinfo[i])],mvinfo[i]),nrow=mvinfo[i])
    tau_temp=matrix(as.vector(abs(tmatrix-t(tmatrix))),ncol=1)
    ee=rbind(ee,ee_temp)
    tau=rbind(tau,tau_temp)
    count=count+mvinfo[i]
  }
  tttau=cbind(rep(1,length(tau)),tau,tau^2)
  phi=solve(t(tttau)%*%tttau)%*%t(tttau)%*%ee
  return(phi)
}

rho.t0.poly.hat=function(t0,phi)
{
  v=matrix(c(1,t0,t0^2),ncol=1)
  rhohat=t(phi)%*%v
  return(rhohat)
}



#    The following function is used to calculate the covariance or 
#    estimated covariance of local linear kernel estimator betahat(t0).
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
cov.hat=function(t0,h,X,Vhat_L,t,n,p,mvinfo,mv)
#    Vhat_L is global information has nothing to do with t0. Thus we need to 
#    calculate outside in order to save computation. And it's about covariance of 
#    error but together with some covariates information.  
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
  count=0
  count1=0
  for(k in 1:n)
  {
    for(i in 1:p)
    {
      for( j in 1:p)
      {
        tk=t[(count+1):(count+mvinfo[k])]    
        Vhat=matrix(Vhat_L[[j+(i-1)*p]][(count1+1):(count1+mvinfo[k]^2)],nrow=mvinfo[k])
        
        L[i,j]=L[i,j]+subject.terms(1,Vhat,tk,t0,h,mvinfo[k])
        L[i,j+p]=L[i,j+p]+subject.terms(2,Vhat,tk,t0,h,mvinfo[k])
        L[i+p,j]=L[i+p,j]+subject.terms(3,Vhat,tk,t0,h,mvinfo[k])
        L[i+p,j+p]=L[i+p,j+p]+subject.terms(4,Vhat,tk,t0,h,mvinfo[k])
      }
    }
    
    count=count+mvinfo[k]
    count1=count1+mvinfo[k]^2
  }
  
  #  Get the estimated covariance of betahat(t0)
  Kr=cbind(diag(rep(1,p)),diag(rep(0,p)))
  cov=Kr%*%solve(R)%*%L%*%solve(R)%*%t(Kr)
  return(cov)
}


mset0.hat=function(t0,h,X,y,t,Vhat_L,n,p,mvinfo,mv,h_beta_2)
{
  b=bias.hat(t0,h,X,t,y,mv,h_beta_2)
  c=cov.hat(t0,h,X,Vhat_L,t,n,p,mvinfo,mv)
  mse=t(b)%*%b+sum(diag(c))
  return(mse)
}

mse.hat=function(h,X,y,t,n,p,mvinfo,mv,EV,h_beta_2)
{
  nmm=sum(mvinfo^2)
  Vhat_L=replicate(p^2, rep(0,nmm), simplify=F)
  count=0
  count1=0
  for(k in 1:n)
  {
    for(i in 1:p)
    {
      for( j in 1:p)
      {
        Xi=X[(count+1):(count+mvinfo[k]),i]
        Xj=X[(count+1):(count+mvinfo[k]),j]
        XX=(matrix(rep(Xi,mvinfo[k]),nrow=mvinfo[k]))*t(matrix(rep(Xj,mvinfo[k]),nrow=mvinfo[k]))
        Vhat_L[[j+(i-1)*p]][(count1+1):(count1+mvinfo[k]^2)]=as.vector(XX*matrix(EV[(count1+1):(count1+mvinfo[k]^2)],nrow=mvinfo[k]))
      }
    }
    
    count=count+mvinfo[k]
    count1=count1+mvinfo[k]^2
  }
  
  integrand=apply(as.matrix(seq(0.01,1,0.01)),1,function(x) mset0.hat(x,h,X,y,t,Vhat_L,n,p,mvinfo,mv,h_beta_2))
  0.01*sum(integrand)
}

bandwidth_est=function(X,y,t,n,p,mvinfo,mv,EV,h_beta_2)
{
  h_temp=as.matrix(seq(1.1,6.1,by=0.5))
  MSE.hat=apply(h_temp,1,function(x) mse.hat(x,X,y,t,n,p,mvinfo,mv,EV,h_beta_2))
  h=h_temp[which(MSE.hat==min(MSE.hat))]
  #h_est=function(h)
  #{
  #  mse.hat(h,X,t,n,p,mvinfo,mv)
  #}
  #h=nlm(h_est,0.3)$estimate
  return(h)
}


####################The followings will not be used here, just standby##########
omiga.hat=function(t0,h,X,t)
{
  K=ker(t-t0,h)
  o11=sum((X[,1]*X[,1])*K)/sum(K)
  o12=sum((X[,1]*X[,2])*K)/sum(K)
  o21=sum((X[,2]*X[,1])*K)/sum(K)
  o22=sum((X[,2]*X[,2])*K)/sum(K)
  omiga=matrix(c(o11,o21,o12,o22),nrow=2)
  return(omiga)
}

mset0.hat.with.X=function(t0,h,X,t,n,p,mvinfo,mv)
{
  b=bias.hat(t0,h,X,t,y,mv)
  o=omiga.hat(t0,h_gamma=0.5,X,t)
  c=cov.hat(t0,h,X,t,n,p,mvinfo,mv)
  mse=t(b)%*%o%*%b+sum(diag(o%*%c))
  return(mse)
}
