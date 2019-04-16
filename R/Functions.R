#   Square root of a positive definite matrix S
S_half=function(S)
{
  e=eigen(S)
  V=e$vectors
  V %*% diag(sqrt(pmax(e$values,0))) %*% t(V)
}

#    Multivariate Normal with mean 0 and given covariance matrix
gauss_vec_0_mean_given_cov=function(B,m,cov)
  #    B is the number of iid gaussian vectors we want to return;
  #    m is the size of each gaussian vector;
  #    cov is the given covariance matrix  
{
  sigihalf=S_half(cov)
  z_B=matrix(rnorm(B*m,0,1),nrow=m)
  GV_in_cols=sigihalf%*%z_B
  return(GV_in_cols)
}


#   Gaussian Vector with mean 0 and AR1 covariance.
gauss_vec_0_mean_ar1_cov=function(t,sigma2,rho,c)
#   X(t)~N(0,sigma2*rho^{ct}) i.e. the variance component is invariant of t
#   and the correlation function is just rho to the power ct.
#   The returned vector is X(t_1), X(t_2), ..., X(t_pp)
{
  pp=length(t)
  TT=abs(matrix(rep(t,pp),nrow=pp)-t(matrix(rep(t,pp),nrow=pp)))*c
#   If you want to make the variance component as an arbitrary positive function of t,
#   you can specify your function here, e.g. sigma^2(t)=(sin(2*pi*t))^2 in the following.
#   SD=sqrt((sin(2*pi*t))^2)
#   SDSD=matrix(rep(SD,pp),nrow=pp)*t(matrix(rep(SD,pp),nrow=pp))
#   Sig=rho^(TT)*SDSD
  Sig=rho^(TT)
  Sig.half=S_half(Sig)
  Z=rnorm(pp)
  GV=Sig.half%*%Z
}



#   Epanechnikov Kernel function: K_h(u)=K(u/h)/h where K(u) = \frac{3}{4}(1-u^2) \,\mathbf{1}_{\{|u|\leq1\}}
ker=function(x,h)
{
  ans=x
  lo=(abs(x)<h)
  ans[lo]=(3/4)*(1-(x[lo]/h)^2)/h
  ans[!lo]=0
  return(ans)
}

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


#   The estimating equations evaluated at beta(t0).
g_beta_t0=function(betat0,t0_info,betahatt,X,y,t,n,p,h,mvinfo,mv,blanced)
#    betat0: is the argument for the estimating function;
#    t0_info: is the kernel estimation of beta(t) at t0, which is used for the bias correction
#    betahatt: in the kernel estimation for beta(t) at all observed t for all individuals;
#    X, y, t: are the data; n is the sample size; p is dimension;
#    h: selected bandwidth; 
#    mvinfo: is a vector (n) telling us how many repeated measurements for each individual   
#    mv: is a vector (N) telling us how many repeated measurements for each individual  
#    blanced: 0 means repeated measurements are unbalanced and 1 means they are balanced;
#    We are returned with a n*p estimating equations matrix.  
{
  t0=t0_info[1]
  betahat0=as.matrix(t0_info[-1])
  z=y-X%*%as.matrix(betat0,ncol=1)
#   #Residual-adjusted  
#   z=y-X%*%as.matrix(betat0,ncol=1)+X%*%betahat0
#   for(j in 1:p){
#     z=z-X[,j]*betahatt[,j]
#   }
  g_betaijt0=X*(as.vector(z*ker(t-t0,h)))
  g_betat0=matrix(NA,ncol=p,nrow=n);
  #  blanced version is more compuational efficient
  if(blanced==1)
  {
    m=mv[1]
    for(i in 1:p)
    {
      g_betat0[,i]=apply(matrix(as.vector(g_betaijt0[,i]),nrow=m,ncol=n),2,mean)
    }
  }
  #  unbalanced version needs more computation due to more for loops
  if(blanced==0)
  {
    count=0
    for(i in 1:n)
    {
      for(j in 1:p)
      {
        g_betat0[i,j]=mean(g_betaijt0[(count+1):(count+mvinfo[i]),j])
      }
      count=count+mvinfo[i]
    }
  }
  return(g_betat0)
}

#   # The recent empirical likelihood R package could be used 
#   library("el.convex")

#   Empirical Likelihood evaluated at beta(t0).
l_beta_t0=function(betat0,t0_info,betahatt,X,y,t,n,p,h,mvinfo,mv,blanced)
{
  Z=g_beta_t0(betat0,t0_info,betahatt,X,y,t,n,p,h,mvinfo,mv,blanced)
  plot(Z[,1])
  points(Z[,2],col="red")
  points(Z[,3],col="green")
  mu=rep(0,p)
  l_owen=-elm(Z,mu)$logelr
#   The following 4 other ways for calculating EL are from the library("el.convex").  
#   l_newton=0.5*(el.test.newton(Z,mu)$"-2LLR")
#   l_bfgs=0.5*(el.test.bfgs(Z,mu)$"-2LLR")
#   l_damped=0.5*(el.test.damped(Z,mu)$"-2LLR")
#   l_dfp=0.5*(el.test.dfp(Z,mu)$"-2LLR")  
  return(l_owen)
#   return(c(l_owen,l_newton,l_bfgs,l_damped,l_dfp))
  
}
##################################################################################
##################The only function needs to be modified for different setup######
#########################Here is the collection###################################
#1)#########The one used in the example Main_EL.R########
#    Empirical Likelihood evaluated at betatilde(t0).
#    The function is the only function we need to modify for different hypothesis.
# l_beta_tilde=function(t0_info,betahatt,X,y,t,n,p,h,mvinfo,mv,blanced)
# {
#   #  The following function has to be modified for different hypothesis: this is 
#   #  a function at the null hypothesis.
#   lH=function(BETA)
#   {
#     lh=l_beta_t0(c(BETA,BETA),t0_info,betahatt,X,y,t,n,p,h,mvinfo,mv,blanced)
#     return(lh[1])
#   }
#   
#   # The folloiwng optimizaiton is also under the null hypothesis. And for some hypothesis
#   # we may not need to do optimization.
#   t0=t0_info[1]
#   l_beta_tilde_t0=golden(lH, 0.5*sin(t0)-0.3,0.5*sin(t0)+0.3,10^(-3))
#   # l_beta_tilde_t0=nlm(lH,0.5*sin(t0)+0.05)
#   return(l_beta_tilde_t0)
# }

l_beta_tilde=function(t0_info,betahatt,X,y,t,n,p,h,mvinfo,mv,blanced)
{
  #  The following function has to be modified for different hypothesis: this is 
  #  a function at the null hypothesis.
  lH=function(BETA)
  {
    lh=l_beta_t0(c(BETA,0,0),t0_info,betahatt,X,y,t,n,p,h,mvinfo,mv,blanced)
    return(lh[1])
  }
  
  # The folloiwng optimizaiton is also under the null hypothesis. And for some hypothesis
  # we may not need to do optimization.
  #t0=t0_info[1]
  #l_beta_tilde_t0=golden(lH, 0, 1, 10^(-3))
  #l_beta_tilde_t0=lH(t0_info[2])
  l_beta_tilde_t0=nlm(lH,0)
  l_beta_tilde_t0=c(l_beta_tilde_t0[[2]],l_beta_tilde_t0[[1]])#now the first one is the arg and the second one is the value
  return(l_beta_tilde_t0)
}


####################################################################################
####################################################################################

Error_info=function(X,y,t,n,p,mv,mvinfo,h_special,h_sig,h_rho)
{
  ###################################################################################
  #####Error information which will be used several times and computational cost#####
  ###################################################################################
  #  The bandwidth h_speical passed in needs to pay attention to. When the error infortion
  #  is used in bandwidth selection, this h_special just h_beta. While when the error information
  #  is used in bootstrapping, this h_speical needs to be the selected bandwidth.
  betahatt=t(apply(as.matrix(t),1,function(x) betahat(x,X,y,t,h_special,mv)))
  r.hat=y
  for(i in 1:p)
  {
    r.hat=r.hat-as.matrix(X[,i]*betahatt[,i])
  }
  # The following variance estimation could be obtained by several methods: kernel or polynomial
  #  Kernel: 
  sigma2.hat=t(apply(as.matrix(t),1,function(x) sigma.t0.hat(r.hat,x,t,h=h_sig)))
  #   #  Constant: 
  #   psi=sigma.t0.hat.constant.par(r.hat)
  #   sigma2.hat=t(apply(as.matrix(t),1,function(x) sigma.t0.hat.constant(x,psi)))
  #   # Quadratic:
  #   psi=sigma.t0.hat.quadratic.par(t,r.hat)
  #   sigma2.hat=t(apply(as.matrix(t),1,function(x) sigma.t0.hat.quadratic(x,psi)))  
  e.hat=as.vector(r.hat)/sqrt(sigma2.hat)
  
  theta0=theta_haat(e.hat,t,n,mvinfo)
  #   phi=rho.t0.poly.hat.par(e.hat,t,n,mvinfo)
  
  NMM=sum(mvinfo^2)
  EV=rep(0,NMM)
  count=0
  count1=0
  for(k in 1:n)
  {
    tk=t[(count+1):(count+mvinfo[k])]
    tt=(matrix(rep(tk,mvinfo[k]),nrow=mvinfo[k]))-t(matrix(rep(tk,mvinfo[k]),nrow=mvinfo[k]))
    SDSD=sqrt((matrix(rep(sigma2.hat[(count+1):(count+mvinfo[k])],mvinfo[k]),nrow=mvinfo[k]))*t(matrix(rep(sigma2.hat[(count+1):(count+mvinfo[k])],mvinfo[k]),nrow=mvinfo[k])))
    #  The estimation of the correlation can be obtained by: kernel, exp, poly
    #     # Kernel:
    #     Rhott=matrix(apply(matrix(as.vector(abs(tt))),1,function(x) rho.t0.hat(e.hat,t,n,x,h=h_rho,mvinfo)),nrow=mvinfo[k])
    # Exp:            
    Rhott=matrix(apply(matrix(as.vector(abs(tt))),1,function(x) rho.t0.exp.hat(x,theta0)),nrow=mvinfo[k])
    #         # Poly:   
    #         Rhott=matrix(apply(matrix(as.vector(abs(tt))),1,function(x) rho.t0.poly.hat(x,phi)),nrow=mvinfo[k])
    Vhat=Rhott*SDSD
    EV[(count1+1):(count1+mvinfo[k]^2)]=as.vector(Vhat)
    count=count+mvinfo[k]
    count1=count1+mvinfo[k]^2
  }
  ####################################################################################
  ####################################################################################
  return(EV)
}
