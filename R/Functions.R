#' This function calculate the sqrt root of a pd matrix
#'
#' @param S PD matrix
#'
#' @return Returns a aqrt root of the matrix
#'
#' @export


S_half=function(S)
{
  e=eigen(S)
  V=e$vectors
  V %*% diag(sqrt(pmax(e$values,0))) %*% t(V)
}

#' This function generate Multivariate Normal with mean 0 and given covariance matrix
#'
#' @param B is the number of iid gaussian vectors we want to return;
#' @param m is the size of each gaussian vector;
#' @param cov is the given covariance matrix
#' @import stats
#'
#' @return Returns a matrix with multivariate normal with mean 0 and given covariance matrix
#'
#' @export


gauss_vec_0_mean_given_cov=function(B,m,cov)
{
  sigihalf=S_half(cov)
  z_B=matrix(rnorm(B*m,0,1),nrow=m)
  GV_in_cols=sigihalf%*%z_B
  return(GV_in_cols)
}

#' This function generate Gaussian Vector with mean 0 and AR1 covariance.
#' X(t)~N(0,sigma2*rho^{ct}) i.e. the variance component is invariant of t
#'
#' @param t grid points for evaluation
#' @param sigma2 marginal variance
#' @param rho correlation coefficient
#' @param c correlation coefficient parameter
#' @import stats
#'
#' @return Returns The returned vector is X(t_1), X(t_2), ..., X(t_pp)
#'
#' @export


gauss_vec_0_mean_ar1_cov=function(t,sigma2,rho,c)
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


#' This function Epanechnikov Kernel function
#'
#' @param x vector for evaluation
#' @param h bandwidth
#'
#' @return Returns the kernel function values
#'
#' @export


ker=function(x,h)
{
  ans=x
  lo=(abs(x)<h)
  ans[lo]=(3/4)*(1-(x[lo]/h)^2)/h
  ans[!lo]=0
  return(ans)
}


#' This function Local Linear Kernel Estimator vector of beta(t0) (p).
#'
#' @param t0 time point for evaluation
#' @param X is a covariates matrix (N*p),
#' @param y is the response vector (N),
#' @param t is the time vector (N)
#'          For one variable, the longitudinals are vetorized together.
#' @param h is the seleted bandwidth.
#' @param mv is a vector (N) telling us how many repeated measurements for each individual
#'
#' @return Returns betahat(t0)
#'
#' @export


betahat=function(t0,X,y,t,h,mv)
{
  D=cbind(X,as.vector(((t-t0)/h))*X)
  p=ncol(X)
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



