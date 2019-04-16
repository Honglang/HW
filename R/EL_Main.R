source("//sttnas/StudentReds/wangho16/My Documents/Research Pool/EL and Functional Data/ELFCODE/R Version/Updated Version/EL_Mean.R")
source("//sttnas/StudentReds/wangho16/My Documents/Research Pool/EL and Functional Data/ELFCODE/R Version/Updated Version/U_Gaussian.R")
source("//sttnas/StudentReds/wangho16/My Documents/Research Pool/EL and Functional Data/ELFCODE/R Version/Updated Version/Bandwidth_Est.R")
source("//sttnas/StudentReds/wangho16/My Documents/Research Pool/EL and Functional Data/ELFCODE/R Version/Updated Version/Theoretical_Bandwidth.R")
source("//sttnas/StudentReds/wangho16/My Documents/Research Pool/EL and Functional Data/ELFCODE/R Version/Updated Version/Golden_Search.R")
source("//sttnas/StudentReds/wangho16/My Documents/Research Pool/EL and Functional Data/ELFCODE/R Version/Updated Version/Functions.R")
source("//sttnas/StudentReds/wangho16/My Documents/Research Pool/EL and Functional Data/ELFCODE/R Version/Updated Version/Bootstrapping.R")



n=100
m=5
p=2
a=0
b=0
rho=0.2
mvinfo=rep(m,n)
N=m*n
mv=rep(m,N)
h_beta=0.4
h_beta_2=0.4
h_sig=0.4
h_rho=0.4

set.seed(1984)

#repnum: number of repeat sampling
repnum=1
#B: number of bootstrapping
B=2
#T, T1: final test statistics, \int_0^1 2l(t)w(t)dt, one is using weight w(t)=1[0,1](t) and the other is using weight w(t)=1.25*1[0.1,0.9](t)
T=rep(0,repnum)
T1=rep(0,repnum)
#TB,T1B: Bootstrapped version of the above test statistics
TB=matrix(NA,nrow=repnum,ncol=B)
T1B=matrix(NA,nrow=repnum,ncol=B)

#logelr: 2l(t) at t=c(0.01,0.02,\cdots, 1)
logelr=matrix(NA,nrow=repnum,ncol=100)
#h_rec: a recod of selected bandwidth for each sample
h_rec=rep(0,repnum)



t=runif(m*n,min=0, max=1)
X1=rnorm(m*n,0,1)
X1=as.matrix(1+2*exp(t)+X1)
X2=rnorm(m*n,0,1)
X2=as.matrix(3-4*t^2+X2)
X=cbind(X1,X2)
beta1=0.5*sin(t)
beta2=(0.5+a)*sin(t+b)

# #Theoretical bandwidth
# h=h_theoretical(X,t,n,p,mvinfo,mv) 

for(rep in 1:repnum)
{
  eps=as.vector(apply(matrix(t,nrow=n,ncol=m,byrow="T"),1,function(x) gauss_vec_0_mean_ar1_cov(x,sigma2=1,rho,c=100)))
  y=as.matrix(X1*beta1+X2*beta2+eps) 
  #plot(t,eps)
  #plot(t,y)
  
  #EV contains the covariance information about the error
  EV_for_h=Error_info(X,y,t,n,p,mv,mvinfo,h_beta,h_sig,h_rho)
  
  h=bandwidth_est(X,y,t,n,p,mvinfo,mv,EV_for_h,h_beta_2)
  #h=0.5
  betahatti=t(apply(matrix(seq(0.01,1,by=0.01),ncol=1),1,function(x) betahat(x,X,y,t,h,mv)))
  betahatt=t(apply(matrix(t,ncol=1),1,function(x) betahat(x,X,y,t,h,mv))) 
  
  l_plus_betatildei=t(apply(cbind(seq(0.01,1,0.01),betahatti),1,function(x) l_beta_tilde(x,betahatt,X,y,t,n,p,h,mvinfo,mv,blanced=1)))
  integrand=2*(l_plus_betatildei[,2])
  
  T[rep]=0.01*sum(integrand)
  T1[rep]=0.01*sum(integrand[11:90])*5/4  
  
  h_rec[rep]=h
  logelr[rep,]=integrand
  
  EV_for_boots=Error_info(X,y,t,n,p,mv,mvinfo,h,h_sig,h_rho)
 
  epi_B=Bootstrapped_vec(B,n,mvinfo,EV_for_boots)
  
  l_plus_betatildet=t(apply(cbind(t,betahatt),1,function(x) l_beta_tilde(x,betahatt,X,y,t,n,p,h,mvinfo,mv,blanced=1)))
  betatildet=l_plus_betatildet[,1]
  
  for(bb in 1:B)
  {
    epi_b=epi_B[,bb]
    y_b=as.matrix(X1*betatildet+X2*betatildet+epi_b)
    betahatti_b=t(apply(as.matrix(seq(0.01,1,by=0.01)),1,function(x) betahat(x,X,y_b,t,h,mv)))
    betahatt_b=t(apply(as.matrix(t),1,function(x) betahat(x,X,y_b,t,h,mv))) 
    l_plus_betatildei_b=t(apply(cbind(seq(0.01,1,0.01),betahatti_b),1,function(x) l_beta_tilde(x,betahatt_b,X,y_b,t,n,p,h,mvinfo,mv,blanced=1)))
    integrand_b=2*(l_plus_betatildei_b[,2])
    TB[rep,bb]=0.01*sum(integrand_b)
    T1B[rep,bb]=0.01*sum(integrand_b[11:90])*5/4  
  }
}

write.table(cbind(T,T1,h_rec,logelr), file = "EL_Fun_Result.csv", sep = ",", col.names = NA)
