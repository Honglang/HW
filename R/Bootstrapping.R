#######################Wild Bootstrapping (6/22/2014)###############################
#  Lots of functions are called from the following source file
#  There are several pre bandwidths to be chosen, h_beta, h_rho, h_sig
####################################################################################
source("//sttnas/StudentReds/wangho16/My Documents/Research Pool/EL and Functional Data/ELFCODE/R Version/Updated Version/Bandwidth_Est.R")


Bootstrapped_vec=function(B,n,mvinfo,EV)
{
  N=sum(mvinfo)
  BV=matrix(0,N,B)
  count=0
  count1=0
  for(k in 1:n)
  {
    Vhat=matrix(EV[(count1+1):(count1+mvinfo[k]^2)],nrow=mvinfo[k])    
    BV[(count+1):(count+mvinfo[k]),]=gauss_vec_0_mean_given_cov(B,mvinfo[k],Vhat)
    count=count+mvinfo[k]
    count1=count1+mvinfo[k]^2
  }
  
  return(BV)
}