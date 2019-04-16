golden=function(f,a,b,tol)
{
  z=seq(a,b,length=20)
  v=apply(as.matrix(z),1,f)
  arg=which(v==min(v))
  if(arg==1) 
  {
    a=z[1]
    b=z[arg+1]
  }
    
  if(arg==20) 
  {
    a=z[arg-1]
    b=z[20]
  }  
  if(arg>1 && arg <20)
  {
    a=z[arg-1]
    b=z[arg+1]
  }
  
  tau=(sqrt(5)-1)/2
  x1=b-tau*(b-a)
  x2=a+tau*(b-a)
  f1=f(x1)
  f2=f(x2)
  a0=a
  b0=b  
  while(b0-a0>tol)
{
  if(f1<f2)
  {
    b0=x2
    x2=x1
    x1=b0-tau*(b0-a0)
    f2=f1
    f1=f(x1)  
  }
  else
  {
    a0=x1
    x1=x2
    x2=a0+tau*(b0-a0)
    f1=f2
    f2=f(x2)
  }
}
  return(c((a0+b0)/2,f((a0+b0)/2)))
  
}

#f1=function(x)
#{
#  (x-1)^2+2
#}

#golden(f1,-1,2,10^(-7))