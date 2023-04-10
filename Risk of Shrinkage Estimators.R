rm(list=ls())
library(MASS)

set.seed(4)

n=100
p = 10
mu = runif(p,0,2)
sigma = diag(1,ncol=p,nrow=p)
x = mvrnorm(n,mu,sigma)

f = function(alpha,x,p,n)
{
  s=0
  for(i in 1:n)
  {
    xx = as.matrix(x[i,])
    s = s+(1/(t(xx)%*%xx))
  }
  s = s/n
  res = p-(alpha*(2*(p-2) - alpha))*(s)
  return(res)
}

pl = array(0)
alpha = seq(0,16,0.01)
for(i in 1:length(alpha))
{
  pl[i] = f(alpha[i],x,p,n)
}

plot(alpha,pl,type="l",col="red",xlab=bquote(alpha),ylab=bquote(R[hat(mu)[alpha]](mu)))
abline(h=p,col="black")
abline(v=p-2,col="blue",lty=2)
abline(v=0,col="brown",lty=4)
abline(v=2*(p-2),col="brown",lty=4)