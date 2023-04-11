rm(list=ls())

#install.packages("genridge")
# install.packages("glmnet")
library(genridge)
library(glmnet)

d_prostate = prostate
x = d_prostate[,1:8] # matrix of predictors
y = d_prostate[,9] # set of response
x = as.matrix(x)

LASSO_model = glmnet(x,y,alpha=1)

plot(LASSO_model,xvar="lambda",label=T,lwd=2)
#legend("topright",legend=c("lcavol (1)","lweight (2)","age (3)","lbph (4)","svi (5)","lcp (6)","gleason (7)",
 #                          "pgg45 (8)"),pch=19,col=c("black","red","green","skyblue",
 #                                                           "aquamarine","violet","black","red"))
                                                            
# cross-validation to find the optimized complexity parameter, lambda
cross_vad = cv.glmnet(x,y,nfolds=20)
cross_vad$lambda.min

lambda = LASSO_model$lambda
eff_lambda = array(0) # effective degrees of freedom for the ridge fits
sv = svd(x) # invoking SVD of the design matrix X, for singular values

for(i in 1:length(lambda))
{
  eff_lambda[i] = sum(sv$d^2/(sv$d^2+lambda[i]))
}

df_opt = sum(sv$d^2/(sv$d^2+cross_vad$lambda.min)) # optimum value for eff. deg. of freedom

beta_coeff = as.matrix(LASSO_model$beta) # ridge coefficients

plot(eff_lambda,beta_coeff[1,],type="l",col="black",ylim=c(-0.2,0.7),lwd=2,
     xlab=bquote("df("~lambda~")"),ylab="Ridge Coefficients")
lines(eff_lambda,beta_coeff[2,],type="l",col="red",lwd=2)
lines(eff_lambda,beta_coeff[3,],type="l",col="green",lwd=2)
lines(eff_lambda,beta_coeff[4,],type="l",col="blue",lwd=2)
lines(eff_lambda,beta_coeff[5,],type="l",col="brown",lwd=2)
lines(eff_lambda,beta_coeff[6,],type="l",col="violet",lwd=2)
lines(eff_lambda,beta_coeff[7,],type="l",col="aquamarine",lwd=2)
lines(eff_lambda,beta_coeff[8,],type="l",col="orange",lwd=2)
abline(h=0,lty=2)
abline(v=df_opt,lty=2)
