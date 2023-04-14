rm(list=ls())


# James-Stein Estimation using the prostate data --------------------------


# Singh et al. (2002) and B. Efron ----------------------------------------
 
# install.packages("sda")
library(sda) # library to access the prostate data

prostate = singh2002$x # rows: samples; columns: genes

# The microarray 'prostate' data comprises of 6033 (columns) gene expressions for
# a total of 102 (rows) samples, where the first 50 samples (rows) are healthy (controls) 
# and the remaining 52 (rows) samples are prostate cancer patients (cases).

# The main idea is to apply James-Stein Estimation to obtain an estimator of the 
# location parameter


# Few Explorations --------------------------------------------------------

# 6033 z-values from the prostate cancer microarray data on 50 non-cancer and 52
# prostate cancer patients

prostz = read.csv("/Users/somjit/Desktop/prostz.csv")

h = hist(prostz$x,breaks=60,col="#19D5FF",ylab="",xlab="z-values",main="")
x_values = seq(min(prostz$x), max(prostz$x), length = 100)
y_values = dnorm(x_values, mean = mean(prostz$x), sd = sd(prostz$x)) 
y_values = y_values * diff(h$mids[1:2]) * length(prostz$x) 
lines(x_values, y_values, lwd = 1)
abline(v=0,col="red")


# JS Estimation with identity Var-Cov Matrix ------------------------------

# Construction (formulation) of the data that is being observed: We consider the
# gene expression observations for each of the 6033 gene expressions only for the
# first control and case respectively.

X = array(dim=0) # p-dimensional vector; p = 6033*2 = 12066
p = 6033*2
p1 = 6033
k = 1 # controls the dimension of the data
X_controls = prostate[1,] # first control gene expressions
X_cases = prostate[51,] # first case gene expressions

for(i in 1:p1)
{
  X[k] = X_controls[i]
  X[k+1] = X_cases[i]
  k=k+2
}

# grand mean (common mean)

grand.mean = mean(X)

# James-Stein Estimator of Location Parameter

mu.js = grand.mean + (1-((p-3)/sum((X-grand.mean)*(X-grand.mean))))*(X-grand.mean)

# plotting the James-Stein Estimates

x_axis = 1:50
s=1:50
X_plot = X[1:50]
mu.js_plot = mu.js[1:50]

plot(1:50,X_plot,type="p",pch=19,cex=0.5,xlab="Plotting the first 50 location estimates",
     ylab="Estimates/Observations",ylim=c(-2,6))
points(1:50,mu.js_plot,col="yellow",pch=19,cex=.5)
arrows(x_axis[s],X_plot[s],x_axis[s],mu.js_plot[s],code=0,col="brown")
abline(h=grand.mean,col="red",lty=4)
legend("topright",legend=c("obs. gene expressions","shrunken JSE's","grand mean"),
       pch=19,col=c("black","yellow","red"))


