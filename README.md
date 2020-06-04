# mcfda
Mean and covariance estimation for functional data analysis

## Install
devtools::install_github("linulysses/mcfda")


## Demonstration
library(mcfda)

set.seed(1)

### set up mean and cov functions and synthesize data functional dataset
mu <- function(s) sin(2*pi*s)

cov <- synfd::matern

D <- synfd::irreg.fd(mu=mu, X=synfd::gaussian.process(cov), n=100, m=5)

### estimate mean by 'FOURIER' (local polynomial) method
mu.obj <- meanfunc(D$t,D$y,newt=NULL,method='FOURIER',
                   tuning='cv',weig=NULL, domain=c(0,1))

### plot the object
plot(mu.obj)

lines(regular.grid(),mu(regular.grid()))

### estimate mean by 'PACE' (local polynomial) method
mu.obj <- meanfunc(D$t,D$y,newt=NULL,method='PACE',
                   tuning='cv',weig=NULL,kernel='epanechnikov',deg=1)

plot(mu.obj)

lines(regular.grid(),mu(regular.grid()))

### estimate covariance by 'PACE' 
cov.obj <- covfunc(D$t,D$y,newt=NULL,mu=mu.obj,method='PACE',
                   tuning='GCV',weig=NULL,kernel='epanechnikov',delta=Inf)


grid <- regular.grid()

### evaluate the estimated cov at grid
cov.hat <- predict(cov.obj,grid)

### the true covariance at grid
cov0 <- cov(grid) 

### relative error in MISE
mean((cov.hat-cov0)^2) / mean(cov0^2) 


### estimate covariance by 'FOURIER' 
cov.obj <- covfunc(D$t,D$y,newt=NULL,mu=mu.obj,method='FOURIER',
                   tuning='lle',weig=NULL)

### evaluate the estimated cov at grid
cov.hat <- predict(cov.obj,grid)

### relative error in MISE
mean((cov.hat-cov0)^2) / mean(cov0^2) 


### estimate covariance by 'SP' 
cov.obj <- covfunc(D$t,D$y,newt=NULL,method='SP',weig=NULL)

### evaluate the estimated cov at grid
cov.hat <- predict(cov.obj,grid)

### relative error in MISE
mean((cov.hat-cov0)^2) / mean(cov0^2) 


### estimate variance of the measurement noise
sig2 <- sigma2(D$t,D$y)

sig2

### estimation error for the measurement variance
abs(sig2 - attr(D,'sig')^2)


## FPCA via [fdapace](https://cran.r-project.org/web/packages/fdapace/index.html) package
Once the mean and covariance functions are obtained, the fdapace package can used with the options userCov, userMu and userSigma2 to perform FPCA:

error <- ifelse(sig2==0,yes=FALSE,no=TRUE)

grid <- regular.grid(h=0)

optns <- list(userCov=list(t=grid,cov=predict(cov.obj,grid)),
              userMu=list(t=grid,mu=predict(mu.obj,grid)), 
              userSigma2 = sig2,
              error=error,
              verbose=T)

R <- fdapace::FPCA(D$y,D$t,optns)

### plot the first FPC, etc
plot(R$phi[,1])

