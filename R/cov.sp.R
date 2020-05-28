cov.sp <- function(Lt,Ly,newt=NULL,
                   domain=NULL,
                   weig=NULL,
                   corf=NULL, # correlation function(theta,x,y)
                   mu=NULL,
                   sig2e=NULL,
                   sig2x=NULL,
                   pfunc=NULL,
                   theta0=NULL,
                   lb=NULL,
                   ub=NULL,
                   D=NULL)
{
    
    if(is.null(corf)){
        corf <- function(x,y,theta) matern(x,y,nu=theta)
        D <- 1
    }
    else
    {
        if(is.null(theta0) && is.null(D)) stop('The dimension D must be specified')
    }
    
    if(is.null(mu))
    {
        mu <- meanfunc(Lt,Ly,method='PACE')
    }
    
    if(is.function(mu)) mu.hat <- lapply(Lt,mu)
    else mu.hat <- predict(mu,Lt)
    
    if(is.null(sig2e)) sig2e <- sigma2(Lt,Ly)
    
    if(is.null(sig2x))
    {
        sig2x <- varfunc(Lt,Ly,mu=mu,sig2=sig2e)
    }
    
    if(is.function(sig2x))
    {
        var.hat <- lapply(Lt, sig2x)
    }
    else
    {
        var.hat <- predict(sig2x,Lt)
    }
    
    if(is.null(domain))
    {
        t.vec <- unlist(Lt)
        domain <- c(min(t.vec),max(t.vec))
    }
    
    th.est <- estimate.theta(Lt,Ly,
                              D=D,
                              var.hat=var.hat,
                              mu.hat=mu.hat,
                              method='LS',
                              rho=corf,
                              weig=weig,
                              pfunc=pfunc,
                              theta.lb=lb,
                              theta.ub=ub,
                              theta0=theta0,
                              domain=domain)$LS
    
    rslt <- list(sig2e=sig2e,
                 theta=th.est,
                 mu.hat=mu.hat,
                 domain=domain,
                 mu=mu,
                 sig2x=sig2x,
                 rho=function(x,y) corf(x,y,th.est),
                 method='SP')
    class(rslt) <- 'covfunc'
    
    if(!is.null(newt))
        rslt$fitted <- predict(rslt,newt)
    
    return(rslt)
}

# matern correlation
matern <- function(x,y=NULL,nu=1,rho=1)
{
    if(is.null(y)) y <- x
    
    G <- expand.grid(x,x)
    
    S <- apply(G,1,function(z){
        delta <- abs(z[1]-z[2])/rho
        if(delta == 0)
            1
        else
            (sqrt(2*nu)*delta)^nu * besselK(sqrt(2*nu)*delta,nu=nu) / (2^(nu-1)*gamma(nu))
    })
    
    C <- matrix(S,nrow=length(x),ncol=length(y),byrow=F)
    return(C)
}


#' Estimate the parameters of the correlation structure
#' @param weig a vector of length(Lt)
#' @param D the number of parameters
#' @param sig0.hat the estimate of variance of measure error
#' @keywords internal
estimate.theta <- function(Lt,Ly,
                           D,
                           var.hat,
                           mu.hat,
                           method,
                           rho,
                           weig=NULL,
                           theta.lb=NULL,
                           theta.ub=NULL,
                           h=0.1,
                           domain=c(0,1),
                           pfunc=NULL,
                           theta0=NULL)
{
    if(is.null(weig))   weig <- rep(1,length(Lt))
    if(is.null(pfunc)) pfunc <- function(th){return(0)}
    
    sig <- lapply(var.hat,function(x){
        x <- ifelse(x <=0,yes=1e-8,no=x)
        sqrt(x)
    })
    
    n <- length(Lt)
    result <- list()
    
    if('LS' %in% method)
    {
        R <- sapply(1:n,function(i){
            tobs <- Lt[[i]]
            idx <- (tobs >= domain[1]+h) & (tobs <= domain[2]-h)
            yy <- Ly[[i]]
            yy <- yy[idx]
            mu <- mu.hat[[i]]
            mu <- mu[idx]
            resid <- yy-mu
            return(list(cbind(resid) %*% rbind(resid)))
        })
        
        TT <- sapply(1:n,function(i){
            tobs <- Lt[[i]]
            idx <- (tobs >= domain[1]+h) & (tobs <= domain[2]-h)
            tobs[idx]
        })
        
        SS <- sapply(1:n,function(i){
            tobs <- Lt[[i]]
            idx <- (tobs >= domain[1]+h) & (tobs <= domain[2]-h)
            sig.t <- sig[[i]]
            sig.t <- sig.t[idx]
            diag(sig.t)
        })
        
        Q <- function(theta){
            if(any(theta < theta.lb) || any(theta > theta.ub)) return(1e300)
            v <- sapply(1:n,function(i){
                tobs <- TT[[i]]
                if(length(tobs)==1) return(0)
                
                rho.mat <- rho(tobs,tobs,theta)
                if(any(is.infinite(rho.mat)) || any(is.nan(rho.mat)))
                {
                    theta
                    tobs
                }
                S <- SS[[i]]
                # print(tobs)
                # print(dim(S))
                # print(dim(rho.mat))
                resid <- (S %*% rho.mat %*% S - R[[i]])^2
                diag(resid) <- 0 # remove points on diagonal
                return(sum(resid))
            })
            
            if(!is.finite(sum(v*weig)))    return(1e300)
            return(sum(v*weig) + pfunc(theta)) # regularization
        }
        
        # set up initial value for optimization
        if(is.null(theta.lb)) theta.lb <- rep(-Inf,D)
        if(is.null(theta.ub)) theta.ub <- rep(Inf,D)
        
        if(is.null(theta0))
            v <- sapply(1:D,function(d){
                lb <- theta.lb[d]
                ub <- theta.ub[d]
                if(lb > -Inf && ub < Inf) return(runif(1)*(ub-lb)+lb)
                else if(lb > -Inf) return(runif(1)+lb)
                else if(ub < Inf) return(ub-runif(1))
                else return(runif(1))
            })
        else v <- theta0
        res <- optim(v,Q,lower=theta.lb,upper=theta.ub,method='L-BFGS-B')
        result$LS <- res$par
        result$Q <- Q(res$par)
    }
    
    if('QMLE' %in% method)
    {
        R <- sapply(1:n,function(i){
            tobs <- Lt[[i]]
            idx <- (tobs >= domain[1]+h) & (tobs <= domain[2]-h)
            yy <- Ly[[i]]
            yy <- yy[idx]
            mu <- mu.hat[[i]]
            mu <- mu[idx]
            resid <- yy-mu
            return(list(cbind(resid) %*% rbind(resid)))
        })
        
        TT <- sapply(1:n,function(i){
            tobs <- Lt[[i]]
            idx <- (tobs >= domain[1]+h) & (tobs <= domain[2]-h)
            tobs[idx]
        })
        
        SS <- sapply(1:n,function(i){
            tobs <- Lt[[i]]
            idx <- (tobs >= domain[1]+h) & (tobs <= domain[2]-h)
            sig.t <- sig[[i]]
            sig.t <- sig.t[idx]
            diag(1/sig.t)
        })
        
        
        # set up initial value for optimization
        if(is.null(theta.lb)) theta.lb <- rep(-Inf,D)
        if(is.null(theta.ub)) theta.ub <- rep(Inf,D)
        
        Q <- function(theta){
            if(any(theta < theta.lb) || any(theta > theta.ub)) return(1e300)
            v <- sapply(1:n,function(i){
                tobs <- TT[[i]]
                m <- length(tobs)
                rho.mat <- rho(theta,tobs,tobs)
                if(any(is.infinite(rho.mat)) || any(is.nan(rho.mat)))
                {
                    theta
                    tobs
                }
                if(cond(rho.mat) > 1e2)
                {
                    rho.inv <- inv(rho.mat+0.01*diag(rep(1,m)))
                }
                else rho.inv <- inv(rho.mat)
                S <- SS[[i]]
                resid <- (t(R[[i]]) %*% S %*% rho.inv %*% S %*% R[[i]])
                resid <- resid + log(abs(det(rho.mat)))
                return(resid/2)
            })
            rs <- sum(v*weig)
            if(is.nan(rs) || is.infinite(rs))
            {
                rs <- 1e100
            }
            return(rs)
        }
        
        
        v <- sapply(1:D,function(d){
            lb <- theta.lb[d]
            ub <- theta.ub[d]
            if(lb > -Inf && ub < Inf) return(runif(1)*(ub-lb)+lb)
            else if(lb > -Inf) return(runif(1)+lb)
            else if(ub < Inf) return(ub-runif(1))
            else return(runif(1))
        })
        res <- optim(v,Q,lower=theta.lb,upper=theta.ub,method='L-BFGS-B')
        result$QMLE <- res$par
    }
    
    return(result)
}
