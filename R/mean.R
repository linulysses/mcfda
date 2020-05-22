#' Estimate the mean function from functional data/snippets
#' @param t a list of vectors (for irregular design) or a vector (for regular design) containing time points of observations for each individual. Each vector should be in ascending order
#' @param y a list of vectors (for irregular design) or a matrix (for regular design) containing the observed values at \code{t}. If it is a matrix, the columns correspond to the time points in the vector \code{t}
#' @param newt  a list of vectors or a vector containing time points of observations to be evaluated. If NULL, then newt is treated as t
#' @param method estimation method, 'PACE' or 'FOURIER'
#' @param tuning tuning method to select possible tuning parameters
#' @param ... other parameters required depending on the \code{method} and \code{tuning}; see details
#' @details
#'     \itemize{
#'         \item{When \code{method='PACE'}, additional parameters \code{kernel} and \code{deg} can be provided. \code{bw} as a scalar is optional. When \code{bw} is provided, the bandwidth is set to \code{bw}}
#'         \item{When \code{method='FOURIER'}, additional parameters \code{q},\code{rho},\code{ext} and \code{domain} are optional. If they are not provided, then they will be deduced from data or selected by the specified \code{tuning} method.}
#'     }
#'
#' @return an object of the class 'meanfunc' containing necessary information to predict the mean function
#' 
#' @importFrom Rdpack reprompt
#' @references 
#' \insertRef{Lin2020}{mcfda}
#' 
#' \insertRef{Yao2005}{mcfda}
#' 
#' @examples
#' mu <- function(s) sin(2*pi*s)
#' D <- synfd::sparse.fd(mu=mu, X=synfd::gaussian.process(), n=100, m=5)
#' mu.obj <- meanfunc(D$t,D$y,newt=NULL,method='PACE',
#'                 tuning='cv',weig=NULL,kernel='gauss',deg=1)
#' # equivalent to
#' # mu.obj <- meanfunc(D$t,D$y)
#'
#' # plot the object
#' plot(mu.obj)
#' @export meanfunc
meanfunc <- function(t,y,newt=NULL,method=c('PACE','FOURIER'),
                      tuning='cv',weig=NULL,...)
{

    method <- match.arg(method)
    R <- NULL
    if(is.list(t)) # irregular design
    {
        if(!is.list(y)) stop('y must be list when t is list')
        if(length(y) != length(t))
            stop('the length of y must match the length of t')
        if(method=='PACE')
            R <- mean.pace(t,y,tuning,weig,...)
        else
            R <- mean.basis(t,y,tuning,weig,...)
    }
    else if(is.vector(t)) # regular design
    {
        if(!is.matrix(y)) stop('y must be a matrix when t is vector')
        if(ncol(y) != length(t))
            stop('the number of columns must match the length of t')
        if(method=='PACE')
            R <- mean.pace(t,y,tuning,weig,...)
        else
            R <- mean.basis(t,y,tuning,weig,...)
    }
    else stop('t must be a list of vectors (for irregular design) or a vector (for regular design')

    class(R) <- 'meanfunc'
    return(R)
}


mean.pace <- function(t,y,tuning,weig,...)
{
    others <- list(...)

    if(is.null(weig)) weig <- get.weig.mean(t,y,'OBS')
    else if(is.character(weig)) weig <- get.weig.mean(t,y,weig)
    else stop('weig is not recognized.')

    # estimate the domain, used when method == 'FOURIER' and also plot.meanfunc
    domain <- get.optional.param('domain',others,NULL)
    if(is.null(domain))
    {
        domain <- c(0,0)
        domain[1] <- min(unlist(t))
        domain[2] <- max(unlist(t))
        domain <- c(domain[1]-0.01*(domain[2]-domain[1]),
                    domain[2]+0.01*(domain[2]-domain[1]))
    }
    
    if(is.list(t))
    {
        n <- length(t)
        x <- unlist(t)
        y <- unlist(y)
        
        # expand weig into the same format of t
        mi <- sapply(t,length)
        weig <- lapply(1:length(t),function(i) rep(weig[i],mi[i]))
        
        weig <- unlist(weig)
        
        datatype <- 'irregular'
    }
    else
    {
        n <- length(t)
        weig <- rep(1/n,n)
        datatype <- 'regular'
        
        x <- t
        y <- apply(y,2,mean)
    }


    kernel <- toupper(get.optional.param('kernel',others,'EPANECHNIKOV'))
    kernel <- match.arg(kernel,c('GAUSSIAN','EPANECHNIKOV'))
    if(kernel == 'GAUSSIAN') kernel <- locpol::gaussK
    else kernel <- locpol::EpaK

    deg <- get.optional.param('deg',others,1)
    bw <- get.optional.param('bw',others,NULL)

    if(is.null(bw))
        bw <- compute.bw.1D(x,y,tuning,weig,kernel,deg)

    R <- list(bw=bw,x=x,y=y,n=n,method='PACE',domain=domain,
              weig=weig,kernel=kernel,deg=deg,yend=c(NULL,NULL))
    class(R) <- 'meanfunc'

    L0 <- domain[2]-domain[1]
    yend <- predict(R,c(domain[1]+L0/100,domain[2]-L0/100))
    R$yend <- yend


    return(R)
}

mean.basis <- function(t,y,tuning,weig,...)
{

    others <- list(...)

    if(is.list(t))
    {
        if(is.null(weig)) weig <- get.weig.mean(t,y,'OBS')
        else if(is.character(weig)) weig <- get.weig.mean(t,y,weig)
        else stop('weig is not recognized.')

        datatype <- 'irregular'
    }
    else
    {
        n <- nrow(y)
        weig <- rep(1/n,n)
        datatype <- 'regular'
    }


    # estimate the domain, used when method == 'FOURIER' and also plot.meanfunc
    domain <- get.optional.param('domain',others,NULL)
    if(is.null(domain))
    {
        domain <- c(0,0)
        domain[1] <- min(unlist(t))
        domain[2] <- max(unlist(t))
        domain <- c(domain[1]-0.01*(domain[2]-domain[1]),
                    domain[2]+0.01*(domain[2]-domain[1]))
    }

    # helper functions

    # compute the auxillary matrices
    compute.aux.matrices <- function(q,t,y,weig,domain)
    {
        n <- length(t)
        B <- lapply(t,function(s){
            evaluate.basis(K=q,domain=domain,grid=s,type='FOURIER')
        })
        H <- Reduce('+',
                    lapply(1:n,
                           function(i) weig[i]*(base::t(B[[i]]) %*% B[[i]])))
        G <- do.call(cbind,
                     lapply(1:n,
                            function(i) weig[i]*base::t(B[[i]])))
        Y <- do.call(rbind,
                     lapply(y,
                            function(yi) matrix(yi,length(yi),1)))
        M <- 1000
        pts <- regular.grid(m=1000,domain=domain)
        W <- deriv.fourier(K=q,grid=pts,r=2,domain=domain)
        W <- base::t(W) %*% W / M

        return(list(W=W,B=B,G=G,Y=Y,H=H))
    }

    # compute the auxillary matrices
    compute.aux.matrices.reg <- function(q,t,y,weig,domain)
    {
        n <- nrow(y)
        B <- evaluate.basis(K=q,domain=domain,grid=t,type='FOURIER')
        H <- base::t(B) %*% B
        G <- base::t(B)
        Y <- apply(y,2,mean)

        M <- 1000
        pts <- regular.grid(m=1000,domain=domain)
        W <- deriv.fourier(K=q,grid=pts,r=2,domain=domain)
        W <- base::t(W) %*% W / M

        return(list(W=W,B=B,G=G,Y=Y,H=H))
    }

    pred <- function(tobs,bhat,ext,domain)
    {

        if(ext > 0) D <- c(domain[1]-(domain[2]-domain[1])*ext,
                           domain[2]+(domain[2]-domain[1])*ext)
        else D <- domain
        B <- evaluate.basis(K=length(bhat),grid=tobs,domain=D,type='FOURIER')
        return(c(B %*% bhat))
    }

    extend <- function(domain,ext)
    {
        return(c(domain[1]-(domain[2]-domain[1])*ext,
                 domain[2]+(domain[2]-domain[1])*ext))
    }

    sse <- function(bhat,ext,domain,t,y,datatype)
    {
        if(datatype=='irregular')
        {
            E <- sapply(1:length(t),function(i){
                tmp <- pred(t[[i]],bhat,ext,domain)
                sum((tmp-y[[i]])^2)
            })
            return(sum(E))
        }
        else
        {
            tmp <- pred(t,bhat,ext,domain)
            y <- apply(y,2,mean)
            return(sum((tmp-y)^2))
        }

    }

    tune <- function(t,y,tuning,weig,q,rho,ext,domain,datatype)
    {
        if(tolower(tuning) == 'cv')
        {
            Kfold <- ifelse(length(t)<20,yes=length(t),no=5)
            if(datatype=='irregular')
                cvo <- cv.partition(length(t),Kfold)
            else
                cvo <- cv.partition(nrow(y),Kfold)
            # compute cv error for each fold
            err <- lapply(1:cvo$num.test.sets,
                          function(k){
                              trIdx <- cvo$training[[k]]
                              teIdx <- cvo$test[[k]]
                              # for each candidate q value (in columns)
                              sapply(q,function(qi){
                                  if(datatype=='irregular')
                                    aux.mat <- compute.aux.matrices(qi,t[trIdx],y[trIdx],weig[trIdx],domain)
                                  else
                                    aux.mat <- compute.aux.matrices.reg(qi,t,y[trIdx,],weig[trIdx],domain)
                                  # for each candidate rho value (in rows)
                                  sapply(rho,function(ri){
                                      bhat <- solve(ri*aux.mat$W+aux.mat$H,aux.mat$G %*% aux.mat$Y)
                                      if(datatype=='irregular')
                                          sse(bhat,0,domain,t[teIdx],y[teIdx],datatype)
                                      else
                                          sse(bhat,0,domain,t,y[teIdx,],datatype)
                                  })

                              },simplify=TRUE)
                          })
            err <- as.matrix(Reduce('+',err))
            I <- which(err==min(err),arr.ind=T)

            tmp <- list(q=q[I[1,2]],rho=rho[I[nrow(I),1]])

            if(length(ext) > 1)
            {
                err <- lapply(1:cvo$num.test.sets,
                              function(k){
                                  trIdx <- cvo$training[[k]]
                                  teIdx <- cvo$test[[k]]
                                  # for each candidate ext value (in columns)
                                  sapply(ext,function(ei){
                                      if(datatype == 'irregular')
                                        aux.mat <- compute.aux.matrices(tmp$q,t[trIdx],y[trIdx],weig[trIdx],extend(domain,ei))
                                      else
                                        aux.mat <- compute.aux.matrices.reg(tmp$q,t,y[trIdx,],weig[trIdx],extend(domain,ei))
                                      bhat <- solve(tmp$rho*aux.mat$W+aux.mat$H,aux.mat$G %*% aux.mat$Y)
                                      if(datatype == 'irregular')
                                          sse(bhat,ei,domain,t[teIdx],y[teIdx],datatype)
                                      else
                                          sse(bhat,ei,domain,t,y[teIdx,],datatype)
                                  })
                              })
                err <- Reduce('+',err)
                I <- which(err==min(err))
                tmp$ext <- ext[I]
            }

            return(tmp)
        }
        else stop(paste0('The tuning method ', tuning, ' is not supported yet'))

    }

    q <- get.optional.param('q',others,seq(3,15,by=2))
    ext <- get.optional.param('ext',others,seq(0,0.2,by=0.05))
    rho <- get.optional.param('rho',others,NULL)

    if(is.null(rho)) # get a rough idea of rho
    {
        tmp <- tune(t,y,tuning,weig,9,10^seq(-8,2,length.out=21),0.1,domain,datatype)
        rho <- seq(tmp$rho/50,tmp$rho*50,length.out=21)
    }

    if(length(q) > 1 || length(rho) > 1 || length(ext) > 1)
    {
        tmp <- tune(t,y,tuning,weig,q,rho,ext,domain,datatype)
        q <- tmp$q
        rho <- tmp$rho
        ext <- tmp$ext
    }


    if(datatype=='irregular')
        aux.mat <- compute.aux.matrices(q,t,y,weig,extend(domain,ext))
    else
        aux.mat <- compute.aux.matrices.reg(q,t,y,weig,extend(domain,ext))

    bhat <- solve(rho*aux.mat$W+aux.mat$H,aux.mat$G %*% aux.mat$Y)

    R <- list(q=q,rho=rho,ext=ext,weig=weig,bhat=bhat,domain=domain,method='FOURIER')
    class(R) <- 'meanfunc'
    return(R)
}


#' predict mean functions at new locations
#' @param meanfunc.obj the object obtained by calling \code{mean.func}
#' @param newt a vector or a list of vectors of real numbers
#' @return the estimated mean function evaluated at \code{newt}. It has the same format of \code{newt}
#' @export
predict.meanfunc <- function(meanfunc.obj,newt)
{
    pred <- function(newt) # newt must be a vector
    {
        idxl <- newt < meanfunc.obj$domain[1]
        idxu <- newt > meanfunc.obj$domain[2]
        idx <- (!idxl) & (!idxu)

        tmp <- loclin1D(x=meanfunc.obj$x,
                               y=meanfunc.obj$y,
                               newx=newt[idx],
                               bw=meanfunc.obj$bw,
                               weig=meanfunc.obj$weig,
                               kernel=meanfunc.obj$kernel,
                               deg=meanfunc.obj$deg)$fitted

        yhat <- rep(0,length(newt))
        yhat[idx] <- tmp
        yhat[idxl] <- meanfunc.obj$yend[1]
        yhat[idxu] <- meanfunc.obj$yend[2]
        return(yhat)
    }

    if(is.list(newt))
    {
        if(toupper(meanfunc.obj$method) == 'PACE')
        {
            mi <- lapply(newt,length)
            newt <- unlist(newt)
            fitted <- pred(newt)

            cm <- c(0,cumsum(mi))

            R <- sapply(1:length(mi),function(i){
                res <- list()
                res[[1]] <- fitted[(cm[i]+1):cm[i+1]]
                res
            })
            return(R)
        }
        else if(toupper(meanfunc.obj$method) == 'FOURIER')
        {
            domain <- meanfunc.obj$domain
            ext <- meanfunc.obj$ext
            if(ext > 0)
                D <- c(domain[1]-(domain[2]-domain[1])*ext,
                       domain[2]+(domain[2]-domain[1])*ext)
            else D <- domain
            ret <- lapply(newt,function(tobs){
                B <- evaluate.basis(K=length(meanfunc.obj$bhat),
                                    grid=tobs,domain=D,type='FOURIER')
                return(c(B %*% meanfunc.obj$bhat))
            })
            return(ret)
        }
        else stop(paste0('Method ', meanfunc.obj$method, ' is not recognized.'))

    }
    else if(is.vector(newt))
    {
        if(toupper(meanfunc.obj$method) == 'PACE')
        {
            return(pred(newt))
        }
        else if(toupper(meanfunc.obj$method) == 'FOURIER')
        {
            domain <- meanfunc.obj$domain
            ext <- meanfunc.obj$ext
            if(ext > 0)
                D <- c(domain[1]-(domain[2]-domain[1])*ext,
                       domain[2]+(domain[2]-domain[1])*ext)
            else D <- domain
            B <- evaluate.basis(K=length(meanfunc.obj$bhat),
                                    grid=newt,domain=D,type='FOURIER')
            return(c(B %*% meanfunc.obj$bhat))

        }
        else stop(paste0('Method ', meanfunc.obj$method, ' is not recognized.'))
    }
    else stop('newt must be a vector or a list of vectors of real numbers')
}


#' plot the estimated mean function
#' @param meanfunc.obj the object obtained by calling \code{meanfunc}
#' @param ... other parameters passed to \code{plot}
#' @export
plot.meanfunc <- function(meanfunc.obj,...)
{
    plot(regular.grid(50,domain=meanfunc.obj$domain),
         predict(meanfunc.obj,regular.grid(50,domain=meanfunc.obj$domain)),
         ...)
}
