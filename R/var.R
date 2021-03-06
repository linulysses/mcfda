#' Estimate the variance function
#' @param Lt a list (for irregular design) or a vector (for regular design) of observations in time domain
#' @param Ly a list (for irregular design) or a matrix (for regular design) of observations corresponding to \code{Lt}
#' @param sig2 variance of noise
#' @param newt a set of locations where the variance function is estimated. If NULL, then set to \code{Lt}
#' @param method estimation method, 'PACE' or 'FOURIER'
#' @param weig a vector of length(Lt). If NULL, then determined by the algorithm.
#' @param mu an object generated by \code{meanfunc}.
#' @param ... other parameters required depending on the \code{method} and \code{tuning}; see details
#' @details
#'     \itemize{
#'         \item{When \code{method='PACE'}, additional parameters \code{kernel} and \code{deg} can be provided. \code{bw} as a scalar is optional. When \code{bw} is provided, the bandwidth is set to \code{bw}}
#'         \item{When \code{method='FOURIER'}, additional parameters \code{q},\code{rho},\code{ext} and \code{domain} are optional. If they are not provided, then they will be deduced from data or selected by the specified \code{tuning} method.}
#'     }
#' @export
varfunc <- function(Lt,Ly,newt=NULL,sig2=NULL,method=c('PACE','FOURIER'),mu=NULL,weig=NULL,...)
{
    method <- match.arg(method)
    if(is.list(Lt) && is.list(Ly))
    {
        n <- length(Lt)
        datatype <- 'irregular'
    }
    else if(is.vector(Lt) && is.matrix(Ly))
    {
        n <- nrow(Ly)
        datatype <- 'regular'
    }
    else stop('unsupported data type')

    if(is.null(sig2)) sig2 <- sigma2(Lt,Ly)


    if(is.null(mu)) mu <- meanfunc(Lt,Ly,method=method)

    if(datatype == 'irregular')
    {
        vLy <- lapply(1:n,function(i){
            if(is.function(mu)) mui <- mu(Lt[[i]])
            else mui <- predict(mu,Lt[[i]])
            yi <- Ly[[i]] - mui
            yi <- yi^2
            yi
        })
    }
    else
    {
        if(is.function(mu)) vmui <- rep.row(mu(Lt),n)
        else mui <- rep.row(predict(mu,Lt),n)
        vLy <- (Ly - mui)^2
    }

    R <- list(obj=meanfunc(Lt,vLy,weig=weig,method=method,...),sig2=sig2)
    class(R) <- 'varfunc'

    if(is.null(newt))
    {
        newt <- Lt
    }

    R$fitted <- predict(R,newt)

    return(R)
}

#' @export
predict.varfunc <- function(R,newt)
{
    tmp <- predict(R$obj,newt)
    if(is.list(tmp))
    {
        lapply(tmp,function(x){
            raw <- x - R$sig2
            raw[raw < 0] <- 0
            return(raw)
        } )
    }
    else if(is.vector(tmp)){
        raw <- (tmp-R$sig2)
        raw[raw < 0] <- 0
        return(raw)
    }
    else stop('unrecognized/unsupported data type')
}
