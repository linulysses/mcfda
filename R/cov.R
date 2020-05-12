#' Estimate the cov function from functional data/snippets
#' @param t a list of vectors (for irregular design) or a vector (for regular design) containing time points of observations for each individual. Each vector should be in ascending order
#' @param y a list of vectors (for irregular design) or a matrix (for regular design) containing the observed values at \code{t}. If it is a matrix, the columns correspond to the time points in the vector \code{t}
#' @param newt  a list of vectors or a vector containing time points of observations to be evaluated. If NULL, then newt is treated as t
#' @param mu the known or estimated mean function object; it must be a scalar (viewed as a constant function), a function handle, or an object obtained by calling \code{meanfunc}
#' @param method estimation method, 'PACE' or 'FOURIER'
#' @param tuning tuning method to select possible tuning parameters
#' @param ... other parameters required depending on the \code{method} and \code{tuning}; see details
#' @details
#'     \itemize{
#'         \item{When \code{method='PACE'}, additional parameters \code{kernel} and \code{deg} are required. \code{bw} as a scalar is optional. When \code{bw} is provided, the bandwidth is set to \code{bw}}
#'         \item{When \code{method='FOURIER'}, additional parameters \code{q},\code{rho},\code{ext} and \code{domain} are optional. If they are not provided, then they will be deduced from data or selected by the specified \code{tuning} method.}
#'     }
#'
#' @return an object of the class 'covfunc' containing necessary information to predict/evaluate the estimated covariance function
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
#'                    tuning='cv',weig=NULL,kernel='gauss',deg=1)
#' cov.obj <- covfunc(D$t,D$y,newt=NULL,mu=mu.obj,method='FOURIER',
#'                    tuning='cv',weig=NULL,domain=c(0,1))
#' cov.hat <- predict(cov.obj,regular.grid())
#' @export covfunc
#' @useDynLib mcfda
covfunc <- function( t,y,newt=NULL,mu=NULL,
                      method=c('FOURIER','PACE'),...)
{
    method <- match.arg(method)
    R <- NULL
    if(is.list(t)) # irregular design
    {
        if(!is.list(y)) stop('y must be list when t is list')
        if(length(y) != length(t))
            stop('the length of y must match the length of t')
        if(method=='PACE')
            return(cov.pace(t,y,mu=mu,newt=newt,...))
        else
            return(cov.basis(t,y,mu=mu,newt=newt,...))
    }
    else if(is.vector(t)) # regular design
    {
        if(method=='PACE')
            stop('method=PACE for regular design is not supported yet')
        else
            return(cov.basis(t,y,mu=mu,newt=newt,...))
    }
    else stop('t must be a list of vectors (for irregular design) or a vector (for regular design')

    class(R) <- 'covfunc'
    return(R)
}



#' predict cov functions at new locations
#' @param covobj the object obtained by calling \code{covfunc}
#' @param newt a vector or a list of vectors of real numbers
#' @return the estimated mean function evaluated at \code{newt}. It has the same format of \code{newt}
#' @export
predict.covfunc <- function(covobj,newt)
{
    if(toupper(covobj$method) == 'PACE')
    {

        pred <- function(newx)
        {
            R <- Lwls2D(bw=covobj$bw,
                                     kern=covobj$kernel,
                                     xin=covobj$x,
                                     yin=covobj$y,
                                     xout1=newx,
                                     xout2=newx,
                                     win=covobj$weight,
                            delta=covobj$delta)

        }

    }
    else if(toupper(covobj$method) == 'FOURIER')
    {
        pred <- function(newt)
        {
            stopifnot(is.vector(newt))
            domain <- covobj$domain
            if(domain[1]!=0 || domain[2]!=1)
            {
                newt <- (newt-domain[1])/(domain[2]-domain[1])
            }

            p <- dim(covobj$C)[1]
            if(covobj$ext > 0) newt <- shrink(newt,covobj$ext)
            B <- evaluate.basis(p,grid=newt)
            return(B %*% covobj$C %*% t(B))
        }

    }
    else if(covobj$method == 'SP')
    {
        pred <- function(newt) # newt is a vector
        {
            stopifnot(is.vector(newt))

            corr.est <- covobj$rho(newt,newt)
            var.hat <- predict(covobj$varf,newt)
            sig.hat <- sqrt(var.hat)
            cov.fitted <- corr.est * (sig.hat %*% t(sig.hat))
            return(cov.fitted)
        }


    }
    else stop(paste0('Method ',covobj$method, ' is not recognized.'))


    if(is.list(newt))
    {
        return(lapply(newt, function(newx) pred(newx)))
    }
    else if(is.vector(newt))
    {
        return(pred(newt))
    }
    else stop('newt must be a vector or a list of vectors of real numbers')
}



#' estimate the window width of snippets
#' @export
estimate.delta <- function(Lt)
{

    if(is.list(Lt))
    {
        tmp <- lapply(Lt, function(v) max(v)-min(v))
        return(max(unlist(tmp)))
    }
    else if(is.vector(Lt)) return(max(Lt)-min(Lt))
    else stop('unsupported type of Lt')
}

