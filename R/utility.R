get.required.param <- function(name,params)
{
    i <- match(name, names(params))
    if (is.na(i)) {
        stop(paste0('The parameter ', name, ' is required.'))
    }
    else
    {
        params[[i]]
    }
}

get.optional.param <- function(name,params,default)
{
    i <- match(name, names(params))
    if (is.na(i)) {
        default
    }
    else
    {
        params[[i]]
    }
}

#' generate a regular grid on a domain
#' @param m the number of points
#' @param domain the domain
#' @param h the margin of the first and last points to the boundaries of the domain
#' @return a vector representing the grid
#' @examples
#' grid <- regular.grid()
#' grid <- regular.grid(50)
#' @export regular.grid
regular.grid <- function(m=100,domain=c(0,1),h=1/(2*m))
{
    seq(domain[1]+h,domain[2]-h,length.out=m)
}

int <- function(x,y)
{
    pracma::trapz(x,y)
}

norm.L2 <- function(x,t=NULL,domain=NULL)
{
    if(is.null(t) && is.null(domain)) stop('one of t and domain must be supplied.')
    if(is.null(t)) t <- regular.grid(domain,length(x))
    sqrt(int(t,x^2))
}

#' replicate a row vector into a matrix
#' @param x the row vector to be replicated
#' @param n replicate the row vector n times
#' @return a matrix of the dimension of \code{n} by \code{length(x)}
#' @export rep.row
rep.row <- function(x,n){
    matrix(rep(x,each=n),nrow=n)
}

#' replicate a column vector into a matrix
#' @param x the column vector to be replicated
#' @param n replicate the column vector n times
#' @return a matrix of the dimension of \code{length(x)} by \code{n}
#' @export rep.col
rep.col <- function(x,n){
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
